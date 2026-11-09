#!/usr/bin/env python3
"""S6 质量门禁：对关联结果做 go / no-go 检查，输出 Markdown 报告。

用法：
  python3 05_checks.py --input assoc/height_lmm.assoc.tsv \
      --pheno plant_height --config project_config.yaml --out gwas_report.md

检查项：必需列、p 值合法性、有效标记数、λ、显著位点数、NMISS 一致性、
        样本量/功效提示。退出码 1 表示存在阻断项（BLOCK），0 表示通过或仅有警告。
"""
import argparse
import math
import os
from statistics import NormalDist

import numpy as np

from assoc_io import read_association
from thresholds import resolve_threshold

ND = NormalDist()
CHI2_1DF_MEDIAN = 0.454936423119572
REQUIRED_COLS = ["CHR", "BP", "P"]


def read_assoc(path):
    rows, columns = read_association(path)
    return rows, columns, columns.chrom, columns.pos, columns.p, columns.n, columns.effect


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--pheno", default="gwas")
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default="gwas_report.md")
    ap.add_argument("--expected-n", type=int, default=0)
    ap.add_argument("--threshold", type=float, default=None)
    ap.add_argument("--positive-control", default="")
    ap.add_argument("--data-mode", choices=("full", "truncated"), required=True,
                    help="Explicitly declare whether the association table is complete or threshold-truncated")
    a = ap.parse_args()

    cfg = {}
    if a.config:
        try:
            import yaml
            cfg = yaml.safe_load(open(a.config, encoding="utf-8")) or {}
        except Exception as e:
            print("配置读取失败：", e)
    dg = cfg.get("diagnostics", {})
    lam_warn = float(dg.get("lambda_warn", 1.2))
    expected_n = a.expected_n or int(cfg.get("population", {}).get("n_samples_expected", 0) or 0)
    threshold_decision = resolve_threshold(a.config, a.threshold)
    a.threshold = threshold_decision.value

    rows, cols, cc, bc, pc, nc, ec = read_assoc(a.input)
    checks = []

    def add(name, status, detail):
        checks.append((name, status, detail))

    # 1 必需列
    missing = [c for c, v in zip(REQUIRED_COLS, [cc, bc, pc]) if not v]
    add("必需列 CHR/BP/P", "PASS" if not missing else "BLOCK",
        "齐全" if not missing else "缺少：" + ", ".join(missing))
    add("效应量列 BETA/OR", "PASS" if ec else "WARN",
        "存在" if ec else "缺失：必须同时报告效应量与 p 值（书中 p101）")
    if not ec or missing:
        # 仍继续做能做的检查
        pass

    # 2 p 值
    ps = []
    bad = 0
    for r in rows:
        try:
            v = float(r[pc])
            if 0 < v <= 1:
                ps.append(v)
            else:
                bad += 1
        except (TypeError, ValueError):
            bad += 1
    ps = np.array(ps)
    add("p 值合法且非空", "PASS" if len(ps) > 0 and bad == 0 else
        ("WARN" if len(ps) > 0 else "BLOCK"),
        f"有效 {len(ps)} 个，异常/缺失 {bad} 个")

    # 3 数据模式与 λ
    add("数据模式", "PASS", f"data_mode={a.data_mode}")
    if a.data_mode == "full" and len(ps):
        z = np.array([ND.inv_cdf(1 - x / 2) for x in np.clip(ps, 1e-12, 1)])
        lam = float(np.median(z ** 2) / CHI2_1DF_MEDIAN)
        if lam > lam_warn:
            add("基因组膨胀因子 λ", "WARN",
                f"λ={lam:.3f} > {lam_warn}：回退 S3 检查 PC 数 / K 矩阵 / 批次协变量；"
                "注意高度多基因性状本身也会抬高 λ（书中 p340）")
        elif lam < 0.85:
            add("基因组膨胀因子 λ", "WARN",
                f"λ={lam:.3f} 偏低：可能过度校正或样本量过小")
        else:
            add("基因组膨胀因子 λ", "PASS", f"λ={lam:.3f}，处于合理区间")
    elif a.data_mode == "truncated":
        lam = float("nan")
        add("基因组膨胀因子 λ", "PASS", "not evaluated for truncated input")
    else:
        lam = float("nan")

    # 4 显著位点
    nsig = int((ps < a.threshold).sum()) if len(ps) else 0
    add("显著位点", "PASS" if nsig > 0 else "WARN",
        f"p<{a.threshold:g} 的位点 {nsig} 个；"
        "显著位点需成簇且与 LD 一致，孤立单点应回查分型质量")

    # 5 NMISS 一致性
    if nc:
        vals = []
        for r in rows:
            try:
                vals.append(int(float(r[nc])))
            except (TypeError, ValueError):
                pass
        if vals:
            spread = max(vals) - min(vals)
            add("NMISS 一致性", "PASS" if spread <= max(2, 0.05 * max(vals)) else "WARN",
                f"范围 {min(vals)}–{max(vals)}，跨度 {spread}；跨度过大提示表型缺失或样本错配")
    else:
        add("NMISS 一致性", "WARN", "无 NMISS 列，无法检查样本一致性")

    # 6 功效提示
    if expected_n and expected_n < 200:
        add("样本量 / 功效", "WARN",
            f"N={expected_n} 偏小：结论应定为候选区间 + 可验证假设，"
            "不宜宣称全基因组显著（书中 p95、p105 强调微效需大样本）")
    elif expected_n:
        add("样本量 / 功效", "PASS", f"N={expected_n}（仍需结合效应量评估功效）")

    # 7 阳性对照
    add("阳性对照复现", "WARN" if not a.positive_control else "PASS",
        a.positive_control or "未提供已知 QTL/基因复现结果；建议补充")

    # 8 结论边界声明
    add("结论边界", "PASS", "报告须写明：关联 ≠ 因果；效应不可跨群体/环境直接外推"
        "（书中 p59、p127）；显著 SNP ≠ 候选基因")

    blocked = [c for c in checks if c[1] == "BLOCK"]
    warned = [c for c in checks if c[1] == "WARN"]
    verdict = "NO-GO（存在阻断项，需修复后重跑）" if blocked else \
              ("GO（有警告，需在报告中说明）" if warned else "GO")

    lines = [
        f"# GWAS 质量门禁报告 · {a.pheno}",
        "",
        f"**结论：{verdict}**",
        "",
        f"- 输入：`{a.input}`",
        f"- 配置：`{a.config or '未提供'}`",
        f"- data_mode={a.data_mode}",
        f"- 有效标记数：{len(ps)}",
        f"- λ = {lam:.4f}" if a.data_mode == "full" and not math.isnan(lam) else "- λ: not evaluated for truncated input",
        f"- 阈值：p < {a.threshold:g}（method={threshold_decision.method}; threshold={a.threshold:.6g}; "
        f"meff={threshold_decision.meff or 'NA'}; source={threshold_decision.source}）",
        f"- 显著位点数：{nsig}",
        "",
        "| 检查项 | 状态 | 说明 |",
        "| --- | --- | --- |",
    ]
    lines += [f"| {n} | {s} | {d} |" for n, s, d in checks]
    lines += [
        "",
        "## 后续动作",
        "",
        "- λ 异常或 Q-Q 图整体膨胀 → 回退 S3（PC 数、K 矩阵、批次协变量），不得放宽阈值掩盖。",
        "- 显著位点按「领先 SNP + 局部 LD 块」定义候选区间，宽度取物种 LD 衰减距离。",
        "- 完成注释后按：独立群体复现 → 已知基因/突变体 → 表达证据 → 基因编辑 的优先序验证。",
        "",
        "> 依据：Mills, Barban & Tropf (2020), *An Introduction to Statistical Genetic Data Analysis*, MIT Press。",
    ]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print(f"门禁结论：{verdict}  |  报告：{a.out}")
    for n, s, d in checks:
        print(f"  [{s}] {n}: {d}")
    raise SystemExit(1 if blocked else 0)


if __name__ == "__main__":
    main()
