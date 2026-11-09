#!/usr/bin/env python3
"""S6 绘图：曼哈顿图 + Q-Q 图 + 基因组膨胀因子 λ。

用法：
  python3 04_plots.py --input assoc/height_lmm.assoc.tsv --pheno plant_height \
      --out figures/ [--config project_config.yaml] [--threshold 5e-8]

输入为制表符或逗号分隔的关联结果，列名自动兼容：
  CHR / chrom / chromosome, BP / pos / position, P / p / PVALUE / pvalue
仅依赖 numpy + matplotlib（λ 计算用标准库，不需要 scipy）。
中文标签会自动探测系统 CJK 字体；若不可用则回退为英文，避免输出方框乱码。
"""
import argparse
import math
import os
from statistics import NormalDist

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

from assoc_io import read_association
from thresholds import resolve_threshold

ND = NormalDist()
CHI2_1DF_MEDIAN = 0.454936423119572  # chi2(1) 的中位数，用于 genomic control


def setup_font():
    """探测常见中文字体；返回是否可用。跨平台零配置。"""
    candidates = [
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJKsc-Regular.otf",
        "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
        "/usr/share/fonts/truetype/arphic/uming.ttc",
        "/System/Library/Fonts/PingFang.ttc",
        "/System/Library/Fonts/STHeiti Light.ttc",
        "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/simsun.ttc",
    ]
    for c in candidates:
        if os.path.exists(c):
            try:
                font_manager.fontManager.addfont(c)
                plt.rcParams["font.family"] = font_manager.FontProperties(fname=c).get_name()
                return True
            except Exception:
                continue
    return False


HAS_CJK = setup_font()
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"


def t(zh, en):
    return zh if HAS_CJK else en


def read_assoc(path):
    rows, columns = read_association(path)
    chrs, bps, ps = [], [], []
    for r in rows:
        try:
            p = float(r[columns.p])
        except (TypeError, ValueError):
            continue
        if not (0 < p <= 1):
            continue
        try:
            bp = int(float(r[columns.pos]))
        except (TypeError, ValueError):
            raise SystemExit(f"invalid BP/POS value: {r.get(columns.pos)!r}")
        chrs.append(str(r[columns.chrom]).replace("chr", ""))
        bps.append(bp)
        ps.append(p)
    return np.array(chrs), np.array(bps), np.array(ps)


def lambda_gc(pvals):
    """λ = median(chi2_obs) / median(chi2_1df)，chi2 由双尾 p 值反推。"""
    # inv_cdf 在参数极接近 1 时会溢出，clip 到 1e-12 足够（λ 取中位数，不受极端 p 影响）
    p = np.clip(np.asarray(pvals, float), 1e-12, 1.0)
    z = np.array([ND.inv_cdf(1 - x / 2) for x in p])
    return float(np.median(z ** 2) / CHI2_1DF_MEDIAN)


def manhattan(chrs, bps, ps, out, threshold):
    order = sorted(set(chrs.tolist()), key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
    offset, ticks = 0, []
    colors = ["#2E5C8A", "#7FA6C9"]
    fig, ax = plt.subplots(figsize=(11, 4.2))
    for i, c in enumerate(order):
        m = chrs == c
        if not m.any():
            continue
        x = bps[m] - bps[m].min() + offset
        y = -np.log10(np.clip(ps[m], 1e-300, 1))
        ax.scatter(x, y, s=3.0, c=colors[i % 2], linewidths=0, rasterized=True)
        ticks.append(offset + (x.max() - x.min()) / 2)
        offset = x.max() + 1
    ax.axhline(-math.log10(threshold), color="#C0392B", ls="--", lw=1.4,
               label=f"{t('显著性阈值', 'Significance threshold')}  p < {threshold:g}")
    ax.set_xticks(ticks); ax.set_xticklabels(order, fontsize=8)
    ax.set_xlabel(t("基因组位置（染色体）", "Genomic position (chromosome)"))
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(t("曼哈顿图", "Manhattan plot"))
    ax.legend(fontsize=9); ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(out, "manhattan.png"), dpi=140); plt.close(fig)


def qq(ps, out, lam):
    n = len(ps)
    obs = -np.log10(np.sort(np.clip(ps, 1e-300, 1)))[::-1]
    exp = (-np.log10(np.arange(1, n + 1) / (n + 1)))[::-1]
    fig, ax = plt.subplots(figsize=(5.8, 5.4))
    hi = min(float(exp.max()), 8.0)
    ax.plot([0, hi], [0, hi], color="#C0392B", ls="--", lw=1.4,
            label=t("理论期望线", "Expected under null"))
    ax.plot(exp, obs, color="#2E5C8A", lw=1.8,
            label=f"{t('观测', 'Observed')}  λ = {lam:.3f}")
    ax.set_xlim(0, hi); ax.set_ylim(0, max(float(obs.max()) * 1.05, hi))
    ax.set_xlabel(t("期望  ", "Expected  ") + r"$-\log_{10}(p)$")
    ax.set_ylabel(t("观测  ", "Observed  ") + r"$-\log_{10}(p)$")
    ax.set_title(t("Q-Q 图", "Q-Q plot"))
    ax.legend(fontsize=9); ax.grid(alpha=0.25, lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(out, "qq.png"), dpi=140); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--pheno", default="gwas")
    ap.add_argument("--out", default="./figures")
    ap.add_argument("--config", default=None)
    ap.add_argument("--threshold", type=float, default=None)
    ap.add_argument("--data-mode", choices=("full", "truncated"), required=True,
                    help="Explicitly declare whether the association table is complete or threshold-truncated")
    a = ap.parse_args()

    decision = resolve_threshold(a.config, a.threshold)
    a.threshold = decision.value

    os.makedirs(a.out, exist_ok=True)
    chrs, bps, ps = read_assoc(a.input)
    manhattan(chrs, bps, ps, a.out, a.threshold)
    if a.data_mode == "full":
        lam = lambda_gc(ps)
        qq(ps, a.out, lam)
    else:
        lam = float("nan")
        stale_qq = os.path.join(a.out, "qq.png")
        if os.path.exists(stale_qq):
            os.remove(stale_qq)
    nsig = int((ps < a.threshold).sum())
    lam_text = f"{lam:.4f}" if a.data_mode == "full" else "NA"
    print(f"[{a.pheno}] 标记数={len(ps)}  λ={lam_text}  显著位点(p<{a.threshold:g})={nsig}")
    print(f"[{a.pheno}] data_mode={a.data_mode} lambda={lam_text}")
    print(f"[{a.pheno}] threshold={a.threshold:.6g} method={decision.method} source={decision.source} meff={decision.meff or 'NA'}")
    print(f"[{a.pheno}] 图已保存：{a.out}/manhattan.png, {a.out}/qq.png")
    print(f"[{a.pheno}] 图内标签：{'中文' if HAS_CJK else '英文（未检测到 CJK 字体）'}")
    if lam > 1.2:
        print(f"[警告] λ={lam:.3f} > 1.2：请回退 S3 检查 PC 数、K 矩阵与批次协变量；"
              "同时注意高度多基因性状本身也会抬高 λ。")


if __name__ == "__main__":
    main()
