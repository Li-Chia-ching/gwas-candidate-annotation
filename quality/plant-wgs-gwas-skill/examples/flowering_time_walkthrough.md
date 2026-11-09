# 迁移范例：从「株高」切到「开花时间」

> 场景：课题组 A 同学已经用本 Skill 跑完株高 GWAS；B 同学研究水稻开花时间，想复用同一套流程，但**不希望改动任何分析逻辑**。本示例展示最短迁移路径。

## 1. 只改配置，不改流程

B 同学复制已冻结的配置，把 `phenotypes` 换成开花时间，其余（QC 后基因型、PC、K 矩阵）全部沿用 A 同学的产物。

```yaml
project:
  name: "rice_flowering_time_gwas"
  owner: "B 同学"
  note: "复用 A 同学株高 GWAS 的 clean 基因型 / PC / K 矩阵"

# —— 以下整段沿用 A 同学，不动 ——
species:      { ... }        # 同一物种、同一参考基因组
population:   { ... }
paths:
  workdir: "./gwas_run"      # ★ 指向 A 同学同一工作目录，复用 qc/ 与 structure/
genotype:     { ... }
qc:           { ... }
ld_pruning:   { ... }
structure:    { ... }
model:        { ... }

# —— 只改这里 ——
phenotypes:
  - name: "flowering_days"
    file: "./data/phenotypes/flowering_blup.tsv"
    column: "BLUP_days"
    type: "quantitative"
    transform: "log"          # 开花时间常右偏；先画分布再决定是否保留
    blup: true                # 多环境 BLUP，已剥离年份/区组效应
    covariates: ["sowing_date", "vernalization", "photoperiod"]
    model: "LMM"              # 自然关联群体，结构强 → 混合线性模型
    note: "强 G×E；另做分环境 GWAS 比较位点一致性"

threshold:
  method: "bonferroni_meff"   # WGS 标记密，用 0.05/M_eff；M_eff 取 LD pruning 后独立标记数
  suggestive: 1.0e-5

annotation:
  candidate_window_kb: 200    # 取该物种 LD 衰减距离
  use_known_qtl_overlap: true
```

**为什么不能重新跑质控**：不同表型若各自质控，基因型矩阵、PC、K 都会变，两个方向的显著位点就不可比了。这是 `SKILL.md` 的硬约束。

## 2. 表型文件

`data/phenotypes/flowering_blup.tsv`（`IID` 必须与 VCF 样本名一致）：

```
FID	IID	BLUP_days	sowing_date	vernalization	photoperiod
S001	S001	78.2	2024-05-10	0	long
S002	S002	81.0	2024-05-10	0	long
S003	S003	NA	2024-05-12	1	short
```

多环境 BLUP 建议用 `lme4` / `sommer` 拟合 `y ~ G + E + G×E + (1|block%E)` 后取基因型 BLUP；把 BLUP 可靠性写进报告。

## 3. 运行

```bash
# 配置校验
python3 scripts/config.py my_config.yaml validate

# 复用已有 qc/ 与 structure/，因此只跑 S5 及之后
bash scripts/03_gwas.sh my_config.yaml
python3 scripts/04_plots.py --config my_config.yaml \
    --input gwas_run/assoc/flowering_days_lmm.assoc.tsv \
    --pheno flowering_days --out gwas_run/figures/flowering_days --data-mode full
python3 scripts/05_checks.py --config my_config.yaml \
    --input gwas_run/assoc/flowering_days_lmm.assoc.tsv \
    --pheno flowering_days --out gwas_run/gwas_report_flowering.md --data-mode full
```

若 B 同学的材料与 A 同学不同（换了群体/测序批次），则必须完整重跑 `01_qc.sh` → `02_structure.sh`，不能复用。

## 4. 这个方向特有的注意点

| 问题 | 处理 |
| --- | --- |
| 光温敏感、强 G×E | 除汇总 BLUP 外，再按环境各跑一次 GWAS，比较领先位点是否一致；报告写明 |
| 分布右偏或双峰 | 先画直方图与 Q-Q 图；`log` / Box-Cox 后复查；双峰可能提示主效基因，可考虑阈值化为早/迟二分类做 logistic |
| 已知 main-effect 基因 | 把已知 FT / Hd 类基因位置作为**阳性对照**：若 GWAS 完全不复现，先查坐标版本与表型定义，而不是换模型 |
| 与播期共变 | 播期必须进协变量；若播期与亚群/批次相关，检查共线性 |
| 多表型相关 | 开花时间常与株高、产量遗传相关，可用 MTAG 联合分析或 LDSC 估计遗传相关（书中第 12 章） |

## 5. 如何解读输出

1. `qc_report.md` / `structure_report.md`：确认沿用的是同一套 clean/PC/K。
2. `manhattan.png`：看是否有越过阈值的峰；峰是否成片（LD 块），而非孤立单点。
3. `qq.png` 与 λ：λ 在 1.0–1.2 且仅尾部上翘较理想；若整体膨胀，回退 S3。
4. `gwas_report_flowering.md`：门禁结论、警告项、后续动作。
5. 候选区间：以「领先 SNP ± species.ld_decay_kb」定义，与已知 QTL、注释 GFF、表达数据交叉。

## 6. 报告里必须写的边界

- 关联 ≠ 因果；显著 SNP 不等于候选基因。
- 本群体、本环境组合的结论不可直接外推到其他亚群或其他生态区。
- N 偏小时明确功效限制，改述为「候选区间 + 可验证假设」。
- 效应量与 p 值同时给出，阈值写明依据（本例为 Bonferroni 0.05/M_eff）。

## 7. 再换一个方向要改什么

| 新方向 | 主要改动 |
| --- | --- |
| 分支数 / 分蘖数 | `type: count`、`transform: sqrt`、协变量加密度/养分 |
| 抗病等级 | `type: binary`（阈值化）、`model: LMM`、协变量加接种批次/菌系、结果报 OR |
| 代谢物 mGWAS | 每个代谢物一条 `phenotypes` 记录、`transform: quantile_normal`、`threshold` 更严 |
| 换物种 | 改 `species`、`ld_pruning`、`annotation`、完整重跑 S1–S3 |

**一句话**：流程、阈值逻辑、门禁都不动；只改 `phenotypes` 段和（必要时）`threshold`、`annotation`。
