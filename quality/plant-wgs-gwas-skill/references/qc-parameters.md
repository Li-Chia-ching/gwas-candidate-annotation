# 质控、多重检验与结果诊断参数

> 书中默认值来自 Mills, Barban & Tropf (2020) 第 4、8、9、12 章；「植物适配」为迁移补充。命令以 PLINK 2.0 为例，PLINK 1.9 需微调。

## 1. 为什么 QC 是推断前提

书中指出 GWAS 需要"海量统计推断 + 大量质量检查（quality control）"（PDF p107）：个体或标记缺失、call rate 组间差异、分型错误都会制造系统假关联。第 8 章把 QC 分成**个体层**与**标记层**（PDF p219–231）。

## 2. 参数速查表

| 检查项 | 命令 | 书中默认 | 植物 WGS 适配 |
| --- | --- | --- | --- |
| 个体缺失率 | `--mind` | 0.05（PDF p220） | 0.05–0.10；低深度样本 0.05 |
| 个体杂合度 | `--het` | 均值 ±3 SD 离群剔除（PDF p220–221） | 自交系高杂合常为污染/混样，须查；不要一刀切剔除 F1/异交材料 |
| 性别校验 | `--check-sex` | X 染色体杂合判缺失（PDF p222） | 植物多无性染色体 → 改查倍性一致性、细胞器污染、样本混样 |
| 重复/近亲 | `--distance` / IBS | IBS > 0.1875 剔除其一（PDF p223） | 同名不同批次要合并验证；已知家系**保留并用 K 矩阵建模** |
| 标记缺失率 | `--geno` | 0.05，call rate > 95%（PDF p98、p224） | 0.05–0.10；低覆盖区域放宽需记录 |
| MAF | `--maf` | N=100,000→0.01；N=10,000→0.05（PDF p225） | N<300 建议 ≥0.05；N>1000 可 0.01；稀有变异另做 burden/SKAT |
| HWE | `--hwe` | 数量性状 p<1e-6（PDF p225）；练习 <0.0001（p231） | **自交物种/强结构群体默认关闭**；优先用 DP/GQ/等位基因平衡 |
| WGS 分型质量 | DP / GQ / AB | 芯片数据无此项 | 设 min DP、min GQ、杂合 AB 0.25–0.75 |
| 差异缺失 | 分层 call rate | 病例/对照间差异（PDF p224） | 按**测序批次 / 亚群**分层检查 |
| 填补质量 | INFO / R² | 剔除低填补质量位点（PDF p108、p173–175） | WGS 一般不需填补；若做了填补须过滤低质量位点 |

**书中逻辑提醒**：HWE 偏离在人类 GWAS 中被假定主要来自分型错误（PDF p225），但这一假定在自交植物中不成立。HWE 应作为**诊断**而非**硬过滤**。

## 3. 参考命令序列（PLINK2）

```bash
# 0) VCF → binary；保留 ID、允许 scaffold/细胞器 contig
plink2 --vcf "$VCF" --double-id --allow-extra-chr \
       --snps-only just-acgt --max-alleles 2 --make-bed --out raw

# 1) 个体层：缺失率 + 杂合度
plink2 --bfile raw --mind 0.05 --make-bed --out step1_ind
plink2 --bfile step1_ind --het --out het
#   在 R/Python 中按均值 ±3 SD 标记离群 → --remove <ids> --make-bed

# 2) 重复 / 近亲（只用 LD-pruned 独立标记，见下）
plink2 --bfile step1_ind --extract pruned.prune.in --distance square --out ibs
#   IBS > 0.1875 的成对样本保留其一（家系结构除外）

# 3) 标记层：call rate + MAF；HWE 视群体类型决定开关
plink2 --bfile step1_ind --geno 0.05 --maf 0.05 --make-bed --out clean
#   异交/自然群体且确认分型错误为主时，可追加：--hwe 1e-6

# 4) WGS 质量过滤建议在 VCF 阶段用 bcftools 完成（PLINK 不读 DP/GQ）
#   bcftools filter -e 'DP<8 || GQ<20' -Oz -o filtered.vcf.gz genotypes.vcf.gz
```

每一步都要把「过滤前个体数 / 标记数 → 过滤后 → 剔除原因」写入 `qc_report.md`。原书还建议 P–Z 图与 Q–Q 图作为额外诊断（PDF p228）。

## 4. LD pruning 与独立标记集

LD 使相邻标记高度相关；若用全部标记做 PCA/GRM，LD 块会主导结果（PDF p241–243）。

```bash
# 常用：窗口 50 SNP、步长 5、r² 阈值 0.2
plink2 --bfile clean --indep-pairwise 50 5 0.2 --out pruned
# 后续 PCA / K 矩阵只用 pruned.prune.in
```

- 自交作物 LD 衰减慢（可达数百 kb），窗口应放大；异交作物衰减快，可用更严的 r²。
- 同一套独立标记集要同时供 PCA 与 GRM 使用。
- LD clumping（按 p 值选每个 LD 块的代表 SNP，PDF p243）用于结果层面的信号归并，不要与 pruning 混用。

## 5. PCA 与亲缘矩阵

```bash
plink2 --bfile clean --extract pruned.prune.in --pca 10 --out pca        # 书中常用 10 或 20 PC (p245)
plink2 --bfile clean --extract pruned.prune.in --make-king square --out king
# 或 GEMMA：gemma -bfile clean -gk 1 -o kin
# 或 GCTA：gcta64 --bfile clean --make-grm --out grm   (p255–256)
```

PC 数不要固定照抄：结合碎石图、PC 与已知亚群/测序批次的相关性决定（Figure 9.2，PDF p247）。

## 6. 多重检验

| 方法 | 取值 | 依据 |
| --- | --- | --- |
| 固定阈值 | `p < 5e-8` 全基因组显著；`5e-6`（或 `1e-5`）提示性 | 书中 p101–102、p239 |
| Bonferroni（WGS 推荐） | `0.05 / M_eff`，M_eff 取 LD pruning 后独立标记数 | 由书中"100 万 SNP × 0.05 → 5 万假阳性"外推 |
| FDR | q=0.05 | 书中 p106 提及 FDR |
| 置换检验 | 经验阈值 | 书中 p116、p78 提及经验 p 值 |

必须在报告中写明用的是哪种、为什么。不要用 p<0.05。

## 7. 结果诊断

- **曼哈顿图**（PDF p102、p333–336）：x=染色体位置，y=−log₁₀(p)，加显著线；显著峰成片是 LD 所致，非多个独立信号。
- **Q-Q 图与 λ**（PDF p337–340）：
  - λ ≈ 1：分层已充分校正；
  - λ > 1.2：提示残余群体结构（书中案例 λ=1.239）；
  - 但高度多基因性状（如株高、开花时间）本身会让 λ>1.2，**须结合 QQ 图形态与 PC-表型相关性判断**；
  - 可用 genomic control 把检验统计量除以 λ 校正，但不能替代重新检查模型。
- **区域关联图**（LocusZoom，PDF p337）：看显著位点周边的 LD 与重组率，定候选区间。
- **遗传力参照**：GCTA GREML（p256–257）或 LDSC（p341）估计 SNP 遗传力，用于判断"该性状在本群体是否值得做 GWAS"及结果合理性。

## 8. 门禁不通过时的回退

| 现象 | 回退 | 动作 |
| --- | --- | --- |
| λ>1.2 且 QQ 整体上翘 | S3 | 调整 PC 数、换/重建 K 矩阵、补批次协变量 |
| 全基因组均匀抬高 | S2/S3 | 查差异缺失、批次、样本混样 |
| K 矩阵非正定 / 不收敛 | S3 | 去近亲对、增加标记数、检查重复样本 |
| 阳性对照不复现 | S1/S4 | 查坐标版本、表型列、样本 ID 对齐 |
| 孤立显著单点无 LD 支持 | S2 | 查该位点 DP/GQ/等位基因平衡 |

**禁止**：通过放宽显著性阈值来"解决"膨胀问题。
