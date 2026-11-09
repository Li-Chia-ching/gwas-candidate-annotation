# 原书依据索引与迁移说明

## 来源

Mills, M. C., Barban, N., & Tropf, F. C. (2020). *An Introduction to Statistical Genetic Data Analysis*. The MIT Press. ISBN 9780262538381.

扫描版共 433 页，页面为 4200×5400 图像，无文本层；本 Skill 中的引用页码为该扫描版的 **PDF 页码**（书内印刷页码 ≈ PDF 页码 − 17）。页码经 OCR 全文核对。

## 章节 → 方法要点 → PDF 页码

| 章节 | 方法要点 | PDF 页 | 本 Skill 采用位置 |
| --- | --- | --- | --- |
| 第 1 章 | SNP、MAF、单基因/多基因/全基因效应；遗传力与遗传架构 | p20–31、p43 | 立项可行性、功效评估 |
| 第 2 章 | 回归、固定/随机/混合模型；复现与过拟合 | p50–51 | `model-selection.md` |
| 第 2 章 | Q-Q 图原理（Figure 2.2） | p53–54 | `qc-parameters.md`、图文简介 |
| 第 3 章 | 群体结构与分层、PCA（Figure 3.2）；HWE；LD 与单体型块 | p75–90 | `qc-parameters.md` §4–5 |
| 第 4 章 | GWAS 研究设计与 meta 分析流程（Figure 4.1） | p97 | 技术路线主框架 |
| 第 4 章 | 可行性、call rate >95%、分析计划、协变量 | p96–99 | `technical-route.md` §S0/S1 |
| 第 4 章 | p 值、Z 值、多重检验校正、`p<5e-8` | p100–102、p239 | `qc-parameters.md` §6 |
| 第 4 章 | 曼哈顿图（Figure 4.2）；二分类 vs 数量性状 | p102–104 | 图文简介、`qc-parameters.md` |
| 第 4 章 | 固定 vs 随机效应；FDR、加权、填补 | p105–106 | `model-selection.md` |
| 第 4 章 | QC 的重要性与差异 call rate | p107–108 | `qc-parameters.md` §1 |
| 第 7 章 | 基因型填补原理与质量 | p172–176 | 背景（WGS 一般不需填补） |
| 第 7 章 | VCF 格式；PLINK 二进制格式 | p191、p170–171、p205–208 | `technical-route.md` §S1 |
| 第 8 章 | 个体层 QC：`--mind`、`--het`、`--check-sex`、IBS>0.1875 | p219–223 | `qc-parameters.md` §2–3 |
| 第 8 章 | 标记层 QC：`--geno`、`--maf`、`--hwe`；P–Z 与 Q-Q 诊断 | p224–231 | `qc-parameters.md` §2–3 |
| 第 9 章 | 关联分析 `--assoc`/`--linear`/`--logistic`；加性/显性/隐性；`--covar` | p236–239 | `model-selection.md` |
| 第 9 章 | 家系内关联 `--qfam`；多重检验 5e-8 | p239–240 | `model-selection.md` |
| 第 9 章 | LD pruning / clumping | p241–243 | `qc-parameters.md` §4 |
| 第 9 章 | PCA `--pca 10`、前 10/20 PC；Figure 9.2 | p243–248 | `qc-parameters.md` §5 |
| 第 9 章 | LMM 取代单纯 PCA 校正 | p248 | `model-selection.md` 决策树 |
| 第 9 章 | IBS `--distance`、`--make-rel`、GRM `--make-grm`、`--grm-cutoff` | p253–256 | `qc-parameters.md` §5 |
| 第 9 章 | GREML 估计 SNP 遗传力 | p256–257 | `technical-route.md` §S3 |
| 第 12 章 | 曼哈顿图与高亮（Figure 12.1/12.2）；区域关联图 | p332–337 | 图文简介、`qc-parameters.md` §7 |
| 第 12 章 | Q-Q 图与 λ（genomic control），案例 λ=1.239，λ>1.2 警示 | p337–340 | `qc-parameters.md` §7 |
| 第 12 章 | LDSC 遗传力、遗传相关、MTAG 多性状 | p341 起 | `model-selection.md`、表型手册 |
| 第 5、14 章 | 遗传分数不可跨祖先群体迁移；多样性与伦理 | p59、p110、p127、p367 | 结论边界与硬约束 |

## 关键原句（OCR 摘录，供核对）

- 「the agreed upon genome-wide significant threshold is p<5×10⁻⁸」（p101）
- 「p<5×10⁻⁸ ... standard for genome-wide statistical significance and p<5×10⁻⁶ often used to represent suggestive」（p102）
- 「If we take the normal p-value threshold of 0.05, we expect that ... a p-value of 0.05 implies that 50,000 SNPs are false positives ... we adopt a much stricter p-value threshold (5×10⁻⁸)」（p239）
- 「It is common to use the first 10 or 20 principal components from a genetic dataset」（p245）
- 「an approach that has increasingly been adopted in association studies is to use linear mixed models instead of PCA to control for population stratification」（p248）
- 「A value close to 1 suggests that data have been properly adjusted ... If λ>1.2, this suggests the presence of stratification ... a λ of 1.239」（p340）
- 「Typically, pairs of individuals with IBS > 0.1875 are removed from the analysis」（p223）
- 「N=100,000, MAF threshold 0.01, with a moderate sample N=10,000, MAF threshold 0.05」（p225）

## 从人类 GWAS 迁移到植物 WGS 的五处替换

| 维度 | 原书（人类） | 本 Skill（植物 WGS） |
| --- | --- | --- |
| 标记来源 | 芯片 + 基因型填补 | WGS 直接分型；通常不需填补，但需过滤低 DP/GQ 位点 |
| 群体结构 | 祖先/地理群体 | 亚群、生态型、自交系、育种群体，分层更强 → LMM 默认 |
| 亲缘 | 剔除隐性相关 | 已知谱系/家系结构保留并用 K 矩阵建模 |
| 表型 | 单点测量 | 多环境 → BLUP/BLUE；计数/等级/偏态需变换；G×E 需报告 |
| 阈值 | 芯片量级固定 5e-8 | WGS 标记更密 → Bonferroni 0.05/M_eff、FDR 或经验阈值 |
| 伦理 | 隐私、知情同意、数据共享受限 | 无人类受试者问题，但需记录田间设计与种质来源 |

## 插图使用说明

配套《GWAS 图文入门简介》中：

- 原书 Figure 4.1（GWAS 研究设计与 meta 分析流程图，PDF p97）为**直接截图**，清晰可用。
- 其余散点类插图（Figure 4.2 曼哈顿图、Figure 9.2 PCA、Figure 12.1/12.5 曼哈顿与 Q-Q 图）在该扫描件中**数据点丢失**，只剩坐标轴与图注；因此简介中的曼哈顿图、Q-Q 图、PCA 图、LD 图按上述章节原理用示意数据重绘，并在页面注明。
- 所有插图**仅限课题组内部教学使用**，不对外发布。
