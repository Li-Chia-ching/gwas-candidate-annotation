# 完整工作流与案例（实测提炼）

> ⚠️ 本案例已脱敏：基因以 Candidate-A/B/C 代称，坐标与统计量为真实模式但已去标识化。
> 案例的教学价值在于**数值模式与判读逻辑**，不在于具体基因——复用时只需替换为你自己的数据。

## 0. 输入前提（典型场景）

- GWAS 供应商交付：多性状 × 3 模型（GLM/MLM/FarmCPU）的 `*_signals.csv`（P<1e-5 截断）、`*.ann.vcf`（仅显著位点）、`*.annotation.xls`；**全基因组 VCF 未交付**是常见卡点。
- 参考基因组为分型到单倍型的多倍体组装（本案例：32 条单倍型染色体、约 19 万个基因），GFF 可 gz 压缩。

## 1. Stage 1 实测：双模型筛选

```bash
python scripts/screen_dual_model_leads.py \
    --signals-a trait1.MLM_signals.csv --signals-b trait1.FarmCPU_signals.csv \
    --model-a MLM --model-b FarmCPU \
    --gff reference.gff.gz --trait Trait1 --outdir out/Trait1/
```

实测产出模式（三个数量性状）：双模型独立 lead 分别为 8、3、5 个；检出一个跨性状共定位位点（两性状共享同一位置）。

**判读要点**：
- `n_dual_snps` >1 的位点比单 SNP 位点置信更高；
- 窗口基因表按 `dist_to_lead_bp` 排序，直接注释基因（`overlaps_lead=1`）优先级最高；
- "信号最多的染色体"与"最强单位点"常不在同一染色体——数量密度与主效强度必须分开讲。

## 2. Stage 2 实测：CDS/蛋白提取

```bash
python scripts/extract_cds_protein.py --gff reference.gff.gz \
    --genome genome.fasta.gz \
    --genes CandidateA CandidateB --outdir out/annotation/ --proteome
```

sanity 检查全通过的典型形态：

| 基因 | 链 | CDS 数 | 蛋白长 | 起始/终止 | 内部终止 |
|:--|:--|--:|--:|:--|--:|
| Candidate-A | + | 5 | ~450 aa | ATG/TGA | 0 |
| Candidate-B | − | 2 | ~350 aa | ATG/TAA | 0 |

任一检查失败（无 ATG、内部终止子、长度非 3 倍数）= 链向/坐标系错误，先修再 BLAST。

## 3. Stage 3 实测：四库 BLASTp

`dbs.tsv` 分层（同种参考 → 近缘种 → 模式植物 → 自蛋白组）：

```
SameSpecies	/path/to/same_species.pep.fa
CloseRelative	/path/to/close_relative.pep.fa
ModelPlant	/path/to/arabidopsis.pep.fa
SelfProteome	out/annotation/db/<genome>_proteome.fasta
```

实测 Top-hit 模式（教学重点：两种截然不同的证据形态）：

| 基因 | 同种参考 | 近缘种 | 模式植物 | 判定 |
|:--|:--|:--|:--|:--|
| Candidate-A | 100% / 99% | 93% / 100% | 63% / 90% | **保守代谢酶（线粒体定位），可克隆** |
| Candidate-B | 84% / 84%（直系为零表达纯预测） | 碎片化 4 命中 31–39% | 32% / 53% | **核酸酶/逆转录酶相关，疑似 TE 衍生，不克隆** |

**功能描述补全**：本地 fasta 头通常只有 ID，Top hit 的功能必须查 TAIR/Phytozome/STRING（WebSearch 位点号即可）。结构域结论可继承直系同源的 InterPro 注释（Candidate-A → FAD 依赖氧化还原酶域；Candidate-B → 内切/外切核酸酶超家族 + RT 相关域）。

**错义位点评估**：把变异氨基酸位置映射到直系同源的结构域坐标上——落在催化域内的小→大带电替换（如 Gly→Arg）破坏性强，与大效应 GWAS 信号互相印证。

**同源拷贝盘点（自比对，多倍体必做）**：
- Candidate-A：4 个单倍型拷贝完整（96–98% identity，仅 11–16 aa 差异）→ 拷贝判别 SNP 充足，可设计特异引物；
- Candidate-B：一个单倍型拷贝缺失、另一拷贝覆盖度截短、其余高趋异 → 拷贝数异常，进一步支持降级。

## 4. Stage 4 实测：实验优先级输出

本案例最终分级（模式示范）：
1. **某已克隆的激素通路负调控因子**（身份明确、克隆已完成）→ 过表达 + 对应激素梯度回复，管线首发；
2. **Candidate-A** → 先 KASP 位点验证，再克隆；功能验证走其代谢通路的分子指标（**不是**激素回复）；
3. **Candidate-B** → 仅 KASP，不克隆；
4. **Candidate-C**（跨性状最强共享峰，但含逆转录转座子结构域）→ 不克隆，只做位点验证。

**核心教训**：不能按"性状 → 激素通路"默认映射设计激素回复（如"株高 → GA/BR"）。实测两个株高错义候选均非激素通路基因。激素回复的门控条件只有一个：**基因注释身份明确属于该激素通路**。

## 5. 报告骨架（可复用）

每基因一节：同源比对表（库/Top hit/identity/cover/E-value）→ 功能域与定位 → 变异位点在结构域中的位置评估 → 同源拷贝表 → 实验建议（克隆/KASP/激素回复/降级，含理由）。开头放总览表（基因 | 身份 | 通路 | 克隆建议），结尾放数据边界与产物清单。
