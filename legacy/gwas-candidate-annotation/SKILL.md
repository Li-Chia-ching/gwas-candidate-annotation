---
name: gwas-candidate-annotation
description: 从 GWAS 显著信号（含供应商截断交付的 P 值截断数据）到候选基因同源注释的端到端流程：双模型 lead SNP 交叉筛选、GFF 窗口候选基因提取、CDS/蛋白提取与翻译、多库 BLASTp 同源注释、多倍体同源拷贝盘点、转座子(TE)风险判别与实验优先级建议（克隆/KASP/激素回复门控）。当用户拿到重测序 GWAS 结果（rMVP/GAPIT/GLM/MLM/FarmCPU 的 signals 或 summary statistics）需要筛选高置信位点、注释候选基因身份、评估克隆价值或设计过表达/KASP/激素回复验证实验时使用。适用于多倍体（含同源四倍体）作物，同样适用于二倍体。
agent_created: true
---

# GWAS 候选基因筛选与同源注释

## 用途

将 GWAS 关联信号转化为**可实验验证的候选基因清单**，并给出每个基因的克隆/验证优先级。覆盖四个阶段：双模型位点筛选 → 候选基因提取 → 同源注释（BLASTp）→ 实验路线建议。

## 何时使用

- 拿到 GWAS 结果（signals.csv / summary statistics / ann.vcf），需要找主效位点与候选基因；
- 需要回答"这个候选基因是什么、值不值得克隆"（TE 关联？拷贝数？功能域？）；
- 需要为多倍体候选基因评估同源拷贝并设计拷贝特异克隆引物；
- 需要判断候选基因是否适用于激素回复实验（必须属于对应激素通路）。

## 执行流程（四阶段）

按序执行，每阶段产物作为下一阶段输入。脚本均在 `scripts/` 下，纯 Python 3 标准库（仅 BLAST 需本地 BLAST+）。

### Stage 1 — 双模型 lead SNP 筛选

```bash
python scripts/screen_dual_model_leads.py \
    --signals-a <trait>.MLM_signals.csv --signals-b <trait>.FarmCPU_signals.csv \
    --model-a MLM --model-b FarmCPU \
    --gff <reference>.gff.gz --trait <TraitName> --outdir out/<trait>/
```

- 默认 gap=100 kb 聚类、±100 kb 窗口取候选基因；窗口是**操作窗口**，不得表述为实测 LD 区块。
- 产出：`*_dual_model_leads.csv`（每独立位点一行 lead SNP）、`*_window_genes.csv`。

### Stage 2 — CDS 提取与蛋白翻译

```bash
python scripts/extract_cds_protein.py --gff <reference>.gff.gz --genome <genome>.fasta.gz \
    --genes <gene1> <gene2> ... --outdir out/annotation/ --proteome
```

- `--proteome` 生成全蛋白组（每基因最长转录本），供 Stage 3 自比对盘点同源拷贝。
- **必须检查 extraction_report.tsv 四项 sanity**：ATG 起始、终止子、无内部终止子、长度 %3==0；任一失败先查链向/坐标系，不得直接进 BLAST。

### Stage 3 — 多库 BLASTp 同源注释

先建 `dbs.tsv`（`名称<TAB>蛋白fasta路径`），推荐分层：同种参考基因组 → 近缘种 → 模式植物 → 自蛋白组（拷贝盘点）。

```bash
python scripts/run_blastp.py --query out/annotation/candidate_proteins.fasta \
    --dbs dbs.tsv --outdir out/annotation/
python scripts/summarize_annotation.py --blastdir out/annotation/blast_out \
    --outdir out/annotation/ --self-db <自蛋白组库名>
```

- 拟南芥等模式植物 Top hit 的功能描述用 WebSearch（TAIR/Phytozome/STRING）补全，本地 fasta 头通常无注释。
- Pfam/CDD 结构域扫描优先 hmmscan；不可用时继承直系同源的 InterPro 注释并在报告中注明。

### Stage 4 — 判读与实验建议

按下表给出逐基因结论，写入汇总报告：

| 证据模式 | 判定 | 建议 |
|:--|:--|:--|
| 近缘种 identity ≥85% + 全长覆盖 + 单基因家族 | 可信直系同源 | 可克隆；按注释通路设计验证 |
| identity 30–60% + 命中碎片化（多个不同基因各覆盖一部分） | 疑似 TE/重复衍生 | **不克隆**，降级为 KASP 位点验证 |
| PFAM/InterPro 含 RVT_*/rve/Retrotrans_gag | TE 关联 | **不克隆**；统计信号再强也不等于可克隆 |
| 自比对显示 1:4（四倍体）完整拷贝 | 拷贝完整 | 克隆引物必须压在拷贝判别 SNP 上，阳性克隆测序确认 |
| 自比对拷贝缺失/截短 | CNV 或注释缺口 | 先查拷贝数，再定策略 |
| 注释为激素通路组分（如 DELLA→GA） | 通路明确 | 可设计激素回复（过表达→表型→外施激素回复） |
| 注释为非激素类型（代谢酶/转运蛋白/未知） | 通路不明 | **禁止默认激素回复**；改做该基因类型对应的分子指标 |

## 数据边界（必须遵守，详见 references/data_boundary.md）

- 只有 P 值截断信号（如 P<1e-5）时：不可重跑 GWAS、λgc 不可估、LD 仅限显著位点间、双模型显著 ≠ 独立 QTL ≠ 因果基因。
- 候选基因的坐标、结构、引物设计必须使用**同一参考基因组版本**，不得混用两个版本坐标系。

## 参考文件

- `references/workflow.md` — 完整工作流细节 + 株高性状实测案例（已脱敏，含数值模式与判读示范）
- `references/data_boundary.md` — 截断数据的解释边界与汇报口径
- `references/pitfalls.md` — 多倍体克隆、TE 判别、工具环境等已验证的坑
