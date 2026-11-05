# gwas-candidate-annotation

从 GWAS 显著信号到候选基因实验优先级的端到端 Skill。面向作物重测序 GWAS 下游分析，同源四倍体作物实测验证，二倍体同样适用。

## 这个 Skill 做什么

输入 GWAS 结果（rMVP/GAPIT 风格的 signals.csv / summary statistics / ann.vcf），输出：

1. **双模型 lead SNP 清单**：两个模型（如 MLM+FarmCPU）同位置显著 → 聚类成独立位点 → 每位点取 lead SNP；
2. **候选基因表**：lead SNP ±100 kb 窗口基因（GFF）；
3. **候选基因身份**：CDS 提取→蛋白翻译→多库 BLASTp（同种/近缘/模式植物/自蛋白组）→ 功能域与通路注释；
4. **克隆可行性判定**：TE 关联判别（RVT/rve 结构域、碎片化命中）、多倍体同源拷贝盘点（拷贝特异引物前提）、激素回复门控（基因身份属于激素通路才允许）。

## 快速上手（10 分钟跑通）

```bash
# 1. 双模型 lead 筛选（每性状一次）
python scripts/screen_dual_model_leads.py \
    --signals-a 性状.MLM_signals.csv --signals-b 性状.FarmCPU_signals.csv \
    --model-a MLM --model-b FarmCPU \
    --gff 参考基因组.gff.gz --trait 性状名 --outdir out/性状名/

# 2. 提取候选基因 CDS + 蛋白（--proteome 建全蛋白组供拷贝盘点）
python scripts/extract_cds_protein.py --gff 参考基因组.gff.gz \
    --genome 基因组.fasta.gz --genes 基因1 基因2 --outdir out/anno/ --proteome

# 3. 多库 BLASTp（需本地 NCBI BLAST+；dbs.tsv 格式：库名<TAB>蛋白fasta路径）
python scripts/run_blastp.py --query out/anno/candidate_proteins.fasta \
    --dbs dbs.tsv --outdir out/anno/
python scripts/summarize_annotation.py --blastdir out/anno/blast_out \
    --outdir out/anno/ --self-db <自蛋白组库名>

# 4. 按 SKILL.md Stage 4 判定表逐基因写结论
```

## 目录结构

```
gwas-candidate-annotation/
├── SKILL.md                      主文档（YAML 头 + 四阶段流程 + 判定表）
├── README.md                     本文件
├── scripts/
│   ├── screen_dual_model_leads.py   双模型 lead SNP 筛选 + 窗口基因
│   ├── extract_cds_protein.py       CDS/蛋白提取 + 全蛋白组构建
│   ├── run_blastp.py                多库 BLASTp（自动建库）
│   └── summarize_annotation.py      比对汇总 + 同源拷贝表
└── references/
    ├── workflow.md               完整流程细节 + 株高性状实测案例（已脱敏）
    ├── data_boundary.md          截断 GWAS 数据的解释边界与汇报口径
    └── pitfalls.md               多倍体克隆/TE 判别/工具环境的实证坑
```

## 依赖

- Python 3.9+（仅标准库）
- NCBI BLAST+（blastp / makeblastdb；Windows 原生或 Linux）
- 可选：HMMER（hmmscan）用于 Pfam 结构域扫描；不可用时继承直系同源 InterPro 注释

## 三条铁律（详版见 references/）

1. **坐标系隔离**：一个分析只用一个参考基因组版本，禁止混用；
2. **统计显著 ≠ 可克隆**：TE 关联（RVT/rve）、碎片化命中、零表达直系同源 → 一律降级为 KASP 位点验证；
3. **激素回复门控**：基因注释身份明确属于该激素通路才可设计激素回复实验，禁止按性状默认（"株高→GA"是错误默认）。

## 实测来源

流程提炼自一套同源四倍体作物的真实重测序 GWAS 下游分析（3 个数量性状 × 3 种关联模型，群体规模 200+，供应商截断交付场景），包含两个错义候选基因的完整注释与分级判读（一个高度保守的代谢酶 → 可克隆；一个逆转录酶相关疑似 TE 衍生 → 降级）。具体基因身份与项目数据未包含在本 Skill 中，详见使用者所在项目的内部记录。
