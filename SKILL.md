---
name: gwas-candidate-annotation-v3
description: 从植物 GWAS 关联结果生成候选位点、候选基因、同源候选、证据等级、实验优先级和可视化报告。适用于完整或供应商截断的 summary statistics、GFF、基因组 FASTA、BLAST 蛋白库、结构域表和 VCF LD 数据。只输出 Priority A/B/C/D 与 Evidence L0-L4，不宣称因果基因确认。
---

# gwas-candidate-annotation V3

## 使用边界

1. 必须显式声明 `data_mode: full` 或 `data_mode: truncated`。
2. 完整 summary 必须在双模型交集前按各自显著性阈值过滤。
3. 距离聚类只能称 `distance_cluster`；未经 LD 或条件分析，不得称独立位点或独立 QTL。
4. 主模型 lead SNP 与跨模型支持分开报告。
5. REF/ALT、多等位冲突、效应方向冲突和参考基因组版本不一致均为阻断错误。
6. BLAST Top hit 只能写 `homology_candidate`；self-BLAST 只能写 `homoeolog_candidate`。不得写直系同源或真实拷贝数。
7. TE/repeat 风险触发克隆降级；激素回复只在注释通路与目标激素明确匹配时允许。
8. 输出只能是候选优先级和证据等级，不写“确认因果基因”。

## 输入契约

YAML 配置示例字段：

```yaml
scenario: demo
base_dir: ./inputs
trait: plant_height
data_mode: full
reference_build: ZM4_V2_hap1
random_seed: 20261109
model_a:
  name: MLM
  path: height.MLM.csv
  p_threshold: 5.0e-8
model_b:
  name: FarmCPU
  path: height.FarmCPU.csv
  p_threshold: 5.0e-8
gff: reference.gff3
genome: genome.fa
db_config: dbs.tsv
self_db: Self
domains: domains.tsv
vcf: variants.vcf
gap: 100000
flank: 100000
region_flank: 100000
```

关联表必须含 SNP/ID、CHROM/CHR、POS/BP、REF、ALT、Effect/BETA、P 值列。`ann.vcf` 的变异后果自动解析仍为待实现。

## 已实现能力

- 双模型显著性过滤、等位基因交集、效应方向检查、主模型 lead 与跨模型支持。
- GFF 固定窗口、直接注释基因、基因体重叠三类候选证据。
- 严格 CDS/蛋白 QC；失败默认阻断候选提取。
- 多库 BLASTp、Top-hit 汇总和 self-BLAST 同源候选。
- TE/repeat 降级、激素回复门控、Priority A-D 与 Evidence L0-L4。
- 二等位 VCF 剂量 LD r² 热图；二倍体和更高倍性 GT 均可按剂量计算。
- `report.html`、`report.md`、标准化 TSV 和独立 PDF 图。
- full 模式输出 Q-Q；truncated 模式禁用 Q-Q/λ。
- 运行清单包含输入 SHA-256、配置快照、软件版本、命令、随机种子和参考版本。

## 待实现

- statistical fine-mapping、PIP/credible set、条件分析。
- 表达、eQTL、单倍型和独立群体复现数据适配器。
- Pfam/InterPro/hmmscan 自动结构域扫描；当前需要用户提供 `domains.tsv`。
- reciprocal best hit、共线性或基因树支持的直系同源推断。
- CNV/真实同源拷贝数估计。
- MLMM/FarmCPU/BLINK 的直接执行；当前仅支持其结果作为下游输入。
- `ann.vcf` 的 ANN/HGVS 自动解析。
- 上游 PLINK2/GEMMA 实机数据验证；当前只有字段适配与构造输入集成测试。

## 运行

```bash
pip install -r requirements.txt
python scripts/run_pipeline.py --config config.yaml --workdir out/run1
```

主输出：`out/run1/report/report.html`、`report.md`、`tables/` 和 `figures/`。运行审计：`out/run1/run_manifest.json`。
