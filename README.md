# gwas-candidate-annotation V3

从植物 GWAS 关联结果生成候选位点、候选基因、同源候选、证据等级、实验优先级和可视化报告。当前版本只输出 `Priority A/B/C/D` 与 `Evidence L0-L4`，不宣称因果基因确认。

## 5 分钟 quickstart

依赖：Python 3.11+、NumPy 2.5.4、Matplotlib 3.11.2、PyYAML 6.0.3、NCBI BLAST+ 2.17.0。

```bash
pip install -r requirements.txt
python scripts/run_pipeline.py --config examples/minimal/config.yaml --workdir out/minimal
```

预期输出：

- `out/minimal/report/report.html`
- `out/minimal/report/report.md`
- `out/minimal/report/tables/*.tsv`
- `out/minimal/report/figures/*.pdf`
- `out/minimal/run_manifest.json`
- `out/minimal/config_snapshot.yaml`

`examples/minimal/` 是小型合成数据，可直接随包分发。大型 E2E 产物不进发布包；运行测试可重新生成。

## 测试

```bash
python -m unittest discover -s tests -v
# 或使用 Makefile
make test
```

当前验证：14 项 P0 回归、10 项迁移回归、2 项报告测试、6 类端到端构造数据、12 项 v3.1 LD 合并/排序回归，共 44/44 通过。

## 已实现能力

- 双模型 summary 显著性过滤。
- REF/ALT、多等位冲突、效应方向和参考版本硬检查。
- 主模型 lead SNP 与跨模型支持分离。
- 距离聚类显式标记为 `distance_cluster`，不冒充 LD 位点。
- lead-lead 剂量 LD r² 单链合并（`merge_clusters_by_ld`）：配置 `ld_merge: true`（需 `vcf`）与 `ld_merge_r2`（默认 0.5）开启，合并后位点标记 `ld_merged_cluster`/`merged_by_ld`。
- 位点排序（`rank_loci`）：每次运行按 nominal lead P 升序输出 `locus_rank`，rank 1 标记 `priority_peak`。
- GFF 固定窗口、直接注释基因和基因体重叠三类候选证据。
- 严格 CDS/蛋白 QC，失败阻断候选提取。
- 多库 BLASTp、Top-hit 同源候选和 self-BLAST 同源拷贝候选。
- TE/repeat 降级和激素回复门控。
- 二等位 VCF 剂量 LD r² 热图，支持二倍体和更高倍性 GT。
- full/truncated 显式数据模式；truncated 禁用 Q-Q/λ。
- `report.html`、`report.md`、标准化 TSV、独立 600-dpi PDF 图。
- 输入 SHA-256、配置快照、软件版本、命令、随机种子和参考版本清单。

详细能力表：`docs/CAPABILITY_MATRIX.md`。

## 待实现

- statistical fine-mapping、PIP/credible set、条件分析。
- 表达、eQTL、单倍型和独立群体复现数据适配器。
- Pfam/InterPro/hmmscan 自动结构域扫描；当前需要 `domains.tsv`。
- reciprocal best hit、共线性或基因树直系同源推断。
- CNV/真实同源拷贝数估计。
- MLMM/FarmCPU/BLINK 直接执行。
- `ann.vcf` 的 ANN/HGVS 自动解析。

PLINK2/GEMMA 上游实机验证是外部依赖，未在本仓库 CI 覆盖；当前只有字段适配与构造输入集成测试。

## 输出语言约束

BLAST Top hit 只写 `homology_candidate`。self-BLAST 结果只写 `homoeolog_candidate`。lead SNP 一律表述为 nominal lead（位点内最小 P），位点内多独立信号的条件分析（COJO/stepwise）未实现，报告中显式标注。禁止写“直系同源”“真实拷贝数”“确认因果基因”。

## 许可与第三方

- 代码：MIT，见 `LICENSE`。
- 第三方依赖和字体说明：`THIRD_PARTY.md`。
- Arial 字体不随包分发；报告使用用户系统中的授权字体。

## 目录

```text
├── LICENSE / VERSION / CHANGELOG.md / THIRD_PARTY.md
├── SKILL.md / README.md / pyproject.toml / requirements.txt
├── scripts/run_pipeline.py
├── src/gwas_candidate/
├── tests/
├── examples/minimal/
├── docs/
├── legacy/gwas-candidate-annotation/
└── quality/plant-wgs-gwas-skill/
```

`legacy/` 保留审计来源和历史接口。`quality/` 是质量框架来源。新工作流以 `src/gwas_candidate/` 和 `scripts/run_pipeline.py` 为准。
