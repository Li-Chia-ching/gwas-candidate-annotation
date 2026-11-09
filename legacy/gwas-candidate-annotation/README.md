# gwas-candidate-annotation（历史版）

本目录是 V3 重构的审计来源和历史接口，仅供回溯。当前能力、输入契约和输出语言以仓库根目录 `README.md` 与 `SKILL.md` 为准。

## 历史能力边界

- `scripts/screen_dual_model_leads.py`：双模型信号筛选和 GFF 固定窗口候选基因。
- `scripts/extract_cds_protein.py`：CDS/蛋白提取与结构检查。
- `scripts/run_blastp.py`：多库 BLASTp。
- `scripts/summarize_annotation.py`：BLAST 命中长表。
- `references/`：截断数据边界、实验门控和历史案例规则。

历史文档中的“end-to-end”“validated”“independent loci”“ortholog”“copy inventory”等完成式表述不再作为 V3 能力声明。V3 只承认已测试实现：BLAST Top hit 是 `homology_candidate`，self-BLAST 结果是 `homoeolog_candidate`，距离聚类是 `distance_cluster`。

## 三条保留规则

1. 同一分析只使用一个参考基因组版本。
2. 统计显著不等于可克隆；TE/repeat 风险触发降级。
3. 激素回复只在候选基因注释明确属于对应激素通路时允许。

## 未实现或已重写

- `ann.vcf` 自动解析：待实现。
- 自动结构域扫描：待实现。
- 自动实验报告：已由 V3 `src/gwas_candidate/` 重写。
- 真实直系同源或真实拷贝数判定：待实现。
