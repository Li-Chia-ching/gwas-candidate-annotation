# plant-wgs-gwas Skill · 课题组使用与迁移说明

这是植物 WGS × 表型 GWAS 的**质量框架来源**，也是 V3 的上游参考实现。当前只把已通过测试的字段适配、绘图、门禁和配置解析并入 V3；上游 PLINK2/GEMMA 实机流程仍标为待验证，不作为 V3 已交付能力。

---

## 1. 它能帮你做什么

- 从 VCF 开始，提供**部分实现的上游 GWAS 路线**：质控 → 群体结构 → 表型建模 → 关联分析 → 多重检验 → 结果诊断 → 候选区间。
- 把每一步的**参数、命令、依据、自检项**写清楚，避免"换个表型就重跑一遍还换了阈值"。
- 作为 V3 质量框架来源使用时，已测试能力限于字段适配、Manhattan/Q-Q、λ、阈值解析和门禁报告。
- 上游 PLINK2/GEMMA 命令脚本已修复构造输入衔接，但**待实机数据验证**。

## 2. 安装（两种方式）

**方式 A｜作为 WorkBuddy / CodeBuddy Skill 加载（推荐）**

把整个目录复制到技能目录即可：

```bash
cp -r plant-wgs-gwas-skill ~/.codebuddy/skills/plant-wgs-gwas
# 或放到当前项目：cp -r plant-wgs-gwas-skill ./.codebuddy/skills/plant-wgs-gwas
```

之后在对话里说「用 plant-wgs-gwas 帮我做开花时间的 GWAS 流程」即可触发。若分发给同事，直接把 `plant-wgs-gwas-skill.zip` 发过去解压到同一位置。

**方式 B｜不当 Skill，只当流程手册用**

直接读 `SKILL.md` + `references/`，按 `scripts/` 里的命令手动执行，效果相同。

## 3. 依赖

| 用途 | 依赖 | 说明 |
| --- | --- | --- |
| 配置读取 | Python 3.8+、PyYAML | `pip install pyyaml` |
| 基因型质控 / PCA | PLINK 2.0（`plink2`） | 或 PLINK 1.9，命令需微调 |
| 混合线性模型 | GEMMA（推荐）/ EMMAX / GCTA / GAPIT(R) / FarmCPU(R) | 按 `model` 配置选择 |
| 绘图与门禁 | numpy、matplotlib | `pip install numpy matplotlib` |

脚本不自动安装软件；缺少工具时会在命令处明确报错并提示所需依赖。所有脚本都从 `project_config.yaml` 读参数，不在脚本里写死物种或表型。

## 4. 换方向 / 换表型：只改配置

复制模板后，**绝大多数情况只需要改 `phenotypes` 段**：

```bash
cp assets/project_config.yaml ./my_config.yaml
# 编辑 my_config.yaml 的 phenotypes / paths / threshold
bash scripts/01_qc.sh        ./my_config.yaml
bash scripts/02_structure.sh ./my_config.yaml
bash scripts/03_gwas.sh      ./my_config.yaml
python3 scripts/04_plots.py  --config ./my_config.yaml --pheno flowering_days --data-mode full
python3 scripts/05_checks.py --config ./my_config.yaml --data-mode full
```

各字段含义与调参建议见 `assets/project_config.yaml` 内注释和 `references/qc-parameters.md`。完整填写范例见 `examples/flowering_time_walkthrough.md`。

**必须保持一致的三个东西**（这是不同方向结果可比的前提）：同一份 QC 后基因型、同一套 PC、同一套 K 矩阵。换表型时不要重新质控。

## 5. 目录结构

```
plant-wgs-gwas-skill/
├── SKILL.md                         # Skill 入口：启动协议、工作流、决策树、门禁、硬约束
├── README.md                        # 本文件
├── references/
│   ├── technical-route.md           # 8 阶段技术路线（含 mermaid 流程图与质量门禁）
│   ├── qc-parameters.md             # QC 参数表、命令、植物适配注意、多重检验与 λ
│   ├── model-selection.md           # 模型决策树、各模型适用条件与命令示例
│   ├── phenotype-cookbook.md        # 表型处理手册：BLUP / 变换 / 协变量 / G×E
│   └── source-map.md                # 原书章节与方法要点映射、迁移差异说明
├── assets/
│   ├── project_config.yaml          # ★ 唯一需要改的参数文件
│   ├── phenotype_template.tsv       # 表型表列格式模板
│   └── sample_sheet_template.tsv    # 样本信息表模板
├── scripts/
│   ├── config.py                    # 配置读取工具（shell 脚本通过它取参数）
│   ├── 01_qc.sh                     # S1–S2 数据准备与两层质控
│   ├── 02_structure.sh              # S3 LD pruning / PCA / K 矩阵
│   ├── 03_gwas.sh                   # S5 逐表型关联分析
│   ├── 04_plots.py                  # S6 曼哈顿图 / Q-Q 图 / λ
│   └── 05_checks.py                 # S6 质量门禁检查，输出 go/no-go
└── examples/
    └── flowering_time_walkthrough.md  # 开花时间方向的完整填写与运行范例
```

## 6. 常见问题

**Q：我们是自交作物，HWE 过滤要照做吗？**
A：不要照搬。人类 GWAS 假定 HWE 偏离主要来自分型错误（书中 p225），但自交物种/强结构群体的偏离常是真实生物学。优先用 WGS 的 DP/GQ/等位基因平衡过滤，HWE 只作诊断，见 `references/qc-parameters.md`。

**Q：样本只有 200 份，能发 GWAS 吗？**
A：功效有限。应把目标定为「候选区间 + 可验证假设」，MAF 阈值提高到 ≥0.05，考虑多位点模型与经验阈值，并在报告中明确功效限制。不要强行用 p<5e-8 得出"全基因组显著"。

**Q：曼哈顿图一堆峰，是不是找到很多基因？**
A：显著位点常因 LD 成片出现，一个 LD 块通常只对应一个信号。按「领先 SNP + 局部 LD 块」定义候选区间，见 `references/technical-route.md` §S7。

**Q：λ>1.2 一定是群体分层吗？**
A：不一定。高度多基因性状本身会抬高 λ（书中 p340）。要看 Q-Q 图整体形态、PC 与表型的相关性、以及是否用了 K 矩阵综合判断。

## 7. 来源与内部使用说明

方法学骨架来自 Mills, Barban & Tropf (2020), MIT Press。配套的《GWAS 图文入门简介》中原书 Figure 4.1 流程图为直接截图，其余统计图按书中原理用示意数据重绘，**仅用于课题组内部教学**，不对外发布。

---

**维护建议**：若课题组换了参考基因组或新增了表型，只更新 `assets/` 模板与 `references/source-map.md`；`SKILL.md` 的工作流与硬约束不要随单个课题改动，否则"可迁移"就失效了。
