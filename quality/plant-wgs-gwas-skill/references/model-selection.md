# 关联模型选择

> 依据 Mills, Barban & Tropf (2020) 第 2、4、9 章；`PDF pNNN` 为扫描版页码。

## 决策树

```
START：群体类型？structure_expectation？
│
├─ weak（同质育种群体 / 单一 RIL 家系，PCA 无成簇）
│    → GLM + 少量 PC
│      固定效应：标记 + PC + 环境协变量；快、易解释
│
├─ moderate / strong（自然关联群体、自交系面板、多亚群、有谱系）
│    → LMM / MLM   ← 植物场景默认
│      y = Xβ(标记 + PC + 协变量) + Zu(K 矩阵随机效应) + e
│      书中 p248：混合线性模型已逐渐取代单纯 PCA 校正
│
├─ 复杂性状 / 大效应位点掩盖小效应 / 多位点架构
│    → MLMM（逐步加入显著位点作协变量）
│    → FarmCPU / BLINK（多位点，兼顾功效与假阳性）
│
└─ 有明确家系结构（RIL / NAM / MAGIC / F2）
     → 家系固定效应或家系内关联（PLINK --qfam，书中 p240）
     → 或把家系作为随机效应放进 LMM
```

## 四类模型对照

| 模型 | 校正对象 | 优点 | 局限 | 常用工具 |
| --- | --- | --- | --- | --- |
| GLM + PC | 群体分层 | 快、结果直观、易复现 | 对隐性亲缘/家系结构无能为力；PC 数敏感 | `plink2 --linear --covar` |
| LMM / MLM (K+PC) | 分层 + 亲缘 + 家系 | 植物关联群体默认；假阳性控制好 | 计算量大；可能"过度校正"掉真实信号 | GEMMA、EMMAX、FAST-LMM、GAPIT、GCTA |
| MLMM | 分层 + 亲缘 + 多位点 | 逐步纳入大效应位点，释放小效应 | 逐步选择带来稳定性问题 | GAPIT |
| FarmCPU / BLINK | 分层 + 多位点 | 功效高、假阳性控制较好，适合 N 偏小 | 需调参；结果需与其他模型交叉验证 | GAPIT / R 包 |
| 家系内检验 | 家系结构 | 免疫群体分层 | 只利用家系内分离，功效较低 | `plink --qfam`（p240） |

## 基因型编码与性状类型

书中第 9 章（PDF p236–238）：

| 性状 | 模型 | 编码/选项 | 输出要点 |
| --- | --- | --- | --- |
| 数量性状（株高、开花时间、BLUP 值） | `--linear` | 加性 `AA/AB/BB` = 0/1/2 | BETA、SE、STAT、P、NMISS |
| 二分类（存活、抗病/感病阈值化） | `--logistic` | 个体编码 2=case、1=control | OR、STAT、P（p237） |
| 显性假设 | `--linear dominant` | `AA+AB vs BB` | 与加性模型比较 |
| 隐性假设 | `--linear recessive` | `AA vs AB+BB` | 与加性模型比较 |
| 家系内 | `--qfam` | 家系固定效应 | 免疫分层，功效低 |

输出字段至少要有：CHR、BP、A1（效应等位）、NMISS、BETA/OR、SE、P（PDF p239）。**效应量必须与 p 值一起报告**（PDF p101 指出 p 值不能给出效应量）。

## 协变量

- 用 `--covar` 加载 tab 分隔的协变量文件（PDF p238）。
- 只放**先于表型、不由基因型决定**的变量：播期、区组、密度、养分、接种批次、测序批次、PC。
- 不要把"由基因型派生的变量"（如亚群标签本身、其他表型）当普通协变量。
- PC 与批次同时纳入时检查共线性；若批次与 PC 高度相关，优先用 K 矩阵 + 少量 PC，并在报告中说明。

## 参考命令

```bash
# GLM + PC（结构弱时）
plink2 --bfile clean --linear --covar covariates.tsv \
       --covar-variance-standardize --out assoc/height_glm

# LMM（推荐）：GEMMA，K 矩阵由 02_structure.sh 生成
gemma -bfile clean -k output/kin.sXX.txt -c covariates.tsv \
      -lmm 1 -n 1 -o height_lmm

# LMM：GAPIT（R），可切换 MLM / MLMM / FarmCPU / BLINK
# GAPIT(Y = pheno, G = geno, PCA.total = 10, model = "MLMM", ...)

# 家系内关联
plink --bfile clean --qfam --mperm 1000 --out assoc/family_based
```

## 逐表型运行约定

- 每个表型单独跑，但**共享同一份 `clean` 基因型、同一套 PC、同一套 K 矩阵**。
- 模型、变换、协变量从 `project_config.yaml` 的 `phenotypes[]` 读取，禁止在命令行里临时改参数。
- 每个表型输出 `<phenotype>.assoc.tsv` 与 `<phenotype>.log`（含软件版本、命令、收敛信息）。
- 多表型遗传相关可用 MTAG 联合分析或 LDSC（第 12 章，PDF p332 起），能提升功效，但要注明联合分析的假设。

## 模型诊断

| 现象 | 可能原因 | 处理 |
| --- | --- | --- |
| λ ≫ 1.2，QQ 整体膨胀 | 结构/亲缘未校正、批次效应 | 换 LMM、加 K 矩阵、补批次协变量、调整 PC 数 |
| λ ≪ 1（过度保守） | 过度校正、样本量小、标记数少 | 减少 PC、检查 K 矩阵是否含近亲重复 |
| 模型不收敛 | K 矩阵非正定、表型常数、协变量共线 | 回 S3/S4 |
| 无任何信号，但阳性对照也不显著 | 表型定义/坐标版本/ID 对齐错误 | 回 S1/S4，先修数据再谈模型 |
| GLM 有峰、LMM 无峰 | GLM 的峰是结构假阳性 | 以 LMM 为准，报告差异 |

## 交叉验证建议

至少用**两种模型**（如 LMM + FarmCPU）跑同一表型，比较领先位点是否一致；一致则增强可信度，不一致则回到 QQ 图与 PC-表型相关性排查。书中也强调结果复现与过拟合问题（第 2 章 p50–51）。
