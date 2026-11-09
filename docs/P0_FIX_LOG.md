# P0 修复证据日志

每项均遵循：先运行回归测试并确认基线失败，再修复，修复后重跑通过，最后单独提交。

| P0 | 修复前结果 | 修复后结果 | 提交 |
|---|---|---|---|
| 完整 summary 显著性过滤 | `test_p0_01_summary_filter.py` 失败：基线无 `--p-threshold-a/b`，无法保证交集前过滤 | 同一测试通过；只保留两模型均达到阈值的位点 | `c48ef95` |
| 距离簇与 LD 位点分离 | `test_p0_02_cluster_not_ld_locus.py` 失败：stdout 仍输出 `independent loci` | 通过；输出 `distance_cluster` 和 `ld_locus_status=not_evaluated` | `1e7a5d4` |
| 主模型 lead 与跨模型支持分离 | `test_p0_03_primary_lead.py` 失败：lead 取两模型最小 P，得到 100 而非主模型 150 | 通过；`lead_model=MLM`，`support_models=MLM;FarmCPU` | `309a454` |
| REF/ALT、效应方向、参考版本硬检查 | `test_p0_04_allele_effect_reference.py` 4 项失败：基线无硬检查参数 | 4 项通过；等位冲突、反向效应、版本不匹配均阻断 | `af853c0` |
| PLINK2/GEMMA 字段适配 | `test_p0_05_native_adapters.py` 2 项失败：PLINK2 缺 `CHR`，GEMMA 缺 `BP/P` | 2 项通过；两类原生结果均可绘图和门禁 | `c6c0454` |
| GEMMA 表型、协变量、K 矩阵衔接 | `test_p0_06_gemma_inputs.py` 失败：GEMMA 收到 `-n BLUP_days` 而非数值表型列 | 通过；逐表型 `.fam`、无表头数值协变量和显式 GEMMA K 路径生效 | `4afc917` |
| M_eff 阈值一致性 | `test_p0_07_meff_threshold.py` 失败：绘图仍用默认 `5e-8`，未按 M_eff=4 得到 0.0125 | 通过；04/05 共用 `bonferroni_meff`，阈值均为 0.0125 | `2f3a239` |
| full/truncated 识别 | `test_p0_08_data_mode.py` 3 项失败：基线无 `--data-mode` | 3 项通过；truncated 禁用 Q-Q/λ，full 正常输出 | `b219aa8` |

## 汇总测试

```bash
python -m unittest discover -s tests -p "test_p0_*.py" -v
```

结果：14 项测试全部通过。

## 后续全量验证

```bash
python -m unittest discover -s tests -v
```

结果：32 项测试全部通过，包括迁移回归、报告生成和六类端到端构造数据。
