#!/usr/bin/env bash
# =============================================================================
# S5 逐表型关联分析（共享同一套 clean 基因型、PC、K 矩阵）
# 用法：bash 03_gwas.sh [project_config.yaml]
# 依赖：plink2、gemma（LMM 路线）、python3 + pyyaml
# =============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
CONFIG=${1:-"$HERE/../assets/project_config.yaml"}
CONFIG=$(realpath "$CONFIG")
BASE=$(dirname "$CONFIG")

cfg() { python3 "$HERE/config.py" "$CONFIG" get "$1"; }
p() { case "$1" in /*|[A-Za-z]:/*|[A-Za-z]:\\*) echo "$1" ;; *) echo "$BASE/$1" ;; esac; }

WORK=$(p "$(cfg paths.workdir)")
PLINK=$(cfg software.plink); PLINK=${PLINK:-plink2}
GEMMA=$(cfg software.gemma); GEMMA=${GEMMA:-gemma}
CLEAN="$WORK/qc/clean"
NPC=$(cfg structure.n_pc)
NPHENO=$(python3 "$HERE/config.py" "$CONFIG" count phenotypes)

mkdir -p "$WORK/assoc" "$WORK/logs"

if [ ! -f "$WORK/structure/pca.eigenvec" ]; then
  echo "[S5] 未找到 PCA 结果，请先运行 02_structure.sh" >&2
  exit 2
fi

for ((i = 0; i < NPHENO; i++)); do
  NAME=$(cfg "phenotypes.$i.name")
  PFILE=$(p "$(cfg "phenotypes.$i.file")")
  COL=$(cfg "phenotypes.$i.column")
  TYPE=$(cfg "phenotypes.$i.type")
  MODEL=$(cfg "phenotypes.$i.model")
  MODEL=${MODEL:-$(cfg model.default)}

  echo "------------------------------------------------------------"
  echo "[S5] 表型 $i: $NAME ($TYPE)  模型=$MODEL  表型列=$COL"

  PHENO="$WORK/assoc/${NAME}.pheno"
  PLINK_COVAR="$WORK/assoc/${NAME}.plink.covar.tsv"
  GEMMA_PREFIX="$WORK/assoc/${NAME}_gemma"
  GEMMA_COVAR="$WORK/assoc/${NAME}.gemma.covar.txt"
  EXTRA_ARGS=()
  if [ -f "$WORK/structure/covariates.tsv" ]; then
    EXTRA_ARGS+=(--extra-covariates "$WORK/structure/covariates.tsv")
  fi

  python3 "$HERE/prepare_gwas_inputs.py" \
    --fam "$CLEAN.fam" \
    --phenotype-file "$PFILE" \
    --phenotype-column "$COL" \
    --pca "$WORK/structure/pca.eigenvec" \
    --n-pc "$NPC" \
    "${EXTRA_ARGS[@]}" \
    --plink-pheno "$PHENO" \
    --plink-covar "$PLINK_COVAR" \
    --gemma-fam "$GEMMA_PREFIX.fam" \
    --gemma-covar "$GEMMA_COVAR"
  cp "$CLEAN.bed" "$GEMMA_PREFIX.bed"
  cp "$CLEAN.bim" "$GEMMA_PREFIX.bim"

  case "$MODEL" in
    GLM)
      "$PLINK" --bfile "$CLEAN" --pheno "$PHENO" --pheno-name "$COL" \
        --linear --covar "$PLINK_COVAR" --covar-variance-standardize \
        --out "$WORK/assoc/${NAME}_glm" 2>&1 | tee "$WORK/logs/gwas_${NAME}.log"
      ;;
    LMM|MLM)
      KPATH_FILE="$WORK/structure/kinship_gemma.txt"
      if [ ! -f "$KPATH_FILE" ]; then
        echo "LMM 需要 02_structure.sh 生成的 GEMMA K 路径记录：$KPATH_FILE" >&2
        exit 2
      fi
      KFILE=$(tr -d '\r\n' < "$KPATH_FILE")
      if [ -z "$KFILE" ] || [ ! -f "$KFILE" ]; then
        echo "GEMMA K 矩阵不存在：$KFILE" >&2
        exit 2
      fi
      "$GEMMA" -bfile "$GEMMA_PREFIX" -k "$KFILE" -c "$GEMMA_COVAR" \
        -n 1 -lmm 1 -outdir "$WORK/assoc" -o "${NAME}_lmm" \
        2>&1 | tee "$WORK/logs/gwas_${NAME}.log"
      ;;
    MLMM|FarmCPU|BLINK)
      echo "[S5] $MODEL 为 R/GAPIT 路线，请按 references/model-selection.md 运行 GAPIT；"
      echo "     输入：$CLEAN.*、$PLINK_COVAR、$PHENO，输出写入 $WORK/assoc/${NAME}.assoc.tsv"
      ;;
    *)
      echo "未知模型: $MODEL" >&2; exit 2
      ;;
  esac
done

{
  echo "# 关联分析运行记录"
  echo
  echo "配置快照：\`$WORK/config_snapshot.yaml\`"
  echo "表型数：$NPHENO；默认模型：$(cfg model.default)；PC 数：$NPC"
  echo "共享输入：\`$CLEAN.*\`；协变量与表型按样本 ID 显式对齐"
  echo
  echo "每个表型的输出必须至少包含 CHR、BP、A1、NMISS、BETA/OR、SE、P；"
  echo "效应量与 p 值须同时报告（书中 p101）。若使用 R 路线，请把结果另存为"
  echo "\`$WORK/assoc/<name>.assoc.tsv\` 后运行 04_plots.py 与 05_checks.py。"
} > "$WORK/assoc/run_summary.md"

echo "[完成] 结果目录：$WORK/assoc/"
