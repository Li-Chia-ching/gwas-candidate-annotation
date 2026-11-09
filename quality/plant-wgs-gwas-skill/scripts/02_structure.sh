#!/usr/bin/env bash
# =============================================================================
# S3 LD pruning -> PCA -> K / GRM，并输出结构诊断材料
# 用法：bash 02_structure.sh [project_config.yaml]
# 依赖：plink2、gemma 或 gcta64（按 structure.kinship_method 选择）、python3 + pyyaml
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
GCTA=$(cfg software.gcta); GCTA=${GCTA:-gcta64}
WIN=$(cfg ld_pruning.window)
STEP=$(cfg ld_pruning.step)
R2=$(cfg ld_pruning.r2)
NPC=$(cfg structure.n_pc)
KMETHOD=$(cfg structure.kinship_method)
CUTOFF=$(cfg structure.grm_cutoff)
CLEAN="$WORK/qc/clean"

mkdir -p "$WORK/structure" "$WORK/logs"

echo "[S3] LD pruning: --indep-pairwise $WIN $STEP $R2"
"$PLINK" --bfile "$CLEAN" --indep-pairwise "$WIN" "$STEP" "$R2" \
  --out "$WORK/structure/pruned" 2>&1 | tee "$WORK/logs/05_pruning.log"

echo "[S3] PCA（只使用独立标记）: --pca $NPC"
"$PLINK" --bfile "$CLEAN" --extract "$WORK/structure/pruned.prune.in" \
  --pca "$NPC" --out "$WORK/structure/pca" 2>&1 | tee "$WORK/logs/06_pca.log"

echo "[S3] 重复 / 近亲诊断（IBS，阈值 $(cfg qc.ibs_cutoff)）"
"$PLINK" --bfile "$CLEAN" --extract "$WORK/structure/pruned.prune.in" \
  --make-king square --out "$WORK/structure/king" 2>&1 | tee "$WORK/logs/07_king.log"

echo "[S3] 亲缘矩阵 K：method=$KMETHOD"
rm -f "$WORK/structure/kinship_gemma.txt"
case "$KMETHOD" in
  gemma)
    "$GEMMA" -bfile "$CLEAN" -gk 1 -outdir "$WORK/structure" -o kin \
      2>&1 | tee "$WORK/logs/08_kin.log"
    mapfile -t KFILES < <(find "$WORK/structure" -maxdepth 1 -name "kin*.txt" ! -name "kinship_gemma.txt")
    if [ "${#KFILES[@]}" -ne 1 ]; then
      echo "GEMMA K 矩阵结果不唯一或缺失：${KFILES[*]}" >&2
      exit 2
    fi
    KFILE="${KFILES[0]}"
    printf '%s\n' "$KFILE" > "$WORK/structure/kinship_gemma.txt"
    ;;
  gcta_grm)
    "$GCTA" --bfile "$CLEAN" --autosome --maf "$(cfg qc.maf)" --make-grm \
      --out "$WORK/structure/grm" 2>&1 | tee "$WORK/logs/08_kin.log"
    KFILE="$WORK/structure/grm.grm.bin"
    ;;
  plink2_king)
    KFILE="$WORK/structure/king.king.square"
    ;;
  *)
    echo "未知 kinship_method: $KMETHOD" >&2; exit 2
    ;;
esac

{
  echo "# 群体结构与亲缘报告"
  echo
  echo "- LD pruning：窗口 $WIN SNP / 步长 $STEP / r² $R2（按物种 LD 衰减调整，书中 p241–243）"
  echo "- 独立标记数：$(wc -l < "$WORK/structure/pruned.prune.in")（同时作为 M_eff 候选）"
  echo "- 主成分数：$NPC（书中常用 10 或 20，p245；**必须结合碎石图与亚群/批次相关性复核**）"
  echo "- 亲缘矩阵：$KMETHOD → \`$KFILE\`（grm_cutoff=$CUTOFF）"
  echo
  echo "## 必做检查"
  echo "1. 画 PC1×PC2 散点并按已知亚群着色，确认是否有成簇结构。"
  echo "2. 计算各 PC 与表型、测序批次、地理来源的相关系数；若显著相关，需纳入协变量或增加 PC。"
  echo "3. 检查 king 矩阵中过近的个体对，按 qc.ibs_cutoff 决定是否去重（家系结构除外）。"
  echo "4. 确认 K 矩阵正定、无 NaN；否则回退 S3 重建。"
} > "$WORK/structure_report.md"

echo "[完成] 结构输出：$WORK/structure/  |  报告：$WORK/structure_report.md"
echo "        K 矩阵：$KFILE"
