#!/usr/bin/env bash
# =============================================================================
# S1–S2 数据准备与两层质控（个体层 + 标记层）
# 用法：bash 01_qc.sh [project_config.yaml]
# 依赖：plink2、bcftools（可选，用于 WGS 质量过滤）、python3 + pyyaml
# 所有参数从配置文件读取；相对路径以配置文件所在目录为基准。
# =============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
CONFIG=${1:-"$HERE/../assets/project_config.yaml"}
CONFIG=$(realpath "$CONFIG")
BASE=$(dirname "$CONFIG")
PLINK=$(python3 "$HERE/config.py" "$CONFIG" get software.plink)
PLINK=${PLINK:-plink2}

cfg() { python3 "$HERE/config.py" "$CONFIG" get "$1"; }
p() { case "$1" in /*) echo "$1" ;; *) echo "$BASE/$1" ;; esac; }

WORK=$(p "$(cfg paths.workdir)")
VCF=$(p "$(cfg genotype.vcf)")
MIND=$(cfg qc.mind)
GENO=$(cfg qc.geno)
MAF=$(cfg qc.maf)
HWE_ON=$(cfg qc.hwe_enabled)
HWE_P=$(cfg qc.hwe_p)
EXTRA=$(cfg genotype.allow_extra_chr)
MIN_DP=$(cfg genotype.min_dp)
MIN_GQ=$(cfg genotype.min_gq)

mkdir -p "$WORK/qc" "$WORK/logs"
cp "$CONFIG" "$WORK/config_snapshot.yaml"

echo "[S1] VCF -> PLINK binary : $VCF"
EXTRA_FLAG=""
[ "$EXTRA" = "true" ] && EXTRA_FLAG="--allow-extra-chr"

# 若配置了 min_dp / min_gq，先在 VCF 层面用 bcftools 过滤（PLINK 不读 FORMAT/DP、GQ）
if [ -n "$MIN_DP" ] && [ -n "$MIN_GQ" ] && command -v bcftools >/dev/null 2>&1; then
  FILTERED="$WORK/qc/genotypes.qc.vcf.gz"
  echo "[S1] WGS 质量过滤: DP>=$MIN_DP, GQ>=$MIN_GQ"
  bcftools filter -e "DP<$MIN_DP || GQ<$MIN_GQ" -Oz -o "$FILTERED" "$VCF" 2> "$WORK/logs/bcftools.log"
  VCF="$FILTERED"
else
  echo "[S1] 跳过 VCF 层 DP/GQ 过滤（bcftools 不可用或未配置）；请在报告中注明"
fi

"$PLINK" --vcf "$VCF" --double-id $EXTRA_FLAG \
  --snps-only just-acgt --max-alleles 2 --make-bed \
  --out "$WORK/qc/raw" 2>&1 | tee "$WORK/logs/01_make_bed.log"

echo "[S2] 个体层 QC: --mind $MIND"
"$PLINK" --bfile "$WORK/qc/raw" --mind "$MIND" --make-bed \
  --out "$WORK/qc/step1_ind" 2>&1 | tee "$WORK/logs/02_mind.log"

echo "[S2] 个体杂合度（离群判定请在 R/Python 中按 ±$(cfg qc.het_sd) SD 完成）"
"$PLINK" --bfile "$WORK/qc/step1_ind" --het \
  --out "$WORK/qc/het" 2>&1 | tee "$WORK/logs/03_het.log"

echo "[S2] 标记层 QC: --geno $GENO --maf $MAF"
HWE_FLAG=""
if [ "$HWE_ON" = "true" ]; then
  HWE_FLAG="--hwe $HWE_P"
  echo "[S2] HWE 过滤已开启：$HWE_FLAG（仅限确认分型错误为主的异交/自然群体）"
else
  echo "[S2] HWE 过滤关闭（自交物种 / 强结构群体的默认安全设置）"
fi
"$PLINK" --bfile "$WORK/qc/step1_ind" --geno "$GENO" --maf "$MAF" $HWE_FLAG \
  --make-bed --out "$WORK/qc/clean" 2>&1 | tee "$WORK/logs/04_marker_qc.log"

# 汇总每一步的个体 / 标记数，写入 qc_report.md
{
  echo "# QC 报告"
  echo
  echo "配置快照：\`$WORK/config_snapshot.yaml\`"
  echo
  echo "| 步骤 | 个体数 | 标记数 |"
  echo "| --- | --- | --- |"
  for tag in raw step1_ind clean; do
    if [ -f "$WORK/qc/$tag.fam" ]; then
      n=$(wc -l < "$WORK/qc/$tag.fam")
      m=$(wc -l < "$WORK/qc/$tag.bim")
      echo "| $tag | $n | $m |"
    fi
  done
  echo
  echo "参数：mind=$MIND, geno=$GENO, maf=$MAF, hwe_enabled=$HWE_ON, min_dp=$MIN_DP, min_gq=$MIN_GQ"
  echo
  echo "待人工确认：杂合度离群样本、重复/近亲样本（见 02_structure.sh 的 IBS 步骤）、"
  echo "按批次/亚群分层的差异缺失检查。书中依据见 references/qc-parameters.md。"
} > "$WORK/qc_report.md"

echo "[完成] QC 输出：$WORK/qc/clean.*  |  报告：$WORK/qc_report.md"
