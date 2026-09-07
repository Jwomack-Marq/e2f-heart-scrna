#!/usr/bin/env bash
# ============================================================================
# Re-count each library with BOTH flow-cell lanes combined (the standard fix).
#
# WHY: lane1 and lane6 are the SAME library sequenced twice for depth
# (~99% shared cell barcodes). They were originally PIPseeker'd SEPARATELY and
# then merged in Seurat, which double-counts nearly every cell at ~half depth.
# The correct workflow is to give PIPseeker BOTH lanes' FASTQs as ONE sample so
# it concatenates reads and dedupes UMIs across lanes -> one matrix per library,
# each cell counted ONCE at its full combined depth.
#
# RUN ON THE CLUSTER (the FASTQs are there, not in the OneDrive project folder).
# This is also the template for processing the future replicated cohort.
# ============================================================================
set -euo pipefail

PIPSEEKER=~/PIPseeker/pipseeker-v3.3.0-linux/pipseeker
IDX=/user/jzhou54/Desktop/tliu4/db/genome/Refdata_scRNA_MAESTRO_GRCm38_1.2.2/GRCm38_STAR_2.7.6a
RAW1=/projects/rpci/tliu4/working_projects/2025_Han_scRNA/Raw_data/01.RawData          # lane1
RAW6=/projects/rpci/tliu4/working_projects/2025_Han_scRNA/Raw_data_Lane6/01.RawData    # lane6
WORK=/projects/rpci/tliu4/jiaojiaozhou/Projects/2025_Han_scRNA                          # output root
COMB="$WORK/combined_fastq"   # per-library dirs of symlinks to BOTH lanes' fastqs

# sample -> Novogene library ID (identical across the two lanes; see run_config.csv)
declare -A CKDL=( [P0WT]=CKDL250003755 [P0KO]=CKDL250003756 [P7WT]=CKDL250003753 [P7KO]=CKDL250003754 )

for S in P0WT P0KO P7WT P7KO; do
  LIB="${S}_${CKDL[$S]}"
  L1="$RAW1/$S/$LIB"
  L6="$RAW6/$S/$LIB"
  DST="$COMB/$S"
  echo "=== $S : combining $L1  +  $L6 ==="
  mkdir -p "$DST"

  # Symlink both lanes' FASTQs into one directory. The two lanes are on different
  # flow-cell lanes, so filenames carry distinct _L00X_ tags and do not collide.
  # (If any filename DID collide, prefix it -- check before running.)
  ln -sf "$L1"/*.fastq.gz "$DST"/ 2>/dev/null || ln -sf "$L1"/*.fq.gz "$DST"/
  ln -sf "$L6"/*.fastq.gz "$DST"/ 2>/dev/null || ln -sf "$L6"/*.fq.gz "$DST"/
  echo "  $(ls "$DST"/*_R1* 2>/dev/null | wc -l) R1 files combined (expect lane1 + lane6 sets)"

  # One PIPseeker run over BOTH lanes -> one matrix/library, UMI-deduped, full depth.
  # Sensitivity 4-5 kept for consistency; REVIEW knee plots since depth ~doubled.
  "$PIPSEEKER" full \
    --fastq "$DST" \
    --star-index-path "$IDX" \
    --chemistry V \
    --min-sensitivity 4 --max-sensitivity 5 \
    --output-path "$WORK/${S}_combined"
done

# ALTERNATIVE (if your PIPseeker build accepts multiple --fastq inputs; check
#  `pipseeker full --help`): skip the symlink dir and pass a comma-separated list:
#   --fastq "$L1,$L6"

echo "=== DONE. New per-library matrices: <SAMPLE>_combined/filtered_matrix/sensitivity_5/ ==="
echo "Expected sanity check: each combined library has ~the SAME cell count as a single"
echo "old lane (not the sum of both), and ~2x the median UMIs/genes per cell."
