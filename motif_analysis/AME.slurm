#!/bin/bash
#SBATCH --job-name=AME
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

MODULE_NAME=$1
INPUT_FASTA=$2
OUT_DIR="${4}/Module${MODULE_NAME}"
EVALUE_THRESH=$3
mkdir -p "${OUT_DIR}"

# ==========================================================
# JASPAR（human+mouse）
JASPAR_VERTEBRATES="/export/home/lisijie/motif_db/JASPAR2026/JASPAR2026_CORE_vertebrates.meme.txt"

# HOCOMOCO（human+mouse）
HOCOMOCO="/export/home/lisijie/motif_db/HOCOMOCOv14/H14CORE_meme_format.meme"

# CIS-BP (find from /export/home/lisijie/motif_db/MEME_motif_databases/motif_databases)
CISBP_HUMAN="/export/home/lisijie/motif_db/MEME_motif_databases/motif_databases/CIS-BP_2.00/Homo_sapiens.meme"
CISBP_MOUSE="/export/home/lisijie/motif_db/MEME_motif_databases/motif_databases/CIS-BP_2.00/Mus_musculus.meme"
# ==========================================================

echo "Running AME for Module ${MODULE_NAME}"
echo "Input: ${INPUT_FASTA}"
echo "Output: ${OUT_DIR}"

time ame \
  --oc "${OUT_DIR}" \
  --control --shuffle-- --kmer 2 \
  --method fisher \
  --scoring totalhits \
  --evalue-report-threshold "${EVALUE_THRESH}" --noseq \
  "${INPUT_FASTA}" "${JASPAR_VERTEBRATES}"
# "${JASPAR_VERTEBRATES}" "${CISBP_MOUSE}" "${HOCOMOCO}"
echo "Done: Module ${MODULE_NAME}"