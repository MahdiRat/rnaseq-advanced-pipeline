#!/usr/bin/env bash
set -euo pipefail

# Usage:
# bash scripts/01_qc_trim_align.sh \
#   /path/to/fastqs \
#   /path/to/output \
#   /path/to/hisat2_index_prefix \
#   /path/to/trimmomatic.jar \
#   /path/to/adapters.fa \
#   /path/to/annotation.gtf \
#   8
#
# Example:
# bash scripts/01_qc_trim_align.sh \
#   data/raw results \
#   data/reference/hisat2_index/chrX_tran \
#   tools/Trimmomatic/trimmomatic.jar \
#   tools/Trimmomatic/adapters/combined_adapters.fa \
#   data/reference/annotation.gtf \
#   8

FASTQ_DIR="${1:?ERROR: FASTQ_DIR missing}"
OUT_DIR="${2:?ERROR: OUT_DIR missing}"
HISAT2_INDEX="${3:?ERROR: HISAT2_INDEX prefix missing}"
TRIMMOMATIC_JAR="${4:?ERROR: TRIMMOMATIC_JAR missing}"
ADAPTERS_FA="${5:?ERROR: ADAPTERS_FA missing}"
GTF="${6:?ERROR: GTF missing}"
THREADS="${7:-4}"

mkdir -p "${OUT_DIR}"/{fastqc_raw,fastqc_trimmed,trimmed_fastqs,alignments,sorted_bams,counts,logs}

log(){ echo "[$(date '+%F %T')] $*" | tee -a "${OUT_DIR}/logs/pipeline.log"; }

need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing tool: $1"; exit 1; }; }
need fastqc
need hisat2
need samtools
need java
need featureCounts

[[ -d "${FASTQ_DIR}" ]] || { echo "ERROR: FASTQ_DIR not found: ${FASTQ_DIR}"; exit 1; }
[[ -f "${TRIMMOMATIC_JAR}" ]] || { echo "ERROR: trimmomatic jar not found: ${TRIMMOMATIC_JAR}"; exit 1; }
[[ -f "${ADAPTERS_FA}" ]] || { echo "ERROR: adapters file not found: ${ADAPTERS_FA}"; exit 1; }
[[ -f "${GTF}" ]] || { echo "ERROR: GTF not found: ${GTF}"; exit 1; }

# detect fastq extension
shopt -s nullglob
R1_LIST=( "${FASTQ_DIR}"/*_1.fastq "${FASTQ_DIR}"/*_1.fastq.gz )
[[ ${#R1_LIST[@]} -gt 0 ]] || { echo "ERROR: No *_1.fastq or *_1.fastq.gz found in ${FASTQ_DIR}"; exit 1; }

log "FASTQ_DIR: ${FASTQ_DIR}"
log "OUT_DIR:   ${OUT_DIR}"
log "INDEX:     ${HISAT2_INDEX}"
log "THREADS:   ${THREADS}"

# ---- FastQC raw ----
log "Running FastQC on raw reads..."
for R1 in "${R1_LIST[@]}"; do
  bn="$(basename "${R1}")"
  SAMPLE="${bn%_1.fastq}"
  SAMPLE="${SAMPLE%_1.fastq.gz}"
  R2="${FASTQ_DIR}/${SAMPLE}_2.fastq"
  [[ -f "${R2}" ]] || R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
  [[ -f "${R2}" ]] || { log "WARNING: Missing R2 for ${SAMPLE}, skipping"; continue; }

  log "FastQC raw: ${SAMPLE}"
  fastqc "${R1}" "${R2}" -o "${OUT_DIR}/fastqc_raw" >> "${OUT_DIR}/logs/fastqc_raw.log" 2>&1
done
log "FastQC raw done."

# ---- Trimmomatic ----
log "Running Trimmomatic..."
for R1 in "${R1_LIST[@]}"; do
  bn="$(basename "${R1}")"
  SAMPLE="${bn%_1.fastq}"
  SAMPLE="${SAMPLE%_1.fastq.gz}"
  R2="${FASTQ_DIR}/${SAMPLE}_2.fastq"
  [[ -f "${R2}" ]] || R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
  [[ -f "${R2}" ]] || { log "WARNING: Missing R2 for ${SAMPLE}, skipping"; continue; }

  log "Trim: ${SAMPLE}"
  java -jar "${TRIMMOMATIC_JAR}" PE \
    "${R1}" "${R2}" \
    "${OUT_DIR}/trimmed_fastqs/${SAMPLE}_1_trimmed.fastq.gz" \
    "${OUT_DIR}/trimmed_fastqs/${SAMPLE}_1_unpaired.fastq.gz" \
    "${OUT_DIR}/trimmed_fastqs/${SAMPLE}_2_trimmed.fastq.gz" \
    "${OUT_DIR}/trimmed_fastqs/${SAMPLE}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS_FA}":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    >> "${OUT_DIR}/logs/trimmomatic.log" 2>&1
done
log "Trimmomatic done."

# ---- FastQC trimmed ----
log "Running FastQC on trimmed reads..."
for R1 in "${OUT_DIR}/trimmed_fastqs"/*_1_trimmed.fastq.gz; do
  [[ -f "${R1}" ]] || continue
  SAMPLE="$(basename "${R1}" _1_trimmed.fastq.gz)"
  R2="${OUT_DIR}/trimmed_fastqs/${SAMPLE}_2_trimmed.fastq.gz"
  [[ -f "${R2}" ]] || { log "WARNING: Missing trimmed R2 for ${SAMPLE}, skipping"; continue; }

  log "FastQC trimmed: ${SAMPLE}"
  fastqc "${R1}" "${R2}" -o "${OUT_DIR}/fastqc_trimmed" >> "${OUT_DIR}/logs/fastqc_trimmed.log" 2>&1
done
log "FastQC trimmed done."

# ---- Alignment + sort/index ----
log "Aligning with HISAT2..."
for R1 in "${OUT_DIR}/trimmed_fastqs"/*_1_trimmed.fastq.gz; do
  [[ -f "${R1}" ]] || continue
  SAMPLE="$(basename "${R1}" _1_trimmed.fastq.gz)"
  R2="${OUT_DIR}/trimmed_fastqs/${SAMPLE}_2_trimmed.fastq.gz"
  [[ -f "${R2}" ]] || { log "WARNING: Missing trimmed R2 for ${SAMPLE}, skipping"; continue; }

  log "Align+Sort: ${SAMPLE}"
  hisat2 -p "${THREADS}" -x "${HISAT2_INDEX}" \
    -1 "${R1}" -2 "${R2}" \
    2>> "${OUT_DIR}/logs/hisat2.log" \
  | samtools view -bS - \
  | samtools sort -o "${OUT_DIR}/sorted_bams/${SAMPLE}_sorted.bam" -

  samtools index "${OUT_DIR}/sorted_bams/${SAMPLE}_sorted.bam"
done
log "Alignment done."

# ---- featureCounts ----
log "Running featureCounts..."
featureCounts -T "${THREADS}" -p -B -C \
  -a "${GTF}" \
  -o "${OUT_DIR}/counts/gene_counts.txt" \
  "${OUT_DIR}"/sorted_bams/*_sorted.bam \
  >> "${OUT_DIR}/logs/featurecounts.log" 2>&1

log "Counts: ${OUT_DIR}/counts/gene_counts.txt"
log "Pipeline finished ✅"

