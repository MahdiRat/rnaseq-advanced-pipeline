#!/bin/bash

# --- Configuration ---
FASTQ_DIR="/home/mahdiel/rnaseq/samples"        # Directory containing raw FASTQ files
OUTPUT_DIR="/home/mahdiel/rnaseq/output"        # Main output directory
TRIMMOMATIC_PATH="/home/mahdiel/rnaseq/Trimmomatic-0.36/trimmomatic-0.36.jar"  # Path to Trimmomatic
GENOME_DIR="/home/mahdiel/rnaseq/indexes/genome_index/chrX_tran"   # HISAT2 genome index directory

# Create necessary output directories if they don't exist
mkdir -p "${OUTPUT_DIR}/fastqc_raw"
mkdir -p "${OUTPUT_DIR}/fastqc_trimmed"
mkdir -p "${OUTPUT_DIR}/trimmed_fastqs"
mkdir -p "${OUTPUT_DIR}/alignments"
mkdir -p "${OUTPUT_DIR}/sorted_bams"

# --- FastQC on raw reads ---
echo "Running FastQC on raw reads..."
for R1_FILE in "${FASTQ_DIR}"/*_1.fastq; do
    SAMPLE_NAME=$(basename "${R1_FILE}" _1.fastq)
    R2_FILE="${FASTQ_DIR}/${SAMPLE_NAME}_2.fastq"
    
    echo "Running FastQC on raw reads for sample: ${SAMPLE_NAME}"
    fastqc "${R1_FILE}" "${R2_FILE}" -o "${OUTPUT_DIR}/fastqc_raw"
done

echo "FastQC on raw reads completed."

# --- Trimming with Trimmomatic ---
echo "Running Trimmomatic to trim adapters and filter low-quality reads..."
for R1_FILE in "${FASTQ_DIR}"/*_1.fastq; do
    SAMPLE_NAME=$(basename "${R1_FILE}" _1.fastq)
    R2_FILE="${FASTQ_DIR}/${SAMPLE_NAME}_2.fastq"
    
    echo "Trimming adapters for sample: ${SAMPLE_NAME}"
    java -jar "${TRIMMOMATIC_PATH}" PE \
        "${R1_FILE}" "${R2_FILE}" \
        "${OUTPUT_DIR}/trimmed_fastqs/${SAMPLE_NAME}_1_trimmed.fastq" \
        "${OUTPUT_DIR}/trimmed_fastqs/${SAMPLE_NAME}_1_unpaired.fastq" \
        "${OUTPUT_DIR}/trimmed_fastqs/${SAMPLE_NAME}_2_trimmed.fastq" \
        "${OUTPUT_DIR}/trimmed_fastqs/${SAMPLE_NAME}_2_unpaired.fastq" \
        ILLUMINACLIP:/home/mahdiel/rnaseq/Trimmomatic-0.36/adapters/combined_adapters.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "Trimming completed."

# --- FastQC on trimmed reads ---
echo "Running FastQC on trimmed reads..."
for R1_FILE in "${OUTPUT_DIR}/trimmed_fastqs"/*_1_trimmed.fastq; do
    SAMPLE_NAME=$(basename "${R1_FILE}" _1_trimmed.fastq)
    R2_FILE="${OUTPUT_DIR}/trimmed_fastqs/${SAMPLE_NAME}_2_trimmed.fastq"
    
    echo "Running FastQC on trimmed reads for sample: ${SAMPLE_NAME}"
    fastqc "${R1_FILE}" "${R2_FILE}" -o "${OUTPUT_DIR}/fastqc_trimmed"
done

echo "FastQC on trimmed reads completed."

# --- Alignment with HISAT2 ---
echo "Aligning reads to reference genome with HISAT2..."
for R1_FILE in "${OUTPUT_DIR}/trimmed_fastqs"/*_1_trimmed.fastq; do
    SAMPLE_NAME=$(basename "${R1_FILE}" _1_trimmed.fastq)
    R2_FILE="${OUTPUT_DIR}/trimmed_fastqs/${SAMPLE_NAME}_2_trimmed.fastq"

    # Check if R2 file exists
    if [[ ! -f "${R2_FILE}" ]]; then
        echo "Warning: ${SAMPLE_NAME}_2_trimmed.fastq not found, skipping this sample."
        continue
    fi

    echo "Processing sample: ${SAMPLE_NAME}"

    # HISAT2 alignment
    hisat2 -p 4 -x "${GENOME_DIR}" \
           -1 "${R1_FILE}" \
           -2 "${R2_FILE}" \
           -S "${OUTPUT_DIR}/alignments/${SAMPLE_NAME}.sam"

    # Convert SAM to BAM and sort
    samtools view -bS "${OUTPUT_DIR}/alignments/${SAMPLE_NAME}.sam" > "${OUTPUT_DIR}/alignments/${SAMPLE_NAME}.bam"
    samtools sort "${OUTPUT_DIR}/alignments/${SAMPLE_NAME}.bam" -o "${OUTPUT_DIR}/sorted_bams/${SAMPLE_NAME}_sorted.bam"

    # Index BAM file
    samtools index "${OUTPUT_DIR}/sorted_bams/${SAMPLE_NAME}_sorted.bam"

    # Remove the SAM file to save space
    rm "${OUTPUT_DIR}/alignments/${SAMPLE_NAME}.sam"
done

# 5) featureCounts (gene-level counts)
GTF="${7:-data/reference/annotation.gtf}"
mkdir -p "${OUT_DIR}/counts"

if [[ -f "${GTF}" ]]; then
  log "Running featureCounts"
  featureCounts -T "${THREADS}" -p -B -C \
    -a "${GTF}" \
    -o "${OUT_DIR}/counts/gene_counts.txt" \
    "${OUT_DIR}"/sorted_bams/*_sorted.bam \
    >> "${OUT_DIR}/logs/featurecounts.log" 2>&1
  log "Counts written to ${OUT_DIR}/counts/gene_counts.txt"
else
  log "Skipping featureCounts: GTF not found at ${GTF}"
fi

echo "Alignment and conversion completed."


