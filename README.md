# RNA-seq Analysis Pipeline

## Overview

This repository contains a complete RNA-seq analysis pipeline including:

- Quality control (FastQC)
- Adapter trimming (Trimmomatic)
- Read alignment (HISAT2)
- BAM sorting and indexing (SAMtools)
- Gene-level quantification (featureCounts)
- Downstream QC analysis in R (PCA, heatmap, exploratory analysis)

The pipeline is designed to be reproducible and portable across systems.

---

## Pipeline Workflow

FASTQ → QC → Trimming → QC → Alignment → Sorted BAM → Gene Counts → R Analysis

---

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/rnaseq-advanced-pipeline.git
cd rnaseq-advanced-pipeline


Create conda environment:

conda env create -f envs/rnaseq.yml
conda activate rnaseq


Required Input Files

You must provide:

1- Paired-end FASTQ files
Example:

sample1_1.fastq.gz
sample1_2.fastq.gz

2- HISAT2 genome index (prefix)

3- Trimmomatic jar file

3- Adapter FASTA file

4- Gene annotation file (GTF)


bash scripts/01_qc_trim_align.sh \
  <FASTQ_DIR> \
  <OUTPUT_DIR> \
  <HISAT2_INDEX_PREFIX> \
  <TRIMMOMATIC_JAR> \
  <ADAPTERS_FA> \
  <ANNOTATION_GTF> \
  <THREADS>


bash scripts/01_qc_trim_align.sh \
  data/raw \
  results \
  data/reference/hisat2_index/genome \
  tools/trimmomatic.jar \
  tools/adapters.fa \
  data/reference/annotation.gtf \
  4


Output

The pipeline generates:

FastQC reports

Sorted and indexed BAM files

Gene count matrix


results/counts/gene_counts.txt


Downstream R Analysis

After generating gene counts:

Rscript src/r/01_qc_pca_heatmap.R \
  --counts results/counts/gene_counts.txt \
  --out results/qc


Outputs:

PCA plot

Heatmap of top variable genes

Optional DESeq2 results (if metadata provided)


Data Source

## Data Source

RNA-seq data were derived from:
Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown

https://doi.org/10.1038/nprot.2016.095
