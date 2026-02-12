#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

counts_file <- get_arg("--counts", "results/counts/gene_counts.txt")
meta_file   <- get_arg("--meta",   NA)
out_dir     <- get_arg("--out",    "results/deseq2")
top_n       <- as.integer(get_arg("--top_n", "50"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Packages (prefer installing via conda/renv; this just loads)
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(matrixStats)
})

message("Counts: ", counts_file)
message("Meta:   ", meta_file)
message("Out:    ", out_dir)

# Read featureCounts output
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
count_mat <- as.matrix(fc[, 7:ncol(fc)])
rownames(count_mat) <- fc$Geneid

# Clean column names (remove paths + bam suffix)
colnames(count_mat) <- basename(colnames(count_mat))
colnames(count_mat) <- sub("_sorted\\.bam$", "", colnames(count_mat))
colnames(count_mat) <- sub("\\.bam$", "", colnames(count_mat))

# Filter low-count genes
keep <- rowSums(count_mat >= 10) >= 2
count_mat_f <- count_mat[keep, , drop = FALSE]

# Metadata
if (!is.na(meta_file) && file.exists(meta_file)) {
  coldata <- read.delim(meta_file, check.names = FALSE)
  stopifnot(all(c("sample") %in% colnames(coldata)))
  rownames(coldata) <- coldata$sample
  coldata <- coldata[colnames(count_mat_f), , drop = FALSE]
} else {
  coldata <- data.frame(sample = colnames(count_mat_f))
  rownames(coldata) <- coldata$sample
}

# Decide design: if condition exists with >=2 levels -> run DESeq2; else QC only
run_de <- ("condition" %in% colnames(coldata)) && (length(unique(coldata$condition)) >= 2)

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat_f),
  colData = coldata,
  design = if (run_de) ~ condition else ~ 1
)

# VST
vsd <- varianceStabilizingTransformation(dds, blind = !run_de)
mat_vst <- assay(vsd)

# Heatmap: top variable genes
vars <- rowVars(mat_vst)
ord <- order(vars, decreasing = TRUE)
top_genes <- rownames(mat_vst)[ord][1:min(top_n, nrow(mat_vst))]
mat_top <- mat_vst[top_genes, , drop = FALSE]
mat_top_z <- t(scale(t(mat_top)))

png(file.path(out_dir, "heatmap_top_variable_genes.png"), width = 1400, height = 1200, res = 150)
pheatmap(
  mat_top_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = paste0("Top ", min(top_n, nrow(mat_vst)), " variable genes (VST, z-score)")
)
dev.off()

# PCA
pca <- prcomp(t(mat_vst), scale. = FALSE)
pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2) * 100

if ("condition" %in% colnames(coldata)) {
  pca_df$condition <- coldata[pca_df$sample, "condition"]
}

p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 3, aes(color = if ("condition" %in% colnames(pca_df)) condition else NULL)) +
  geom_text(vjust = -0.6, size = 3) +
  xlab(sprintf("PC1 (%.1f%%)", var_expl[1])) +
  ylab(sprintf("PC2 (%.1f%%)", var_expl[2])) +
  ggtitle("PCA (VST normalized counts)") +
  theme_minimal()

ggsave(filename = file.path(out_dir, "pca_vst.png"), plot = p, width = 8, height = 6, dpi = 160)

# Optional: Differential expression if possible
if (run_de) {
  dds <- DESeq(dds)
  res <- results(dds)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$pvalue), ]
  write.table(res_df, file.path(out_dir, "deseq2_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("DESeq2 results written: ", file.path(out_dir, "deseq2_results.tsv"))
} else {
  message("No valid 'condition' found (or only one level). Skipping DESeq2 differential expression; QC plots saved.")
}

message("Done.")
