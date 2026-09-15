library(DESeq2)
library(ggplot2)
library(pheatmap)

# ── Capture all output in Snakemake's log file when run via the pipeline ──
if (exists("snakemake") && length(snakemake@log) > 0) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con)
  sink(log_con, type = "message")
}

counts_file   <- "results/counts/gene_counts.txt"
  gene_map_file <- "data/reference/gene_id2name.tsv"
  all_file      <- "results/de/all_results.csv"
  sig_file      <- "results/de/significant_genes.csv"

# ── Load counts ──
raw <- read.delim(counts_file, comment.char = '#')
counts <- raw[, 7:ncol(raw)]
rownames(counts) <- raw$Geneid

# Clean column names
colnames(counts) <- gsub(".*align\\.|.sorted.bam", "", colnames(counts))

# ── Sample metadata ──
sample_condition <- c(
  SRR1039508 = "control", 
  SRR1039512 = "control",
  SRR1039516 = "control", 
  SRR1039520 = "control",
  SRR1039509 = "treated", 
  SRR1039513 = "treated",
  SRR1039517 = "treated", 
  SRR1039521 = "treated"
)

coldata <- data.frame(
  condition = factor(sample_condition[colnames(counts)], levels = c("control", "treated")),
  row.names = colnames(counts)
)

stopifnot(!any(is.na(coldata$condition)))

# ── DESeq2 ──
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","treated","control"), alpha=0.05)

# ── Annotate with gene symbols (rownames are Ensembl IDs) ──
id2name <- read.delim(gene_map_file, header=FALSE,
                      col.names=c("gene_id","symbol"), stringsAsFactors=FALSE)
res$symbol <- id2name$symbol[match(rownames(res), id2name$gene_id)]

res <- res[order(res$padj), ]

# ── Summary ──
summary(res)
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
cat("\nSignificant DE genes:", nrow(sig), "\n")
cat("Upregulated:", sum(sig$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig$log2FoldChange < 0), "\n")

# ── Validate: check for known dex-responsive genes ──
known_genes <- c("CRISPLD2","DUSP1","KLF15","PER1","TSC22D3")
cat("\nValidation - known dex-responsive genes:\n")
for (g in known_genes) {
  idx <- which(res$symbol == g)
  if (length(idx) > 0) {
    r <- res[idx[1], ]
    cat(sprintf("  %s: log2FC=%.2f, padj=%.2e %s\n",
                g, r$log2FoldChange, r$padj,
                ifelse(!is.na(r$padj) && r$padj < 0.05, "SIGNIFICANT", "not significant")))
  } else {
    cat(sprintf("  %s: not found\n", g))
  }
}

# ── Save results ──
dir.create(dirname(all_file), showWarnings = FALSE, recursive = TRUE)
write.csv(as.data.frame(res), all_file, row.names = TRUE)
write.csv(as.data.frame(sig), sig_file, row.names = TRUE)