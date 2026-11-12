#!/usr/bin/env Rscript
# PCA scores: VST counts, 500 genes with the highest variance.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "results/gene_expression"
out_dir <- if (length(args) >= 2) args[2] else file.path(data_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(file.path(data_dir, "qc/all_tissues_pca_scores.csv"))
variance <- read.csv(file.path(data_dir, "qc/all_tissues_pca_variance.csv"))$variance_percent
tissues <- c("VAF", "LIV", "GM", "BM")
groups <- paste(rep(c("CB", "TS"), each = 4), rep(tissues, 2), sep = "_")
colors <- rep(c("#F67088", "#96A331", "#35ACA4", "#A38CF4"), 2)
shapes <- c(3, 5, 10, 12, 15, 16, 23, 25)
index <- match(paste(scores$breed, scores$tissue, sep = "_"), groups)

pdf(file.path(out_dir, "PCA.pdf"), width = 8, height = 6, useDingbats = FALSE)
layout(matrix(1:2, nrow = 1), widths = c(1, 0.35))
par(mai = c(0.7, 0.7, 0.1, 0.1), mgp = c(2, 0.5, 0), tcl = -0.3)
plot(scores$PC1, scores$PC2,
     xlab = sprintf("PC1 (%.1f%%)", variance[1]),
     ylab = sprintf("PC2 (%.1f%%)", variance[2]),
     col = colors[index], pch = shapes[index],
     cex = 1.5, lwd = 1.2, cex.lab = 1.2, cex.axis = 1.2)
par(mai = c(0.7, 0, 0.1, 0.1))
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("topleft", legend = groups, col = colors, pch = shapes, bty = "n")
dev.off()
