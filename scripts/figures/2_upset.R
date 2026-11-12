#!/usr/bin/env Rscript
# Exact intersections of within-tissue DEGs (BH-adjusted P < 0.05).
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "results/gene_expression"
out_dir <- if (length(args) >= 2) args[2] else file.path(data_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(file.path(data_dir, "tables/DEG_overlap_intersections.csv"))
dat <- dat[order(dat$genes, decreasing = TRUE), ]
tissues <- c("VAF", "LIV", "GM", "BM")

pdf(file.path(out_dir, "UpSet.pdf"), width = 8.3, height = 5.4, useDingbats = FALSE)
layout(matrix(c(1, 2), 2, 1), heights = c(3.5, 1.6))
par(mar = c(0.3, 4.4, 2, 1))
x <- barplot(dat$genes, col = "#474D58", border = NA, axes = FALSE,
             ylim = c(0, max(dat$genes) * 1.16))
axis(2, las = 1)
mtext("Genes in each exact intersection", 2, 2.8)
text(x, dat$genes, labels = dat$genes, pos = 3, cex = 0.7)

par(mar = c(2, 4.4, 0.2, 1))
plot(NA, xlim = c(min(x) - 0.6, max(x) + 0.6), ylim = c(0.5, 4.5),
     axes = FALSE, xlab = "", ylab = "")
axis(2, at = 4:1, labels = tissues, las = 1, tick = FALSE)
for (i in seq_len(nrow(dat))) {
    included <- vapply(dat[i, tissues], function(v) {
        tolower(as.character(v)) %in% c("true", "1")
    }, logical(1))
    points(rep(x[i], 4), 4:1, pch = 16, col = "#DEDEDE", cex = 1.05)
    y <- (4:1)[included]
    if (length(y) > 1) segments(x[i], min(y), x[i], max(y), lwd = 1.7, col = "#474D58")
    points(rep(x[i], length(y)), y, pch = 16, col = "#474D58", cex = 1.05)
}
dev.off()
