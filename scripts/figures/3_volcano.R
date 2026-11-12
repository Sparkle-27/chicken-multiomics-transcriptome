#!/usr/bin/env Rscript
# TS/CB contrast; BH-adjusted P < 0.05; top 10 labels by absolute log2FC.
options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggrepel)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "results/gene_expression"
out_dir <- if (length(args) >= 2) args[2] else file.path(data_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tissues <- c("VAF", "LIV", "BM", "GM")
alpha <- 0.05
top_n <- 10
colors <- c(Down = "#5CA7BA", Stable = "grey", Up = "#FE4365")
annotation <- read.delim(file.path(data_dir, "inputs/gene_annotation.tsv"))
gene_names <- if ("gene_name" %in% names(annotation)) annotation$gene_name else annotation$gene_symbol
missing_name <- is.na(gene_names) | gene_names == ""
gene_names[missing_name] <- annotation$gene_id[missing_name]
labels <- list()

for (tissue in tissues) {
    dat <- read.csv(file.path(data_dir, "results", tissue, "primary/results_all.csv"))
    dat <- dat[is.finite(dat$log2FoldChange) & is.finite(dat$pvalue), ]
    significant <- !is.na(dat$padj) & dat$padj < alpha
    dat$Significant <- factor(ifelse(significant,
        ifelse(dat$log2FoldChange > 0, "Up", "Down"), "Stable"),
        levels = c("Down", "Stable", "Up"))
    dat$logp <- -log10(pmax(dat$pvalue, 1e-300))
    dat$gene_name <- gene_names[match(dat$gene_id, annotation$gene_id)]
    missing_name <- is.na(dat$gene_name) | dat$gene_name == ""
    dat$gene_name[missing_name] <- dat$gene_id[missing_name]

    selected <- dat[significant, ]
    selected <- selected[order(-abs(selected$log2FoldChange), selected$padj,
                               selected$pvalue, selected$gene_id), ]
    selected <- head(selected, top_n)
    labels[[tissue]] <- data.frame(tissue = tissue, rank = seq_len(nrow(selected)),
        selected[c("gene_id", "gene_name", "log2FoldChange", "baseMean", "pvalue", "padj")])
    dat$label <- ifelse(dat$gene_id %in% selected$gene_id, dat$gene_name, "")
    dat <- dat[order(dat$Significant != "Stable"), ]

    p <- ggplot(dat, aes(log2FoldChange, logp)) +
        geom_point(aes(color = Significant), size = 2) +
        geom_text_repel(aes(label = label), size = 4, fontface = "italic",
            seed = 123, max.overlaps = Inf, max.time = 10,
            box.padding = 0.4, point.padding = 0.2,
            min.segment.length = 0, segment.color = "grey40", segment.size = 0.3) +
        scale_color_manual(values = colors, drop = FALSE) +
        scale_x_continuous(expand = expansion(mult = c(0.12, 0.16))) +
        labs(x = "log2(fold change)", y = "-log10 (p-value)", title = tissue) +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.text = element_text(face = "bold", color = "black", size = 18),
            axis.title = element_text(face = "bold", color = "black", size = 24),
            plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
            legend.position = "bottom",
            legend.key.width = grid::unit(0.8, "cm"),
            legend.key.height = grid::unit(0.8, "cm"),
            legend.text = element_text(face = "bold", size = 18),
            legend.title = element_text(face = "bold", size = 18))
    ggsave(file.path(out_dir, paste0(tissue, "_volcano.pdf")), p,
           width = 8, height = 8, device = "pdf", useDingbats = FALSE)
}
write.csv(do.call(rbind, labels), file.path(out_dir, "volcano_top10_genes.csv"), row.names = FALSE)
