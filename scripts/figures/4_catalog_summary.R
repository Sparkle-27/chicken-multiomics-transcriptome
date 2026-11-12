#!/usr/bin/env Rscript
# Transcript classification, exon counts and transcript lengths.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Usage: Rscript 4_catalog_summary.R catalog_manifest.tsv reference_features.tsv output_dir")
}
manifest <- read.delim(args[1L], colClasses = "character", check.names = FALSE)
required <- c("group", "panel", "classification_tsv")
if (!all(required %in% names(manifest)) || anyNA(manifest[required]) ||
    !all(manifest$panel %in% c("breed_tissue", "merged_tissue"))) {
    stop("Manifest needs group, panel (breed_tissue or merged_tissue), classification_tsv.")
}
if (anyDuplicated(manifest[c("group", "panel")]) ||
    !all(c("breed_tissue", "merged_tissue") %in% manifest$panel)) {
    stop("Supply unique group/panel pairs and both types of catalog.")
}

columns <- c("isoform", "exons", "length", "structural_category", "subcategory")
catalogs <- lapply(seq_len(nrow(manifest)), function(i) {
    tab <- read.delim(manifest$classification_tsv[i], check.names = FALSE)
    if (!all(columns %in% names(tab)) || !nrow(tab) || anyDuplicated(tab$isoform) ||
        anyNA(tab$isoform) || anyNA(tab$structural_category)) {
        stop("Invalid final-catalog classification table: ", manifest$classification_tsv[i])
    }
    data.frame(tab[columns], group = manifest$group[i], panel = manifest$panel[i])
})
catalogs <- do.call(rbind, catalogs)
breed_tissue <- catalogs[catalogs$panel == "breed_tissue", ]
merged_tissue <- catalogs[catalogs$panel == "merged_tissue", ]
categories <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog",
                "novel_not_in_catalog", "Other")
category <- factor(ifelse(breed_tissue$structural_category %in% categories,
                          breed_tissue$structural_category, "Other"), levels = categories)
counts <- table(structural_category = category, group = breed_tissue$group)
proportions <- prop.table(counts, margin = 2L)
subcategory <- ifelse(is.na(merged_tissue$subcategory) | !nzchar(merged_tissue$subcategory),
                      "Unspecified", merged_tissue$subcategory)
sub_counts <- table(subcategory = subcategory, group = merged_tissue$group)

reference <- read.delim(args[2L], check.names = FALSE)
features <- c("exons", "length")
if (!all(features %in% names(reference)) || !nrow(reference)) {
    stop("Reference features need exons and length columns.")
}
shape <- rbind(data.frame(reference[features], group = "GRCg7b Ensembl 110"),
               breed_tissue[c(features, "group")])
if (any(!is.finite(shape$exons) | shape$exons < 1) ||
    any(!is.finite(shape$length) | shape$length < 1)) {
    stop("Exon counts and transcript lengths must be positive and finite.")
}
out <- args[3L]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(counts), file.path(out, "Fig3a_structural_class_counts.csv"), row.names = FALSE)
write.csv(shape, file.path(out, "Fig3bc_exons_and_lengths.csv"), row.names = FALSE)
write.csv(as.data.frame(sub_counts), file.path(out, "Fig3d_structural_subcategories.csv"), row.names = FALSE)

pdf(file.path(out, "Fig3_catalog_summary.pdf"), width = 12, height = 7)
par(mar = c(9, 4, 2, 2))
barplot(proportions, col = hcl.colors(nrow(proportions), "Dark 3"), las = 2,
        ylab = "Proportion", legend.text = rownames(proportions), args.legend = list(bty = "n"))
boxplot(exons ~ group, data = shape, las = 2, xlab = "", ylab = "Exons per isoform")
boxplot(length ~ group, data = shape, las = 2, xlab = "", ylab = "Isoform length (nt)")
barplot(t(sub_counts), beside = TRUE, col = hcl.colors(ncol(sub_counts), "Dark 3"),
        las = 2, cex.names = 0.68, ylab = "Isoforms by structural subcategory",
        legend.text = colnames(sub_counts), args.legend = list(bty = "n"))
invisible(dev.off())
