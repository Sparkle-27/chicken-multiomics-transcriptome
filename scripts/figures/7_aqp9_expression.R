#!/usr/bin/env Rscript
# Descriptive template for the retained AQP9 transcript example.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: Rscript 7_aqp9_expression.R AQP9_sample_TPM.tsv output_prefix")
}
data <- read.delim(args[1L], check.names = FALSE)
required <- c("sample_id", "tissue", "breed", "transcript_id", "TPM")
if (!all(required %in% names(data)) || !nrow(data) || anyNA(data[required]) ||
    anyDuplicated(data[c("sample_id", "tissue", "transcript_id")]) ||
    any(!is.finite(data$TPM) | data$TPM < 0)) {
    stop("Provide one finite nonnegative TPM per sample, tissue and transcript ID.")
}
out <- args[2L]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(data[required], paste0(out, ".csv"), row.names = FALSE)
pdf(paste0(out, ".pdf"), width = 7, height = 5)
for (tissue in unique(data$tissue)) {
    selected <- data[data$tissue == tissue, ]
    groups <- interaction(selected$transcript_id, selected$breed, drop = TRUE, sep = " / ")
    boxplot(selected$TPM ~ groups, main = paste("AQP9", tissue), ylab = "TPM",
            xlab = "", cex.axis = 0.8, outline = FALSE)
    stripchart(selected$TPM ~ groups, vertical = TRUE, method = "stack",
               pch = 19, cex = 0.7, add = TRUE)
}
invisible(dev.off())
