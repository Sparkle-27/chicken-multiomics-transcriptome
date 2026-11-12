#!/usr/bin/env Rscript
# Adapted from script/22-lengthPlot-v2.r and script/11-lengthStatCCSbam.r.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 5_length_summary.R one_length_per_line.txt output_prefix")
lengths <- scan(args[1L], quiet = TRUE)
if (!length(lengths) || any(!is.finite(lengths) | lengths < 1 | lengths != floor(lengths))) {
    stop("Supply positive integer sequence lengths.")
}
descending <- sort(lengths, decreasing = TRUE)
quartiles <- quantile(lengths, c(0.25, 0.75), names = FALSE)
stats <- data.frame(sequences = length(lengths), mean_nt = mean(lengths), median_nt = median(lengths),
                    q1_nt = quartiles[1L], q3_nt = quartiles[2L], iqr_nt = diff(quartiles),
                    n50_nt = descending[which(cumsum(descending) >= sum(descending) / 2)[1L]])
dir.create(dirname(args[2L]), recursive = TRUE, showWarnings = FALSE)
write.csv(stats, paste0(args[2L], ".csv"), row.names = FALSE)
pdf(paste0(args[2L], ".pdf"), width = 6, height = 5)
hist(lengths, breaks = "FD", main = "Iso-Seq sequence lengths", xlab = "Length (nt)", col = "steelblue3")
abline(v = stats$median_nt, lty = 2, lwd = 2)
invisible(dev.off())
