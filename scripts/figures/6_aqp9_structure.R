#!/usr/bin/env Rscript
# AQP9 exon structures from reference and Iso-Seq annotations.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop("Usage: Rscript 6_aqp9_structure.R Ensembl110.gtf PacBio_AQP9.gtf gene_symbol output_prefix")
}
attribute <- function(text, key) {
    pattern <- paste0('(?:^|;)[[:space:]]*', key, '[[:space:]]+"([^"]+)"')
    hits <- regmatches(text, regexec(pattern, text, perl = TRUE))
    vapply(hits, function(hit) if (length(hit)) hit[2L] else NA_character_, character(1L))
}
read_exons <- function(path, source, gene = NULL) {
    gtf <- read.delim(path, header = FALSE, comment.char = "#", quote = "",
                      colClasses = "character")
    if (ncol(gtf) != 9L) stop("Expected a 9-column GTF: ", path)
    gtf <- gtf[gtf[[3L]] == "exon", ]
    if (!is.null(gene)) {
        gtf <- gtf[which(attribute(gtf[[9L]], "gene_name") == gene |
                          attribute(gtf[[9L]], "gene_id") == gene), ]
    }
    exons <- data.frame(source = rep(source, nrow(gtf)), seqid = gtf[[1L]], start = as.integer(gtf[[4L]]),
                        end = as.integer(gtf[[5L]]), strand = gtf[[7L]],
                        transcript_id = attribute(gtf[[9L]], "transcript_id"))
    if (!nrow(exons) || anyNA(exons) || any(exons$end < exons$start) ||
        !all(exons$strand %in% c("+", "-"))) stop("Missing or invalid exon records: ", path)
    exons
}
exons <- rbind(read_exons(args[1L], "Reference", args[3L]), read_exons(args[2L], "PacBio"))
if (length(unique(exons$seqid)) != 1L || length(unique(exons$strand)) != 1L) {
    stop("Reference and PacBio exons must use the same reference sequence and strand.")
}
if (max(min(exons$start[exons$source == "Reference"]), min(exons$start[exons$source == "PacBio"])) >
    min(max(exons$end[exons$source == "Reference"]), max(exons$end[exons$source == "PacBio"]))) {
    stop("The supplied reference and PacBio GTF intervals do not overlap.")
}
label <- paste(exons$source, exons$transcript_id)
transcripts <- unique(label)
exons$track <- match(label, transcripts)
dir.create(dirname(args[4L]), recursive = TRUE, showWarnings = FALSE)
write.csv(exons, paste0(args[4L], ".csv"), row.names = FALSE)
pdf(paste0(args[4L], ".pdf"), width = 10, height = max(4, 0.32 * length(transcripts) + 2))
par(mar = c(4, 12, 3, 1))
plot(NA, xlim = range(c(exons$start, exons$end)), ylim = c(length(transcripts) + 0.6, 0.4),
     yaxt = "n", ylab = "", xlab = paste(exons$seqid[1L], "genomic position (bp)"),
     main = paste(args[3L], "transcripts", paste0("(", exons$strand[1L], " strand)")))
axis(2, seq_along(transcripts), transcripts, las = 1, cex.axis = 0.78)
for (i in seq_along(transcripts)) {
    tx <- exons[exons$track == i, ]
    segments(min(tx$start), i, max(tx$end), i, col = "grey35")
    rect(tx$start, i - 0.25, tx$end + 1L, i + 0.25, border = NA,
         col = if (tx$source[1L] == "Reference") "#43796C" else "#596CA4")
}
invisible(dev.off())
