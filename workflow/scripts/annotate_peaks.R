# Peak annotation with ChIPseeker (as in the mcf7 analyses 02, 05, 14, 15),
# but with a TxDb built from reference.gtf so it works for any genome, and
# gene names taken from the GTF instead of an OrgDb.
#   annotatePeak(tssRegion = +/- annotate.tss_region) per peak set
#   -> <set>.annotation.tsv, a feature-distribution summary, and plots.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(rtracklayer)
  library(ggplot2)
})

names <- unlist(snakemake@params[["names"]])
beds <- unlist(snakemake@input[["beds"]])
tables <- unlist(snakemake@output[["tables"]])
tss <- as.numeric(snakemake@params[["tss"]])

gtf <- snakemake@input[["gtf"]]
txdb <- suppressWarnings(makeTxDbFromGFF(gtf, format = "gtf"))
genes <- import(gtf)
genes <- genes[genes$type == "gene"]
gene_names <- if ("gene_name" %in% colnames(mcols(genes))) {
  setNames(genes$gene_name, genes$gene_id)
} else setNames(genes$gene_id, genes$gene_id)

summaries <- list(); anno_list <- list()
for (i in seq_along(names)) {
  nm <- names[i]
  bed <- if (file.size(beds[i]) > 0) {
    read.delim(beds[i], header = FALSE, comment.char = "#", colClasses = "character")
  } else data.frame()
  if (nrow(bed) == 0) {
    cat(nm, ": no peaks\n")
    write.table(data.frame(), tables[i], sep = "\t", quote = FALSE, row.names = FALSE)
    next
  }
  gr <- GRanges(bed[[1]], IRanges(as.integer(bed[[2]]) + 1, as.integer(bed[[3]])))
  anno <- annotatePeak(gr, TxDb = txdb, tssRegion = c(-tss, tss), verbose = FALSE)
  df <- as.data.frame(anno)
  df$gene_name <- unname(gene_names[df$geneId])
  write.table(df, tables[i], sep = "\t", quote = FALSE, row.names = FALSE)
  stat <- anno@annoStat
  summaries[[nm]] <- data.frame(Set = nm, Peaks = length(gr), Feature = stat$Feature,
                                Percent = round(stat$Frequency, 2))
  anno_list[[nm]] <- anno
  cat(nm, ":", length(gr), "peaks annotated\n")
}

summary_df <- if (length(summaries)) do.call(rbind, summaries) else
  data.frame(Set = character(), Peaks = integer(), Feature = character(), Percent = numeric())
write.table(summary_df, snakemake@output[["summary"]], sep = "\t", quote = FALSE,
            row.names = FALSE)

pdf(snakemake@output[["plots"]], width = 9, height = 1.5 + 0.5 * max(1, length(anno_list)))
if (length(anno_list)) {
  print(plotAnnoBar(anno_list) + ggtitle("Genomic feature distribution"))
  print(plotDistToTSS(anno_list, title = "Distance to nearest TSS"))
} else {
  plot.new(); text(0.5, 0.5, "no peaks to annotate")
}
dev.off()
cat("\n"); print(sessionInfo())
