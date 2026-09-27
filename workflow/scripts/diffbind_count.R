# DiffBind: build the DBA object and count reads for one target (port of the
# mcf7 DiffBind notebooks, e.g. 06_ap1_dynamics 03_diffbind_fos_pjun.Rmd).
#   dba(sampleSheet) -> dba.count(minOverlap = 2)   (diff.peaks: individual)
#   or dba.count over a fixed peak set             (any other diff.peaks)
# The counted object is saved once and reused under every normalization
# (diff.normalization and the normcheck methods).

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(DiffBind))

sheet <- read.csv(snakemake@input[["sheet"]], stringsAsFactors = FALSE)
cat("DiffBind sample sheet:", nrow(sheet), "samples\n")
print(table(sheet$Condition))

dba_obj <- dba(sampleSheet = sheet)
dba_obj$config$cores <- snakemake@threads

count_args <- list(DBA = dba_obj, minOverlap = snakemake@params[["min_overlap"]])
if (!is.null(snakemake@params[["summits"]])) {
  count_args$summits <- snakemake@params[["summits"]]
}
if (!identical(snakemake@params[["peakset"]], "individual")) {
  peaks <- read.table(snakemake@input[["consensus"]], sep = "\t", header = FALSE,
                      comment.char = "#")
  count_args$peaks <- GenomicRanges::GRanges(
    peaks$V1, IRanges::IRanges(peaks$V2 + 1, peaks$V3))
  cat("Counting over peak set", snakemake@params[["peakset"]], ":", nrow(peaks), "peaks\n")
}
dba_obj <- do.call(dba.count, count_args)
print(dba_obj)

saveRDS(dba_obj, snakemake@output[[1]])
cat("\n"); print(sessionInfo())
