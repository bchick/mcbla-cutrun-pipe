# DiffBind differential binding for one target under one normalization
# (port of the mcf7 DiffBind notebooks and cutrun_greenlist_helpers.R):
#   depth     : dba.normalize() (DiffBind default, library size)
#   greenlist : dba.normalize(library = size factors, normalize = DBA_NORM_LIB)
#   spikein   : as greenlist, with spike-in size factors
#   csaw      : dba.normalize(DBA_NORM_NATIVE, background = TRUE), csaw-style
#               15 kb background bins
# -> contrasts (~Condition = samplesheet group) -> dba.analyze -> dba.report.

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages({
  library(DiffBind)
  library(ggplot2)
})

method <- snakemake@params[["method"]]
fdr <- as.numeric(snakemake@params[["fdr"]])
lfc <- as.numeric(snakemake@params[["lfc"]])
labels <- unlist(snakemake@params[["labels"]])
tab_dir <- snakemake@output[["tables"]]
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

dba_obj <- readRDS(snakemake@input[["dba"]])
dba_obj$config$cores <- snakemake@threads
ct <- read_contrasts(snakemake@input[["contrasts"]])
ct <- ct[ct$label %in% labels, , drop = FALSE]
stopifnot(nrow(ct) > 0)

samples <- dba_obj$samples$SampleID
if (method %in% c("greenlist", "spikein")) {
  sf <- read.delim(snakemake@input[["sf"]], stringsAsFactors = FALSE)
  idx <- match(samples, sf$SampleID)
  if (anyNA(idx)) stop("no ", method, " size factor for: ",
                       paste(samples[is.na(idx)], collapse = ", "))
  dba_obj <- dba.normalize(dba_obj, method = DBA_DESEQ2,
                           library = sf$size_factor[idx], normalize = DBA_NORM_LIB)
} else if (method == "csaw") {
  dba_obj <- dba.normalize(dba_obj, method = DBA_DESEQ2,
                           normalize = DBA_NORM_NATIVE, background = TRUE)
} else {
  dba_obj <- dba.normalize(dba_obj)
}
norm <- dba.normalize(dba_obj, bRetrieve = TRUE)
print(norm)
write_tsv(data.frame(SampleID = samples,
                     Condition = dba_obj$samples$Condition,
                     lib_size = norm$lib.sizes,
                     norm_factor = norm$norm.factors,
                     method = method),
          snakemake@output[["sizefactors"]])

dba_obj <- add_contrasts(dba_obj, ct, snakemake@params[["batch"]])
dba_obj <- analyze(dba_obj)
check_size_factors(dba_obj)
print(dba.show(dba_obj, bContrasts = TRUE))

summaries <- list()
for (i in seq_len(nrow(ct))) {
  label <- ct$label[i]
  df <- report_df(dba_obj, i)
  write_tsv(df, file.path(tab_dir, paste0(label, "_all.tsv")))
  write_tsv(df[!is.na(df$FDR) & df$FDR < fdr, ], file.path(tab_dir, paste0(label, "_sig.tsv")))
  summaries[[label]] <- summarise_contrast(df, label, fdr, lfc)
}
summary_df <- do.call(rbind, summaries)
summary_df$method <- method
write_tsv(summary_df, snakemake@output[["summary"]])
print(summary_df)

pdf(snakemake@output[["plots"]], width = 9, height = 7)
try(dba.plotHeatmap(dba_obj, main = sprintf("Sample correlation (%s)", method)))
try(dba.plotPCA(dba_obj, label = DBA_CONDITION, attributes = DBA_CONDITION))
for (i in seq_len(nrow(ct))) try(dba.plotMA(dba_obj, contrast = i, sub = ct$label[i]))
long <- rbind(
  data.frame(Contrast = summary_df$Contrast, Direction = "Gained", N = summary_df$Gained),
  data.frame(Contrast = summary_df$Contrast, Direction = "Lost", N = -summary_df$Lost)
)
long$Contrast <- factor(long$Contrast, levels = rev(ct$label))
print(ggplot(long, aes(Contrast, N, fill = Direction)) +
  geom_col() + geom_hline(yintercept = 0) + coord_flip() +
  scale_fill_manual(values = c(Gained = "#1B7837", Lost = "#762A83")) +
  labs(x = NULL, y = sprintf("Differential peaks (FDR < %g)", fdr),
       title = sprintf("Differential binding (%s normalization)", method)) +
  theme_bw())
dev.off()

saveRDS(dba_obj, snakemake@output[["rds"]])
cat("\n"); print(sessionInfo())
