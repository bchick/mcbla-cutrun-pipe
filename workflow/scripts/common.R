# Shared helpers for the R steps (sourced via snakemake@params$helpers).

start_log <- function() {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log)
  sink(log, type = "message")
  invisible(log)
}

read_contrasts <- function(path) {
  ct <- read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  stopifnot(all(c("group1", "group2", "label") %in% colnames(ct)))
  ct
}

# Add explicit two-group contrasts (group1 vs group2 by Condition, which holds
# the samplesheet group), in the order of contrasts.tsv, under a ~Condition
# design, or ~Factor + Condition
# when `batch` is TRUE (the sample sheet's Factor column holds the batch).
#
# The design matters: contrasts given as group masks (as 08 NB01 did) put
# DiffBind 3.x in its legacy per-contrast mode, where the DESeq2 test sets its
# own size factors (RLE on the contrast's peak counts, or full library size)
# and ignores dba.normalize(). Only the reported Conc/Fold used the stored
# normalization, so every normalization gave identical p-values and FDRs.
add_contrasts <- function(dba_obj, ct, batch = FALSE) {
  design <- if (isTRUE(batch)) "~Factor + Condition" else "~Condition"
  cat("DiffBind design:", design, "\n")
  dba_obj$contrasts <- NULL
  for (i in seq_len(nrow(ct))) {
    g1 <- DiffBind::dba.mask(dba_obj, DiffBind::DBA_CONDITION, ct$group1[i])
    g2 <- DiffBind::dba.mask(dba_obj, DiffBind::DBA_CONDITION, ct$group2[i])
    if (sum(g1) == 0 || sum(g2) == 0) {
      stop(sprintf("Contrast %s: no samples for %s or %s",
                   ct$label[i], ct$group1[i], ct$group2[i]))
    }
    args <- list(dba_obj, contrast = c("Condition", ct$group1[i], ct$group2[i]))
    if (i == 1) args$design <- design
    dba_obj <- do.call(DiffBind::dba.contrast, args)
  }
  dba_obj
}

# dba.analyze without DiffBind's own blacklist/greylist step: reads and peaks
# are already filtered against reference.blacklist (or deliberately not, when
# it is empty), and IgG is handled by the pipeline's IgG gate, not a greylist.
analyze <- function(dba_obj) {
  DiffBind::dba.analyze(dba_obj, bBlacklist = FALSE, bGreylist = FALSE)
}

# Stop if the DESeq2 fit did not use the size factors dba.normalize() stored.
# Guards against DiffBind silently falling back to its legacy per-contrast
# normalization (see add_contrasts).
check_size_factors <- function(dba_obj) {
  stored <- dba_obj$norm$DESeq2$norm.facs
  used <- DESeq2::sizeFactors(dba_obj$DESeq2$DEdata)
  if (is.null(used) || length(used) != length(stored) ||
      !isTRUE(all.equal(unname(stored), unname(used), tolerance = 1e-6))) {
    stop("DESeq2 size factors differ from dba.normalize():\n  stored: ",
         paste(signif(stored, 4), collapse = " "), "\n  used:   ",
         paste(signif(used, 4), collapse = " "))
  }
  cat("DESeq2 size factors match dba.normalize():",
      paste(signif(used, 4), collapse = " "), "\n")
  invisible(used)
}

# Full (th = 1) report for contrast i as a data.frame with stable columns.
report_df <- function(dba_obj, i) {
  res <- DiffBind::dba.report(dba_obj, contrast = i, th = 1, bCounts = FALSE)
  df <- as.data.frame(res)
  conc_cols <- grep("^Conc_", colnames(df))
  if (length(conc_cols) == 2) colnames(df)[conc_cols] <- c("Conc_group1", "Conc_group2")
  df$peak_id <- sprintf("%s:%d-%d", df$seqnames, df$start, df$end)
  df
}

summarise_contrast <- function(df, label, fdr, lfc) {
  sig <- !is.na(df$FDR) & df$FDR < fdr
  data.frame(
    Contrast   = label,
    Total_peaks = nrow(df),
    Gained     = sum(sig & df$Fold > 0),
    Lost       = sum(sig & df$Fold < 0),
    Sig_total  = sum(sig),
    Gained_lfc = sum(sig & df$Fold >= lfc),
    Lost_lfc   = sum(sig & df$Fold <= -lfc),
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
