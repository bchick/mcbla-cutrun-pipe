# Compare the diff.normalization results with each other normalization and
# flag normalization-sensitive contrasts (mcbla-bulkatac-pipe normcheck, with
# a >20 % rule).
#
# A contrast is normalization-sensitive under a method when its gained OR
# lost count changes by more than `threshold` (fraction, default 0.20) of
# the primary count AND by at least `min_abs` peaks, or when the direction of
# the net change flips (more gained than lost under one method, more lost
# under the other).

source(snakemake@params[["helpers"]])
start_log()
suppressPackageStartupMessages(library(ggplot2))

primary <- snakemake@params[["primary"]]
methods <- unlist(snakemake@params[["methods"]])
labels <- unlist(snakemake@params[["labels"]])
fdr <- as.numeric(snakemake@params[["fdr"]])
thr <- as.numeric(snakemake@params[["threshold"]])
min_abs <- as.numeric(snakemake@params[["min_abs"]])

read_tab <- function(p) read.delim(p, stringsAsFactors = FALSE)
tab_path <- function(summary, lab) file.path(dirname(summary), "tables", paste0(lab, "_all.tsv"))
base <- read_tab(snakemake@input[["primary"]])
sig_ids <- function(df) df$peak_id[!is.na(df$FDR) & df$FDR < fdr]
moved <- function(new, old) {
  d <- abs(new - old)
  d >= min_abs & (old == 0 | d > thr * old)
}

rows <- list()
for (k in seq_along(methods)) {
  m <- methods[k]
  other <- snakemake@input[["others"]][k]
  summ <- read_tab(other)
  for (lab in labels) {
    d <- base[base$Contrast == lab, ]
    s <- summ[summ$Contrast == lab, ]
    dt <- read_tab(tab_path(snakemake@input[["primary"]], lab))
    mt <- read_tab(tab_path(other, lab))
    a <- sig_ids(dt); b <- sig_ids(mt)
    jac <- if (length(union(a, b)) == 0) NA else length(intersect(a, b)) / length(union(a, b))
    shared <- intersect(dt$peak_id, mt$peak_id)
    r <- if (length(shared) > 2) {
      cor(dt$Fold[match(shared, dt$peak_id)], mt$Fold[match(shared, mt$peak_id)])
    } else NA
    flip <- sign(s$Gained - s$Lost) * sign(d$Gained - d$Lost) < 0 &
      abs((s$Gained - s$Lost) - (d$Gained - d$Lost)) >= min_abs
    rows[[length(rows) + 1]] <- data.frame(
      Contrast = lab, method = m, primary = primary,
      Gained_primary = d$Gained, Lost_primary = d$Lost,
      Gained = s$Gained, Lost = s$Lost,
      dGained = s$Gained - d$Gained, dLost = s$Lost - d$Lost,
      direction_flip = flip,
      jaccard_sig = round(jac, 3), lfc_pearson = round(r, 4),
      sensitive = moved(s$Gained, d$Gained) | moved(s$Lost, d$Lost) | flip,
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
write_tsv(tab, snakemake@output[["table"]])
print(tab)

verdict <- do.call(rbind, lapply(labels, function(lab) {
  x <- tab[tab$Contrast == lab, ]
  sens <- x$method[x$sensitive]
  anchored <- intersect(sens, c("greenlist", "spikein", "csaw"))
  rec <- if (any(x$direction_flip)) {
    "direction flips between normalizations: a global shift is likely; do not report the depth result alone"
  } else if (length(anchored) > 0 && primary == "depth") {
    "sensitive to an anchor-based normalization: global shift suspected; report the greenlist/spike-in result"
  } else if (length(sens) > 0) {
    "normalization-sensitive: report the methods side by side"
  } else {
    "robust: normalizations agree"
  }
  data.frame(Contrast = lab, primary = primary,
             normalization_sensitive = length(sens) > 0,
             sensitive_methods = if (length(sens)) paste(sens, collapse = ",") else "",
             recommendation = rec, stringsAsFactors = FALSE)
}))
write_tsv(verdict, snakemake@output[["verdict"]])
print(verdict)

long <- rbind(
  data.frame(Contrast = base$Contrast, method = primary, Direction = "Gained", N = base$Gained),
  data.frame(Contrast = base$Contrast, method = primary, Direction = "Lost", N = -base$Lost),
  data.frame(Contrast = tab$Contrast, method = tab$method, Direction = "Gained", N = tab$Gained),
  data.frame(Contrast = tab$Contrast, method = tab$method, Direction = "Lost", N = -tab$Lost)
)
long <- long[long$Contrast %in% labels, ]
long$Contrast <- factor(long$Contrast, levels = rev(labels))
flag <- verdict$Contrast[verdict$normalization_sensitive]
p <- ggplot(long, aes(Contrast, N, fill = method, alpha = Direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0) + coord_flip() +
  scale_alpha_manual(values = c(Gained = 1, Lost = 0.5)) +
  labs(x = NULL, y = sprintf("Differential peaks (FDR < %g); lost plotted negative", fdr),
       title = "Normalization sensitivity",
       subtitle = if (length(flag)) paste("Sensitive:", paste(flag, collapse = ", ")) else
         "No contrast is normalization-sensitive") +
  theme_bw()
ggsave(snakemake@output[["barplot"]], p, width = 10, height = 2 + 0.45 * length(labels))
