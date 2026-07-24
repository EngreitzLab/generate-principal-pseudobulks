shhh <- suppressPackageStartupMessages
shhh({
  library(tidyverse)
  library(data.table)
})

# ── usage ──────────────────────────────────────────────────────────────────────
# Single threshold set:
#   Rscript explore_qc_thresholds.R <meta_path> [--flag value ...]
#
# Multiple threshold sets (semicolon-separated, quoted):
#   Rscript explore_qc_thresholds.R <meta_path> \
#     --sets "--tss-min 3; --tss-min 5 --frip-min 0.1; --tss-min 7 --frip-min 0.2"
#
# Show per-subsample breakdown:
#   Rscript explore_qc_thresholds.R <meta_path> [flags] --show-subsamples

# ── defaults ───────────────────────────────────────────────────────────────────

DEFAULTS <- list(
  rna_min        = 1e3,
  rna_max        = Inf,
  gene_min       = 1000,
  gene_max       = Inf,
  pct_mt_max     = 30,
  pct_ribo_max   = 100,
  atac_min       = 1e3,
  atac_max       = Inf,
  tss_enr_min    = 3,
  nuc_signal_max = 1.5,
  pct_dup_max    = 100,
  frip_min       = 0
)

# Thresholds always shown regardless of whether they drop cells.
# All others only shown if overridden from default OR if they drop >0 cells.
ALWAYS_SHOW <- c("rna_min", "atac_min", "tss_enr_min", "nuc_signal_max", "pct_mt_max")

FLAG_MAP <- c(
  "--rna-min"      = "rna_min",
  "--rna-max"      = "rna_max",
  "--gene-min"     = "gene_min",
  "--gene-max"     = "gene_max",
  "--pct-mt-max"   = "pct_mt_max",
  "--pct-ribo-max" = "pct_ribo_max",
  "--atac-min"     = "atac_min",
  "--atac-max"     = "atac_max",
  "--tss-min"      = "tss_enr_min",
  "--nuc-max"      = "nuc_signal_max",
  "--pct-dup-max"  = "pct_dup_max",
  "--frip-min"     = "frip_min"
)

THRESHOLD_LABELS <- c(
  rna_min        = "RNA reads >= ",
  rna_max        = "RNA reads <= ",
  gene_min       = "genes >= ",
  gene_max       = "genes <= ",
  pct_mt_max     = "% mito < ",
  pct_ribo_max   = "% ribo < ",
  atac_min       = "ATAC frags >= ",
  atac_max       = "ATAC frags <= ",
  tss_enr_min    = "TSS enrichment >= ",
  nuc_signal_max = "nucleosomal signal < ",
  pct_dup_max    = "% dup reads < ",
  frip_min       = "FRIP >= "
)

# ── helpers ────────────────────────────────────────────────────────────────────

parse_flag_string <- function(flag_str) {
  p      <- DEFAULTS
  tokens <- strsplit(trimws(flag_str), "\\s+")[[1]]
  tokens <- tokens[tokens != ""]
  if (length(tokens) == 0) return(p)
  i <- 1
  while (i <= length(tokens)) {
    flag <- tokens[i]
    if (!flag %in% names(FLAG_MAP)) stop("Unknown flag: ", flag, call. = FALSE)
    if (i + 1 > length(tokens))     stop("No value for: ", flag, call. = FALSE)
    val <- suppressWarnings(as.numeric(tokens[i + 1]))
    if (is.na(val)) stop("Non-numeric value for ", flag, ": ", tokens[i + 1], call. = FALSE)
    p[[FLAG_MAP[flag]]] <- val
    i <- i + 2
  }
  p
}

label_threshold_set <- function(p) {
  changed <- names(p)[sapply(names(p), function(k) {
    d <- DEFAULTS[[k]]; v <- p[[k]]
    if (is.infinite(d) && is.infinite(v)) return(FALSE)
    !isTRUE(all.equal(d, v))
  })]
  if (length(changed) == 0) return("defaults")
  paste(sapply(changed, function(k) paste0(THRESHOLD_LABELS[k], p[[k]])), collapse = ", ")
}

apply_filters <- function(meta, p) {
  meta %>% filter(
    rna_read_count       >= p$rna_min,       rna_read_count       <= p$rna_max,
    gene_count           >= p$gene_min,       gene_count           <= p$gene_max,
    pct_mito             <  p$pct_mt_max,
    pct_ribo             <  p$pct_ribo_max,
    num_frags            >= p$atac_min,       num_frags            <= p$atac_max,
    nucleosomal_signal   <  p$nuc_signal_max,
    tss_enrichment       >= p$tss_enr_min,
    pct_duplicated_reads <  p$pct_dup_max,
    frip                 >= p$frip_min
  )
}

# Build per-cell failure matrix; return list with:
#   $fail_mat : logical matrix, one col per threshold, TRUE = cell fails
#   $total    : # cells failing each threshold (regardless of other thresholds)
#   $alone    : # cells failing ONLY that threshold (would be saved if it were removed)
build_failure_stats <- function(meta, p) {
  fail_mat <- data.frame(
    rna_min        = meta$rna_read_count       <  p$rna_min,
    rna_max        = meta$rna_read_count       >  p$rna_max,
    gene_min       = meta$gene_count           <  p$gene_min,
    gene_max       = meta$gene_count           >  p$gene_max,
    pct_mt_max     = meta$pct_mito             >= p$pct_mt_max,
    pct_ribo_max   = meta$pct_ribo             >= p$pct_ribo_max,
    atac_min       = meta$num_frags            <  p$atac_min,
    atac_max       = meta$num_frags            >  p$atac_max,
    nuc_signal_max = meta$nucleosomal_signal   >= p$nuc_signal_max,
    tss_enr_min    = meta$tss_enrichment       <  p$tss_enr_min,
    pct_dup_max    = meta$pct_duplicated_reads >= p$pct_dup_max,
    frip_min       = meta$frip                 <  p$frip_min
  )

  any_fail    <- apply(fail_mat, 1, any)
  n_fail_each <- rowSums(fail_mat)  # how many thresholds each cell fails

  total <- colSums(fail_mat & any_fail)  # cells dropped, touching this threshold
  alone <- colSums(fail_mat & n_fail_each == 1)  # cells dropped ONLY by this threshold

  list(fail_mat = fail_mat, total = total, alone = alone)
}

# ── report ─────────────────────────────────────────────────────────────────────

print_report <- function(meta, threshold_sets, show_subsamples) {
  n_total    <- nrow(meta)
  subsamples <- sort(unique(meta$subsample))
  n_subs     <- length(subsamples)

  cat(rep("=", 72), "\n", sep = "")
  cat("QC THRESHOLD EXPLORATION REPORT\n")
  cat(sprintf("Input: %d cells across %d subsamples\n", n_total, n_subs))
  cat(rep("=", 72), "\n\n", sep = "")

  for (s in seq_along(threshold_sets)) {
    p     <- threshold_sets[[s]]
    label <- label_threshold_set(p)

    meta_filt  <- apply_filters(meta, p)
    n_kept     <- nrow(meta_filt)
    n_dropped  <- n_total - n_kept

    # subsample survival
    sub_kept_count <- sum(sapply(subsamples, function(ss)
      any(meta_filt$subsample == ss)))
    sub_summary <- sprintf("%d / %d subsamples with ≥1 cell passing", sub_kept_count, n_subs)

    cat(sprintf("── Set %d: %s\n", s, label))
    cat(rep("-", 72), "\n", sep = "")
    cat(sprintf("  Overall:  %d kept  |  %d dropped (%.1f%% of total)\n",
                n_kept, n_dropped, 100 * n_dropped / n_total))
    cat(sprintf("  Subsamples: %s\n\n", sub_summary))

    # threshold breakdown table
    fs <- build_failure_stats(meta, p)
    active <- names(which(fs$total > 0))  # only show thresholds that drop something

    # determine which thresholds to display:
    #   always show the 5 core thresholds
    #   also show any threshold overridden from its default
    #   also show any threshold (even permissive-default) that drops >0 cells
    overridden <- names(p)[sapply(names(p), function(k) {
      d <- DEFAULTS[[k]]; v <- p[[k]]
      if (is.infinite(d) && is.infinite(v)) return(FALSE)
      !isTRUE(all.equal(d, v))
    })]
    show_thresholds <- union(ALWAYS_SHOW, union(overridden, names(which(fs$total > 0))))

    top_thresh <- if (any(fs$total[show_thresholds] > 0))
      names(which.max(fs$total[show_thresholds])) else NULL

    if (!is.null(top_thresh)) {
      cat(sprintf("  Most stringent threshold (most cells dropped): %s%s  [%d cells]\n\n",
                  THRESHOLD_LABELS[top_thresh], p[[top_thresh]], fs$total[top_thresh]))
    } else {
      cat("  No cells dropped.\n\n")
    }

    cat(sprintf("  %-26s  %8s  %8s\n", "Threshold", "Total", "Alone"))
    cat(sprintf("  %-26s  %8s  %8s\n", strrep("-", 26), strrep("-", 8), strrep("-", 8)))

    # sort: always-show first (in fixed order), then others by total descending
    extra <- setdiff(show_thresholds, ALWAYS_SHOW)
    extra_ord <- extra[order(fs$total[extra], decreasing = TRUE)]
    display_order <- c(ALWAYS_SHOW[ALWAYS_SHOW %in% show_thresholds], extra_ord)

    for (thr in display_order) {
      marker <- if (!is.null(top_thresh) && thr == top_thresh) " ◀" else ""
      cat(sprintf("  %-26s  %8d  %8d%s\n",
                  paste0(THRESHOLD_LABELS[thr], p[[thr]]),
                  fs$total[thr], fs$alone[thr], marker))
    }
    # cat(sprintf("\n  Note: 'Total' = cells failing this threshold (may overlap others).\n"))
    # cat(sprintf("        'Alone' = cells that ONLY fail this threshold.\n\n"))

    # optional per-subsample table
    if (show_subsamples) {
      cat(sprintf("  Per-subsample breakdown:\n"))
      cat(sprintf("  %-24s  %7s  %7s  %7s  %6s  %-30s\n",
                  "Subsample", "Total", "Kept", "Dropped", "% Drop", "Most stringent threshold"))
      cat(sprintf("  %-24s  %7s  %7s  %7s  %6s  %-30s\n",
                  strrep("-", 24), strrep("-", 7), strrep("-", 7),
                  strrep("-", 7), strrep("-", 6), strrep("-", 30)))

      for (ss in subsamples) {
        sub_meta    <- filter(meta, subsample == ss)
        sub_filt    <- apply_filters(sub_meta, p)
        sub_n       <- nrow(sub_meta)
        sub_kept    <- nrow(sub_filt)
        sub_dropped <- sub_n - sub_kept
        sub_pct     <- 100 * sub_dropped / sub_n

        sub_fs      <- build_failure_stats(sub_meta, p)
        sub_active  <- names(which(sub_fs$total > 0))
        if (length(sub_active) > 0) {
          top       <- names(which.max(sub_fs$total[sub_active]))
          top_label <- sprintf("%s%s [%d]", THRESHOLD_LABELS[top], p[[top]], sub_fs$total[top])
        } else {
          top_label <- "—"
        }
        cat(sprintf("  %-24s  %7d  %7d  %7d  %5.1f%%  %s\n",
                    ss, sub_n, sub_kept, sub_dropped, sub_pct, top_label))
      }
      cat("\n")
    }
  }
  cat(rep("=", 72), "\n", sep = "")
}

# ── main ───────────────────────────────────────────────────────────────────────

raw_args <- commandArgs(trailingOnly = TRUE)
if (length(raw_args) < 1) {
  stop("Usage: Rscript explore_qc_thresholds.R <meta_path> [--flag value ...] [--show-subsamples]\n",
       "       Rscript explore_qc_thresholds.R <meta_path> --sets \"set1; set2; ...\"",
       call. = FALSE)
}

meta_path       <- raw_args[1]
rest            <- raw_args[-1]
show_subsamples <- "--show-subsamples" %in% rest
rest            <- rest[rest != "--show-subsamples"]

if (!file.exists(meta_path)) stop("File not found: ", meta_path, call. = FALSE)

meta <- fread(meta_path, sep = "\t", header = TRUE) %>%
  mutate(subsample = as.character(subsample))

sets_idx <- which(rest == "--sets")
if (length(sets_idx) > 0) {
  sets_str       <- rest[sets_idx + 1]
  set_strings    <- strsplit(sets_str, ";")[[1]]
  threshold_sets <- lapply(set_strings, parse_flag_string)
} else {
  threshold_sets <- list(parse_flag_string(paste(rest, collapse = " ")))
}

print_report(meta, threshold_sets, show_subsamples)