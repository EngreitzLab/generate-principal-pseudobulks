shhh <- suppressPackageStartupMessages
shhh({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(scales)
})

# ── helpers ────────────────────────────────────────────────────────────────────

read_meta <- function(meta_path) {
  meta <- fread(meta_path, sep = "\t", header = TRUE) %>%
    mutate(
      nUMI_non_MT = rna_read_count * (100 - pct_mito) / 100,
      subsample   = as.character(subsample)
    ) %>%
    arrange(subsample)
  meta$subsample <- factor(meta$subsample,
                           levels = unique(meta$subsample), ordered = TRUE)
  meta
}

# ── RNA QC plots ───────────────────────────────────────────────────────────────

make_rna_qc_plots <- function(meta_path, cp, greys,
                               rna_min, rna_max,
                               gene_min, gene_max,
                               pct_mt_max, pct_ribo_max,
                               rna_qc_out) {
  meta <- read_meta(meta_path)

  pdf(rna_qc_out, width = 10, height = 12)

  # cells per subsample
  cell_count <- table(meta$subsample) %>% as.data.frame()
  colnames(cell_count) <- c("subsample", "n_cells")

  p1 <- ggplot(cell_count, aes(x = subsample, y = n_cells)) +
    geom_bar(stat = "identity", fill = greys[3]) +
    geom_text(aes(label = n_cells), position = position_dodge(width = 0.9),
              vjust = -0.25, size = 3.5, angle = 45) +
    labs(title = "Cells per subsample", x = "Subsample", y = "# cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # UMI density
  p2 <- ggplot(meta, aes(x = rna_read_count)) +
    geom_density(linewidth = 1, fill = cp[10], color = cp[10], alpha = 0.5) +
    geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
    scale_x_log10(limits = c(100, 70000), n.breaks = 10) +
    labs(title = "RNA read count density", x = "# RNA reads", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # non-MT UMI density
  p3 <- ggplot(meta, aes(x = nUMI_non_MT)) +
    geom_density(linewidth = 1, fill = cp[10], color = cp[10], alpha = 0.5) +
    scale_x_log10(limits = c(100, 70000), n.breaks = 10) +
    labs(title = "Non-mito read count density", x = "# non-mito RNA reads", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # gene count density
  p4 <- ggplot(meta, aes(x = gene_count)) +
    geom_density(linewidth = 1, fill = cp[10], color = cp[10], alpha = 0.5) +
    geom_vline(xintercept = c(gene_min, gene_max), color = greys[5], linetype = "dashed") +
    scale_x_log10(limits = c(100, 25000), n.breaks = 10) +
    labs(title = "# genes/cell distribution", x = "# genes per cell", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "none")

  # log10 genes per subsample
  p5 <- ggplot(meta, aes(x = subsample, y = log10(gene_count))) +
    geom_boxplot(fill = cp[10], outlier.size = 1.5, outlier.shape = 16,
                 outlier.color = greys[3]) +
    geom_hline(yintercept = log10(c(gene_min, gene_max)), color = greys[5], linetype = "dashed") +
    labs(title = "log10(# genes/cell) per subsample",
         y = "log10(# genes per cell)", x = "Subsample") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # genes x reads, colored by % mito
  p6 <- meta %>% arrange(pct_mito) %>%
    ggplot(aes(x = rna_read_count, y = gene_count)) +
    geom_point(aes(color = pct_mito), shape = 16, size = 1) +
    geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
    geom_hline(yintercept = c(gene_min, gene_max), color = greys[5], linetype = "dashed") +
    scale_color_gradient(low = "#d3a9ce", high = "#430b4e") +
    scale_x_log10(limits = c(100, 70000), n.breaks = 10) +
    scale_y_log10(limits = c(100, 25000), n.breaks = 10) +
    labs(title = "Genes x reads, colored by % mito",
         x = "# RNA reads per cell", y = "# genes per cell", color = "% mito") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          aspect.ratio = 1, legend.position = "right")

  # genes x reads, colored by % ribo
  p7 <- meta %>% arrange(pct_ribo) %>%
    ggplot(aes(x = rna_read_count, y = gene_count)) +
    geom_point(aes(color = pct_ribo), shape = 16, size = 1) +
    geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
    geom_hline(yintercept = c(gene_min, gene_max), color = greys[5], linetype = "dashed") +
    scale_color_gradient(low = "#d3a9ce", high = "#430b4e") +
    scale_x_log10(limits = c(100, 70000), n.breaks = 10) +
    scale_y_log10(limits = c(100, 25000), n.breaks = 10) +
    labs(title = "Genes x reads, colored by % ribo",
         x = "# RNA reads per cell", y = "# genes per cell", color = "% ribo") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          aspect.ratio = 1, legend.position = "none")

  # reads x % mito, colored by % ribo
  p7.5 <- meta %>% arrange(pct_ribo) %>%
    ggplot(aes(x = rna_read_count, y = pct_mito)) +
    geom_point(aes(color = pct_ribo), shape = 16, size = 1) +
    geom_vline(xintercept = c(rna_min, rna_max), color = greys[5], linetype = "dashed") +
    geom_hline(yintercept = pct_mt_max, color = greys[5], linetype = "dashed") +
    scale_color_gradient(low = "#d3a9ce", high = "#430b4e") +
    scale_x_log10(limits = c(100, 70000), n.breaks = 10) +
    labs(title = "Reads x % mito, colored by % ribo",
         x = "# RNA reads per cell", y = "% mitochondrial reads", color = "% ribo") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          aspect.ratio = 1, legend.position = "right")

  # % mito density
  p8 <- ggplot(meta, aes(x = pct_mito)) +
    geom_density(linewidth = 1, fill = cp[10], color = cp[10], alpha = 0.5) +
    geom_vline(xintercept = pct_mt_max, color = greys[5], linetype = "dashed") +
    scale_x_log10(n.breaks = 8) +
    labs(title = "% mito density", x = "% mitochondrial reads", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "none")

  # % ribo density
  p9 <- ggplot(meta, aes(x = pct_ribo)) +
    geom_density(linewidth = 1, fill = cp[10], color = cp[10], alpha = 0.5) +
    geom_vline(xintercept = pct_ribo_max, color = greys[5], linetype = "dashed") +
    scale_x_log10(n.breaks = 10) +
    labs(title = "% ribo density", x = "% ribosomal reads", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "none")

  # log10(genes / reads)
  p10 <- meta %>%
    mutate(log10_genes_per_read = log10(gene_count / rna_read_count)) %>%
    ggplot(aes(x = log10_genes_per_read)) +
    geom_density(linewidth = 1, fill = cp[10], color = cp[10], alpha = 0.5) +
    labs(title = "log10(genes per read)", x = "log10(# genes / # reads)", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "none")

  figure <- ggarrange(
    ggarrange(p1, p2, p3, ncol = 3, widths = c(2, 1.5, 1.5)),
    ggarrange(p4, p5, ncol = 2, widths = c(1, 1.5)),
    ggarrange(p6, p7, p7.5, ncol = 3, widths = c(1.3, 1, 1.2)),
    ggarrange(p8, p9, p10, ncol = 3, widths = c(1, 1, 1)),
    nrow = 4, heights = c(1.2, 1, 1.3, 1))
  print(figure)
  dev.off()
}

# ── apply RNA filter ───────────────────────────────────────────────────────────

apply_rna_qc_filter <- function(meta_path, greys,
                                 rna_min, rna_max,
                                 gene_min, gene_max,
                                 pct_mt_max, pct_ribo_max,
                                 rna_filt_out) {
  meta <- read_meta(meta_path)

  meta_filt <- meta %>%
    filter(rna_read_count > rna_min, rna_read_count < rna_max,
           gene_count     > gene_min, gene_count     < gene_max,
           pct_mito       < pct_mt_max,
           pct_ribo       < pct_ribo_max)

  cell_count <- table(meta_filt$subsample) %>% as.data.frame()
  colnames(cell_count) <- c("subsample", "n_cells")

  p1 <- ggplot(cell_count, aes(x = subsample, y = n_cells)) +
    geom_bar(stat = "identity", fill = greys[3]) +
    geom_text(aes(label = n_cells), position = position_dodge(width = 0.9),
              vjust = -0.25, size = 3.5, angle = 45) +
    labs(title = "Cells per subsample\nafter RNA QC filtering",
         x = "Subsample", y = "# cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")
  ggsave(rna_filt_out, p1, width = 4, height = 3.5)
}

# ── ATAC QC plots and filter ───────────────────────────────────────────────────

make_atac_qc_plots_and_filter <- function(meta_path, greys, cp,
                                           atac_min, atac_max,
                                           tss_enr_min, nuc_signal_max,
                                           pct_dup_max, frip_min,
                                           atac_qc_out, atac_filt_out) {
  meta <- read_meta(meta_path)

  pdf(atac_qc_out, width = 12, height = 10)

  cell_count <- table(meta$subsample) %>% as.data.frame()
  colnames(cell_count) <- c("subsample", "n_cells")

  p1 <- ggplot(cell_count, aes(x = subsample, y = n_cells)) +
    geom_bar(stat = "identity", fill = greys[3]) +
    geom_text(aes(label = n_cells), position = position_dodge(width = 0.9),
              vjust = -0.25, size = 3.5, angle = 45) +
    labs(title = "Cells per subsample", x = "Subsample", y = "# cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # fragment count density
  p2 <- ggplot(meta, aes(x = num_frags)) +
    geom_density(linewidth = 1, fill = cp[7], color = cp[7], alpha = 0.5) +
    geom_vline(xintercept = c(atac_min, atac_max), color = greys[5], linetype = "dashed") +
    scale_x_continuous(n.breaks = 15) +
    labs(title = "Fragment density", x = "# fragments", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # log10 fragments per subsample
  p3 <- ggplot(meta, aes(x = subsample, y = log10(num_frags))) +
    geom_boxplot(fill = cp[7], outlier.size = 1.5, outlier.shape = 16,
                 outlier.color = greys[3]) +
    geom_hline(yintercept = log10(c(atac_min, atac_max)), color = greys[5], linetype = "dashed") +
    labs(title = "log10(# fragments/cell) per subsample",
         y = "log10(# fragments per cell)", x = "Subsample") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # nucleosomal signal density
  p4 <- ggplot(meta, aes(x = nucleosomal_signal)) +
    geom_density(linewidth = 1, fill = cp[7], color = cp[7], alpha = 0.5) +
    geom_vline(xintercept = nuc_signal_max, color = greys[5], linetype = "dashed") +
    labs(title = "Nucleosomal signal density", x = "Nucleosomal signal", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "none")

  # TSS enrichment per subsample
  p5 <- ggplot(meta, aes(x = subsample, y = tss_enrichment)) +
    geom_boxplot(fill = cp[7], outlier.size = 1.5, outlier.shape = 16,
                 outlier.color = greys[3]) +
    geom_hline(yintercept = tss_enr_min, color = greys[5], linetype = "dashed") +
    labs(title = "TSS enrichment per subsample",
         x = "Subsample", y = "TSS enrichment") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  # fragments x TSS enrichment (2D density)
  num_bins <- 8
  contour_colors <- colorRampPalette(c("#ffffff", "#c5e5fb", "#002359"))(num_bins)

  p6 <- ggplot(meta, aes(x = num_frags, y = tss_enrichment)) +
    geom_density_2d_filled(contour_var = "density", bins = num_bins) +
    geom_vline(xintercept = c(atac_min, atac_max), color = greys[5], linetype = "dashed") +
    geom_hline(yintercept = tss_enr_min, color = greys[5], linetype = "dashed") +
    scale_fill_manual(values = contour_colors) +
    scale_x_log10(limits = c(100, 1e5), n.breaks = 10) +
    labs(title = "Fragments x TSS enrichment",
         x = "# fragments", y = "TSS enrichment", fill = "Density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          aspect.ratio = 1, legend.position = "right")

  # % duplicated reads density
  p7 <- ggplot(meta, aes(x = pct_duplicated_reads)) +
    geom_density(linewidth = 1, fill = cp[7], color = cp[7], alpha = 0.5) +
    geom_vline(xintercept = pct_dup_max, color = greys[5], linetype = "dashed") +
    labs(title = "% duplicated reads density",
         x = "% duplicated ATAC reads", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "none")

  # FRIP density
  p8 <- ggplot(meta, aes(x = frip)) +
    geom_density(linewidth = 1, fill = cp[7], color = cp[7], alpha = 0.5) +
    geom_vline(xintercept = frip_min, color = greys[5], linetype = "dashed") +
    labs(title = "FRIP density",
         x = "Fraction of reads in peaks (FRIP)", y = "Cell density") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "none")

  # FRIP per subsample
  p9 <- ggplot(meta, aes(x = subsample, y = frip)) +
    geom_boxplot(fill = cp[7], outlier.size = 1.5, outlier.shape = 16,
                 outlier.color = greys[3]) +
    geom_hline(yintercept = frip_min, color = greys[5], linetype = "dashed") +
    labs(title = "FRIP per subsample", x = "Subsample", y = "FRIP") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")

  figure <- ggarrange(
    ggarrange(p1, p2, p3, ncol = 3, widths = c(1.5, 1, 1)),
    ggarrange(p4, p5, p6, ncol = 3, widths = c(1, 1, 1.3)),
    ggarrange(p7, p8, p9, ncol = 3, widths = c(1, 1, 1.3)),
    nrow = 3, heights = c(1, 1, 1))
  print(figure)
  dev.off()

  # cells per subsample after ATAC filter
  meta_filt <- meta %>%
    filter(num_frags            > atac_min,       num_frags            < atac_max,
           nucleosomal_signal   < nuc_signal_max,
           tss_enrichment       > tss_enr_min,
           pct_duplicated_reads < pct_dup_max,
           frip                 > frip_min)

  cell_count <- table(meta_filt$subsample) %>% as.data.frame()
  colnames(cell_count) <- c("subsample", "n_cells")

  p_filt <- ggplot(cell_count, aes(x = subsample, y = n_cells)) +
    geom_bar(stat = "identity", fill = greys[3]) +
    geom_text(aes(label = n_cells), position = position_dodge(width = 0.9),
              vjust = -0.25, size = 3.5, angle = 45) +
    labs(title = "Cells per subsample\nafter ATAC QC filtering",
         x = "Subsample", y = "# cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")
  ggsave(atac_filt_out, p_filt, width = 4, height = 3.5)
}

# ── apply all filters, save barcodes and stats ────────────────────────────────

apply_all_filters <- function(meta_path, greys,
                               rna_min, rna_max,
                               gene_min, gene_max,
                               pct_mt_max, pct_ribo_max,
                               atac_min, atac_max,
                               tss_enr_min, nuc_signal_max,
                               pct_dup_max, frip_min,
                               cluster_stats_out, all_filt_out,
                               cluster_barcodes_out, qc_thresholds_out) {
  meta <- read_meta(meta_path)

  meta_filt <- meta %>%
    filter(rna_read_count      > rna_min,       rna_read_count      < rna_max,
           gene_count          > gene_min,       gene_count          < gene_max,
           pct_mito            < pct_mt_max,
           pct_ribo            < pct_ribo_max,
           num_frags           > atac_min,       num_frags           < atac_max,
           nucleosomal_signal  < nuc_signal_max,
           tss_enrichment      > tss_enr_min,
           pct_duplicated_reads < pct_dup_max,
           frip                > frip_min)

  # summary per subsample
  meta_summ <- meta_filt %>%
    group_by(subsample) %>%
    summarize(n_cells            = n(),
              total_fragments    = sum(num_frags),
              total_RNA_reads    = sum(rna_read_count),
              mean_frag_per_cell = total_fragments / n_cells,
              mean_RNA_per_cell  = total_RNA_reads / n_cells,
              mean_frip          = mean(frip),
              mean_tss           = mean(tss_enrichment),
              .groups = "drop") %>%
    as_tibble()
  n_total   <- nrow(meta)
  n_kept    <- nrow(meta_filt)
  n_dropped <- n_total - n_kept
  cat(sprintf("\nFiltering summary:\n  Total cells:   %d\n  Cells kept:    %d (%.1f%%)\n  Cells dropped: %d (%.1f%%)\n\n",
    n_total, n_kept, 100 * n_kept / n_total, n_dropped, 100 * n_dropped / n_total))

  fwrite(meta_summ, cluster_stats_out, sep = "\t", row.names = FALSE,
         col.names = TRUE, quote = FALSE)

  # plot
  p1 <- ggplot(meta_summ, aes(x = subsample, y = n_cells)) +
    geom_bar(stat = "identity", fill = greys[4]) +
    geom_text(aes(label = n_cells), position = position_dodge(width = 0.9),
              vjust = -0.25, size = 3.5, angle = 45) +
    labs(title = "Cells per subsample\nafter RNA+ATAC QC filtering",
         x = "Subsample", y = "# cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none")
  ggsave(all_filt_out, p1, width = 4, height = 3.5)

  # barcode file
  meta_bc <- meta_filt %>% select(barcode, subsample, analysis_accession)
  fwrite(meta_bc, cluster_barcodes_out, sep = "\t", row.names = FALSE,
         col.names = TRUE, quote = FALSE)

  # thresholds
  thresh <- data.frame(
    threshold = c("Minimum RNA reads per cell", "Maximum RNA reads per cell",
                  "Minimum genes per cell",     "Maximum genes per cell",
                  "Maximum % mitochondrial",    "Maximum % ribosomal",
                  "Minimum ATAC fragments",     "Maximum ATAC fragments",
                  "Minimum TSS enrichment",     "Maximum nucleosomal signal",
                  "Maximum % duplicated reads", "Minimum FRIP"),
    value = c(rna_min, rna_max, gene_min, gene_max,
              pct_mt_max, pct_ribo_max, atac_min, atac_max,
              tss_enr_min, nuc_signal_max, pct_dup_max, frip_min))
  fwrite(thresh, qc_thresholds_out, sep = "\t", row.names = FALSE,
         col.names = TRUE, quote = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  usage <- paste(
    "Usage: Rscript qc_per_cell.R <meta_path> <out_dir> [options]",
    "Required:",
    "  meta_path              Path to per_cell_qc.tsv",
    "  out_dir                Output directory",
    "Options (all optional):",
    "  --rna-min      NUM     Min RNA reads per cell         [default: 1000]",
    "  --rna-max      NUM     Max RNA reads per cell         [default: Inf]",
    "  --gene-min     NUM     Min genes per cell             [default: 0]",
    "  --gene-max     NUM     Max genes per cell             [default: Inf]",
    "  --pct-mt-max   NUM     Max % mitochondrial reads      [default: 30]",
    "  --pct-ribo-max NUM     Max % ribosomal reads          [default: 100]",
    "  --atac-min     NUM     Min ATAC fragments per cell    [default: 1000]",
    "  --atac-max     NUM     Max ATAC fragments per cell    [default: Inf]",
    "  --tss-min      NUM     Min TSS enrichment             [default: 3]",
    "  --nuc-max      NUM     Max nucleosomal signal         [default: 1.5]",
    "  --pct-dup-max  NUM     Max % duplicated ATAC reads    [default: 100]",
    "  --frip-min     NUM     Min FRIP                       [default: 0]",
    sep = "\n")

  if (length(args) < 2) stop(usage, call. = FALSE)

  meta_path <- args[1]
  out_dir   <- args[2]

  # defaults
  params <- list(
    rna_min        = 1e3,
    rna_max        = Inf,
    gene_min       = 0,
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

  # parse optional flags
  flag_map <- c(
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

  extra <- args[-(1:2)]
  i <- 1
  while (i <= length(extra)) {
    flag <- extra[i]
    if (!flag %in% names(flag_map)) stop("Unknown option: ", flag, "\n", usage, call. = FALSE)
    if (i + 1 > length(extra)) stop("No value provided for: ", flag, call. = FALSE)
    val <- suppressWarnings(as.numeric(extra[i + 1]))
    if (is.na(val)) stop("Non-numeric value for ", flag, ": ", extra[i + 1], call. = FALSE)
    params[[flag_map[flag]]] <- val
    i <- i + 2
  }

  if (!file.exists(meta_path)) stop("Input file not found: ", meta_path, call. = FALSE)

  c(list(meta_path = meta_path, out_dir = out_dir), params)
}

args          <- parse_args()
meta_path     <- args$meta_path
out_dir       <- args$out_dir
rna_min       <- args$rna_min
rna_max       <- args$rna_max
gene_min      <- args$gene_min
gene_max      <- args$gene_max
pct_mt_max    <- args$pct_mt_max
pct_ribo_max  <- args$pct_ribo_max
atac_min      <- args$atac_min
atac_max      <- args$atac_max
tss_enr_min   <- args$tss_enr_min
nuc_signal_max <- args$nuc_signal_max
pct_dup_max   <- args$pct_dup_max
frip_min      <- args$frip_min

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rna_qc_out        <- file.path(out_dir, "RNA_QC_plots.pdf")
rna_filt_out      <- file.path(out_dir, "cells_per_subsample_after_RNA_QC.pdf")
atac_qc_out       <- file.path(out_dir, "ATAC_QC_plots.pdf")
atac_filt_out     <- file.path(out_dir, "cells_per_subsample_after_ATAC_QC.pdf")
cluster_stats_out <- file.path(out_dir, "filtered_cell_subsample_metrics.tsv")
all_filt_out      <- file.path(out_dir, "cells_per_subsample_after_all_QC.pdf")
qc_thresholds_out <- file.path(out_dir, "qc_thresholds.tsv")
cluster_barcodes_out <- file.path(out_dir, "filtered_barcodes_with_subsamples.tsv.gz")

## color palettes
cp    <- c("#429130", "#2f9a71", "#159594", "#0096a0", "#0083ab",
           "#0075b3", "#006eae", "#5b5da3", "#8d4b9b", "#a64791",
           "#b03e67", "#c5373d", "#d8571f", "#e96a00", "#ca9b23")
greys <- c("#e5e5e9", "#c5cad7", "#96a0b3", "#6e788d", "#435369", "#1c2a43")

## RUN
make_rna_qc_plots(meta_path, cp, greys,
                  rna_min, rna_max, gene_min, gene_max,
                  pct_mt_max, pct_ribo_max, rna_qc_out)

apply_rna_qc_filter(meta_path, greys,
                    rna_min, rna_max, gene_min, gene_max,
                    pct_mt_max, pct_ribo_max, rna_filt_out)

make_atac_qc_plots_and_filter(meta_path, greys, cp,
                               atac_min, atac_max,
                               tss_enr_min, nuc_signal_max,
                               pct_dup_max, frip_min,
                               atac_qc_out, atac_filt_out)

apply_all_filters(meta_path, greys,
                  rna_min, rna_max, gene_min, gene_max,
                  pct_mt_max, pct_ribo_max,
                  atac_min, atac_max,
                  tss_enr_min, nuc_signal_max,
                  pct_dup_max, frip_min,
                  cluster_stats_out, all_filt_out,
                  cluster_barcodes_out, qc_thresholds_out)