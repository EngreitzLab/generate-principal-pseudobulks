# generate-principal-pseudobulks

A lightweight, self-contained pipeline for QC-filtering IGVF multiome primary
pseudobulks: from raw per-cell QC metrics through a reviewed filtering
decision to a released ATAC fragment file and RNA count matrix per cluster.

This repo is a focused extract of [QC-and-Predictions](https://github.com/kaybrand/QC-and-Predictions),
which additionally wires the output of this pipeline into scE2G and IGVF
Data Portal uploads. It exists so that anyone who wants to reproduce IGVF
principal pseudobulk data or generate new principal pseudobulks without
running scE2G on them too has a focused place to start.

## Pipeline overview

```
primary pseudobulk output                 datatables/{dataset}_data/            plots/{dataset}/{cell_type}/
annotation-{cell_type}-IGVF*/   Step 0     {cell_type}_per_cell_qc.tsv  Step 1    filtered_barcodes_with_subsamples.tsv.gz
  per_cell_qc.tsv             --------->                              --------->  (+ QC plots, thresholds, metrics)
  (per subsample)             concatenate                             manual judgment, using
                                                                       explore_qc_thresholds.R
                                                                       and plot_per_cell_qc.R

                                                                                        |
                                                                                        | Step 2 (Snakefile)
                                                                                        v

                                                    multiome_data/{dataset}/{cell_type}/
                                                      atac_fragments_{dataset}_{cell_type}.tsv.gz(.tbi)
                                                      rna_count_matrix_{dataset}_{cell_type}/{matrix,barcodes,features}.tsv.gz
                                                    config/tables/{dataset}_config.tsv
```

A **cluster** is one `(dataset, cell_type)` pair. `cell_type` is always the
exact string used after `annotation-` in the primary pseudobulk's directory
naming (`annotation-{cell_type}-IGVF*`, or `annotation-{cell_type}` when a
pseudobulk wasn't generated at the subsample level) — not necessarily an
ontology term.

## Environment setup

Two conda environments are used:

| Env file | Used by | Contents |
|---|---|---|
| `workflow/envs/qc_per_cell_env.yaml` (`qc_per_cell`) | Step 1 (R scripts) | r-tidyverse, r-data.table, r-ggplot2, r-ggpubr, r-scales |
| `workflow/envs/filter_multiome_env.yaml` (`filter_multiome`) | Step 2 (Snakefile) | bedops (`sort-bed`), htslib (`bgzip`/`tabix`), anndata, scipy, numpy, pandas, pyyaml, snakemake-minimal |

`scripts/build_per_cell_qc_datatable.py` (Step 0) has no third-party
dependencies — any Python 3 works.

```bash
conda env create -f workflow/envs/qc_per_cell_env.yaml
conda env create -f workflow/envs/filter_multiome_env.yaml
```

---

## Step 0 — Build a cluster-level per-cell QC datatable

Each `annotation-{cell_type}-IGVF*` directory in the primary pseudobulk
output carries its own per-subsample `per_cell_qc.tsv` (same columns
everywhere: `analysis_accession`, `barcode`, `subsample`, `rna_read_count`,
`gene_count`, `pct_mito`, `pct_ribo`, `num_frags`, `pct_duplicated_reads`,
`nucleosomal_signal`, `tss_enrichment`, `frip`). `build_per_cell_qc_datatable.py`
concatenates every subsample's `per_cell_qc.tsv` for one cell type into a
single cluster-level datatable, skipping any 0-byte `per_cell_qc.tsv` files
it encounters.

```bash
python3 scripts/build_per_cell_qc_datatable.py \
    --pseudobulks /path/to/{dataset}/pseudobulks \
    --cell-type   k562 \
    --out         datatables/{dataset}_data/k562_per_cell_qc.tsv
```

This datatable — one row per cell, before any QC filtering — is the input to
Step 1.

## Step 1 — Decide on QC thresholds

This isn't two independent scripts run in sequence — it's one step of
manual judgment, using both scripts together:

1. Use `explore_qc_thresholds.R` to see, for a candidate threshold set, how
   many cells each threshold removes **in total** (cells falling below that
   cutoff anywhere in the dataset) and **alone** (cells that would be kept
   if only that one threshold were dropped). This tells you which
   thresholds are actually doing work and which are redundant with others.
2. Use `plot_per_cell_qc.R` alongside it to look at the underlying QC
   distributions per subsample and see where it makes sense to set QC
   thresholds.
3. Once you've settled on thresholds, run `plot_per_cell_qc.R` **again**
   with your final flags — this run is what writes
   `filtered_barcodes_with_subsamples.tsv.gz` and records your choices in
   `qc_thresholds.tsv`.

While iterating, watch for:
- Is there a bimodal distribution suggesting mixed cell populations of
  different quality?
- Are any subsamples mislabeled or a different cell type?

**Explore thresholds:**
```bash
conda run -n qc_per_cell Rscript scripts/plotting_scripts/explore_qc_thresholds.R \
    datatables/{dataset}_data/{cell_type}_per_cell_qc.tsv \
    --rna-min 1000 --gene-min 1000 --pct-mt-max 30 --tss-min 3

# Compare multiple threshold sets at once:
conda run -n qc_per_cell Rscript scripts/plotting_scripts/explore_qc_thresholds.R \
    datatables/{dataset}_data/{cell_type}_per_cell_qc.tsv \
    --sets "--tss-min 3; --tss-min 5 --pct-mt-max 20; --tss-min 7 --pct-mt-max 15"

# Show per-subsample breakdown:
conda run -n qc_per_cell Rscript scripts/plotting_scripts/explore_qc_thresholds.R \
    datatables/{dataset}_data/{cell_type}_per_cell_qc.tsv \
    --tss-min 5 --show-subsamples
```

**Record the final thresholds:**
```bash
conda run -n qc_per_cell Rscript scripts/plotting_scripts/plot_per_cell_qc.R \
    datatables/{dataset}_data/{cell_type}_per_cell_qc.tsv \
    plots/{dataset}/{cell_type}/ \
    --rna-min 1000 --gene-min 1000 --pct-mt-max 30 \
    --atac-min 1000 --nuc-max 1.5 --tss-min 3
```

Available threshold flags (all optional; see script `--help` output for
current defaults): `--rna-min`, `--rna-max`, `--gene-min`, `--gene-max`,
`--pct-mt-max`, `--pct-ribo-max`, `--atac-min`, `--atac-max`, `--tss-min`,
`--nuc-max`, `--pct-dup-max`, `--frip-min`.

**Caution with `--frip-min`:** FRIP is calculated against peaks called on
primary pseudobulks, not on a cluster level. A primary pseudobulk with very
few cells (often because cell annotation doesn't match the MULTI-seq/
subsample tag) will have a low FRIP, making it a proxy for pseudobulk
concordance rather than a conventional signal about each individual cell's
quality. Treat it as a diagnostic for "does this pseudobulk's annotation
make sense" rather than a clean per-cell threshold to filter on.

`plot_per_cell_qc.R` writes to `plots/{dataset}/{cell_type}/`:

| File | Description |
|---|---|
| `filtered_barcodes_with_subsamples.tsv.gz` | **The QC guide** — see [FILE_SPEC_QC_FILTERED_BARCODE_LIST.md](FILE_SPEC_QC_FILTERED_BARCODE_LIST.md). Sole input to Step 2. |
| `qc_thresholds.tsv` | The exact threshold values used to produce the QC guide above (metadata only — not needed downstream). |
| `RNA_QC_plots.pdf`, `ATAC_QC_plots.pdf` | Diagnostic distributions before filtering. |
| `cells_per_subsample_after_RNA_QC.pdf`, `cells_per_subsample_after_ATAC_QC.pdf`, `cells_per_subsample_after_all_QC.pdf` | Cell counts per subsample at each filtering stage. |
| `filtered_cell_subsample_metrics.tsv` | Per-subsample cell counts after filtering. |

If a filtering strategy needs custom logic beyond metric thresholds (e.g.
retaining only cells with concordant cell-type calls across two annotation
methods, or excluding a specific subsample entirely), construct
`filtered_barcodes_with_subsamples.tsv.gz` by hand from the Step 0 datatable
instead — it just needs the three columns documented in the file spec.

## Step 2 — Filter ATAC fragments and RNA counts

`Snakefile` runs `workflow/scripts/filter_atac_fragments.py` and
`workflow/scripts/filter_rna_counts.py` across every `(dataset, cell_type)`
pair in your config, then writes a combined config table.

Both scripts automatically discover and concatenate *every*
`annotation-{cell_type}-IGVF*` directory matching `--cell-type` under
`--pseudobulks` — so one cluster's output can merge cells across many
separate primary pseudobulk subsamples, not just a single directory.
`filter_rna_counts.py` additionally uses `--gtf` (the IGVF Consortium's
official GENCODE 43 / GRCh38 transcriptome reference, accession
`IGVFFI9573KOZR`) to convert each cell's Ensembl gene IDs to the gene
symbols scE2G expects, summing counts across Ensembl IDs that share a
symbol. That conversion behavior — and whether to convert IDs at all, or
restrict to standard chromosomes — is controlled by the script's own flags
(`--gtf`, `--ensemblIDs-as-genes`, `--standard-chromosomes-only`; see the
argument reference below).

`filter_atac_fragments.py` drops any fragments on a chromosome not listed
in `--chrom-sizes`, with a warning per dropped chromosome — make sure the
sizes file covers every contig you want retained (e.g. `chrEBV`, alt
contigs).

```bash
cp config/config_QC_pseudobulks.example.yaml config/config_QC_pseudobulks.yaml
# edit config/config_QC_pseudobulks.yaml: QC_plots_dir, pseudobulk_dir, out_dir,
# transcriptome, chrom_sizes, and the datasets/cell-types to run.

conda activate filter_multiome
snakemake -n --use-conda   # dry-run first
snakemake --use-conda
```

For each `(dataset, cell_type)` this produces, under `{out_dir}/{dataset}/{cell_type}/`:

| Output | Description |
|---|---|
| `atac_fragments_{dataset}_{cell_type}.tsv.gz(.tbi)` | Filtered, sorted, bgzipped + tabix-indexed ATAC fragments for QC-passing barcodes. |
| `rna_count_matrix_{dataset}_{cell_type}/{matrix.mtx,barcodes.tsv,features.tsv}.gz` | Filtered RNA count matrix (genes x cells, gene symbols). See [FILE_SPEC_RNA_COUNT_MATRIX.md](FILE_SPEC_RNA_COUNT_MATRIX.md). |

and one combined table per dataset at `{out_dir}/config/tables/{dataset}_config.tsv`,
listing every cluster's `rna_matrix_file` / `atac_frag_file` paths (ready to
feed into a downstream pipeline such as scE2G).

Both scripts can also be run directly (outside Snakemake) for a single
cluster:

**`filter_atac_fragments.py`**

| Flag | Required | Description |
|---|---|---|
| `--qc-guide` | yes | Path to the gzipped QC guide (`filtered_barcodes_with_subsamples.tsv.gz` from Step 1). |
| `--pseudobulks` | yes | Path to the dataset's `pseudobulks/` directory. |
| `--cell-type` | yes | Cell type identifier — the string between `annotation-` and `-IGVF` in directory names. |
| `--chrom-sizes` | yes | Chromosome sizes file defining the expected sort order. |
| `--out` | yes | Output path for the filtered, sorted, bgzipped fragment file. |
| `--clean` | no | Write only the 16-bp barcode sequence, stripping the 10x-lane suffix. Default: keep the full barcode. |

**`filter_rna_counts.py`**

| Flag | Required | Description |
|---|---|---|
| `--qc-guide` | yes | Path to the gzipped QC guide. |
| `--pseudobulks` | yes | Path to the dataset's `pseudobulks/` directory. |
| `--cell-type` | yes | Cell type identifier — the string between `annotation-` and `-IGVF` in directory names. |
| `--out` | yes | Output path; format inferred from extension (`.csv`/`.csv.gz`, `.h5ad`/`.h5`, `.mtx`/`.mtx.gz` — the Snakefile uses `.mtx`, writing a matrix directory). |
| `--gtf` | no | GTF for Ensembl ID → gene symbol conversion (the Snakefile passes the IGVF GENCODE 43 reference). |
| `--ensemblIDs-as-genes` | no | Skip gene-symbol conversion and keep Ensembl IDs. Default: convert, summing counts across Ensembl IDs that share a symbol. |
| `--standard-chromosomes-only` | no | Restrict to chr1–22/X/Y/M genes (requires `--gtf`). Always used by the Snakefile. |
| `--log` | no | Path to write a detailed log of GTF mapping decisions. |

---

## Script reference

| Script | Step | Environment |
|---|---|---|
| `scripts/build_per_cell_qc_datatable.py` | 0 | none (stdlib) |
| `scripts/plotting_scripts/explore_qc_thresholds.R` | 1 | `qc_per_cell` |
| `scripts/plotting_scripts/plot_per_cell_qc.R` | 1 | `qc_per_cell` |
| `Snakefile` (→ `workflow/scripts/filter_atac_fragments.py`, `filter_rna_counts.py`) | 2 | `filter_multiome` |

## File specs

- [FILE_SPEC_QC_FILTERED_BARCODE_LIST.md](FILE_SPEC_QC_FILTERED_BARCODE_LIST.md) — the QC guide produced in Step 1.
- [FILE_SPEC_RNA_COUNT_MATRIX.md](FILE_SPEC_RNA_COUNT_MATRIX.md) — the RNA count matrix produced in Step 2.
