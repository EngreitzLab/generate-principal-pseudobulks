# File Specification: QC Filtered Barcode List (TSV)

## Summary

A gzipped, header-ed TSV (`TabularFile`) listing the 16-bp 10x cell barcodes that
passed RNA+ATAC QC and compose one annotated cluster. Produced by
`apply_all_filters()` in `plotting_scripts/plot_per_cell_qc.R`:

```
Rscript plot_per_cell_qc.R path/to/primary/pseudobulk/per_cell_qc.tsv plots_dir \
    --rna-min 1000 --gene-min 1000 --pct-mt-max 30 \
    --atac-min 1000 --nuc-max 1.5 --tss-min 3
```

which writes `filtered_barcodes_with_subsamples.tsv.gz` alongside a **QC threshold
Document** (`qc_thresholds.tsv`) from the same call — the two are a matched pair.
The threshold document records the exact cutoff values applied to `per_cell_qc.tsv`
to produce this barcode list; it's useful metadata for understanding the QC
rationale, but the QC filtered barcode list itself is the sole file necessary to
reproduce the filtering choice — it already contains the exact resulting cell set.

This file is the sole determinant of cluster membership passed downstream: it is
consumed as `--qc-guide` by the principal pseudobulking pipeline
(`filter_rna_counts.py`, `filter_atac_fragments.py`) to subset the unfiltered,
per-subsample RNA count matrices and ATAC fragment files down to the QC-passing
barcode set that becomes the cluster's released RNA count matrix and ATAC
fragments file. No other filtering is applied downstream — every barcode listed
here is included, and no barcode absent from this list is.

Note that the QC metric thresholds in the paired QC threshold Document may not be
the only determinant of this file's barcode content. Where additional filtering
choices were made — such as merging two clusters or removing cell barcodes whose
MULTI-seq tag didn't match the cluster's annotation — the IGVF `submitter_comment`
will describe said choices.

## Format

- Gzip-compressed TSV, one header row, no row index.
- One row per QC-passing cell barcode; rows are not required to be sorted.

| Column | Description |
|---|---|
| `barcode` | 10x Genomics 16-bp cell barcode, suffixed with the IGVF accession of the multiplexed 10x lane it was sequenced in: `{16-bp barcode}_{IGVF 10x lane accession}`. |
| `subsample` | IGVF accession of the MULTI-seq-tagged in-vitro system (subsample) this cell was demultiplexed to. |
| `analysis_accession` | IGVF accession of the intermediate analysis set this cell was derived from. |

### Example

```
barcode                          subsample       analysis_accession
AAGCCTTAGGAACGGT_IGVFSM0539NUPM  IGVFSM0143DYAD  IGVFDS1244UUGQ
```

## Guarantees

1. Every `barcode` value is unique (the 10x-lane suffix disambiguates barcodes
   collected on different lanes; no cell appears twice).
2. Row count exactly matches the `n_cells` total for this cluster, and the
   `barcode` column is the same *set* of barcodes as the downstream paired
   **QC-Filtered RNA Count Matrix**'s `barcodes.tsv` — not necessarily in the
   same row order, since that file's row *n* defines column *n* of the matrix
   independent of this file's row order.
3. Every row satisfies every threshold recorded in the paired QC threshold
   Document — both files are written from the same filtered table in the same
   script run.
4. `subsample` and `analysis_accession` are carried through unmodified from
   `per_cell_qc.tsv`; no cell is reassigned to a different subsample or
   analysis set during filtering.
