# File Specification: QC-Filtered RNA Count Matrix (Matrix Market)

## Summary

A `.tar.gz` archive of a sparse RNA count matrix in Matrix Market format, restricted
to QC-passing cells for one cluster. Produced by
`workflow/scripts/filter_rna_counts.py`:

```
python filter_rna_counts.py \
    --qc-guide    {QC filtered barcode list} \
    --pseudobulks {path to pseudobulk directory} \
    --cell-type   {cluster annotation in primary pseudobulk} \
    --out         {out}.mtx \
    --gtf         {IGVFFI9573KOZR} \
    --standard-chromosomes-only
```

`--cell-type` is the cluster annotation string used in the primary pseudobulk's
directory naming — not "cell type" in the IGVF Data Portal sense (Portal "cell type"
usually means a CL ontology term ID).

`--qc-guide` is a **QC filtered barcode list** — a separate, defined file type that
ships alongside this RNA count matrix in the same file set (together with the matching
filtered ATAC fragments file). It defines the exact set of QC-passing cell barcodes
that compose a cluster; it is not derived from the ATAC data.

## Archive contents

`filter_rna_counts.py` itself writes gzipped members (`matrix.mtx.gz`, `barcodes.tsv.gz`,
`features.tsv.gz`); the packaging step that produces the distributed `.tar.gz` decompresses
each one first, so the archive actually contains:

| File | Description |
|---|---|
| `matrix.mtx` | Sparse count matrix, **genes (rows) x cells (columns)** |
| `barcodes.tsv` | One cell barcode per line — defines the column order of the matrix |
| `features.tsv` | One gene symbol per line — defines the row order of the matrix |

Both `.tsv` files are plain, single-column, header-less lists: line *n* of
`barcodes.tsv` labels column *n* of the matrix, and line *n* of `features.tsv`
labels row *n*.

## matrix.mtx

- Standard MatrixMarket sparse-matrix format (`%%MatrixMarket matrix coordinate integer general`).
- Values are **raw UMI counts** (no normalization), summed across every Ensembl gene ID
  that maps to the same gene symbol.
- Orientation is genes x cells (the transpose of the usual cells x genes AnnData layout).

## barcodes.tsv

- Each line is a full cell identifier: `{16-bp cell barcode}_{IGVF 10x lane accession}`,
  e.g. `CCATATTTCGATAACC_IGVFSM4662QKFQ`.
- This is exactly the cluster's QC filtered barcode list — the sibling file (in the same
  file set as this matrix) that defines which QC-passing cell barcodes compose the cluster.
- No duplicate barcodes; the number of lines equals the number of cells in the matrix
  and the number of barcodes in the QC filtered barcode list.

## features.tsv

- Each line is a gene symbol (HGNC-style), not an Ensembl ID.
- Gene symbols are pulled directly from **IGVFFI9573KOZR** (the official GRCh38 / GENCODE 43
  transcriptome reference GTF) metadata for each versioned Ensembl gene ID.
- Restricted to genes on standard chromosomes only: `chr1`-`chr22`, `chrX`, `chrY`, `chrM`.
  Genes found only on scaffolds/patches/alt contigs are dropped; genes with Ensembl IDs on
  both a standard and a nonstandard chromosome retain only their standard-chromosome counts.
- Where multiple Ensembl IDs map to the same gene symbol, their counts are summed into a
  single row for that symbol — so each gene symbol appears exactly once (no duplicates).

## Guarantees

Every archive passes these checks before being written:
1. No duplicate cell barcodes.
2. No duplicate gene symbols.
3. Cell count (matrix columns / `barcodes.tsv` lines) exactly matches the QC filtered barcode list.
4. Every Ensembl ID in the source data was uniquely resolved to a standard-chromosome gene
   symbol (a run fails rather than silently dropping unmapped genes).
