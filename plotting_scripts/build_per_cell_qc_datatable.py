#!/usr/bin/env python3
"""
Build a cluster-level per-cell QC datatable by concatenating each subsample's
per_cell_qc.tsv from a cell type's directories in the primary pseudobulk output.

This produces the input to Step 1 of the QC-filtering pipeline (see README.md):
plot_per_cell_qc.R and explore_qc_thresholds.R both operate on the file this
script writes.

Each `annotation-{cell_type}-IGVF*` directory in the primary pseudobulk root
carries its own per-subsample per_cell_qc.tsv with identical columns; this
script concatenates them, keeping a single header, and skips any 0-byte
per_cell_qc.tsv files it encounters.

Usage:
    python build_per_cell_qc_datatable.py \
        --pseudobulks /path/to/{dataset}/pseudobulks \
        --cell-type   k562 \
        --out         datatables/{dataset}_data/k562_per_cell_qc.tsv
"""

import argparse
import glob
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--pseudobulks", required=True,
                    help="Path to the dataset's pseudobulks/ directory.")
    p.add_argument("--cell-type", required=True,
                    help="Cell type identifier -- the string after 'annotation-' "
                         "in the pseudobulk directory name.")
    p.add_argument("--out", required=True,
                    help="Output path for the concatenated per-cell QC datatable.")
    return p.parse_args()


def find_subsample_dirs(pseudobulks_dir, cell_type):
    """Subsample-level pseudobulks (annotation-{cell_type}-IGVF*), falling back
    to a single non-subsampled directory (annotation-{cell_type}) when the
    pseudobulk wasn't generated at the subsample level."""
    dirs = sorted(glob.glob(os.path.join(pseudobulks_dir, f"annotation-{cell_type}-IGVF*")))
    if dirs:
        return dirs
    single = os.path.join(pseudobulks_dir, f"annotation-{cell_type}")
    if os.path.isdir(single):
        return [single]
    return []


def main():
    args = parse_args()

    dirs = find_subsample_dirs(args.pseudobulks, args.cell_type)
    if not dirs:
        sys.exit(f"[error] No annotation-{args.cell_type}[-IGVF*] directories "
                  f"found under {args.pseudobulks}")

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    header = None
    n_rows = 0
    n_used = 0
    with open(args.out, "w") as out_fh:
        for d in dirs:
            qc_path = os.path.join(d, "per_cell_qc.tsv")
            if not os.path.exists(qc_path):
                print(f"[warn] No per_cell_qc.tsv in {d}, skipping", file=sys.stderr)
                continue
            if os.path.getsize(qc_path) == 0:
                print(f"[warn] Empty per_cell_qc.tsv in {d}, skipping", file=sys.stderr)
                continue
            with open(qc_path) as in_fh:
                this_header = in_fh.readline()
                if header is None:
                    header = this_header
                    out_fh.write(header)
                elif this_header != header:
                    sys.exit(f"[error] Column mismatch in {qc_path}:\n"
                              f"  expected: {header.strip()}\n"
                              f"  got:      {this_header.strip()}")
                for line in in_fh:
                    out_fh.write(line)
                    n_rows += 1
            n_used += 1

    if header is None:
        os.remove(args.out)
        sys.exit(f"[error] No usable per_cell_qc.tsv files found for cell type "
                  f"'{args.cell_type}' under {args.pseudobulks}")

    plural = "y" if n_used == 1 else "ies"
    print(f"[info] Wrote {n_rows} cells from {n_used} subsample director{plural} "
          f"to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
