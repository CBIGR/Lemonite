#!/usr/bin/env python3
"""
Prepare phosphoproteomics data as an additional LemonTree regulator layer for the
Wang GBM Lemonite analysis, WITHOUT re-clustering (existing modules are kept).

Input : phosphoproteome_normalized sheet of the Wang2021 supplement (mmc3.xlsx).
        Already log2-transformed, median-polished and ComBat batch-corrected by the
        authors, so NO re-normalisation is performed here. We only:
          1. subset to the samples present in the existing LemonTree run,
          2. set missing values (NA) to 0,
          3. z-score each phosphosite (row-wise scale to mean 0, sd 1) so the layer
             matches the scale of every other feature in LemonPreprocessed_complete.txt.
        This mirrors exactly the pipeline's PTM branch in
        nextflow/scripts/Preprocessing_TFA_Proteomics.R  ( t(scale(t(hPTM))) ).

Output (written to --out-dir):
    phospho.txt                          - regulator list (one phosphosite id per line)
    LemonPreprocessed_phospho.txt        - phospho abundance matrix (Gene_symbol, Protein_id, <samples>)
    LemonPreprocessed_complete_phospho.txt - existing complete matrix + phospho rows appended
                                           (this is the -data_file for the regulators task)

The regulator assignment itself (LemonTree -task regulators) is run separately by
run_phospho_regulators.sh against the EXISTING Lemon_out/tight_clusters.txt.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import openpyxl


def sanitize(name: str) -> str:
    """Match the R preprocessing: replace spaces and dashes in row names with '_'."""
    return name.replace(" ", "_").replace("-", "_")


def read_complete_header(path: Path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    # columns: symbol, ensembl_gene_id (or Gene_symbol, Protein_id), then samples
    return hdr[:2], hdr[2:]


def load_phospho_sheet(xlsx: Path, sheet: str, sample_index: dict):
    """Yield (site_id, values-aligned-to-our-samples) for each phosphosite row.

    values is a float array (len == number of our samples); missing entries are NaN
    for now (NA handling / scaling happens later so we can compute per-row stats first).
    """
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    header = list(next(it))
    # metadata columns: site_id, symbol, phosphosites, peptide  -> samples start at idx 4
    meta_n = 4
    phospho_samples = header[meta_n:]
    # position in the phospho row for each of our samples
    col_for_sample = {}
    for our_sample, our_pos in sample_index.items():
        try:
            col_for_sample[our_sample] = meta_n + phospho_samples.index(our_sample)
        except ValueError:
            pass  # sample absent in phospho data (padded later)

    n_our = len(sample_index)
    for row in it:
        site_id = row[0]
        if site_id is None:
            continue
        vals = np.full(n_our, np.nan, dtype=float)
        for our_sample, our_pos in sample_index.items():
            c = col_for_sample.get(our_sample)
            if c is None:
                continue
            v = row[c]
            if v is None or v == "NA" or v == "":
                continue
            try:
                vals[our_pos] = float(v)
            except (TypeError, ValueError):
                continue
        yield str(site_id), vals


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mmc3", required=True, help="Path to 1-s2.0-...-mmc3.xlsx")
    ap.add_argument("--sheet", default="phosphoproteome_normalized")
    ap.add_argument("--complete", required=True,
                    help="Existing LemonPreprocessed_complete.txt (defines samples and is the base to append to)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-valid-frac", type=float, default=0.5,
                    help="Drop phosphosites whose fraction of non-NA samples is below this (default 0.5). "
                         "Set to 0 to keep every site with any signal.")
    ap.add_argument("--top-var", type=int, default=0,
                    help="If >0, after filtering keep only the N most variable phosphosites "
                         "(variance computed on the NA->0 abundances, before z-scoring). "
                         "0 (default) keeps all sites that pass --min-valid-frac.")
    ap.add_argument("--prefix", default="phospho", help="Row-name prefix / reg list file stem (default: phospho)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    complete_path = Path(args.complete)

    id_cols, samples = read_complete_header(complete_path)
    sample_index = {s: i for i, s in enumerate(samples)}
    print(f"[info] complete matrix samples: {len(samples)}")
    print(f"[info] id columns: {id_cols}")

    # ---- load, filter, NA->0, record variance, z-score per row -------------
    # Variance is computed on the NA->0 abundances (pre-scaling) so it can be used
    # for --top-var selection; z-scoring afterwards would flatten every variance to ~1.
    kept_names = []
    kept_rows = []
    kept_var = []
    n_total = n_all_na = n_low_valid = n_zero_var = 0
    for site_id, vals in load_phospho_sheet(Path(args.mmc3), args.sheet, sample_index):
        n_total += 1
        valid = ~np.isnan(vals)
        n_valid = int(valid.sum())
        if n_valid == 0:
            n_all_na += 1
            continue
        if n_valid / len(samples) < args.min_valid_frac:
            n_low_valid += 1
            continue
        row = np.where(valid, vals, 0.0)          # NA -> 0 (matches R: hPTM[is.na] <- 0)
        var = float(row.var(ddof=1))
        sd = np.sqrt(var)
        if sd == 0 or not np.isfinite(sd):
            n_zero_var += 1
            continue
        zrow = (row - row.mean()) / sd            # z-score per phosphosite (t(scale(t(x))))
        kept_names.append(sanitize(site_id))
        kept_rows.append(zrow)
        kept_var.append(var)

    print(f"[info] phosphosites read: {n_total}")
    print(f"[info]   dropped all-NA: {n_all_na}")
    print(f"[info]   dropped < {args.min_valid_frac:.0%} valid samples: {n_low_valid}")
    print(f"[info]   dropped zero-variance: {n_zero_var}")
    print(f"[info]   passed filtering: {len(kept_names)} phosphosites")

    if not kept_names:
        sys.exit("[error] no phosphosites survived filtering")

    # ---- optional: keep only the N most variable phosphosites --------------
    if args.top_var and args.top_var < len(kept_names):
        order = np.argsort(kept_var)[::-1][:args.top_var]   # highest variance first
        order = sorted(order.tolist())                      # keep original row order
        kept_names = [kept_names[i] for i in order]
        kept_rows = [kept_rows[i] for i in order]
        kept_var = [kept_var[i] for i in order]
        print(f"[info]   --top-var {args.top_var}: KEPT {len(kept_names)} most variable phosphosites")
    else:
        print(f"[info]   KEPT: {len(kept_names)} phosphosites (no top-var cap)")

    # Guard against duplicate sanitized ids (LemonTree keys rows by name)
    seen, dup = set(), 0
    uniq_names, uniq_rows = [], []
    for name, row in zip(kept_names, kept_rows):
        if name in seen:
            dup += 1
            continue
        seen.add(name)
        uniq_names.append(name)
        uniq_rows.append(row)
    if dup:
        print(f"[warn] dropped {dup} duplicate phosphosite ids after sanitising")
    kept_names, kept_rows = uniq_names, uniq_rows

    fmt = lambda x: f"{x:.6g}"

    # ---- write phospho-only matrix ----------------------------------------
    phospho_matrix = out_dir / "LemonPreprocessed_phospho.txt"
    with open(phospho_matrix, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene_symbol", "Protein_id"] + samples)
        for name, row in zip(kept_names, kept_rows):
            w.writerow([name, name] + [fmt(v) for v in row])
    print(f"[ok] wrote {phospho_matrix}")

    # ---- write regulator list ---------------------------------------------
    reg_list = out_dir / f"{args.prefix}.txt"
    with open(reg_list, "w") as fh:
        fh.write("\n".join(kept_names) + "\n")
    print(f"[ok] wrote {reg_list}")

    # ---- append phospho rows to the existing complete matrix --------------
    # Base rows use 2 id columns; our phospho id columns fill both slots.
    complete_phospho = out_dir / "LemonPreprocessed_complete_phospho.txt"
    n_base = 0
    with open(complete_path) as src, open(complete_phospho, "w", newline="") as dst:
        w = csv.writer(dst, delimiter="\t", lineterminator="\n")
        header = src.readline().rstrip("\n").split("\t")
        w.writerow(header)
        for line in src:
            dst.write(line if line.endswith("\n") else line + "\n")
            n_base += 1
        for name, row in zip(kept_names, kept_rows):
            w.writerow([name, name] + [fmt(v) for v in row])
    print(f"[ok] wrote {complete_phospho} ({n_base} base rows + {len(kept_names)} phospho rows)")


if __name__ == "__main__":
    main()
