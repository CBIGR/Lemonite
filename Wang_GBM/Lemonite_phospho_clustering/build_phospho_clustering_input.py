#!/usr/bin/env python3
"""
Build the phosphoproteomics CLUSTERING-INPUT matrix for LemonTree (this variant clusters ON
phospho, unlike Wang_GBM/Lemonite_phospho/ where phospho was a regulator).

Reads mmc3 `phosphoproteome_normalized` (already log2 / median-polished / ComBat), then:
  1. subset to the requested sample set,
  2. drop sites with < min_valid_frac non-NA samples,
  3. NA -> 0,
  4. z-score each phosphosite (row-wise; matches the pipeline's PTM scaling t(scale(t(x)))),
  5. optionally keep the top-N most variable sites (variance on NA->0 abundances, pre-scaling).

Writes LemonTree's primary clustering input `LemonPreprocessed_expression.txt`
(`symbol <TAB> id <TAB> <samples>`), plus a `phospho_sites.txt` list and a `samples.txt`.
Two id columns are the phosphosite id (LemonTree keys rows by the first column).

Usage:
  build_phospho_clustering_input.py --out-dir <dir>                 # all sites
  build_phospho_clustering_input.py --out-dir <dir> --top-var 5000  # top-5000 HVG
"""
import argparse, csv, os
import numpy as np
import openpyxl

MMC3 = "/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/data/1-s2.0-S1535610821000507-mmc3.xlsx"
SHEET = "phosphoproteome_normalized"


def sanitize(name):
    return str(name).replace(" ", "_").replace("\t", "_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--mmc3", default=MMC3)
    ap.add_argument("--min-valid-frac", type=float, default=0.5)
    ap.add_argument("--top-var", type=int, default=0, help="keep top-N most variable sites; 0 = all")
    ap.add_argument("--sample-prefix", default="C3L,C3N",
                    help="comma-separated sample-id prefixes to keep (default CPTAC discovery cohort)")
    args = ap.parse_args()
    prefixes = tuple(p.strip() for p in args.sample_prefix.split(","))
    os.makedirs(args.out_dir, exist_ok=True)

    wb = openpyxl.load_workbook(args.mmc3, read_only=True, data_only=True)
    ws = wb[SHEET]
    it = ws.iter_rows(values_only=True)
    hdr = next(it)
    # cols: site_id, symbol, phosphosites, peptide, <samples...>
    meta_n = 4
    all_samples = [str(c) for c in hdr[meta_n:]]
    # The phospho sheet mixes the CPTAC discovery cohort (C3L-/C3N-) with 10 pediatric CBTTC
    # validation samples (PT-*) that have no matching proteome/metabolome/lipidome. Keep only
    # samples matching --sample-pattern (default CPTAC) so all omics layers are on one cohort.
    keep_idx = [i for i, s in enumerate(all_samples) if s.startswith(prefixes)]
    samples = [all_samples[i] for i in keep_idx]
    print(f"[info] samples: {len(all_samples)} in sheet -> kept {len(samples)} "
          f"matching prefixes {prefixes}")

    names, rows, variances = [], [], []
    n_all_na = n_lowvalid = 0
    for r in it:
        site_id = r[0]
        if site_id is None:
            continue
        allvals = r[meta_n:]
        vals = [allvals[i] for i in keep_idx]
        arr = np.array([np.nan if (v is None or v == "NA" or v == "") else float(v) for v in vals],
                       dtype=float)
        valid = ~np.isnan(arr)
        if valid.sum() == 0:
            n_all_na += 1
            continue
        if valid.mean() < args.min_valid_frac:
            n_lowvalid += 1
            continue
        row = np.where(valid, arr, 0.0)          # NA -> 0
        variances.append(float(np.var(row)))
        names.append(sanitize(site_id))
        rows.append(row)

    print(f"[info] samples: {len(samples)} | sites read (pass valid-frac): {len(names)} "
          f"| dropped all-NA {n_all_na}, low-valid {n_lowvalid}")

    order = list(range(len(names)))
    if args.top_var and args.top_var < len(order):
        order = sorted(order, key=lambda i: variances[i], reverse=True)[:args.top_var]
        order.sort()
        print(f"[info] --top-var {args.top_var}: kept {len(order)} most-variable sites")

    # de-dup by name, z-score per row
    seen = set(); fnames = []; frows = []
    for i in order:
        nm = names[i]
        if nm in seen:
            continue
        seen.add(nm)
        row = rows[i]
        sd = row.std()
        z = (row - row.mean()) / sd if sd > 0 else row * 0.0
        fnames.append(nm); frows.append(z)

    fmt = lambda x: f"{x:.6g}"
    expr = os.path.join(args.out_dir, "LemonPreprocessed_expression.txt")
    with open(expr, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["symbol", "id"] + samples)
        for nm, z in zip(fnames, frows):
            w.writerow([nm, nm] + [fmt(v) for v in z])
    with open(os.path.join(args.out_dir, "phospho_sites.txt"), "w") as fh:
        fh.write("\n".join(fnames) + "\n")
    with open(os.path.join(args.out_dir, "samples.txt"), "w") as fh:
        fh.write("\n".join(samples) + "\n")
    print(f"[ok] wrote {expr}  ({len(fnames)} sites x {len(samples)} samples)")


if __name__ == "__main__":
    main()
