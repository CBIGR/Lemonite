#!/usr/bin/env python3
"""
Build the regulator abundance layers for the phospho-clustering variant and assemble
LemonPreprocessed_complete.txt (the matrix LemonTree reads regulator profiles from).

Regulator layers (each z-scored per feature, aligned to the phospho clustering samples):
  - Proteins     (Proteome_normalized.csv)           -> Proteins  / proteins.txt
  - Metabolites  (metabolome.csv)                     -> Metabolites / metabolites.txt
  - Lipids       (lipidome_pos.csv + lipidome_neg.csv)-> Lipids / lipids.txt
  - TFs          (lovering_TF_list.txt, subset of the protein layer that are TFs)
  - KinaseActivity (from kinase_and_tf_activity.R)    -> already written by that R script

complete.txt = phospho clustering rows + all regulator rows, so `-task regulators` can look
up each regulator's profile. Columns = the phospho clustering samples (missing values padded).
"""
import argparse, csv, os
import numpy as np

DATA = "/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/data"


def sanitize(name):
    """LemonTree splits rows on whitespace, so feature names must contain no spaces/tabs.
    Mirror the phospho builder: spaces and tabs -> underscore."""
    return str(name).replace(" ", "_").replace("\t", "_")


def read_matrix(path, id_cols, sep="\t"):
    """Return (samples, {feature: {sample: value}}) from a wide table; id_cols = # leading id cols."""
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split(sep)
        samples = hdr[id_cols:]
        data = {}
        for line in f:
            p = line.rstrip("\n").split(sep)
            if len(p) <= id_cols:
                continue
            feat = sanitize(p[0])
            vals = {}
            for s, v in zip(samples, p[id_cols:]):
                if v not in ("", "NA", "NaN"):
                    try:
                        vals[s] = float(v)
                    except ValueError:
                        pass
            data[feat] = vals
    return samples, data


def zscore_rows(data, samples):
    """Return {feature: [z per sample]} z-scored per feature over available samples; NA->0."""
    out = {}
    for feat, vals in data.items():
        arr = np.array([vals.get(s, np.nan) for s in samples], dtype=float)
        valid = ~np.isnan(arr)
        if valid.sum() < 3:
            continue
        row = np.where(valid, arr, np.nanmean(arr[valid]))
        sd = row.std()
        out[feat] = ((row - row.mean()) / sd) if sd > 0 else row * 0.0
    return out


def write_reg(out_dir, prefix, zmat, samples):
    fmt = lambda x: f"{x:.6g}"
    path = os.path.join(out_dir, f"LemonPreprocessed_{prefix.lower()}.txt")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene_symbol", "Protein_id"] + samples)
        for feat, z in zmat.items():
            w.writerow([feat, feat] + [fmt(v) for v in z])
    with open(os.path.join(out_dir, f"{prefix.lower()}.txt"), "w") as fh:
        fh.write("\n".join(zmat.keys()) + "\n")
    print(f"[ok] {prefix}: {len(zmat)} features -> {path}")
    return zmat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant-dir", required=True, help="<variant> results dir (has Preprocessing/, Activity/)")
    args = ap.parse_args()
    prep = os.path.join(args.variant_dir, "Preprocessing")
    act = os.path.join(args.variant_dir, "Activity")

    # clustering samples (the phospho matrix columns)
    with open(os.path.join(prep, "samples.txt")) as fh:
        samples = [l.strip() for l in fh if l.strip()]
    print(f"[info] clustering samples: {len(samples)}")

    layers = {}
    # Proteins
    _, prot = read_matrix(os.path.join(DATA, "Proteome_normalized.csv"), id_cols=3)
    layers["Proteins"] = write_reg(prep, "Proteins", zscore_rows(prot, samples), samples)
    # Metabolites
    _, met = read_matrix(os.path.join(DATA, "metabolome.csv"), id_cols=1)
    layers["Metabolites"] = write_reg(prep, "Metabolites", zscore_rows(met, samples), samples)
    # Lipids (pos + neg)
    _, lp = read_matrix(os.path.join(DATA, "lipidome_pos.csv"), id_cols=1)
    _, ln = read_matrix(os.path.join(DATA, "lipidome_neg.csv"), id_cols=1)
    lp.update(ln)
    layers["Lipids"] = write_reg(prep, "Lipids", zscore_rows(lp, samples), samples)
    # TFs = protein-layer features that are in the lovering TF list
    tfset = set()
    with open(os.path.join(DATA, "lovering_TF_list.txt")) as fh:
        next(fh)
        for line in fh:
            tfset.add(line.split("\t")[0].strip())
    tf_z = {f: z for f, z in layers["Proteins"].items() if f in tfset}
    write_reg(prep, "TFs", tf_z, samples)
    # KinaseActivity: written by kinase_and_tf_activity.R as LemonPreprocessed_kinaseactivity.txt
    kin_reg = os.path.join(act, "LemonPreprocessed_kinaseactivity.txt")
    kin_z = {}
    if os.path.exists(kin_reg):
        _, kin = read_matrix(kin_reg, id_cols=2)
        kin_z = zscore_rows(kin, samples)
        # copy list + matrix into Preprocessing for consistency
        write_reg(prep, "KinaseActivity", kin_z, samples)
    else:
        print("[warn] kinase activity regulator file not found; run kinase_and_tf_activity.R first")

    # ---- assemble LemonPreprocessed_complete.txt: phospho rows + all regulator rows ----------
    # LemonTree indexes rows by NAME; duplicate names corrupt its array indexing
    # (ArrayIndexOutOfBoundsException). So the complete matrix must have UNIQUE row names.
    #  - TFs are a value-identical subset of Proteins -> drop duplicates (keep the protein row).
    #  - KinaseActivity is DIFFERENT data (activity, not abundance) but shares gene names with
    #    Proteins -> suffix its names '_kact' so it is a distinct feature. Rewrite the
    #    kinaseactivity.txt reg-list to the suffixed names so `-task regulators` matches.
    if kin_z:
        kin_z = {f"{k}_kact": v for k, v in kin_z.items()}
        with open(os.path.join(prep, "kinaseactivity.txt"), "w") as fh:
            fh.write("\n".join(kin_z.keys()) + "\n")

    complete = os.path.join(prep, "LemonPreprocessed_complete.txt")
    fmt = lambda x: f"{x:.6g}"
    n = 0
    seen = set()
    with open(complete, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["symbol", "id"] + samples)
        # phospho clustering rows (already unique site ids)
        with open(os.path.join(prep, "LemonPreprocessed_expression.txt")) as src:
            next(src)
            for line in src:
                nm = line.split("\t", 1)[0]
                if nm in seen:
                    continue
                seen.add(nm)
                fh.write(line if line.endswith("\n") else line + "\n"); n += 1
        # regulator rows (dedupe by name across layers)
        for prefix, zmat in [("Proteins", layers["Proteins"]), ("Metabolites", layers["Metabolites"]),
                             ("Lipids", layers["Lipids"]), ("TFs", tf_z), ("KinaseActivity", kin_z)]:
            for feat, z in zmat.items():
                if feat in seen:
                    continue
                seen.add(feat)
                w.writerow([feat, feat] + [fmt(v) for v in z]); n += 1
    print(f"[ok] complete matrix: {n} unique rows -> {complete}")


if __name__ == "__main__":
    main()
