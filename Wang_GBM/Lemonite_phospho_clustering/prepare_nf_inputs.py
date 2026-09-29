#!/usr/bin/env python3
"""
Build the shared `data/` directory used to feed the phospho-as-primary matrix through the
REAL nextflow pipeline (main.nf, preprocessing_type=proteomics), instead of hand-rolled
LemonTree calls. One data/ dir is shared by both variants (top5000_hvg / all_sites); they
differ only by --top_n_genes at `nextflow run` time.

Produces, under --out-dir:
  phospho_expression_full.tsv   full z-scored phospho primary matrix (all valid sites, no
                                 top-N truncation -- main.nf's --top_n_genes does that itself)
  metadata.tsv                  Sample_ID + gender + vital_status for the 99 CPTAC samples
  proteins_raw.tsv              Proteome_normalized.csv, symbol + sample columns only
                                 (refseq_prot_id/hgnc_id dropped -- the generic regulator
                                 loader only strips ONE id column)
  lipids_raw.tsv                lipidome_pos.csv + lipidome_neg.csv concatenated
  metabolome.csv                copied as-is (native --metabolomics_file hook handles it)
"""
import argparse, csv, os, shutil

DATA = "/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/data"


def build_metadata(out_dir, samples):
    with open(os.path.join(DATA, "clinical_metadata.csv"), newline="") as f:
        r = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        hdr = next(r)
        idx = {h: i for i, h in enumerate(hdr)}
        rows = {}
        for row in r:
            if len(row) <= max(idx.get("case_submitter_id", 0), idx.get("gender", 0), idx.get("vital_status", 0)):
                continue
            sid = row[idx["case_submitter_id"]]
            if sid in rows:
                continue
            rows[sid] = (row[idx.get("gender", 0)] or "unknown", row[idx.get("vital_status", 0)] or "unknown")
    out = os.path.join(out_dir, "metadata.tsv")
    n_missing = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Sample_ID", "gender", "vital_status"])
        for s in samples:
            gender, vs = rows.get(s, ("unknown", "unknown"))
            if s not in rows:
                n_missing += 1
            w.writerow([s, gender, vs])
    print(f"[ok] metadata.tsv: {len(samples)} samples ({n_missing} without a clinical_metadata match)")


def build_proteins_raw(out_dir):
    src = os.path.join(DATA, "Proteome_normalized.csv")
    out = os.path.join(out_dir, "proteins_raw.tsv")
    with open(src) as f, open(out, "w", newline="") as fh:
        r = csv.reader(f, delimiter="\t")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        hdr = next(r)
        w.writerow([hdr[0]] + hdr[3:])   # symbol + samples (drop refseq_prot_id, hgnc_id)
        n = 0
        for row in r:
            if len(row) < 4 or not row[0]:
                continue
            w.writerow([row[0]] + row[3:])
            n += 1
    print(f"[ok] proteins_raw.tsv: {n} proteins -> {out}")


def build_lipids_raw(out_dir):
    out = os.path.join(out_dir, "lipids_raw.tsv")
    n = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        first = True
        for fname in ("lipidome_pos.csv", "lipidome_neg.csv"):
            with open(os.path.join(DATA, fname)) as f:
                r = csv.reader(f, delimiter="\t")
                hdr = next(r)
                if first:
                    w.writerow(hdr)
                    first = False
                for row in r:
                    w.writerow(row)
                    n += 1
    print(f"[ok] lipids_raw.tsv: {n} lipids (pos+neg) -> {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--phospho-dir", required=True,
                    help="dir containing LemonPreprocessed_expression.txt + samples.txt from build_phospho_clustering_input.py --top-var 0")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    with open(os.path.join(args.phospho_dir, "samples.txt")) as fh:
        samples = [l.strip() for l in fh if l.strip()]

    build_metadata(args.out_dir, samples)
    build_proteins_raw(args.out_dir)
    build_lipids_raw(args.out_dir)
    shutil.copy(os.path.join(DATA, "metabolome.csv"), os.path.join(args.out_dir, "metabolome.csv"))
    print("[ok] metabolome.csv copied")

    shutil.copy(os.path.join(args.phospho_dir, "LemonPreprocessed_expression.txt"),
                os.path.join(args.out_dir, "phospho_expression_full.tsv"))
    print("[ok] phospho_expression_full.tsv copied")


if __name__ == "__main__":
    main()
