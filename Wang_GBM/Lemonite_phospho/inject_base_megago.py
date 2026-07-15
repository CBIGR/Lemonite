#!/usr/bin/env python3
"""
Inject the base GBM run's module-level MegaGO / rrvgo functional-clustering labels
into a phospho-variant Module_Overview_Complete.csv.

The gene co-expression modules are identical between the base run and the phospho
variants (regulator assignment does not change modules), and MegaGO clustering is a
module-level operation. So rather than recompute MegaGO, we reuse the base run's
per-module MegaGO_cluster / MegaGO_label / MegaGO_representative_GO_ID by joining on
the Module id.

Both overview files are tab-separated.
"""
import argparse
import pandas as pd

# base-run (lowercase) MegaGO column  ->  canonical column the variant overview uses
COL_MAP = {
    "MegaGO_cluster": "MegaGO_Cluster",
    "MegaGO_label": "MegaGO_Label",
    "MegaGO_representative_GO_ID": "MegaGO_Representative_GO_ID",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overview", required=True, help="phospho-variant Module_Overview.csv (updated in place)")
    ap.add_argument("--base-overview", required=True, help="base run Module_Overview_Complete.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.overview, sep="\t")
    base = pd.read_csv(args.base_overview, sep="\t")

    if "Module" not in df.columns or "Module" not in base.columns:
        raise SystemExit(f"[error] 'Module' column missing (overview cols={list(df.columns)[:5]}...)")

    src_cols = [c for c in COL_MAP if c in base.columns]
    if not src_cols:
        raise SystemExit(f"[error] no MegaGO columns in base overview; has {list(base.columns)}")

    base_m = base[["Module"] + src_cols].rename(columns=COL_MAP).copy()
    tgt_cols = [COL_MAP[c] for c in src_cols]
    base_m["Module"] = base_m["Module"].astype(str)
    df["Module"] = df["Module"].astype(str)

    # drop the placeholder MegaGO columns from the --n_clusters 1 run before merging
    df = df.drop(columns=[c for c in tgt_cols if c in df.columns], errors="ignore")
    merged = df.merge(base_m, on="Module", how="left")

    lead = COL_MAP["MegaGO_cluster"]
    n_labeled = merged[lead].notna().sum() if lead in merged else 0
    have = tgt_cols
    merged.to_csv(args.overview, sep="\t", index=False)
    # keep an xlsx alongside if the original had one
    try:
        merged.to_excel(args.overview.replace(".csv", ".xlsx"), index=False)
    except Exception as e:
        print(f"[warn] could not write xlsx: {e}")

    print(f"[ok] injected base MegaGO labels for {n_labeled}/{len(merged)} modules "
          f"into {args.overview} (columns: {have})")


if __name__ == "__main__":
    main()
