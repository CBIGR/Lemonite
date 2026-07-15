#!/usr/bin/env python3
"""
Build the STANDALONE phospho prior-knowledge network (phospho_pkn.tsv).

Do NOT modify the shared nextflow/PKN/Lemonite_PKN.tsv. Output schema matches the main
PKN (Node1, Node2, Source, Type) so the evaluation can treat it identically.

Edge types produced here:
  A. kinase-substrate    Node1=kinase(enzyme) gene  -> Node2=substrate gene   (OmniPath KSN)
  B. phospho-TF-target   Node1=TF gene              -> Node2=target gene      (CollecTRI, local)
(C. PPI is read directly from Lemonite_PKN.tsv at evaluation time, not copied here.)

OmniPath KSN is fetched via decoupleR inside the pipeline .sif and cached as TSV so this
is reproducible offline once fetched. See PHOSPHO_PKN_RESOURCES.md for why OmniPath's KSN
is the key add (it aggregates ~18 curated+predicted DBs incl. PhosphoSitePlus/SIGNOR/DEPOD).
"""
import argparse, os, subprocess, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "/home/borisvdm/repo/LemonIte"
SIF = os.path.join(REPO, "nextflow", "lemontree-pipeline_v1.0.0.sif")
COLLECTRI = os.path.join(REPO, "nextflow", "PKN", "CollecTRI_network.txt")
CACHE_DIR = os.path.join(HERE, "phospho_pkn_cache")
KSN_CACHE = os.path.join(CACHE_DIR, "omnipath_ksn.tsv")
OUT = os.path.join(HERE, "phospho_pkn.tsv")

COLS = ["Node1", "Node2", "Source", "Type"]


def fetch_omnipath_ksn(force=False):
    """Fetch OmniPath kinase-substrate network via decoupleR inside the .sif; cache TSV."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    if os.path.exists(KSN_CACHE) and not force:
        print(f"[A] Using cached OmniPath KSN: {KSN_CACHE}")
        return pd.read_csv(KSN_CACHE, sep="\t")
    print("[A] Fetching OmniPath kinase-substrate network via decoupleR (inside .sif)...")
    r_code = (
        'ksn <- decoupleR::get_ksn_omnipath(); '
        f'write.table(ksn, "{KSN_CACHE}", sep="\\t", quote=FALSE, row.names=FALSE); '
        'cat("KSN_ROWS", nrow(ksn), "\\n")'
    )
    res = subprocess.run(
        ["singularity", "exec", SIF, "Rscript", "-e", r_code],
        capture_output=True, text=True,
    )
    sys.stdout.write(res.stdout)
    if res.returncode != 0 or not os.path.exists(KSN_CACHE):
        sys.stderr.write(res.stderr)
        raise RuntimeError("OmniPath KSN fetch failed")
    ksn = pd.read_csv(KSN_CACHE, sep="\t")
    print(f"[A] OmniPath KSN cached: {len(ksn)} edges -> {KSN_CACHE}")
    return ksn


def build_kinase_substrate(ksn):
    """OmniPath KSN -> Node1=kinase, Node2=substrate GENE, Type=kinase-substrate.

    decoupleR's KSN `target` is SITE-level (e.g. 'AP2M1_T154'): GENE_<residue>. Strip the
    residue suffix to get the substrate gene symbol (keep residue as metadata so the edge can
    be site-resolved later). Node1 (kinase 'source') is already a bare gene symbol.
    """
    df = ksn.rename(columns={"source": "Node1", "target": "_site"})[["Node1", "_site"]].copy()
    df["_site"] = df["_site"].astype(str)
    # split GENE_RESIDUE on the LAST underscore (gene symbols can contain underscores rarely,
    # but the residue token is <letter><digits>, so rsplit on '_' once is safe).
    df["Node2"] = df["_site"].str.rsplit("_", n=1).str[0]
    df["Residue"] = df["_site"].str.rsplit("_", n=1).str[1]
    df["Source"] = "OmniPath_KSN"
    df["Type"] = "kinase-substrate"
    df = df.dropna(subset=["Node1", "Node2"])
    df = df[df["Node2"].str.len() > 0]
    df = df.drop_duplicates(subset=["Node1", "Node2"])  # gene-level dedupe
    print(f"[A] kinase-substrate edges: {len(df)} "
          f"({df['Node1'].nunique()} kinases -> {df['Node2'].nunique()} substrate genes)")
    return df[COLS + ["Residue"]]


def build_phospho_tf_target():
    """CollecTRI -> Node1=TF, Node2=target, Type=phospho-TF-target."""
    ct = pd.read_csv(COLLECTRI, sep="\t")
    df = ct.rename(columns={"source": "Node1", "target": "Node2"})[["Node1", "Node2"]].copy()
    df["Source"] = "CollecTRI"
    df["Type"] = "phospho-TF-target"
    df["Residue"] = ""
    df = df.dropna(subset=["Node1", "Node2"]).drop_duplicates(subset=["Node1", "Node2"])
    print(f"[B] phospho-TF-target edges: {len(df)} "
          f"({df['Node1'].nunique()} TFs -> {df['Node2'].nunique()} targets)")
    return df[COLS + ["Residue"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force-fetch", action="store_true", help="re-fetch OmniPath even if cached")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    frames = []
    ksn = fetch_omnipath_ksn(force=args.force_fetch)
    frames.append(build_kinase_substrate(ksn))
    frames.append(build_phospho_tf_target())

    pkn = pd.concat(frames, ignore_index=True)
    pkn = pkn[COLS + ["Residue"]]
    # dedupe within (Node1,Node2,Type) but keep both a KSN and TF edge for the same pair
    pkn = pkn.drop_duplicates(subset=["Node1", "Node2", "Type"])
    pkn.to_csv(args.out, sep="\t", index=False)
    print(f"\n[phospho_pkn] wrote {len(pkn)} edges -> {args.out}")
    print(pkn["Type"].value_counts().to_string())


if __name__ == "__main__":
    main()
