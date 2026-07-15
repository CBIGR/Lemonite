#!/usr/bin/env python3
"""
Evaluate phospho regulator -> module assignments against the standalone phospho-PKN.

Mirrors calculate_metabolite_gene_enrichment (nextflow/scripts/evaluate_against_PKN.py:499):
per module, count observed (phosphoprotein, module-gene) pairs that are KNOWN edges in the
phospho-PKN vs. expected under the global density; hypergeometric SF + BH-FDR.

Edge types (scored separately AND combined):
  A. kinase-substrate   (OmniPath KSN, from phospho_pkn.tsv)
  B. phospho-TF-target  (CollecTRI,   from phospho_pkn.tsv)
  C. PPI                (phosphoprotein-gene physical interaction, from Lemonite_PKN.tsv)

Matching is GENE-LEVEL: phosphosite ID GENE:RefSeq:Residue -> phosphoprotein gene = token
before the first ':'. Node2 side is the module gene symbol.

Outputs (per variant, into <run>/PhosphoPKN_Evaluation/):
  Phospho_PKN_enrichment_<edgeset>.csv   per-module hypergeometric table (A / B / C / combined)
  Phospho_PKN_supported_edges.csv        every supported (module, phosphosite, gene, edge_type)
  phospho_LemoniteKG_interactions.mvf    ModuleViewer file (heatmaps), same format as metabolite mvf
"""
import argparse, os
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "/home/borisvdm/repo/LemonIte"
PHOSPHO_PKN = os.path.join(HERE, "phospho_pkn.tsv")
MAIN_PKN = os.path.join(REPO, "nextflow", "PKN", "Lemonite_PKN.tsv")

MVF_COLOR = "ORANGE"


def phosphoprotein(psite):
    """GENE:RefSeq:Residue -> GENE (token before first ':')."""
    return str(psite).split(":")[0].strip()


def read_module_map(path):
    """module_id \\t a|b|c  ->  {module_id(int): [items]}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            mid, items = line.split("\t", 1)
            try:
                mid = int(float(mid))
            except ValueError:
                continue
            out[mid] = [x for x in items.split("|") if x]
    return out


def load_known_edges(phospho_pkn_path, main_pkn_path, need_ppi_nodes):
    """Return {edge_type: set((node1_gene, node2_gene))}. C(PPI) restricted to phosphoprotein
    nodes for tractability."""
    known = {"kinase-substrate": set(), "phospho-TF-target": set(), "PPI": set()}
    ppkn = pd.read_csv(phospho_pkn_path, sep="\t")
    # phospho-TF-target is directional: phosphoprotein must be the TF (Node1) acting on the
    # module gene (target). kinase-substrate is treated direction-agnostically because our
    # phosphoprotein may be EITHER the kinase (its substrates in the module) OR the substrate
    # (its upstream kinase in the module) — both are meaningful "known connection" evidence.
    tf = ppkn[ppkn["Type"] == "phospho-TF-target"]
    known["phospho-TF-target"] = set(zip(tf["Node1"].astype(str), tf["Node2"].astype(str)))
    ks = ppkn[ppkn["Type"] == "kinase-substrate"]
    for a, b in zip(ks["Node1"].astype(str), ks["Node2"].astype(str)):
        known["kinase-substrate"].add((a, b))
        known["kinase-substrate"].add((b, a))
    # C: PPI edges from main PKN touching any phosphoprotein node (either direction -> store both)
    print("  loading PPI subset from main PKN (phosphoprotein nodes only)...")
    for chunk in pd.read_csv(main_pkn_path, sep="\t", chunksize=500000):
        ppi = chunk[chunk["Type"] == "PPI"]
        n1 = ppi["Node1"].astype(str); n2 = ppi["Node2"].astype(str)
        m = n1.isin(need_ppi_nodes) | n2.isin(need_ppi_nodes)
        for a, b in zip(n1[m], n2[m]):
            known["PPI"].add((a, b)); known["PPI"].add((b, a))
    for t, s in known.items():
        print(f"  known {t} edges loaded: {len(s)}")
    return known


def enrich(module2genes, module2pp, known_edges, edge_types, all_genes, all_pp):
    """Hypergeometric per module over (phosphoprotein x module-gene) pairs that are known."""
    known = set()
    for t in edge_types:
        known |= known_edges[t]

    total_pairs = len(all_pp) * len(all_genes)
    total_known = sum(1 for p in all_pp for g in all_genes if (p, g) in known)
    rows = []
    for m, genes in module2genes.items():
        pps = module2pp.get(m, [])
        if not pps or not genes:
            continue
        n_possible = len(pps) * len(genes)
        n_obs = sum(1 for p in pps for g in genes if (p, g) in known)
        expected = n_possible * (total_known / total_pairs) if total_pairs else 0
        p = (hypergeom.sf(n_obs - 1, M=total_pairs, n=total_known, N=n_possible)
             if total_pairs and n_possible and total_known else 1.0)
        rows.append({"Module": m, "N_phosphoproteins": len(pps), "N_genes": len(genes),
                     "N_observed": n_obs, "N_expected": expected,
                     "Fold_enrichment": (n_obs / expected) if expected else 0.0, "P_value": p})
    df = pd.DataFrame(rows)
    if len(df):
        df["FDR"] = multipletests(df["P_value"].values, method="fdr_bh")[1]
        df = df.sort_values("P_value")
    return df, total_known, total_pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, help="<variant>/run directory")
    ap.add_argument("--out-subdir", default="PhosphoPKN_Evaluation")
    args = ap.parse_args()

    run = args.run_dir
    outdir = os.path.join(run, args.out_subdir)
    os.makedirs(outdir, exist_ok=True)

    clusters = os.path.join(run, "ModuleViewer_files", "clusters_list.txt")
    if not os.path.exists(clusters):
        clusters = os.path.join(run, "clusters_list.txt")
    phos = os.path.join(run, "ModuleViewer_files", "Phospho.selected_regs_list.txt")
    spec = os.path.join(run, "Networks", "specific_modules.txt")

    module2genes = read_module_map(clusters)
    module2psites = read_module_map(phos)
    keep = set(int(float(l)) for l in open(spec) if l.strip())
    module2genes = {m: g for m, g in module2genes.items() if m in keep}
    module2psites = {m: p for m, p in module2psites.items() if m in keep}

    # map phosphosites -> phosphoprotein genes per module (keep psite for the mvf/edge list)
    module2pp = {m: sorted({phosphoprotein(p) for p in ps}) for m, ps in module2psites.items()}
    all_genes = sorted({g for gs in module2genes.values() for g in gs})
    all_pp = sorted({p for ps in module2pp.values() for p in ps})
    print(f"Modules: {len(module2genes)} | genes: {len(all_genes)} | phosphoproteins: {len(all_pp)}")

    known = load_known_edges(PHOSPHO_PKN, MAIN_PKN, set(all_pp))

    edge_sets = {"kinase-substrate": ["kinase-substrate"],
                 "phospho-TF-target": ["phospho-TF-target"],
                 "PPI": ["PPI"],
                 "combined": ["kinase-substrate", "phospho-TF-target", "PPI"]}
    summary = {}
    for name, types in edge_sets.items():
        df, tk, tp = enrich(module2genes, module2pp, known, types, all_genes, all_pp)
        df.to_csv(os.path.join(outdir, f"Phospho_PKN_enrichment_{name}.csv"), index=False)
        nsig = int((df["FDR"] < 0.05).sum()) if len(df) else 0
        summary[name] = (len(df), nsig, tk)
        print(f"[{name}] modules tested={len(df)} FDR<0.05={nsig} total_known_pairs={tk}")

    # supported edge list + mvf (combined known set), per (module, phosphosite, gene, edge_type)
    supported = []
    mvf_lines = {}  # module -> {psite: set(genes)}
    combined_by_type = {t: known[t] for t in ("kinase-substrate", "phospho-TF-target", "PPI")}
    for m, psites in module2psites.items():
        genes = set(module2genes.get(m, []))
        for ps in psites:
            pp = phosphoprotein(ps)
            hit_genes = set()
            for t, edges in combined_by_type.items():
                g_hit = {g for g in genes if (pp, g) in edges}
                for g in g_hit:
                    supported.append({"Module": m, "Phosphosite": ps, "Phosphoprotein": pp,
                                      "Gene": g, "Edge_type": t})
                hit_genes |= g_hit
            if hit_genes:
                mvf_lines.setdefault(m, {}).setdefault(ps, set()).update(hit_genes)

    sup_df = pd.DataFrame(supported)
    sup_df.to_csv(os.path.join(outdir, "Phospho_PKN_supported_edges.csv"), index=False)

    # Cross-module phospho-TF -> target links: a phosphosite on a TF whose CollecTRI targets
    # fall in ANY network module (not necessarily the TF's own). These are the mechanistic
    # phospho-TF connections worth surfacing/visualizing even when same-module support is 0.
    tf_edges = known["phospho-TF-target"]
    gene2modules = {}
    for m, gs in module2genes.items():
        for g in gs:
            gene2modules.setdefault(g, []).append(m)
    xmod = []
    for m, psites in module2psites.items():
        for ps in psites:
            pp = phosphoprotein(ps)
            for g in all_genes:
                if (pp, g) in tf_edges:
                    for tgt_mod in gene2modules.get(g, []):
                        xmod.append({"TF_phosphosite": ps, "TF": pp, "TF_module": m,
                                     "Target_gene": g, "Target_module": tgt_mod,
                                     "Same_module": (tgt_mod == m)})
    xmod_df = pd.DataFrame(xmod)
    xmod_df.to_csv(os.path.join(outdir, "Phospho_TF_target_links.csv"), index=False)
    print(f"Cross-module phospho-TF->target links: {len(xmod_df)} "
          f"(TFs: {sorted(set(xmod_df['TF'])) if len(xmod_df) else []})")

    mvf = os.path.join(run, "ModuleViewer_files", "phospho_LemoniteKG_interactions.mvf")
    with open(mvf, "w") as h:
        h.write("::TYPE=Phospho_KG\n::TITLE:Phospho_KG\n::OBJECT=GENES\n")
        h.write(f"::COLOR={MVF_COLOR}\n")
        n = 0
        for m in sorted(mvf_lines):
            for ps, gset in mvf_lines[m].items():
                h.write(f"{int(m)}\t{'|'.join(sorted(gset))}\t{ps}\n"); n += 1
        if n == 0:
            h.write("# No phospho-gene interactions supported by phospho-PKN\n")
    print(f"\nSupported (module,phosphosite,gene) edges: {len(sup_df)}")
    print(f"mvf written: {mvf}")
    print("Summary (edgeset: modules_tested, FDR<0.05, known_pairs):")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
