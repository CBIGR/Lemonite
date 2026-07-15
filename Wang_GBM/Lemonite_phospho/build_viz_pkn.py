#!/usr/bin/env python3
"""
Build the viz-only combined PKN = main Lemonite_PKN.tsv + standalone phospho_pkn.tsv.

Used ONLY by the subnetwork figures (create_subnetwork_graphs.py via run_downstream.sh
PKN_VIZ) so phospho regulator->target edges render alongside metabolite/PPI edges. The main
PKN is NOT modified. The output (Lemonite_PKN_plus_phospho.tsv, ~90 MB) is a regenerable
derived artifact and is gitignored — rebuild it with this script.
"""
import os
import shutil
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MAIN = "/home/borisvdm/repo/LemonIte/nextflow/PKN/Lemonite_PKN.tsv"
PHOS = os.path.join(HERE, "phospho_pkn.tsv")
OUT = os.path.join(HERE, "Lemonite_PKN_plus_phospho.tsv")


def main():
    pp = pd.read_csv(PHOS, sep="\t")[["Node1", "Node2", "Source", "Type"]]
    shutil.copyfile(MAIN, OUT)                       # main PKN header + edges
    pp.to_csv(OUT, sep="\t", index=False, header=False, mode="a")  # append phospho edges
    n = sum(1 for _ in open(OUT))
    print(f"viz PKN: {n} rows (incl header), {len(pp)} phospho edges appended -> {OUT}")


if __name__ == "__main__":
    main()
