"""
Step-3 analysis: summary statistics for the final combined PKN.

Logs node/edge counts and per-source breakdowns. The reproducible summary figures
(per-source counts, database-overlap UpSet plots, composition) are drawn separately
in visualization.py; the notebook's MetalinksDB/MEBOCOST comparison plots require
extra database files that are not needed to build or validate the network and are
omitted.

Re-derived from Build_final_PKN.ipynb cells 25-30.
"""

import logging

import pandas as pd

import config

logger = logging.getLogger('pkn.analysis')


def analyze():
    """Compute and log summary statistics for LemonIte_PKN.tsv."""
    pkn = pd.read_csv(config.FINAL_PKN_FILE, sep='\t')

    mg = pkn[pkn['Type'] == 'metabolite-gene']
    ppi = pkn[pkn['Type'] == 'PPI']

    metabolites = set(mg['Node1'].unique())
    genes = set(mg['Node2'].unique()) | set(ppi['Node1'].unique()) | set(ppi['Node2'].unique())

    stats = {
        'unique_metabolites': len(metabolites),
        'unique_genes': len(genes),
        'total_nodes': len(metabolites) + len(genes),
        'metabolite_gene_edges': len(mg),
        'ppi_edges': len(ppi),
        'total_edges': len(pkn),
    }
    logger.info("=== Final PKN statistics ===")
    for k, v in stats.items():
        logger.info("  %-24s %d", k, v)

    logger.info("Metabolite-gene edges by source:")
    for src, n in mg['Source'].value_counts().items():
        logger.info("  %-18s %d", src, n)
    logger.info("PPI edges by source:")
    for src, n in ppi['Source'].value_counts().items():
        logger.info("  %-18s %d", src, n)
    return stats
