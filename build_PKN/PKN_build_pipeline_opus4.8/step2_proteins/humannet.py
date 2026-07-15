"""
HumanNet protein-protein interactions (local file, already gene symbols).

Re-derived from Collect_PKNdata_proteins.ipynb cell 14.
"""

import logging

import pandas as pd

import config

logger = logging.getLogger('pkn.humannet')


def run(genes):
    """Return pruned HumanNet GeneA/GeneB interactions for the gene set."""
    humannet = pd.read_csv(config.HUMANNET_LOCATION, sep='\t', header=None,
                           usecols=[0, 1], names=['GeneA', 'GeneB'])
    pruned = humannet[humannet['GeneA'].isin(genes) | humannet['GeneB'].isin(genes)]
    pruned.to_csv(config.HUMANNET_GENES_FILE, index=False)
    logger.info("HumanNet PPI: %d interactions", len(pruned))
    return pruned
