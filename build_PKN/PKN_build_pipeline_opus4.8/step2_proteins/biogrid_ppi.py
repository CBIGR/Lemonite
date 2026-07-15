"""
BioGRID protein-protein interactions (local file).

Filters BioGRID to human physical interactions involving at least one gene from the
metabolite-gene PKN (+ first neighbours).

Re-derived from Collect_PKNdata_proteins.ipynb cells 10-12.
"""

import logging

import pandas as pd

import config

logger = logging.getLogger('pkn.biogrid_ppi')


def run(genes):
    """Return pruned BioGRID GeneA/GeneB interactions for the gene set."""
    biogrid = pd.read_csv(config.BIOGRID_PPI_LOCATION, sep='\t', low_memory=False)
    biogrid = biogrid[(biogrid['Organism Name Interactor A'] == 'Homo sapiens')
                      | (biogrid['Organism Name Interactor B'] == 'Homo sapiens')]
    biogrid = biogrid[['Official Symbol Interactor A', 'Official Symbol Interactor B']]
    biogrid.columns = ['GeneA', 'GeneB']
    pruned = biogrid[biogrid['GeneA'].isin(genes) | biogrid['GeneB'].isin(genes)]
    pruned.to_csv(config.BIOGRID_GENES_FILE, index=False)
    logger.info("BioGRID PPI: %d interactions", len(pruned))
    return pruned
