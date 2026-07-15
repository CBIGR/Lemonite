"""
Step-2 orchestration: build the combined PPI network from STRING, BioGRID, HuRI
and HumanNet for the genes in the metabolite-gene PKN.

Output: PPI_network.tsv with columns GeneA, GeneB, Source, combined_score.

Re-derived from Collect_PKNdata_proteins.ipynb cell 22.
"""

import logging

import pandas as pd

import config
from step2_proteins import string_api, biogrid_ppi, huri, humannet

logger = logging.getLogger('pkn.ppi')


def build_ppi_network(resume=True):
    """Collect, label and concatenate the four PPI sources; write PPI_network.tsv."""
    mg = pd.read_csv(config.METABOLITE_GENE_PKN, sep='\t')
    genes = mg['Gene'].unique()
    logger.info("PPI: seeding from %d genes in metabolite-gene PKN", len(genes))

    string_df = string_api.run(genes, resume=resume)
    biogrid_df = biogrid_ppi.run(genes)
    huri_df = huri.run(genes)
    humannet_df = humannet.run(genes)

    string_df['Source'] = 'STRING'
    biogrid_df['Source'] = 'BioGRID_genes'
    huri_df['Source'] = 'HuRI'
    humannet_df['Source'] = 'HumanNet'

    ppi = pd.concat([string_df, biogrid_df, huri_df, humannet_df], ignore_index=True)
    ppi['GeneA'] = ppi['GeneA'].astype(str)
    ppi['GeneB'] = ppi['GeneB'].astype(str)
    ppi = ppi[~ppi['GeneA'].str.contains('nan') & ~ppi['GeneB'].str.contains('nan')]

    # Keep a combined_score column even where a source lacks one
    if 'combined_score' not in ppi.columns:
        ppi['combined_score'] = pd.NA

    ppi.to_csv(config.PPI_NETWORK_FILE, sep='\t', index=False)
    n_genes = len(set(ppi['GeneA']) | set(ppi['GeneB']))
    logger.info("PPI network: %d interactions, %d unique genes -> %s",
                len(ppi), n_genes, config.PPI_NETWORK_FILE)
    return ppi
