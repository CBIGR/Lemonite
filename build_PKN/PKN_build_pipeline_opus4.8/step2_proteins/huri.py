"""
HuRI (Human Reference Interactome) protein-protein interactions (local file).

HuRI is keyed by Ensembl gene ids; these are mapped to HGNC symbols via the local
Ensembl mapping table before pruning to the gene set of interest.

Re-derived from Collect_PKNdata_proteins.ipynb cells 16-20.
"""

import logging

import pandas as pd

import config

logger = logging.getLogger('pkn.huri')


def run(genes):
    """Return pruned HuRI GeneA/GeneB interactions (mapped to symbols)."""
    huri = pd.read_csv(config.HURI_LOCATION, sep='\t', header=None)
    huri.columns = ['GeneA', 'GeneB']

    mapping = pd.read_csv(config.ENSEMBL_MAPPING_FILE, sep='\t')
    symbol_to_ensembl = dict(zip(mapping['hgnc_symbol'], mapping['ensembl_gene_id']))
    ensembl_to_symbol = dict(zip(mapping['ensembl_gene_id'], mapping['hgnc_symbol']))

    genes_ensembl = [symbol_to_ensembl[g] for g in genes if g in symbol_to_ensembl]
    pruned = huri[huri['GeneA'].isin(genes_ensembl) | huri['GeneB'].isin(genes_ensembl)].copy()

    pruned['GeneA'] = pruned['GeneA'].map(ensembl_to_symbol)
    pruned['GeneB'] = pruned['GeneB'].map(ensembl_to_symbol)
    failed = pruned[['GeneA', 'GeneB']].isna().any(axis=1).sum()
    pruned = pruned.dropna(subset=['GeneA', 'GeneB'])
    pruned.to_csv(config.HURI_GENES_FILE, index=False)
    logger.info("HuRI PPI: %d interactions (%d unmapped Ensembl ids dropped)",
                len(pruned), failed)
    return pruned
