"""
Step-3 combiner: merge the metabolite-gene, PPI and enzyme-histone-mark networks
into the final PKN.

Standardizes each to Node1/Node2/Source, tags its edge Type (metabolite-gene / PPI
/ histone-modification), concatenates, and dedupes on (Node1, Node2). Writes
LemonIte_PKN.tsv.

Re-derived from Build_final_PKN.ipynb cells 2-6; the histone-modification block is
new (step 2b).
"""

import logging
import os

import pandas as pd

import config

logger = logging.getLogger('pkn.combiner')


def combine_networks():
    """Combine metabolite-gene + PPI + hPTM networks; dedupe on node pair; write PKN."""
    mg = pd.read_csv(config.METABOLITE_GENE_PKN, sep='\t')
    ppi = pd.read_csv(config.PPI_NETWORK_FILE, sep='\t')

    # Metabolite-gene -> Node1/Node2/Source (drop the url column for the combined view)
    mg = mg[['Metabolite', 'Gene', 'Source']].copy()
    mg.columns = ['Node1', 'Node2', 'Source']
    mg['Type'] = 'metabolite-gene'

    # PPI -> Node1/Node2/Source (drop combined_score)
    ppi = ppi[['GeneA', 'GeneB', 'Source']].copy()
    ppi.columns = ['Node1', 'Node2', 'Source']
    ppi['Type'] = 'PPI'

    frames = [ppi, mg]

    # Enzyme -> histone-mark (step 2b). Optional: only present if step 2b ran.
    if os.path.exists(config.HPTM_NETWORK_FILE):
        hptm = pd.read_csv(config.HPTM_NETWORK_FILE, sep='\t')
        if not hptm.empty:
            hptm = hptm[['Enzyme', 'Mark', 'Source']].copy()
            hptm.columns = ['Node1', 'Node2', 'Source']
            hptm['Type'] = 'histone-modification'
            frames.append(hptm)

    combined = pd.concat(frames, axis=0, ignore_index=True)
    pruned = combined.drop_duplicates(subset=['Node1', 'Node2'])
    pruned.to_csv(config.FINAL_PKN_FILE, sep='\t', index=False)

    n_mg = (pruned['Type'] == 'metabolite-gene').sum()
    n_ppi = (pruned['Type'] == 'PPI').sum()
    n_hptm = (pruned['Type'] == 'histone-modification').sum()
    logger.info("Combined PKN: %d edges (%d metabolite-gene, %d PPI, %d histone-modification) -> %s",
                len(pruned), n_mg, n_ppi, n_hptm, config.FINAL_PKN_FILE)
    return pruned
