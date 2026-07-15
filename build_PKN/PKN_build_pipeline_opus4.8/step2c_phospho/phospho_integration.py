"""
Step-2c orchestration: build the STANDALONE phospho-regulator network from OmniPath
enzyme-substrate (kinase-substrate) + CollecTRI (phospho-TF-target).

Output: phospho_pkn.tsv with columns
    Node1, Node2, Source, Type, Residue, URL
where Type is 'kinase-substrate' or 'phospho-TF-target'. This layer is deliberately kept
separate from the main LemonIte_PKN.tsv (step 3 does NOT merge it) and is consumed by the
phospho-regulator evaluation.
"""

import logging

import pandas as pd

import config
from step2c_phospho import omnipath_ksn, collectri

logger = logging.getLogger('pkn.phospho')

_COLS = ['Node1', 'Node2', 'Source', 'Type', 'Residue', 'URL']


def build_phospho_network(resume=True):
    """Collect OmniPath KSN + CollecTRI edges; write phospho_pkn.tsv."""
    ksn_df = omnipath_ksn.run(resume=resume)
    tf_df = collectri.run(resume=resume)

    frames = [df for df in (ksn_df, tf_df) if df is not None and not df.empty]
    if not frames:
        phospho = pd.DataFrame(columns=_COLS)
    else:
        phospho = pd.concat(frames, ignore_index=True)[_COLS]
        # keep one row per (Node1, Node2, Type): a KSN and a TF edge for the same pair coexist
        phospho = phospho.drop_duplicates(subset=['Node1', 'Node2', 'Type'])

    phospho.to_csv(config.PHOSPHO_PKN_FILE, sep='\t', index=False)
    if not phospho.empty:
        counts = phospho['Type'].value_counts()
        logger.info("phospho network: %d edges (%d kinase-substrate, %d phospho-TF-target) -> %s",
                    len(phospho), int(counts.get('kinase-substrate', 0)),
                    int(counts.get('phospho-TF-target', 0)), config.PHOSPHO_PKN_FILE)
    else:
        logger.warning("phospho network: 0 edges -> %s", config.PHOSPHO_PKN_FILE)
    return phospho
