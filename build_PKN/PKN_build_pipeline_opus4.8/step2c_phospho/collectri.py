"""
Load the CollecTRI TF -> target network (shipped with the repo) as phospho-TF-target edges.

CollecTRI (source, target, mor) is a curated collection of TF -> transcriptional-target
regulons. For the phospho layer we treat a phospho-event on a TF as potentially modulating
that TF's known targets, so each TF -> target becomes a phospho-TF-target edge. Local file,
no network access needed.
"""

import logging
import os

import pandas as pd

import config

logger = logging.getLogger('pkn.phospho.collectri')

COLLECTRI_URL = 'https://github.com/saezlab/CollecTRI'


def run(resume=True):
    """Return phospho-TF-target edges: Node1, Node2, Source, Type, Residue, URL."""
    path = config.COLLECTRI_FILE
    if not os.path.exists(path):
        logger.warning("CollecTRI file not found: %s -> skipping phospho-TF-target edges", path)
        return pd.DataFrame(columns=['Node1', 'Node2', 'Source', 'Type', 'Residue', 'URL'])

    ct = pd.read_csv(path, sep='\t')
    if 'source' not in ct.columns or 'target' not in ct.columns:
        logger.warning("CollecTRI file missing source/target columns (has %s)", list(ct.columns))
        return pd.DataFrame(columns=['Node1', 'Node2', 'Source', 'Type', 'Residue', 'URL'])

    df = ct.rename(columns={'source': 'Node1', 'target': 'Node2'})[['Node1', 'Node2']].copy()
    df['Node1'] = df['Node1'].astype(str)
    df['Node2'] = df['Node2'].astype(str)
    df['Source'] = 'CollecTRI'
    df['Type'] = 'phospho-TF-target'
    df['Residue'] = ''
    df['URL'] = COLLECTRI_URL
    df = df.dropna(subset=['Node1', 'Node2']).drop_duplicates(subset=['Node1', 'Node2'])
    logger.info("Phospho-TF-target edges: %d (%d TFs -> %d targets) from %s",
                len(df), df['Node1'].nunique(), df['Node2'].nunique(), path)
    return df[['Node1', 'Node2', 'Source', 'Type', 'Residue', 'URL']].reset_index(drop=True)
