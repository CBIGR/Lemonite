"""
MetalinksDB ligand-receptor interactions (local file).

Keyed directly by HMDB id. The per-interaction provenance URL comes from the
``source`` column of metalinks.csv mapped through METALINKS_SOURCE_URLS.

Re-derived from Collect_PKNdata_metabolites.ipynb cells 59-60.
"""

import logging

import numpy as np
import pandas as pd

import config

logger = logging.getLogger('pkn.metalinks')


def run(df):
    """Build MetalinksDB metabolite-gene interactions; write links + processed files."""
    metalinks = pd.read_csv(config.METALINKS_PATH)

    hmdb2gene = (metalinks.groupby('hmdb')['gene_symbol']
                 .apply(lambda x: list(set(x))).to_dict())

    df_proc = df[['HMDB']].copy()
    df_proc['MetalinksDB'] = df_proc['HMDB'].apply(
        lambda x: '|'.join(hmdb2gene[x]) if x in hmdb2gene else np.nan)
    df_proc.to_csv(config.DB_OUTPUT_FILES['MetalinksDB'], index=False, sep='\t')

    # Provenance URL: MetalinksDB itself has no per-interaction web page, and its
    # source DBs (CellPhoneDB/scConnect/...) only expose homepages. The metabolite's
    # HMDB page DOES carry a "Protein Associations" section listing the proteins it
    # interacts with — concrete evidence the user can inspect for the metabolite
    # side of the edge. Link there per HMDB id, and keep the originating source DB
    # name + its homepage as a secondary provenance column.
    rows = []
    for _, m in metalinks.iterrows():
        src = m.get('source')
        hmdb = m['hmdb']
        url = config.URL_TEMPLATES['HMDB'].format(hmdb_id=hmdb) if pd.notna(hmdb) \
            else config.METALINKS_DEFAULT_URL
        src_url = (config.METALINKS_SOURCE_URLS.get(src, config.METALINKS_DEFAULT_URL)
                   if pd.notna(src) else config.METALINKS_DEFAULT_URL)
        rows.append({'HMDB': hmdb, 'Gene': m['gene_symbol'], 'Source': src,
                     'URL': url, 'Source_DB_URL': src_url})
    pd.DataFrame(rows).to_csv(config.URL_FILES['MetalinksDB'], index=False, sep='\t')
    logger.info("MetalinksDB: %d interaction links", len(rows))
    return df_proc
