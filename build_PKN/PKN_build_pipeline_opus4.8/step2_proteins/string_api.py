"""
STRING protein-protein interactions (REST API, chunked gene-set query).

For the gene list derived from step 1, resolves STRING ids then fetches interaction
partners in chunks (threaded). Each chunk request is retry-wrapped. Result is cached
to STRING_INTERACTIONS_FILE so re-runs skip the API.

Re-derived from Collect_PKNdata_proteins.ipynb cell 3.
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

import config
from utils.api_retry import retry_api_call
from utils.http import get_session

logger = logging.getLogger('pkn.string')


@retry_api_call(db_name='STRING')
def _get_chunk(gene_list_chunk):
    """Resolve STRING ids for a gene chunk and fetch their interaction partners."""
    session = get_session()
    base = config.STRING_API_URL
    out_fmt = 'tsv-no-header'

    # 1) gene symbols -> STRING ids
    ids_url = '/'.join([base, out_fmt, 'get_string_ids'])
    params = {'identifiers': '\r'.join(gene_list_chunk), 'species': config.STRING_SPECIES,
              'limit': 1, 'echo_query': 1, 'caller_identity': config.STRING_CALLER_IDENTITY}
    res = session.post(ids_url, data=params, timeout=config.API_RETRY_CONFIG['STRING']['timeout'])
    string_ids = []
    for line in res.text.strip().split('\n'):
        cols = line.split('\t')
        if len(cols) >= 3 and cols[2]:
            string_ids.append(cols[2])
    if not string_ids:
        return pd.DataFrame(columns=['query_name', 'partner_name', 'combined_score'])

    # 2) STRING ids -> interaction partners
    part_url = '/'.join([base, out_fmt, 'interaction_partners'])
    params = {'identifiers': '%0d'.join(string_ids), 'species': config.STRING_SPECIES,
              'caller_identity': config.STRING_CALLER_IDENTITY}
    res = session.post(part_url, data=params, timeout=config.API_RETRY_CONFIG['STRING']['timeout'])
    rows = []
    for line in res.text.strip().split('\n'):
        cols = line.strip().split('\t')
        if len(cols) >= 6:
            rows.append({'query_name': cols[2], 'partner_name': cols[3], 'combined_score': cols[5]})
    return pd.DataFrame(rows)


def run(genes, resume=True):
    """Fetch the STRING network for ``genes``; return GeneA/GeneB/combined_score df."""
    if resume and os.path.exists(config.STRING_INTERACTIONS_FILE):
        df = pd.read_csv(config.STRING_INTERACTIONS_FILE)
        logger.info("STRING: loaded %d cached interactions", len(df))
        return df

    genes = list(genes)
    chunks = [genes[i:i + config.STRING_CHUNK_SIZE]
              for i in range(0, len(genes), config.STRING_CHUNK_SIZE)]
    logger.info("STRING: %d genes in %d chunks", len(genes), len(chunks))

    results = []
    with ThreadPoolExecutor(max_workers=config.STRING_MAX_WORKERS) as executor:
        futures = {executor.submit(_get_chunk, c): i for i, c in enumerate(chunks)}
        for future in as_completed(futures):
            res = future.result()
            if isinstance(res, pd.DataFrame) and not res.empty:
                results.append(res)

    df = pd.concat(results, ignore_index=True) if results else \
        pd.DataFrame(columns=['query_name', 'partner_name', 'combined_score'])
    df.columns = ['GeneA', 'GeneB', 'combined_score']
    df.to_csv(config.STRING_INTERACTIONS_FILE, index=False)
    logger.info("STRING: %d interactions", len(df))
    return df
