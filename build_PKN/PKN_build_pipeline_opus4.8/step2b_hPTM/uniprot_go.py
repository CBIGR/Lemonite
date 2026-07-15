"""
UniProtKB enzyme - histone-mark annotations (REST API, GO-function filter).

A second, independent source for the same enzyme-histone-mark edges as QuickGO,
queried from UniProtKB by GO molecular-function term (config.HISTONE_GO_TERMS).
UniProt's ``go:`` filter matches the term and its descendants over the GO
hierarchy, so each curated parent term captures the leaf enzymes too. Restricted
to reviewed (Swiss-Prot) human entries for high precision.

Running both QuickGO and UniProtKB gives cross-validation: an edge found by both
carries two Source labels in the final network. Whole-file cache like QuickGO /
STRING. Each row deep-links to the enzyme's UniProt function/GO section.

Re-uses the pagination idiom from step1_metabolites/uniprot.py (RFC5988 Link
header cursor via response.links).
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from requests.exceptions import RequestException

import config
from utils.api_retry import retry_api_call
from utils.http import get_session

logger = logging.getLogger('pkn.uniprot_hptm')

_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
_PAGE_SIZE = 500

_OUTPUT_COLS = ['Gene', 'UniProt_ID', 'Mark', 'Activity', 'GO_ID', 'Source', 'URL']


@retry_api_call(db_name='UniProt_hPTM')
def _fetch_term(go_id):
    """Fetch reviewed human proteins annotated with ``go_id`` (self + descendants).

    Returns a list of (gene, accession) tuples. Raises on a bad HTTP status so the
    retry decorator can take over.
    """
    session = get_session()
    timeout = config.API_RETRY_CONFIG['UniProt_hPTM']['timeout']
    query = f"go:{go_id.split(':', 1)[1]} AND organism_id:9606 AND reviewed:true"
    url = (f"{_SEARCH_URL}?query={query}"
           f"&fields=accession,gene_primary&format=tsv&size={_PAGE_SIZE}")
    out = []
    while url:
        resp = session.get(url, timeout=timeout)
        if resp.status_code != 200:
            raise RequestException(f"HTTP {resp.status_code}")
        lines = resp.text.strip().split('\n')
        for line in lines[1:]:  # skip header
            cols = line.split('\t')
            if len(cols) >= 2 and cols[0] and cols[1]:
                out.append((cols[1], cols[0]))  # (gene, accession)
        url = resp.links.get('next', {}).get('url')
    return out


def run(resume=True):
    """Retrieve enzyme-histone-mark annotations from UniProtKB by GO term.

    Returns a dataframe with one row per (gene, mark, activity, GO term) with a
    provenance URL. Cached whole-file to UNIPROT_HPTM_FILE.
    """
    if resume and os.path.exists(config.UNIPROT_HPTM_FILE):
        df = pd.read_csv(config.UNIPROT_HPTM_FILE)
        logger.info("UniProt hPTM: loaded %d cached annotations", len(df))
        return df

    go_terms = config.HISTONE_GO_TERMS
    logger.info("UniProt hPTM: querying %d histone GO terms", len(go_terms))

    max_workers = config.API_RETRY_CONFIG['UniProt_hPTM']['max_workers']
    rows = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_fetch_term, go_id): go_id for go_id in go_terms}
        for future in as_completed(futures):
            go_id = futures[future]
            mark, activity = go_terms[go_id]
            res = future.result()
            if not isinstance(res, list):
                continue  # 'none' sentinel from an exhausted-retry term
            for gene, acc in res:
                rows.append({
                    'Gene': gene,
                    'UniProt_ID': acc,
                    'Mark': mark,
                    'Activity': activity,
                    'GO_ID': go_id,
                    'Source': 'UniProtKB_GO',
                    'URL': config.URL_TEMPLATES['UniProt_entry'].format(accession=acc),
                })

    df = pd.DataFrame(rows, columns=_OUTPUT_COLS)
    if not df.empty:
        df = df.drop_duplicates(subset=['Gene', 'Mark', 'Activity', 'GO_ID'])
    df.to_csv(config.UNIPROT_HPTM_FILE, index=False)
    logger.info("UniProt hPTM: %d enzyme-histone-mark annotations (%d genes)",
                len(df), df['Gene'].nunique() if not df.empty else 0)
    return df
