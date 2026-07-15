"""
QuickGO enzyme - histone-mark annotations (EBI REST API).

Queries the QuickGO annotation-search endpoint for every residue-specific histone
GO term in config.HISTONE_GO_TERMS (batched, human only, protein gene products),
mapping each annotated gene product to the histone mark and activity (writer /
eraser / reader) that term encodes. QuickGO expands each queried term to its
descendants, so a curated set of parent terms captures the leaf enzymes too.

Each row is one (enzyme gene, histone mark, activity) association carrying the GO
term, the GO evidence code, its reference (GO_REF or PMID), and a provenance URL
that deep-links into the QuickGO annotation viewer filtered to exactly that gene
product + GO term — so the user sees the evidence for that specific activity.

Whole-file cache (like step2 STRING): one query pass, cached to QUICKGO_HPTM_FILE.
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from requests.exceptions import RequestException

import config
from utils.api_retry import retry_api_call
from utils.http import get_session

logger = logging.getLogger('pkn.quickgo')

_ANNOTATION_URL = 'https://www.ebi.ac.uk/QuickGO/services/annotation/search'

_OUTPUT_COLS = ['Gene', 'UniProt_ID', 'Mark', 'Activity', 'GO_ID',
                'GO_Evidence', 'Reference', 'Source', 'URL']


@retry_api_call(db_name='QuickGO')
def _fetch_term(go_id):
    """Fetch all human protein annotations for one GO id, self + descendants.

    Returns a list of (gene_symbol, uniprot_accession, go_evidence, reference)
    tuples for the queried term. ``goUsage=descendants`` ensures enzymes annotated
    only to a more specific child of ``go_id`` are still captured under the queried
    (parent) term's mark/activity. Raises on a bad HTTP status so the retry
    decorator can take over.
    """
    session = get_session()
    timeout = config.API_RETRY_CONFIG['QuickGO']['timeout']
    out = []
    page = 1
    while True:
        params = {
            'goId': go_id,
            'goUsage': 'descendants',
            'taxonId': config.QUICKGO_TAXON,
            'geneProductType': 'protein',
            'limit': config.QUICKGO_PAGE_LIMIT,
            'page': page,
        }
        resp = session.get(_ANNOTATION_URL, params=params,
                           headers={'Accept': 'application/json'}, timeout=timeout)
        if resp.status_code != 200:
            raise RequestException(f"HTTP {resp.status_code}")
        data = resp.json()
        for r in data.get('results', []):
            gp = r.get('geneProductId', '')
            acc = gp.split(':', 1)[1] if ':' in gp else gp
            out.append((r.get('symbol'), acc, r.get('goEvidence'), r.get('reference')))
        info = data.get('pageInfo') or {}
        if not info or info.get('current', page) >= info.get('total', page):
            break
        page += 1
    return out


def _reference_url(reference):
    """A PubMed link when the reference is a PMID, else None (QuickGO URL is used)."""
    if isinstance(reference, str) and reference.startswith('PMID:'):
        return config.URL_TEMPLATES['PubMed'].format(pmid=reference.split(':', 1)[1])
    return None


def run(resume=True):
    """Retrieve enzyme-histone-mark annotations from QuickGO.

    Returns a dataframe with one row per (gene, mark, activity, GO term) with a
    provenance URL. Cached whole-file to QUICKGO_HPTM_FILE.
    """
    if resume and os.path.exists(config.QUICKGO_HPTM_FILE):
        df = pd.read_csv(config.QUICKGO_HPTM_FILE)
        logger.info("QuickGO: loaded %d cached annotations", len(df))
        return df

    go_terms = config.HISTONE_GO_TERMS
    logger.info("QuickGO: querying %d histone GO terms", len(go_terms))

    max_workers = config.API_RETRY_CONFIG['QuickGO']['max_workers']
    rows = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_fetch_term, go_id): go_id for go_id in go_terms}
        for future in as_completed(futures):
            go_id = futures[future]
            mark, activity = go_terms[go_id]
            res = future.result()
            if not isinstance(res, list):
                continue  # 'none' sentinel from an exhausted-retry term
            for gene, acc, evidence, reference in res:
                if not gene or not acc:
                    continue
                url = _reference_url(reference) or \
                    config.URL_TEMPLATES['QuickGO_annotation'].format(accession=acc, go_id=go_id)
                rows.append({
                    'Gene': gene,
                    'UniProt_ID': acc,
                    'Mark': mark,
                    'Activity': activity,
                    'GO_ID': go_id,
                    'GO_Evidence': evidence,
                    'Reference': reference,
                    'Source': 'QuickGO',
                    'URL': url,
                })

    df = pd.DataFrame(rows, columns=_OUTPUT_COLS)
    if not df.empty:
        df = df.drop_duplicates(subset=['Gene', 'Mark', 'Activity', 'GO_ID'])
    df.to_csv(config.QUICKGO_HPTM_FILE, index=False)
    logger.info("QuickGO: %d enzyme-histone-mark annotations (%d genes)",
                len(df), df['Gene'].nunique() if not df.empty else 0)
    return df
