"""
IntAct metabolite-protein interactions (EBI REST API, ChEBI query).

Queries the IntAct interaction list endpoint per metabolite ChEBI id, keeps human
(taxId 9606) interactions, and accumulates per-interaction link rows with the
IntAct interaction-details URL. Checkpointed for safe resume.

Re-derived from Collect_PKNdata_metabolites.ipynb cell 44.
"""

import json
import logging

import pandas as pd
from requests.exceptions import RequestException

import config
from utils import cache
from utils.api_retry import retry_api_call
from utils.http import get_session

logger = logging.getLogger('pkn.intact')


@retry_api_call(db_name='IntAct')
def get_intact_interactions(chebi_id, min_mi_score=0):
    """Return pipe-joined gene symbols for a ChEBI id, storing link rows globally."""
    if chebi_id == 'NA' or pd.isna(chebi_id):
        return 'none'
    try:
        chebi_id = int(chebi_id)
    except (ValueError, TypeError):
        return 'none'

    timeout = config.API_RETRY_CONFIG['IntAct']['timeout']
    url = (f'https://www.ebi.ac.uk/intact/ws/interaction/list?draw=1&maxMIScore=1'
           f'&minMIScore={min_mi_score}&negativeFilter=POSITIVE_ONLY&page=0'
           f'&pageSize=10000&query=CHEBI%3A{chebi_id}')

    response = get_session().post(url, timeout=timeout)
    if response.status_code != 200:
        raise RequestException(f"HTTP {response.status_code}")

    data = json.loads(response.text)
    idf = pd.DataFrame(data.get('data', []))
    if idf.empty:
        return 'none'

    idf = idf[idf['taxIdAStyled'].str.contains('9606')
              | idf['taxIdBStyled'].str.contains('9606')]

    store = cache.url_data_storage.setdefault('IntAct', [])
    interactions = []
    for _, row in idf.iterrows():
        ebi = row['ac']
        gene = None
        if row['idA'] == f'CHEBI:{chebi_id} (chebi)':
            gene = row['moleculeB']
        elif row['idB'] == f'CHEBI:{chebi_id} (chebi)':
            gene = row['moleculeA']
        if gene:
            interactions.append(gene)
            store.append({
                'ChEBI_ID': chebi_id, 'Gene': gene, 'Interaction_ID': ebi,
                'URL': config.URL_TEMPLATES['IntAct'].format(ebi_id=ebi),
            })
    return '|'.join(interactions) if interactions else 'none'


def run(df, resume=True):
    """Process IntAct for all metabolites with checkpointing."""
    return cache.process_single_database(
        df, 'IntAct', 'ChEBI', get_intact_interactions, resume=resume)
