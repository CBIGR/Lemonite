"""
ChEMBL bioactivity metabolite-gene interactions (REST API via web client).

Uses the pre-computed ChEMBL id (from preprocessing) to fetch human bioactivities,
keeps records that indicate a genuine interaction, maps targets to gene symbols, and
accumulates per-interaction link rows.

ACTIVITY FILTER: the reference notebook kept only activities whose ``activity_comment``
was exactly ``active`` or ``substrate``. In current ChEMBL this drops many genuine
interactions: substrates and inhibitors whose comment carries an embedded constant
(e.g. ``substrate [Km=330uM]``, ``inhibitor [487uM]``) and the large fraction of
activities that have no comment but a measured potency. The filter is therefore
broadened to accept (a) comments containing an interaction keyword (active, substrate,
inhibitor, agonist, antagonist, binder, inducer) while excluding explicit negatives,
or (b) a measured pChEMBL value at or above ACTIVE_PCHEMBL_THRESHOLD (potency of about
10 micromolar or better, consistent with the LINCS IC50 threshold used elsewhere).

URL FIX: the reference notebook used the deprecated compound_report_card /
target_report_card URLs, which now serve only an empty single-page-app shell. This
module uses the current canonical explore/compound and explore/target browse URLs
(config.URL_TEMPLATES['chEMBL_compound' / 'chEMBL_target']).

Re-derived from Collect_PKNdata_metabolites.ipynb cell 47.
"""

import logging

import pandas as pd

import config
from utils import cache
from utils.api_retry import retry_api_call

logger = logging.getLogger('pkn.chembl')

# Interaction keywords accepted as a substring of activity_comment.
_ACTIVE_KEYWORDS = ('active', 'substrate', 'inhibitor', 'agonist',
                    'antagonist', 'binder', 'inducer')
# Comments that explicitly indicate no interaction (override the keyword match).
_NEGATIVE_KEYWORDS = ('not active', 'inactive', 'not determined', 'inconclusive')
# pChEMBL >= 5 corresponds to a potency of about 10 micromolar or better.
ACTIVE_PCHEMBL_THRESHOLD = 5.0


def _filter_active(dat):
    """Keep activities that indicate a genuine metabolite-target interaction."""
    comment = dat.get('activity_comment')
    comment = comment.fillna('').str.lower() if comment is not None \
        else pd.Series([''] * len(dat), index=dat.index)

    is_negative = comment.apply(lambda c: any(n in c for n in _NEGATIVE_KEYWORDS))
    has_keyword = comment.apply(lambda c: any(k in c for k in _ACTIVE_KEYWORDS))

    pchembl = pd.to_numeric(dat.get('pchembl_value'), errors='coerce')
    is_potent = pchembl >= ACTIVE_PCHEMBL_THRESHOLD

    keep = (has_keyword & ~is_negative) | is_potent.fillna(False)
    return dat[keep]


@retry_api_call(db_name='chEMBL')
def get_chembl_interactions(chembl_id, canonical_smiles=None):
    """Return pipe-joined gene symbols for a ChEMBL id, storing link rows globally."""
    if chembl_id in ('none', '', 'NA') or pd.isna(chembl_id):
        return 'none'

    from chembl_webresource_client.new_client import new_client

    act = new_client.activity.filter(molecule_chembl_id=chembl_id, target_tax_id='9606')
    dat = pd.DataFrame(act)
    if dat.empty:
        return 'none'

    dat = _filter_active(dat)
    if dat.empty:
        return 'none'
    target_ids = dat['target_chembl_id'].dropna().unique()
    if target_ids.size == 0:
        return 'none'

    store = cache.url_data_storage.setdefault('chEMBL', [])
    interactions = []
    for target in target_ids:
        mol = new_client.target.filter(chembl_id=target)
        if not mol:
            continue
        components = mol[0].get('target_components', [])
        if not components or 'target_component_synonyms' not in components[0]:
            continue
        syn = pd.DataFrame(components[0]['target_component_synonyms'])
        if syn.empty or 'syn_type' not in syn.columns:
            continue
        symbols = syn[syn['syn_type'] == 'GENE_SYMBOL']['component_synonym'].values
        if len(symbols) == 0:
            continue
        gene = symbols[0]
        interactions.append(gene)
        store.append({
            'Canonical_SMILES': canonical_smiles, 'Gene': gene,
            'ChEMBL_ID': chembl_id, 'Target_ID': target,
            'URL': config.URL_TEMPLATES['chEMBL_compound'].format(chembl_id=chembl_id),
            'URL_target': config.URL_TEMPLATES['chEMBL_target'].format(target_id=target),
        })
    return '|'.join(interactions) if interactions else 'none'


def run(df, resume=True):
    """Process ChEMBL for all metabolites with checkpointing."""
    return cache.process_single_database(
        df, 'chEMBL', 'ChEMBL_id', get_chembl_interactions, resume=resume)
