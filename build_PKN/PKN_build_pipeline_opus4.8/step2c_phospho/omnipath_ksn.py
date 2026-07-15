"""
Fetch the OmniPath enzyme-substrate (kinase-substrate) network via the OmniPath web
service (HTTP TSV) — no R / decoupleR dependency.

OmniPath's enzsub domain aggregates ~18 curated + predicted resources (PhosphoSitePlus,
SIGNOR, DEPOD, dbPTM, phospho.ELM, HPRD, KEA, PhosphoNetworks, NetworKIN via MIMP, ...), so
this single pull subsumes most site-level kinase-substrate knowledge. We keep phosphorylation
enzsub edges with gene symbols and collapse to gene level (Node1 = enzyme_genesymbol,
Node2 = substrate_genesymbol), retaining residue as metadata.

Resumable: the raw TSV is cached at config.OMNIPATH_KSN_FILE; a resume run reuses it.
"""

import logging
import os

import pandas as pd

import config
from utils import http

logger = logging.getLogger('pkn.phospho.ksn')

# request gene symbols + provenance sources; restrict to phosphorylation
_PARAMS = {
    'types': 'phosphorylation',
    'genesymbols': '1',
    'fields': 'sources,references',
    'format': 'tsv',
}


def _fetch_raw():
    """GET the enzsub TSV from OmniPath into a DataFrame (retry-hardened session)."""
    sess = http.get_session()
    logger.info("Fetching OmniPath enzsub (phosphorylation) from %s", config.OMNIPATH_ENZSUB_URL)
    resp = sess.get(config.OMNIPATH_ENZSUB_URL, params=_PARAMS, timeout=300)
    resp.raise_for_status()
    from io import StringIO
    df = pd.read_csv(StringIO(resp.text), sep='\t')
    logger.info("OmniPath enzsub returned %d rows", len(df))
    return df


def run(resume=True):
    """Return a gene-level kinase-substrate DataFrame.

    Columns: Node1, Node2, Source, Type, Residue, URL
    """
    cache = config.OMNIPATH_KSN_FILE
    if resume and os.path.exists(cache) and os.path.getsize(cache) > 0:
        logger.info("Using cached OmniPath KSN: %s", cache)
        raw = pd.read_csv(cache, sep='\t')
    else:
        raw = _fetch_raw()
        raw.to_csv(cache, sep='\t', index=False)
        logger.info("Cached OmniPath enzsub raw TSV -> %s", cache)

    if raw.empty or 'enzyme_genesymbol' not in raw.columns:
        logger.warning("OmniPath enzsub empty or missing gene symbols")
        return pd.DataFrame(columns=['Node1', 'Node2', 'Source', 'Type', 'Residue', 'URL'])

    df = raw.copy()
    df['Node1'] = df['enzyme_genesymbol'].astype(str)
    df['Node2'] = df['substrate_genesymbol'].astype(str)
    res_type = df.get('residue_type', pd.Series('', index=df.index)).astype(str)
    res_off = df.get('residue_offset', pd.Series('', index=df.index)).astype(str)
    df['Residue'] = (res_type + res_off).str.replace('nan', '', regex=False)
    df['Source'] = 'OmniPath_enzsub'
    df['Type'] = 'kinase-substrate'
    df['URL'] = config.OMNIPATH_ENZSUB_URL

    df = df[(df['Node1'].str.strip() != '') & (df['Node2'].str.strip() != '')]
    df = df[~df['Node1'].isin(['nan', 'None']) & ~df['Node2'].isin(['nan', 'None'])]
    # gene-level dedupe (a kinase-substrate gene pair can have many sites/refs)
    df = df.drop_duplicates(subset=['Node1', 'Node2'])
    out = df[['Node1', 'Node2', 'Source', 'Type', 'Residue', 'URL']].reset_index(drop=True)
    logger.info("Kinase-substrate edges (gene-level): %d (%d kinases -> %d substrates)",
                len(out), out['Node1'].nunique(), out['Node2'].nunique())
    return out
