"""
STITCH protein-chemical links (local file + MyGene.info for symbol mapping).

Matches HMDB metabolites to STITCH chemicals by PubChem CID, maps STITCH Ensembl
protein ids to gene symbols via MyGene.info, emits per-interaction link rows with
the STITCH network URL, then writes the processed cache.

Re-derived from Collect_PKNdata_metabolites.ipynb cells 23-29, 49.
"""

import logging

import pandas as pd

import config

logger = logging.getLogger('pkn.stitch')


def run(df):
    """Build STITCH metabolite-gene interactions; write links + processed files."""
    stitch = pd.read_csv(config.STITCH_LOCATION, sep='\t')
    stitch['chemical'] = stitch['chemical'].str.replace(r'\D', '', regex=True)
    stitch['protein'] = stitch['protein'].str.replace('9606.', '', regex=False)

    pubchem_ids = df['PubChem'].dropna().astype(int).astype(str).tolist()
    stitch = stitch[stitch['chemical'].isin(pubchem_ids)]
    logger.info("STITCH: %d links after PubChem filter (%d chemicals, %d proteins)",
                len(stitch), stitch['chemical'].nunique(), stitch['protein'].nunique())

    # Ensembl protein -> gene symbol via MyGene.info
    import mygene
    mg = mygene.MyGeneInfo()
    proteins = stitch['protein'].unique().tolist()
    result = mg.querymany(proteins, scopes='ensembl.protein', fields='symbol',
                          species='human', verbose=False)
    protein_to_symbol = {r.get('query'): r.get('symbol') for r in result}
    stitch['symbol'] = stitch['protein'].map(protein_to_symbol)

    pubchem_to_hmdb = {}
    for _, r in df.dropna(subset=['PubChem']).iterrows():
        pubchem_to_hmdb[str(int(r['PubChem']))] = r['HMDB']

    rows = []
    for cid, group in stitch.groupby('chemical'):
        if cid not in pubchem_to_hmdb:
            continue
        url = config.URL_TEMPLATES['STITCH'].format(cid=cid)
        for sym in group['symbol']:
            if isinstance(sym, str) and sym:
                rows.append({'HMDB': pubchem_to_hmdb[cid], 'STITCH': sym, 'url': url})

    links = pd.DataFrame(rows).drop_duplicates(subset=['HMDB', 'STITCH'])
    links.to_csv(config.URL_FILES['STITCH'], index=False, sep='\t',
                 columns=['HMDB', 'STITCH', 'url'] if not links.empty else None)
    logger.info("STITCH: %d interaction links", len(links))

    df_proc = df[['HMDB', 'PubChem']].copy()
    if not links.empty:
        grouped = (links.groupby('HMDB')['STITCH']
                   .apply(lambda x: '|'.join(x.astype(str))).reset_index())
        df_proc = df_proc.merge(grouped, on='HMDB', how='left')
    else:
        df_proc['STITCH'] = pd.NA
    df_proc.to_csv(config.DB_OUTPUT_FILES['STITCH'], index=False, sep='\t')
    return df_proc
