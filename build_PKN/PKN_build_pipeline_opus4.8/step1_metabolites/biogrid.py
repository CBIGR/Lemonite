"""
BioGRID chemical-protein interactions (local file).

Matches HMDB metabolites to BioGRID chemicals by InChIKey, emits per-interaction
link rows with the BioGRID chemical page URL + supporting PubMed URL, then writes
the processed (one-row-per-metabolite) cache.

Re-derived from Collect_PKNdata_metabolites.ipynb cells 19-21, 50.
"""

import logging

import pandas as pd

import config

logger = logging.getLogger('pkn.biogrid')


def run(df):
    """Build BioGRID metabolite-gene interactions; write links + processed files."""
    biogrid = pd.read_csv(config.BIOGRID_LOCATION, sep='\t')
    biogrid = biogrid[biogrid['InChIKey'].isin(df['InChIKey'])]
    logger.info("BioGRID: %d chemical-protein rows after InChIKey filter", len(biogrid))

    # InChIKey -> (chemical_id, chemical_name) for URL construction
    chem_meta = (biogrid.drop_duplicates('InChIKey')
                 .set_index('InChIKey')[['BioGRID Chemical ID', 'Chemical Name']]
                 .to_dict('index'))

    rows = []
    key_to_hmdb = df.dropna(subset=['InChIKey']).set_index('InChIKey')['HMDB'].to_dict()
    for key in biogrid['InChIKey'].unique():
        if key not in key_to_hmdb:
            continue
        meta = chem_meta[key]
        url = config.URL_TEMPLATES['BioGRID'].format(
            chemical_id=meta['BioGRID Chemical ID'], chemical_name=meta['Chemical Name'])
        for _, br in biogrid[biogrid['InChIKey'] == key].iterrows():
            pubmed = br.get('Pubmed ID')
            pmid_url = (config.URL_TEMPLATES['PubMed'].format(pmid=int(pubmed))
                        if pd.notna(pubmed) else None)
            rows.append({'HMDB': key_to_hmdb[key], 'BioGRID': br['Official Symbol'],
                         'url': url, 'pubmed_url': pmid_url})

    links = pd.DataFrame(rows).drop_duplicates(subset=['HMDB', 'BioGRID'])
    links.to_csv(config.URL_FILES['BioGRID'], index=False, sep='\t')
    logger.info("BioGRID: %d interaction links", len(links))

    # Processed: one row per metabolite, pipe-joined genes
    df_proc = df[['HMDB', 'InChIKey']].copy()
    if not links.empty:
        grouped = (links.groupby('HMDB')['BioGRID']
                   .apply(lambda x: '|'.join(x.astype(str))).reset_index())
        df_proc = df_proc.merge(grouped, on='HMDB', how='left')
    else:
        df_proc['BioGRID'] = pd.NA
    df_proc.to_csv(config.DB_OUTPUT_FILES['BioGRID'], index=False, sep='\t')
    return df_proc
