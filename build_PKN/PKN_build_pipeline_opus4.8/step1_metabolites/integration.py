"""
Step-1 integration: merge per-database processed files into the long-format
metabolite-gene PKN, attaching provenance URLs from the per-database link files.

Outputs:
  * OUTPUT_FILE_FINAL          wide table (one row per metabolite, one col per DB)
  * METABOLITE_GENE_PKN        long table: Metabolite, Gene, Source, url

Re-derived from Collect_PKNdata_metabolites.ipynb cells 63 & 67.
"""

import logging
import os

import numpy as np
import pandas as pd

import config

logger = logging.getLogger('pkn.integration')

SOURCES = ['UniProtKB', 'STITCH', 'IntAct', 'BioGRID', 'chEMBL', 'LINCS',
           'Human1_GEM_dist1', 'Human1_GEM_dist2', 'MetalinksDB']


def _merge_wide(df_base):
    """Merge each processed DB column onto the base annotation table by HMDB."""
    df = df_base.copy()
    for db_name, path in config.DB_OUTPUT_FILES.items():
        if not os.path.exists(path):
            logger.warning("Missing processed file for %s: %s", db_name, path)
            continue
        df_db = pd.read_csv(path, sep='\t', low_memory=False)
        if db_name not in df_db.columns:
            continue
        df = df.merge(df_db[['HMDB', db_name]], on='HMDB', how='left',
                      suffixes=('', f'_{db_name}'))
        dup = f'{db_name}_{db_name}'
        if dup in df.columns:
            df[db_name] = df[dup].fillna(df[db_name])
            df.drop(columns=[dup], inplace=True)
    df.to_csv(config.OUTPUT_FILE_FINAL, index=False, sep='\t')
    return df


def _build_url_lookup(df_final):
    """Build a (HMDB, Gene, Source) -> URL lookup from the per-database link files."""
    # ID -> HMDB reverse maps (some link files are keyed by non-HMDB ids)
    inchikey_to_hmdb = dict(zip(df_final['InChIKey'].dropna(), df_final['HMDB']))
    chebi_to_hmdb = {}
    for _, r in df_final.dropna(subset=['ChEBI']).iterrows():
        chebi_to_hmdb[int(r['ChEBI'])] = r['HMDB']
    chembl_to_hmdb = {}
    if 'ChEMBL_id' in df_final.columns:
        for _, r in df_final.iterrows():
            if pd.notna(r.get('ChEMBL_id')) and r['ChEMBL_id'] != 'none':
                chembl_to_hmdb[str(r['ChEMBL_id'])] = r['HMDB']

    lookup = {}

    def _add(path, source, hmdb_resolver, gene_col, url_col, row_filter=None):
        if not path or not os.path.exists(path):
            return
        d = pd.read_csv(path, sep='\t', low_memory=False)
        if gene_col not in d.columns or url_col not in d.columns:
            return
        if row_filter is not None:
            d = d[row_filter(d)]
        for _, r in d.iterrows():
            hmdb = hmdb_resolver(r)
            if hmdb is None or pd.isna(hmdb):
                continue
            gene = r[gene_col]
            url = r[url_col]
            if pd.notna(gene) and pd.notna(url):
                # First reaction wins as the representative trace URL for the edge.
                lookup.setdefault((hmdb, str(gene), source), url)

    def _add_gem(path, source, row_filter):
        """GEM link rows carry two reaction URLs; join them for distance-2 edges.

        A distance-1 edge has only URL_1; a distance-2 edge has URL_1 (metabolite ->
        intermediate) and URL_2 (intermediate -> gene), joined so the single lookup
        URL still exposes both connecting reactions.
        """
        if not path or not os.path.exists(path):
            return
        d = pd.read_csv(path, sep='\t', low_memory=False)
        if 'URL_1' not in d.columns or 'Gene' not in d.columns:
            return
        d = d[row_filter(d)]
        for _, r in d.iterrows():
            hmdb, gene, url1 = r.get('HMDB'), r.get('Gene'), r.get('URL_1')
            if pd.isna(hmdb) or pd.isna(gene) or pd.isna(url1):
                continue
            url2 = r.get('URL_2')
            url = f"{url1} -> {url2}" if pd.notna(url2) else url1
            lookup.setdefault((hmdb, str(gene), source), url)

    _add(config.URL_FILES['BioGRID'], 'BioGRID', lambda r: r.get('HMDB'), 'BioGRID', 'url')
    _add(config.URL_FILES['STITCH'], 'STITCH', lambda r: r.get('HMDB'), 'STITCH', 'url')
    _add(config.URL_FILES['LINCS'], 'LINCS', lambda r: r.get('HMDB'), 'Gene', 'URL')
    _add(config.URL_FILES['MetalinksDB'], 'MetalinksDB', lambda r: r.get('HMDB'), 'Gene', 'URL')
    _add(config.URL_FILES['UniProtKB'], 'UniProtKB',
         lambda r: inchikey_to_hmdb.get(r.get('InChIKey')), 'Gene', 'URL')
    _add(config.URL_FILES['IntAct'], 'IntAct',
         lambda r: chebi_to_hmdb.get(int(r['ChEBI_ID'])) if pd.notna(r.get('ChEBI_ID')) else None,
         'Gene', 'URL')
    _add(config.URL_FILES['chEMBL'], 'chEMBL',
         lambda r: chembl_to_hmdb.get(str(r.get('ChEMBL_ID'))), 'Gene', 'URL')
    # GEM dist1/dist2 share one link file keyed by HMDB. Each row carries the
    # Metabolic Atlas URL(s) of the reaction(s) that connect metabolite -> gene:
    # one reaction at distance 1 (URL_1), two at distance 2 (URL_1 -> URL_2). A
    # distance-2 edge cannot be represented by a single URL, so both are joined
    # with ' -> '. Filter by Distance so each source maps to its own reactions.
    _add_gem(config.URL_FILES['Human1_GEM_dist1'], 'Human1_GEM_dist1',
             lambda d: d['Distance'] == 1)
    _add_gem(config.URL_FILES['Human1_GEM_dist2'], 'Human1_GEM_dist2',
             lambda d: d['Distance'] == 2)

    return lookup


def integrate(df_annotated):
    """Build the wide table and the long-format metabolite_gene_PKN.tsv with URLs."""
    df_final = _merge_wide(df_annotated)
    url_lookup = _build_url_lookup(df_final)

    rows = []
    for _, row in df_final.iterrows():
        metabolite = f"{row['Name']}_{row['HMDB']}"
        for source in SOURCES:
            val = row.get(source)
            if pd.isna(val) or val == 'none':
                continue
            for gene in str(val).split('|'):
                gene = gene.strip()
                if gene and gene != 'none':
                    url = url_lookup.get((row['HMDB'], gene, source))
                    rows.append({'Metabolite': metabolite, 'Gene': gene,
                                 'Source': source, 'url': url})

    pkn = pd.DataFrame(rows)
    pkn['Gene'] = pkn['Gene'].replace('', np.nan)
    pkn = pkn.dropna(subset=['Gene']).drop_duplicates(subset=['Metabolite', 'Gene', 'Source'])
    pkn.to_csv(config.METABOLITE_GENE_PKN, index=False, sep='\t')

    covered = pkn['url'].notna().sum()
    logger.info("metabolite_gene_PKN: %d edges (%d with URL, %.1f%%)",
                len(pkn), covered, 100 * covered / max(len(pkn), 1))
    return pkn
