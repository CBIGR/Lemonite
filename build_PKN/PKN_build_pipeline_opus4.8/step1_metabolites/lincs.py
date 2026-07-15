"""
LINCS / LSP biochemical binding metabolite-gene interactions (local files).

No API calls: pure indexed lookups ChEMBL id -> lspci_id -> gene targets below an
IC50 threshold (human targets only). Per-interaction provenance prefers a PubMed
reference; only when none exists does it fall back to the HMS-LINCS small-molecule
browser.

URL FIX: the reference notebook fell back to
labsyspharm.github.io/lsp-cheminformatics/compounds/{id}, which is 404. The fallback
is now config.URL_TEMPLATES['LINCS_fallback'] (HMS-LINCS browser).

Re-derived from Collect_PKNdata_metabolites.ipynb cells 32-38.
"""

import logging

import numpy as np
import pandas as pd

import config

logger = logging.getLogger('pkn.lincs')


def _build_indexes():
    """Build ChEMBL->lspci, lspci->genes, and (lspci,target)->pubmed lookups."""
    compound_mapping = pd.read_csv(config.LINCS_COMPOUND_MAPPING)
    target_dictionary = pd.read_csv(config.LINCS_TARGET_DICTIONARY)
    references = pd.read_csv(config.LINCS_REFERENCES)

    chembl_to_lspci = {}
    chembl_rows = compound_mapping[compound_mapping['source'] == 'chembl']
    for _, r in chembl_rows.iterrows():
        chembl_to_lspci.setdefault(r['external_id'], []).append(r['lspci_id'])

    biochem = pd.read_csv(
        config.LINCS_BIOCHEM_AGG,
        usecols=['lspci_id', 'lspci_target_id', 'value', 'value_unit', 'symbol'])
    biochem['value_unit'] = biochem['value_unit'].str.lower().str.strip()
    biochem.loc[biochem['value_unit'] == 'um', 'value'] *= 1000  # uM -> nM

    human = set(target_dictionary[target_dictionary['tax_id'] == 9606]['lspci_target_id'])
    biochem = biochem[biochem['lspci_target_id'].isin(human)]

    lspci_to_genes = {}
    for _, r in biochem.iterrows():
        lspci_to_genes.setdefault(r['lspci_id'], []).append(
            {'symbol': r['symbol'], 'value': r['value'], 'lspci_target_id': r['lspci_target_id']})

    pubmed_index = {}
    pubmed_refs = references[references['reference_type'] == 'pubmed']
    for _, r in pubmed_refs.iterrows():
        key = (r['lspci_id'], r.get('lspci_target_id'))
        pubmed_index.setdefault(key, str(int(r['reference_value'])))

    return chembl_to_lspci, lspci_to_genes, pubmed_index


def run(df):
    """Build LINCS interactions; write links + processed files (no API)."""
    chembl_to_lspci, lspci_to_genes, pubmed_index = _build_indexes()
    thr = config.LINCS_IC50_THRESHOLD

    lincs_dict, link_rows = {}, []
    for _, row in df.iterrows():
        chembl_id = row.get('ChEMBL_id')
        if pd.isna(chembl_id) or chembl_id == 'none':
            continue
        lspci_ids = chembl_to_lspci.get(chembl_id, [])
        if not lspci_ids:
            continue

        genes = set()
        for lspci_id in lspci_ids:
            for inter in lspci_to_genes.get(lspci_id, []):
                if inter['value'] <= thr and pd.notna(inter['symbol']):
                    genes.add(inter['symbol'])
        if not genes:
            continue
        lincs_dict[row['HMDB']] = genes

        for gene in genes:
            pmid_url = None
            for lspci_id in lspci_ids:
                for inter in lspci_to_genes.get(lspci_id, []):
                    if inter['symbol'] == gene:
                        pmid = pubmed_index.get((lspci_id, inter.get('lspci_target_id')))
                        if pmid:
                            pmid_url = config.URL_TEMPLATES['PubMed'].format(pmid=pmid)
                            break
                if pmid_url:
                    break
            # Fallback when no PubMed reference is indexed: the LINCS biochemical
            # binding values are ChEMBL-derived, so the ChEMBL compound page shows
            # the actual measured bioactivities (real evidence for the binding),
            # which is far better than the generic HMS-LINCS browser landing page.
            if not pmid_url:
                pmid_url = config.URL_TEMPLATES['chEMBL_compound'].format(chembl_id=chembl_id)
            link_rows.append({
                'HMDB': row['HMDB'], 'ChEMBL_ID': chembl_id,
                'lspci_id': lspci_ids[0], 'Gene': gene,
                'IC50_threshold_nM': thr, 'URL': pmid_url,
            })

    df_proc = df[['HMDB', 'ChEMBL_id']].copy()
    df_proc['LINCS'] = df_proc['HMDB'].apply(
        lambda x: '|'.join(sorted(lincs_dict[x])) if x in lincs_dict else np.nan)
    df_proc.to_csv(config.DB_OUTPUT_FILES['LINCS'], index=False, sep='\t')

    if link_rows:
        pd.DataFrame(link_rows).to_csv(config.URL_FILES['LINCS'], index=False, sep='\t')
    logger.info("LINCS: %d metabolites, %d interaction links",
                len(lincs_dict), len(link_rows))
    return df_proc
