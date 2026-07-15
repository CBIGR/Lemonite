"""
Step-1 preprocessing: HMDB load -> canonical SMILES -> ChEMBL ID mapping.

Produces the annotated metabolite dataframe (with Canonical_smiles and ChEMBL_id
columns) used by every downstream retriever. The ChEMBL-ID mapping is cached so it
is computed once and reused across runs (and by both the chEMBL and LINCS sources).

Re-derived from Collect_PKNdata_metabolites.ipynb cells 5, 13, 16.
"""

import logging
import os

import numpy as np
import pandas as pd

import config
from utils import file_io
from utils.api_retry import retry_api_call
from utils.cache import apply_function_with_multithreading

logger = logging.getLogger('pkn.preprocessing')


@retry_api_call(db_name='ChEMBL_Mapping')
def _chembl_id_from_smiles(canonical_smiles):
    """Look up a ChEMBL molecule id from a canonical SMILES (exact structure match)."""
    from chembl_webresource_client.new_client import new_client
    if pd.isna(canonical_smiles) or canonical_smiles in ('NA', ''):
        return 'none'
    mol_data = list(new_client.molecule.filter(
        molecule_structures__canonical_smiles=canonical_smiles))
    if not mol_data:
        return 'none'
    return mol_data[0]['molecule_chembl_id']


def load_or_parse_hmdb(max_metabolites=None):
    """Load the cached annotated HMDB table or parse the XML to create it."""
    if os.path.exists(config.HMDB_ANNOTATED_FILE):
        df = pd.read_csv(config.HMDB_ANNOTATED_FILE, sep='\t')
        logger.info("Loaded %d metabolites from %s", len(df), config.HMDB_ANNOTATED_FILE)
    else:
        df = file_io.parse_hmdb(config.HMDB_METABOLITES_XML, max_metabolites=max_metabolites)
        df.to_csv(config.HMDB_ANNOTATED_FILE, index=False, sep='\t')
        logger.info("Wrote annotated HMDB table -> %s", config.HMDB_ANNOTATED_FILE)
    if max_metabolites is not None:
        df = df.head(max_metabolites).copy()
        logger.info("Limited to first %d metabolites", max_metabolites)
    return df


def preprocess_metabolites(df, resume=True):
    """
    Add Canonical_smiles and ChEMBL_id columns to the metabolite table.

    Canonical SMILES are computed with RDKit. ChEMBL ids are looked up from the
    canonical SMILES via the ChEMBL web client (threaded, retry-wrapped) and cached
    in CHEMBL_MAPPING_FILE so subsequent runs skip the API calls.
    """
    if 'Canonical_smiles' not in df.columns:
        logger.info("Canonicalizing SMILES for %d metabolites...", len(df))
        df['Canonical_smiles'] = df['SMILES'].apply(file_io.smiles_to_canonical)

    # Load existing ChEMBL mapping (resume). Skip the merge if df already carries a
    # ChEMBL_id column (e.g. the annotated table was persisted by a prior run) to
    # keep this function idempotent and avoid a duplicate-column merge error.
    if resume and 'ChEMBL_id' not in df.columns and os.path.exists(config.CHEMBL_MAPPING_FILE):
        mapping = pd.read_csv(config.CHEMBL_MAPPING_FILE, sep='\t')
        df = df.merge(mapping[['HMDB', 'ChEMBL_id']], on='HMDB', how='left')
        logger.info("Loaded %d cached ChEMBL mappings", mapping['ChEMBL_id'].notna().sum())
    if 'ChEMBL_id' not in df.columns:
        df['ChEMBL_id'] = 'none'

    needs = df['ChEMBL_id'].isna() | (df['ChEMBL_id'] == 'none')
    to_map = df[needs & df['Canonical_smiles'].notna()]
    logger.info("ChEMBL mapping: %d already mapped, %d to map",
                (~needs).sum(), len(to_map))

    if len(to_map):
        ids = apply_function_with_multithreading(
            to_map, 'Canonical_smiles', 'ChEMBL_Mapping', _chembl_id_from_smiles)
        for idx, val in ids.items():
            df.loc[idx, 'ChEMBL_id'] = val
        df[['HMDB', 'Canonical_smiles', 'ChEMBL_id']].to_csv(
            config.CHEMBL_MAPPING_FILE, index=False, sep='\t')
        logger.info("Mapped %d metabolites to ChEMBL ids", (df['ChEMBL_id'] != 'none').sum())

    # Persist Canonical_smiles + ChEMBL_id back into the annotated table so that
    # step 3's URL remapping (Canonical_SMILES -> HMDB) works. Defensively drop any
    # transient merge-suffix columns so the persisted table never accumulates
    # ChEMBL_id_x / _y artefacts across re-runs.
    stray = [c for c in df.columns if c.endswith('_x') or c.endswith('_y')]
    if stray:
        df = df.drop(columns=stray)
    df.to_csv(config.HMDB_ANNOTATED_FILE, index=False, sep='\t')
    return df
