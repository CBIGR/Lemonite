"""
Preprocessing module for metabolite data.

Handles:
- Canonical SMILES generation
- ChEMBL ID mapping
- Metabolite dataframe preparation
"""

import pandas as pd
import numpy as np
import logging
from typing import List, Dict
import os
import sys

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.api_retry import retry_api_call
from concurrent.futures import ThreadPoolExecutor, as_completed
import config


def smiles_to_canonical(smiles: str) -> str:
    """
    Convert SMILES to canonical SMILES using RDKit.
    
    Parameters:
    -----------
    smiles : str
        SMILES string
    
    Returns:
    --------
    str
        Canonical SMILES or NaN if conversion fails
    """
    from rdkit import Chem
    
    try:
        if smiles != 'NA' and pd.notna(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                return canonical
    except Exception:
        pass
    
    return np.nan


@retry_api_call(db_name='ChEMBL_Mapping')
def get_chembl_id_from_smiles(canonical_smiles: str) -> str:
    """
    Get ChEMBL ID from canonical SMILES string.
    
    Parameters:
    -----------
    canonical_smiles : str
        Canonical SMILES string
    
    Returns:
    --------
    str
        ChEMBL ID or 'none' if not found
    """
    from chembl_webresource_client.new_client import new_client
    
    if pd.isna(canonical_smiles) or canonical_smiles == 'NA' or canonical_smiles == '':
        return 'none'
    
    molecule = new_client.molecule
    mol_data = molecule.filter(molecule_structures__canonical_smiles=canonical_smiles)
    mol_data = list(mol_data)
    
    if not mol_data:
        return 'none'
    
    return mol_data[0]['molecule_chembl_id']


def add_canonical_smiles(metabolites_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add canonical SMILES column to metabolite dataframe.
    
    Parameters:
    -----------
    metabolites_df : pd.DataFrame
        Must contain 'SMILES' column
    
    Returns:
    --------
    pd.DataFrame
        With added 'Canonical_smiles' column
    """
    logger = logging.getLogger('preprocessing')
    logger.info("Adding canonical SMILES...")
    
    if 'SMILES' not in metabolites_df.columns:
        logger.warning("No SMILES column found, skipping canonical SMILES generation")
        metabolites_df['Canonical_smiles'] = np.nan
        return metabolites_df
    
    metabolites_df['Canonical_smiles'] = metabolites_df['SMILES'].apply(smiles_to_canonical)
    
    count = metabolites_df['Canonical_smiles'].notna().sum()
    logger.info(f"✓ Generated canonical SMILES for {count}/{len(metabolites_df)} metabolites")
    
    return metabolites_df


def map_chembl_ids(metabolites_df: pd.DataFrame, max_workers: int = None) -> pd.DataFrame:
    """
    Map canonical SMILES to ChEMBL IDs using multithreading.
    
    Caches results for resume capability.
    
    Parameters:
    -----------
    metabolites_df : pd.DataFrame
        Must contain 'HMDB', 'Canonical_smiles' columns
    max_workers : int, optional
        Number of concurrent API workers
    
    Returns:
    --------
    pd.DataFrame
        With added 'ChEMBL_id' column
    """
    logger = logging.getLogger('preprocessing')
    
    # Use config for max workers
    if max_workers is None:
        max_workers = config.API_RETRY_CONFIG['ChEMBL_Mapping']['max_workers']
    
    mapping_file = os.path.join(config.OUTPUT_DIR, 'HMDB_metabolites_ChEMBL_mapping.csv')
    
    # Load existing mapping if available
    if os.path.exists(mapping_file):
        logger.info(f"Loading existing ChEMBL mapping from {mapping_file}")
        chembl_mapping = pd.read_csv(mapping_file, sep='\t')
        
        metabolites_df = metabolites_df.merge(
            chembl_mapping[['HMDB', 'ChEMBL_id']],
            on='HMDB',
            how='left',
            suffixes=('', '_mapped')
        )
        
        # Use mapped values if available
        if 'ChEMBL_id_mapped' in metabolites_df.columns:
            if 'ChEMBL_id' not in metabolites_df.columns:
                metabolites_df['ChEMBL_id'] = metabolites_df['ChEMBL_id_mapped']
            else:
                metabolites_df['ChEMBL_id'] = metabolites_df['ChEMBL_id_mapped'].fillna(
                    metabolites_df['ChEMBL_id']
                )
            metabolites_df.drop('ChEMBL_id_mapped', axis=1, inplace=True)
        
        mapped_count = (metabolites_df['ChEMBL_id'] != 'none').sum() if 'ChEMBL_id' in metabolites_df.columns else 0
        logger.info(f"Loaded {mapped_count} existing ChEMBL mappings")
    else:
        logger.info("No existing ChEMBL mapping found, will create new mapping")
        metabolites_df['ChEMBL_id'] = 'none'
    
    # Identify metabolites that need mapping
    needs_mapping = metabolites_df['ChEMBL_id'].isna() | (metabolites_df['ChEMBL_id'] == 'none')
    to_map = metabolites_df[needs_mapping & metabolites_df['Canonical_smiles'].notna()].copy()
    
    logger.info(f"\nChEMBL Mapping Status:")
    logger.info(f"  Total metabolites: {len(metabolites_df)}")
    logger.info(f"  Already mapped: {(~needs_mapping).sum()}")
    logger.info(f"  Need mapping: {len(to_map)}")
    
    if len(to_map) > 0:
        logger.info(f"\nMapping {len(to_map)} SMILES to ChEMBL IDs...")
        logger.info("This is a one-time operation that will speed up both chEMBL and LINCS queries")
        logger.info("="*80)
        
        # Use multithreading for faster mapping
        chembl_ids = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(get_chembl_id_from_smiles, row['Canonical_smiles']): idx
                for idx, row in to_map.iterrows()
            }
            
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    chembl_id = future.result()
                    chembl_ids.append((idx, chembl_id))
                except Exception as e:
                    logger.error(f"Failed to map {idx}: {e}")
                    chembl_ids.append((idx, 'none'))
        
        # Update dataframe
        for idx, chembl_id in chembl_ids:
            metabolites_df.loc[idx, 'ChEMBL_id'] = chembl_id
        
        # Save mapping for future use
        mapping_df = metabolites_df[['HMDB', 'Canonical_smiles', 'ChEMBL_id']].copy()
        mapping_df.to_csv(mapping_file, index=False, sep='\t')
        
        mapped_count = (metabolites_df['ChEMBL_id'] != 'none').sum()
        logger.info(f"\nChEMBL mapping complete!")
        logger.info(f"  Successfully mapped: {mapped_count} metabolites")
        logger.info(f"  Mapping saved to: {mapping_file}")
    else:
        logger.info("\nAll metabolites already mapped, skipping ChEMBL API calls")
    
    logger.info("="*80)
    
    return metabolites_df


def preprocess_metabolites(metabolites: List[Dict]) -> pd.DataFrame:
    """
    Preprocess metabolite data for retriever pipeline.
    
    Steps:
    1. Convert to DataFrame
    2. Add canonical SMILES
    3. Map ChEMBL IDs
    
    Parameters:
    -----------
    metabolites : list of dict
        Raw metabolite records from HMDB
    
    Returns:
    --------
    pd.DataFrame
        Preprocessed metabolite data ready for retrievers
    """
    logger = logging.getLogger('preprocessing')
    logger.info("Preprocessing metabolites...")
    
    # Convert to DataFrame
    df = pd.DataFrame(metabolites)
    logger.info(f"Loaded {len(df)} metabolites")
    
    # Ensure HMDB column is named consistently
    if 'HMDB_ID' in df.columns and 'HMDB' not in df.columns:
        df['HMDB'] = df['HMDB_ID']
    
    # Add canonical SMILES
    df = add_canonical_smiles(df)
    
    # Map ChEMBL IDs
    df = map_chembl_ids(df)
    
    logger.info("Preprocessing complete!")
    
    return df


if __name__ == '__main__':
    # Test preprocessing
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Test with first 100
    test_metabolites = metabolites[:100]
    
    df = preprocess_metabolites(test_metabolites)
    print(f"\nPreprocessed {len(df)} metabolites")
    print(df[['HMDB', 'name', 'Canonical_smiles', 'ChEMBL_id']].head())
