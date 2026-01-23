"""
File I/O utilities for the PKN pipeline.

Handles loading metabolites, saving results, managing progress/resume,
and data format conversions.
"""

import pandas as pd
import os
import json
import logging
from typing import List, Dict, Optional
import xml.etree.ElementTree as ET


def load_hmdb_metabolites(xml_path: str) -> List[Dict]:
    """
    Load metabolite data from HMDB XML file.
    
    Parameters:
    -----------
    xml_path : str
        Path to HMDB metabolites XML file
    
    Returns:
    --------
    list of dict
        List of metabolite records with keys:
        - HMDB_ID
        - name
        - InChIKey
        - InChI
        - SMILES
    """
    logger = logging.getLogger('file_io')
    logger.info(f"Loading HMDB metabolites from {xml_path}")
    
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # Namespace for HMDB XML
    ns = {'hmdb': 'http://www.hmdb.ca'}
    
    metabolites = []
    for metabolite in root.findall('hmdb:metabolite', ns):
        record = {}
        
        # HMDB ID
        accession = metabolite.find('hmdb:accession', ns)
        if accession is not None:
            record['HMDB_ID'] = accession.text
        
        # Name
        name = metabolite.find('hmdb:name', ns)
        if name is not None:
            record['name'] = name.text
        
        # InChIKey
        inchikey = metabolite.find('hmdb:inchikey', ns)
        if inchikey is not None:
            record['InChIKey'] = inchikey.text
        
        # InChI
        inchi = metabolite.find('hmdb:inchi', ns)
        if inchi is not None:
            record['InChI'] = inchi.text
        
        # SMILES
        smiles = metabolite.find('hmdb:smiles', ns)
        if smiles is not None:
            record['SMILES'] = smiles.text
        
        if 'HMDB_ID' in record:
            metabolites.append(record)
    
    logger.info(f"Loaded {len(metabolites)} metabolites from HMDB")
    return metabolites


def save_interactions(df: pd.DataFrame, output_file: str, add_source: bool = True):
    """
    Save interaction dataframe to CSV.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Interaction data
    output_file : str
        Output file path
    add_source : bool
        Whether to add 'Source' column if missing
    """
    if add_source and 'Source' not in df.columns:
        df['Source'] = os.path.basename(output_file).split('_')[0]
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, index=False)
    logging.info(f"Saved {len(df)} interactions to {output_file}")


def load_progress(progress_file: str) -> Dict:
    """
    Load processing progress from JSON file.
    
    Used for resuming long-running operations (e.g., L1000 processing).
    
    Returns:
    --------
    dict
        Progress data with keys like 'last_index', 'processed_count', etc.
    """
    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            return json.load(f)
    return {}


def save_progress(progress_file: str, progress_data: Dict):
    """
    Save processing progress to JSON file.
    
    Parameters:
    -----------
    progress_file : str
        Path to progress file
    progress_data : dict
        Progress information to save
    """
    os.makedirs(os.path.dirname(progress_file), exist_ok=True)
    with open(progress_file, 'w') as f:
        json.dump(progress_data, f, indent=2)
    logging.info(f"Saved progress to {progress_file}")


def combine_database_results(db_files: Dict[str, str]) -> pd.DataFrame:
    """
    Combine results from multiple database files into single dataframe.
    
    Parameters:
    -----------
    db_files : dict
        Dictionary mapping database names to file paths
    
    Returns:
    --------
    pd.DataFrame
        Combined dataframe with all interactions
    """
    logger = logging.getLogger('file_io')
    all_dfs = []
    
    for db_name, file_path in db_files.items():
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            if 'Source' not in df.columns:
                df['Source'] = db_name
            all_dfs.append(df)
            logger.info(f"Loaded {len(df)} interactions from {db_name}")
        else:
            logger.warning(f"File not found: {file_path}")
    
    if not all_dfs:
        logger.warning("No database files found to combine")
        return pd.DataFrame()
    
    combined = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Combined {len(combined)} total interactions from {len(all_dfs)} databases")
    
    return combined


def ensure_output_dir(file_path: str):
    """
    Ensure the directory for a file path exists.
    
    Parameters:
    -----------
    file_path : str
        Path to file (directory will be created if needed)
    """
    directory = os.path.dirname(file_path)
    if directory:
        os.makedirs(directory, exist_ok=True)
