"""
LINCS biochemical binding data retriever.

Processes local LINCS files for metabolite-gene interactions based on IC50 values.
Uses ChEMBL IDs for compound mapping.
"""

import pandas as pd
import numpy as np
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pipeline import LocalFileRetriever
import config


class LINCSRetriever(LocalFileRetriever):
    """
    Retrieve metabolite-gene interactions from LINCS biochemical binding data.
    
    LINCS provides IC50 values for compound-target interactions.
    Requires pre-computed ChEMBL IDs for compound matching.
    """
    
    def __init__(self, ic50_threshold=10000, filter_human=True, cache_file=None):
        """
        Initialize LINCS retriever.
        
        Parameters:
        -----------
        ic50_threshold : float
            IC50 threshold in nM (default: 10000 = 10 μM)
        filter_human : bool
            Filter for human targets only
        """
        super().__init__(
            db_name='LINCS',
            file_path=config.LINCS_BIOCHEM_AGG,
            cache_file=cache_file or config.DB_OUTPUT_FILES['LINCS']
        )
        self.ic50_threshold = ic50_threshold
        self.filter_human = filter_human
        self._initialized = False
        self.chembl_to_lspci = {}
        self.lspci_to_genes = {}
    
    def _initialize_lookups(self):
        """
        Initialize LINCS lookup tables (one-time operation).
        
        Builds mappings:
        - ChEMBL ID → LINCS compound ID (lspci_id)
        - LINCS compound ID → gene symbols with IC50 values
        """
        if self._initialized:
            return
        
        self.logger.info("="*80)
        self.logger.info("Initializing LINCS lookup tables (one-time operation)...")
        self.logger.info("="*80)
        
        # Load compound mapping: ChEMBL → lspci_id
        self.logger.info(f"Loading LINCS compound mapping from {config.LINCS_COMPOUND_MAPPING}")
        lsp_compound_mapping = pd.read_csv(config.LINCS_COMPOUND_MAPPING)
        
        chembl_compounds = lsp_compound_mapping[lsp_compound_mapping['source'] == 'chembl']
        
        for _, row in chembl_compounds.iterrows():
            chembl_id = row['external_id']
            lspci_id = row['lspci_id']
            
            if chembl_id not in self.chembl_to_lspci:
                self.chembl_to_lspci[chembl_id] = []
            self.chembl_to_lspci[chembl_id].append(lspci_id)
        
        self.logger.info(f"✓ Indexed {len(self.chembl_to_lspci)} ChEMBL → LINCS compound mappings")
        
        # Load biochemical data
        self.logger.info(f"Loading LINCS biochemical data from {self.file_path}")
        biochem_data = pd.read_csv(
            self.file_path,
            usecols=['lspci_id', 'lspci_target_id', 'value', 'value_unit', 'symbol']
        )
        
        # Convert IC50 to nM (standardize units)
        biochem_data.loc[biochem_data['value_unit'] == 'uM', 'value'] *= 1000
        self.logger.info(f"✓ Loaded {len(biochem_data)} biochemical interactions")
        
        # Filter for human targets if requested
        if self.filter_human:
            self.logger.info(f"Loading target dictionary from {config.LINCS_TARGET_DICTIONARY}")
            lsp_target_dictionary = pd.read_csv(config.LINCS_TARGET_DICTIONARY)
            
            human_targets = set(
                lsp_target_dictionary[lsp_target_dictionary['tax_id'] == 9606]['lspci_target_id']
            )
            
            biochem_data = biochem_data[biochem_data['lspci_target_id'].isin(human_targets)]
            self.logger.info(f"✓ Filtered to {len(biochem_data)} human target interactions")
        
        # Build lspci_id → genes mapping
        for _, row in biochem_data.iterrows():
            lspci_id = row['lspci_id']
            
            if lspci_id not in self.lspci_to_genes:
                self.lspci_to_genes[lspci_id] = []
            
            self.lspci_to_genes[lspci_id].append({
                'symbol': row['symbol'],
                'value': row['value']
            })
        
        self.logger.info(f"✓ Indexed {len(self.lspci_to_genes)} LINCS compounds with interactions")
        self.logger.info("="*80)
        self.logger.info("LINCS initialization complete!")
        self.logger.info("="*80)
        
        self._initialized = True
    
    def get_genes_for_chembl_id(self, chembl_id: str) -> List[str]:
        """
        Get gene symbols for a ChEMBL ID.
        
        Parameters:
        -----------
        chembl_id : str
            ChEMBL compound ID
        
        Returns:
        --------
        list of str
            Gene symbols passing IC50 threshold
        """
        if pd.isna(chembl_id) or chembl_id == 'none' or chembl_id == '':
            return []
        
        # Ensure lookups are initialized
        self._initialize_lookups()
        
        # Fast lookup: ChEMBL → LINCS IDs
        lspci_ids = self.chembl_to_lspci.get(chembl_id, [])
        if not lspci_ids:
            return []
        
        # Fast lookup: LINCS IDs → genes (filtered by IC50)
        genes = set()
        for lspci_id in lspci_ids:
            interactions = self.lspci_to_genes.get(lspci_id, [])
            for interaction in interactions:
                if (interaction['value'] <= self.ic50_threshold and 
                    pd.notna(interaction['symbol'])):
                    genes.add(interaction['symbol'])
        
        return list(genes)
    
    def parse_file(self) -> pd.DataFrame:
        """
        Parse LINCS files (delegated to get_interactions since ChEMBL IDs needed).
        """
        raise NotImplementedError(
            "LINCS requires ChEMBL IDs from metabolites. "
            "Use get_interactions() instead of parse_file()."
        )
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get LINCS interactions for metabolites.
        
        Parameters:
        -----------
        metabolites : list of dict
            Metabolite records with 'HMDB_ID' and 'ChEMBL_id'
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'Source']
        """
        # Try cache first
        cached = self.load_cache()
        if cached is not None:
            return cached
        
        # Initialize lookup tables
        self._initialize_lookups()
        
        # Process metabolites
        self.logger.info(f"Processing {len(metabolites)} metabolites for LINCS interactions...")
        
        interactions = []
        processed_count = 0
        
        for metabolite in metabolites:
            hmdb_id = metabolite.get('HMDB_ID')
            chembl_id = metabolite.get('ChEMBL_id')
            
            if not hmdb_id or not chembl_id or chembl_id == 'none':
                continue
            
            genes = self.get_genes_for_chembl_id(chembl_id)
            
            for gene in genes:
                interactions.append({
                    'HMDB_ID': hmdb_id,
                    'Gene': gene,
                    'Source': 'LINCS'
                })
            
            processed_count += 1
            
            if processed_count % 1000 == 0:
                self.logger.info(f"Processed {processed_count}/{len(metabolites)} metabolites...")
        
        df = pd.DataFrame(interactions)
        df = df.drop_duplicates(subset=['HMDB_ID', 'Gene'])
        
        self.logger.info(
            f"Found {len(df)} LINCS interactions for "
            f"{df['HMDB_ID'].nunique() if len(df) > 0 else 0} metabolites "
            f"(IC50 ≤ {self.ic50_threshold} nM)"
        )
        
        # Save cache
        self.save_cache(df)
        
        return df
    
    def generate_urls(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Generate LINCS URLs for metabolites.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'URL', 'ChEMBL_ID']
        """
        # Get interactions
        df = self.get_interactions(metabolites)
        
        # Add ChEMBL IDs
        chembl_mapping = {m['HMDB_ID']: m.get('ChEMBL_id') 
                         for m in metabolites if 'ChEMBL_id' in m}
        
        df['ChEMBL_ID'] = df['HMDB_ID'].map(chembl_mapping)
        
        # Generate URLs (LINCS drug repurposing hub)
        df['URL'] = df['ChEMBL_ID'].apply(
            lambda cid: f"https://clue.io/repurposing#query/{cid}" 
            if pd.notna(cid) and cid != 'none' 
            else ''
        )
        
        result = df[['HMDB_ID', 'Gene', 'URL', 'ChEMBL_ID']].copy()
        
        return result


if __name__ == '__main__':
    # Test the retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Note: LINCS requires ChEMBL IDs which need to be pre-computed
    # For testing, we'll just show the structure
    print("\nLINCS retriever requires ChEMBL IDs to be added to metabolites first")
    print("Run preprocessing.py to add ChEMBL mappings before using LINCS")
    
    # Test with empty ChEMBL IDs (will return empty results)
    retriever = LINCSRetriever(ic50_threshold=10000)
    print(f"Initialized LINCS retriever with IC50 threshold: {retriever.ic50_threshold} nM")
