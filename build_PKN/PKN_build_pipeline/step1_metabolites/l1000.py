"""
L1000 gene expression signature retriever.

Processes L1000 GMT file (2.1GB) for gene expression signatures.
Single-threaded to avoid loading massive file multiple times.
"""

import pandas as pd
import numpy as np
import re
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pipeline import DatabaseRetriever
import config


class L1000Retriever(DatabaseRetriever):
    """
    Retrieve metabolite-gene interactions from L1000 gene expression data.
    
    L1000 provides gene expression signatures (up/down-regulated genes)
    for compounds. Uses InChIKey for compound matching.
    
    IMPORTANT: Single-threaded to avoid loading 2.1GB GMT file multiple times.
    """
    
    def __init__(self, direction='both', cache_file=None):
        """
        Initialize L1000 retriever.
        
        Parameters:
        -----------
        direction : str
            Gene regulation direction: 'up', 'down', or 'both'
        """
        super().__init__(
            db_name='L1000',
            cache_file=cache_file or config.DB_OUTPUT_FILES['L1000']
        )
        self.direction = direction
        self._initialized = False
        self.inchikey_to_compounds = {}
        self.compound_to_genes = {}
    
    def _initialize_lookups(self):
        """
        Initialize L1000 lookup tables (one-time operation).
        
        Parses:
        1. Compound metadata (InChIKey → compound name mapping)
        2. GMT file (compound name → gene signatures)
        
        This takes several minutes due to the 2.1GB file size.
        """
        if self._initialized:
            return
        
        self.logger.info("="*80)
        self.logger.info("Initializing L1000 lookup tables (one-time operation)...")
        self.logger.info("This may take a few minutes - loading 2.1GB GMT file...")
        self.logger.info("="*80)
        
        # Load compound metadata
        self.logger.info(f"Loading L1000 compound metadata from {config.L1000_COMPOUNDS_LOCATION}")
        compounds_df = pd.read_csv(config.L1000_COMPOUNDS_LOCATION, sep='\t')
        
        # Build InChIKey → compound name(s) mapping
        for _, row in compounds_df.iterrows():
            inchi = row.get('inchi_key')
            compound_name = row.get('pert_name')
            
            if pd.notna(inchi) and pd.notna(compound_name):
                if inchi not in self.inchikey_to_compounds:
                    self.inchikey_to_compounds[inchi] = set()
                self.inchikey_to_compounds[inchi].add(compound_name.lower())
        
        self.logger.info(f"✓ Indexed {len(self.inchikey_to_compounds)} InChIKeys → compound names")
        
        # Parse GMT file
        self.logger.info(f"Parsing GMT file from {config.L1000_GMT_LOCATION}")
        self.logger.info("Format: 'ID_compound_concentration direction' followed by gene symbols")
        
        line_count = 0
        
        with open(config.L1000_GMT_LOCATION, 'r') as f:
            for line in f:
                line_count += 1
                
                if line_count % 100000 == 0:
                    self.logger.info(f"  Processed {line_count:,} signatures...")
                
                parts = line.strip().split('\t')
                
                if len(parts) < 3:
                    continue
                
                # First column: "ExperimentID_CellLine_..._compound_concentration direction"
                identifier = parts[0]
                genes = parts[2:]  # Remaining columns are gene symbols
                
                # Extract compound name and direction using regex
                # Pattern: anything_COMPOUNDNAME_NUMBER[uM/nM/etc] [up/down]
                match = re.search(r'_([a-zA-Z0-9\-]+)_[\d\.]+[a-zA-Z]+\s+(up|down)$', identifier)
                
                if not match:
                    continue
                
                compound_name = match.group(1).lower()
                sig_direction = match.group(2)
                
                # Create key: (compound_name, direction)
                key = (compound_name, sig_direction)
                
                if key not in self.compound_to_genes:
                    self.compound_to_genes[key] = []
                
                self.compound_to_genes[key].extend(genes)
        
        self.logger.info(f"✓ Parsed {line_count:,} total signatures")
        self.logger.info(f"✓ Indexed {len(self.compound_to_genes)} unique (compound, direction) combinations")
        self.logger.info("="*80)
        self.logger.info("L1000 initialization complete!")
        self.logger.info("="*80)
        
        self._initialized = True
    
    def get_genes_for_inchikey(self, inchikey: str) -> List[str]:
        """
        Get gene symbols for an InChIKey.
        
        Parameters:
        -----------
        inchikey : str
            InChIKey identifier
        
        Returns:
        --------
        list of str
            Gene symbols (up/down-regulated based on direction setting)
        """
        if pd.isna(inchikey) or inchikey == 'none' or inchikey == '':
            return []
        
        # Ensure lookups are initialized
        self._initialize_lookups()
        
        # Lookup: InChIKey → compound names
        compound_names = self.inchikey_to_compounds.get(inchikey, set())
        
        if not compound_names:
            return []
        
        # Lookup: compound names → genes (filtered by direction)
        all_genes = set()
        
        for compound_name in compound_names:
            if self.direction in ['up', 'both']:
                key_up = (compound_name, 'up')
                genes_up = self.compound_to_genes.get(key_up, [])
                all_genes.update(genes_up)
            
            if self.direction in ['down', 'both']:
                key_down = (compound_name, 'down')
                genes_down = self.compound_to_genes.get(key_down, [])
                all_genes.update(genes_down)
        
        return list(all_genes)
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get L1000 interactions for metabolites.
        
        SINGLE-THREADED to avoid loading 2.1GB file in multiple workers.
        Includes resume capability with progress saving.
        
        Parameters:
        -----------
        metabolites : list of dict
            Metabolite records with 'HMDB_ID' and 'InChIKey'
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'Source']
        """
        output_file = self.cache_file
        
        # Check if we can resume from existing file
        if os.path.exists(output_file):
            self.logger.info(f"Found existing file: {output_file}")
            processed_df = pd.read_csv(output_file, sep='\t')
            self.logger.info(f"Resuming from {len(processed_df)} rows")
            
            # Convert back to interaction format
            interactions = []
            for _, row in processed_df.iterrows():
                genes_str = row.get('L1000', '')
                if pd.notna(genes_str) and genes_str != 'none' and genes_str != '':
                    genes = genes_str.split('|')
                    for gene in genes:
                        interactions.append({
                            'HMDB_ID': row['HMDB'],
                            'Gene': gene,
                            'Source': 'L1000'
                        })
            
            result = pd.DataFrame(interactions)
            self.logger.info(f"Loaded {len(result)} L1000 interactions from cache")
            return result
        
        # Initialize lookups
        self._initialize_lookups()
        
        # Process metabolites one by one (single-threaded)
        self.logger.info(f"Processing {len(metabolites)} metabolites for L1000 interactions...")
        self.logger.info("Using SINGLE-THREADED processing to save RAM")
        
        interactions = []
        processed_count = 0
        interaction_count = 0
        
        # Create temporary storage for processed format
        temp_results = []
        
        for metabolite in metabolites:
            hmdb_id = metabolite.get('HMDB_ID')
            inchikey = metabolite.get('InChIKey')
            
            if pd.isna(inchikey):
                temp_results.append({
                    'HMDB': hmdb_id,
                    'InChIKey': inchikey,
                    'L1000': 'none'
                })
                continue
            
            # Get interactions
            genes = self.get_genes_for_inchikey(inchikey)
            
            if genes:
                temp_results.append({
                    'HMDB': hmdb_id,
                    'InChIKey': inchikey,
                    'L1000': '|'.join(genes)
                })
                
                for gene in genes:
                    interactions.append({
                        'HMDB_ID': hmdb_id,
                        'Gene': gene,
                        'Source': 'L1000'
                    })
                
                interaction_count += 1
            else:
                temp_results.append({
                    'HMDB': hmdb_id,
                    'InChIKey': inchikey,
                    'L1000': 'none'
                })
            
            processed_count += 1
            
            # Save progress every 1000 metabolites
            if processed_count % config.RESUME_SAVE_INTERVAL == 0:
                temp_df = pd.DataFrame(temp_results)
                temp_df.to_csv(output_file, index=False, sep='\t')
                self.logger.info(f"  Progress saved: {processed_count}/{len(metabolites)} metabolites")
        
        # Final save
        temp_df = pd.DataFrame(temp_results)
        temp_df.to_csv(output_file, index=False, sep='\t')
        
        result = pd.DataFrame(interactions)
        result = result.drop_duplicates(subset=['HMDB_ID', 'Gene'])
        
        self.logger.info(
            f"Found {len(result)} L1000 interactions for "
            f"{interaction_count} metabolites "
            f"(direction: {self.direction})"
        )
        
        return result


if __name__ == '__main__':
    # Test the retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Test with first 100 metabolites (GMT file is large)
    print("\nTesting L1000 retriever with first 100 metabolites...")
    print("Note: First run will take several minutes to parse 2.1GB GMT file")
    
    test_metabolites = metabolites[:100]
    
    retriever = L1000Retriever(direction='both')
    interactions = retriever.get_interactions(test_metabolites)
    print(f"\nL1000 interactions (test): {len(interactions)}")
    print(interactions.head())
