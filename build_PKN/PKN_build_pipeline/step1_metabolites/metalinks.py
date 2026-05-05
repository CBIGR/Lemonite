"""
MetalinksDB metabolite-gene interaction retriever.

Processes local MetalinksDB CSV file from liana+ package.
"""

import pandas as pd
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pipeline import LocalFileRetriever
import config


class MetalinksRetriever(LocalFileRetriever):
    """
    Retrieve metabolite-gene interactions from MetalinksDB.
    
    MetalinksDB is a curated database from the liana+ package containing
    metabolite-gene interactions with HMDB IDs and gene symbols.
    """
    
    def __init__(self, cache_file=None):
        """Initialize MetalinksDB retriever."""
        super().__init__(
            db_name='MetalinksDB',
            file_path=config.METALINKS_PATH,
            cache_file=cache_file or config.DB_OUTPUT_FILES['MetalinksDB']
        )
    
    def parse_file(self) -> pd.DataFrame:
        """
        Parse MetalinksDB CSV file.
        
        Expected columns: 'hmdb', 'gene_symbol'
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'Source']
        """
        self.logger.info(f"Loading MetalinksDB from {self.file_path}")
        
        # Load file
        metalinks = pd.read_csv(self.file_path)
        self.logger.info(f"Loaded {len(metalinks)} total interactions from MetalinksDB")
        
        # Check required columns
        if 'hmdb' not in metalinks.columns or 'gene_symbol' not in metalinks.columns:
            raise ValueError(
                f"MetalinksDB file missing required columns. "
                f"Expected: ['hmdb', 'gene_symbol'], "
                f"Found: {metalinks.columns.tolist()}"
            )
        
        # Prepare output format
        interactions = []
        for _, row in metalinks.iterrows():
            hmdb_id = row['hmdb']
            gene = row['gene_symbol']
            
            if pd.notna(hmdb_id) and pd.notna(gene):
                interactions.append({
                    'HMDB_ID': hmdb_id,
                    'Gene': gene,
                    'Source': 'MetalinksDB'
                })
        
        df = pd.DataFrame(interactions)
        
        # Remove duplicates
        df = df.drop_duplicates(subset=['HMDB_ID', 'Gene'])
        
        self.logger.info(
            f"Parsed {len(df)} unique interactions from MetalinksDB: "
            f"{df['HMDB_ID'].nunique()} metabolites × {df['Gene'].nunique()} genes"
        )
        
        return df
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get MetalinksDB interactions for metabolites.
        
        Parameters:
        -----------
        metabolites : list of dict
            Metabolite records with 'HMDB_ID'
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'Source']
        """
        # Try cache first
        cached = self.load_cache()
        if cached is not None:
            return cached
        
        # Parse file
        df = self.parse_file()
        
        # Filter to metabolites of interest
        metabolite_ids = {m.get('HMDB_ID') for m in metabolites if m.get('HMDB_ID')}
        df = df[df['HMDB_ID'].isin(metabolite_ids)]
        
        self.logger.info(
            f"Found {len(df)} MetalinksDB interactions for "
            f"{df['HMDB_ID'].nunique()} metabolites"
        )
        
        # Save cache
        self.save_cache(df)
        
        return df


if __name__ == '__main__':
    # Test the retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Test Metalinks retriever
    retriever = MetalinksRetriever()
    interactions = retriever.get_interactions(metabolites)
    print(f"\nMetalinksDB interactions: {len(interactions)}")
    print(f"Metabolites: {interactions['HMDB_ID'].nunique()}")
    print(f"Genes: {interactions['Gene'].nunique()}")
    print(interactions.head())
