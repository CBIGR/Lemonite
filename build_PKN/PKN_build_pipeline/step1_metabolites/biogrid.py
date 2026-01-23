"""
BioGRID chemical-protein interaction retriever.

Processes local BioGRID CHEMICALS file to extract metabolite-gene interactions.
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


class BioGRIDRetriever(LocalFileRetriever):
    """
    Retrieve metabolite-gene interactions from BioGRID CHEMICALS database.
    
    BioGRID provides chemical-protein interactions with InChIKey identifiers.
    This retriever filters to metabolites of interest and extracts gene symbols.
    """
    
    def __init__(self, cache_file=None):
        """Initialize BioGRID retriever."""
        super().__init__(
            db_name='BioGRID',
            file_path=config.BIOGRID_LOCATION,
            cache_file=cache_file or config.DB_OUTPUT_FILES['BioGRID']
        )
    
    def parse_file(self) -> pd.DataFrame:
        """
        Parse BioGRID CHEMICALS file.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'Source', 'InChIKey', 'BioGRID_ID', 'Chemical_Name']
        """
        self.logger.info(f"Loading BioGRID CHEMICALS from {self.file_path}")
        
        # Load file
        biogrid = pd.read_csv(self.file_path, sep='\t')
        self.logger.info(f"Loaded {len(biogrid)} total interactions from BioGRID")
        
        # Extract relevant columns
        # Expected columns: 'InChIKey', 'Official Symbol', 'BioGRID Chemical ID', 'Chemical Name'
        interactions = []
        
        for _, row in biogrid.iterrows():
            inchikey = row.get('InChIKey')
            gene = row.get('Official Symbol')
            chem_id = row.get('BioGRID Chemical ID')
            chem_name = row.get('Chemical Name')
            
            if pd.notna(inchikey) and pd.notna(gene):
                interactions.append({
                    'InChIKey': inchikey,
                    'Gene': gene,
                    'Source': 'BioGRID',
                    'BioGRID_ID': chem_id,
                    'Chemical_Name': chem_name
                })
        
        df = pd.DataFrame(interactions)
        
        # Remove duplicates
        df = df.drop_duplicates(subset=['InChIKey', 'Gene'])
        
        self.logger.info(f"Parsed {len(df)} unique interactions from BioGRID")
        
        return df
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get BioGRID interactions for metabolites.
        
        Parameters:
        -----------
        metabolites : list of dict
            Metabolite records with 'HMDB_ID' and 'InChIKey'
        
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
        
        # Create mapping of InChIKey to HMDB_ID
        inchikey_to_hmdb = {
            m['InChIKey']: m['HMDB_ID'] 
            for m in metabolites 
            if 'InChIKey' in m and 'HMDB_ID' in m
        }
        
        # Filter to metabolites of interest
        df = df[df['InChIKey'].isin(inchikey_to_hmdb.keys())]
        
        # Add HMDB_ID
        df['HMDB_ID'] = df['InChIKey'].map(inchikey_to_hmdb)
        
        # Select final columns
        result = df[['HMDB_ID', 'Gene', 'Source']].copy()
        
        self.logger.info(f"Found {len(result)} BioGRID interactions for {result['HMDB_ID'].nunique()} metabolites")
        
        # Save cache
        self.save_cache(result)
        
        return result
    
    def generate_urls(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Generate BioGRID URLs for metabolites.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'URL', 'BioGRID_ID', 'Chemical_Name']
        """
        # Parse file to get chemical IDs and names
        df = self.parse_file()
        
        # Create mapping of InChIKey to HMDB_ID
        inchikey_to_hmdb = {
            m['InChIKey']: m['HMDB_ID'] 
            for m in metabolites 
            if 'InChIKey' in m and 'HMDB_ID' in m
        }
        
        # Filter to metabolites of interest
        df = df[df['InChIKey'].isin(inchikey_to_hmdb.keys())]
        df['HMDB_ID'] = df['InChIKey'].map(inchikey_to_hmdb)
        
        # Generate URLs
        df['URL'] = df.apply(
            lambda row: f"https://thebiogrid.org/chemical/{row['BioGRID_ID']}/{row['Chemical_Name']}.html"
            if pd.notna(row['BioGRID_ID']) and pd.notna(row['Chemical_Name'])
            else '',
            axis=1
        )
        
        result = df[['HMDB_ID', 'Gene', 'URL', 'BioGRID_ID', 'Chemical_Name']].copy()
        
        return result


if __name__ == '__main__':
    # Test the retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Test BioGRID retriever
    retriever = BioGRIDRetriever()
    interactions = retriever.get_interactions(metabolites)
    print(f"\nBioGRID interactions: {len(interactions)}")
    print(interactions.head())
    
    # Test URL generation
    urls = retriever.generate_urls(metabolites)
    print(f"\nBioGRID URLs generated: {len(urls)}")
    print(urls.head())
