"""
STITCH protein-chemical interaction retriever.

Processes local STITCH database file with MyGene.info for Ensembl→Gene symbol mapping.
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


class STITCHRetriever(LocalFileRetriever):
    """
    Retrieve metabolite-gene interactions from STITCH database.
    
    STITCH provides protein-chemical links using PubChem IDs.
    Requires MyGene.info API to map Ensembl protein IDs to gene symbols.
    """
    
    def __init__(self, cache_file=None):
        """Initialize STITCH retriever."""
        super().__init__(
            db_name='STITCH',
            file_path=config.STITCH_LOCATION,
            cache_file=cache_file or config.DB_OUTPUT_FILES['STITCH']
        )
        self.protein_to_symbol = None  # Cached mapping
    
    def _map_proteins_to_symbols(self, protein_ids: List[str]) -> Dict[str, str]:
        """
        Map Ensembl protein IDs to gene symbols using MyGene.info.
        
        Parameters:
        -----------
        protein_ids : list of str
            Ensembl protein IDs (e.g., 'ENSP00000000233')
        
        Returns:
        --------
        dict
            Mapping of protein_id → gene_symbol
        """
        import mygene
        
        self.logger.info(f"Mapping {len(protein_ids)} Ensembl protein IDs to gene symbols...")
        
        mg = mygene.MyGeneInfo()
        result = mg.querymany(
            protein_ids, 
            scopes='ensembl.protein', 
            fields='symbol', 
            species='human'
        )
        
        protein_to_symbol = {}
        for item in result:
            protein = item.get('query')
            symbol = item.get('symbol', np.nan)
            protein_to_symbol[protein] = symbol
        
        mapped_count = sum(1 for v in protein_to_symbol.values() if pd.notna(v))
        self.logger.info(f"Mapped {mapped_count}/{len(protein_ids)} proteins to gene symbols")
        
        return protein_to_symbol
    
    def parse_file(self) -> pd.DataFrame:
        """
        Parse STITCH protein-chemical links file.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['PubChem_ID', 'Gene', 'Source', 'Protein_ID']
        """
        self.logger.info(f"Loading STITCH from {self.file_path}")
        
        # Load file
        stitch = pd.read_csv(self.file_path, sep='\t')
        self.logger.info(f"Loaded {len(stitch)} total interactions from STITCH")
        
        # Clean chemical column (extract PubChem CID)
        stitch['chemical'] = stitch['chemical'].str.replace(r'\D', '', regex=True)
        
        # Clean protein column (remove '9606.' prefix for human)
        stitch['protein'] = stitch['protein'].str.replace('9606.', '', regex=False)
        
        self.logger.info(
            f"STITCH contains {len(stitch)} interactions: "
            f"{stitch['chemical'].nunique()} chemicals × "
            f"{stitch['protein'].nunique()} proteins"
        )
        
        # Get unique protein IDs and map to gene symbols
        unique_proteins = stitch['protein'].unique().tolist()
        self.protein_to_symbol = self._map_proteins_to_symbols(unique_proteins)
        
        # Map proteins to symbols
        stitch['symbol'] = stitch['protein'].map(self.protein_to_symbol)
        
        # Filter out unmapped proteins
        stitch_filtered = stitch[stitch['symbol'].notna()].copy()
        
        self.logger.info(
            f"After mapping: {len(stitch_filtered)} interactions with valid gene symbols"
        )
        
        # Prepare output format
        interactions = []
        for _, row in stitch_filtered.iterrows():
            interactions.append({
                'PubChem_ID': row['chemical'],
                'Gene': row['symbol'],
                'Source': 'STITCH',
                'Protein_ID': row['protein']
            })
        
        df = pd.DataFrame(interactions)
        df = df.drop_duplicates(subset=['PubChem_ID', 'Gene'])
        
        return df
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get STITCH interactions for metabolites.
        
        Parameters:
        -----------
        metabolites : list of dict
            Metabolite records with 'HMDB_ID' and 'PubChem'
        
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
        
        # Create mapping of PubChem ID to HMDB_ID
        pubchem_to_hmdb = {}
        for m in metabolites:
            if 'PubChem' in m and 'HMDB_ID' in m and pd.notna(m['PubChem']):
                # Convert to string, handle floats
                pubchem_id = str(int(float(m['PubChem'])))
                pubchem_to_hmdb[pubchem_id] = m['HMDB_ID']
        
        # Filter to metabolites of interest
        df = df[df['PubChem_ID'].isin(pubchem_to_hmdb.keys())]
        
        # Add HMDB_ID
        df['HMDB_ID'] = df['PubChem_ID'].map(pubchem_to_hmdb)
        
        # Select final columns
        result = df[['HMDB_ID', 'Gene', 'Source']].copy()
        
        self.logger.info(
            f"Found {len(result)} STITCH interactions for "
            f"{result['HMDB_ID'].nunique()} metabolites"
        )
        
        # Save cache
        self.save_cache(result)
        
        return result
    
    def generate_urls(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Generate STITCH URLs for metabolites.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'URL', 'PubChem_ID']
        """
        # Parse file to get PubChem IDs
        df = self.parse_file()
        
        # Create mapping
        pubchem_to_hmdb = {}
        for m in metabolites:
            if 'PubChem' in m and 'HMDB_ID' in m and pd.notna(m['PubChem']):
                pubchem_id = str(int(float(m['PubChem'])))
                pubchem_to_hmdb[pubchem_id] = m['HMDB_ID']
        
        # Filter to metabolites of interest
        df = df[df['PubChem_ID'].isin(pubchem_to_hmdb.keys())]
        df['HMDB_ID'] = df['PubChem_ID'].map(pubchem_to_hmdb)
        
        # Generate URLs
        df['URL'] = df['PubChem_ID'].apply(
            lambda cid: f"http://stitch.embl.de/cgi/show_network_section.pl?identifier=CID{cid}&species=9606"
        )
        
        result = df[['HMDB_ID', 'Gene', 'URL', 'PubChem_ID']].copy()
        
        return result


if __name__ == '__main__':
    # Test the retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Test STITCH retriever
    retriever = STITCHRetriever()
    interactions = retriever.get_interactions(metabolites)
    print(f"\nSTITCH interactions: {len(interactions)}")
    print(interactions.head())
    
    # Test URL generation
    urls = retriever.generate_urls(metabolites)
    print(f"\nSTITCH URLs generated: {len(urls)}")
    print(urls.head())
