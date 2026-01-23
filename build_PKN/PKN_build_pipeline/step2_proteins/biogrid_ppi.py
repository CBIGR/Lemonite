"""
BioGRID PPI retriever for protein-protein interactions.

Parses local BioGRID database file to extract human PPI interactions.
"""

import pandas as pd
import logging
from typing import List
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config
from utils.pipeline import LocalFileRetriever


class BioGRIDPPIRetriever(LocalFileRetriever):
    """
    Retrieve protein-protein interactions from BioGRID database.
    
    BioGRID provides manually curated PPI data from high-throughput studies.
    
    Database: https://thebiogrid.org/
    """
    
    def __init__(self):
        """Initialize BioGRID PPI retriever."""
        super().__init__()
        self.logger = logging.getLogger('biogrid_ppi')
        self.file_path = config.BIOGRID_PPI_FILE
    
    def parse_file(self, genes: List[str]) -> pd.DataFrame:
        """
        Parse BioGRID file and extract PPIs involving input genes.
        
        Parameters:
        -----------
        genes : list
            List of gene symbols to filter for
        
        Returns:
        --------
        pd.DataFrame
            PPIs with columns: ['GeneA', 'GeneB', 'Source']
        """
        self.logger.info(f"Reading BioGRID PPI file: {self.file_path}")
        
        # Read BioGRID file
        # BioGRID is tab-separated with many columns
        biogrid = pd.read_csv(self.file_path, sep='\t', low_memory=False)
        
        self.logger.info(f"Loaded {len(biogrid)} total interactions")
        
        # Filter for human interactions only
        human_mask = (
            (biogrid['Organism Name Interactor A'] == 'Homo sapiens') | 
            (biogrid['Organism Name Interactor B'] == 'Homo sapiens')
        )
        biogrid = biogrid[human_mask]
        self.logger.info(f"Filtered to {len(biogrid)} human interactions")
        
        # Extract gene symbols
        biogrid = biogrid[[
            'Official Symbol Interactor A', 
            'Official Symbol Interactor B'
        ]].copy()
        
        # Rename columns
        biogrid.columns = ['GeneA', 'GeneB']
        
        # Filter for interactions involving input genes
        # Include interactions where either partner is in the input gene list
        mask = (
            biogrid['GeneA'].isin(genes) | 
            biogrid['GeneB'].isin(genes)
        )
        biogrid_filtered = biogrid[mask].copy()
        
        self.logger.info(f"Filtered to {len(biogrid_filtered)} interactions involving input genes")
        
        # Add source column
        biogrid_filtered['Source'] = 'BioGRID_PPI'
        
        # Remove duplicates
        biogrid_filtered = biogrid_filtered.drop_duplicates(subset=['GeneA', 'GeneB'])
        
        # Get unique genes
        all_genes = set(biogrid_filtered['GeneA']) | set(biogrid_filtered['GeneB'])
        input_genes_found = len(all_genes & set(genes))
        self.logger.info(f"Found {input_genes_found} of {len(genes)} input genes in network")
        
        return biogrid_filtered
    
    def get_interactions(self, genes: List[str]) -> pd.DataFrame:
        """
        Get protein-protein interactions from BioGRID.
        
        Parameters:
        -----------
        genes : list
            List of gene symbols
        
        Returns:
        --------
        pd.DataFrame
            PPIs with columns: ['GeneA', 'GeneB', 'Source']
        """
        self.logger.info(f"Querying BioGRID PPI for {len(genes)} genes...")
        
        # Parse file
        interactions = self.parse_file(genes)
        
        if len(interactions) == 0:
            self.logger.warning("No interactions found!")
            return pd.DataFrame(columns=['GeneA', 'GeneB', 'Source'])
        
        # Save to file
        output_file = config.DB_OUTPUT_FILES.get('BioGRID_PPI',
                                                 os.path.join(config.OUTPUT_DIR, 'BioGRID_PPI_interactions.csv'))
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        interactions.to_csv(output_file, index=False)
        self.logger.info(f"Saved to: {output_file}")
        
        return interactions


if __name__ == '__main__':
    # Test BioGRID PPI retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # Test with a small set of genes
    test_genes = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS']
    
    retriever = BioGRIDPPIRetriever()
    interactions = retriever.get_interactions(test_genes)
    
    print(f"\nTest results:")
    print(f"  Total interactions: {len(interactions)}")
    print(f"  Unique genes: {len(set(interactions['GeneA']) | set(interactions['GeneB']))}")
    print(f"\nFirst 10 interactions:")
    print(interactions.head(10))
