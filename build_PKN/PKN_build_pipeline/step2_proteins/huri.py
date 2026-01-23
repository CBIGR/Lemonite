"""
HuRI (Human Reference Interactome) retriever for protein-protein interactions.

Parses local HuRI database file and maps Ensembl IDs to gene symbols.
"""

import pandas as pd
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config
from utils.pipeline import LocalFileRetriever


class HuRIRetriever(LocalFileRetriever):
    """
    Retrieve protein-protein interactions from HuRI database.
    
    HuRI provides high-quality binary protein-protein interactions
    from systematic yeast two-hybrid screens.
    
    Database: http://www.interactome-atlas.org/
    """
    
    def __init__(self):
        """Initialize HuRI retriever."""
        super().__init__()
        self.logger = logging.getLogger('huri')
        self.file_path = config.HURI_FILE
        self.ensembl_mapping_file = config.ENSEMBL_MAPPING_FILE
        self.ensembl_to_symbol = None
    
    def _load_ensembl_mapping(self) -> Dict[str, str]:
        """
        Load Ensembl ID to gene symbol mapping.
        
        Returns:
        --------
        dict
            Mapping of Ensembl ID → gene symbol
        """
        if self.ensembl_to_symbol is not None:
            return self.ensembl_to_symbol
        
        self.logger.info(f"Loading Ensembl mapping from: {self.ensembl_mapping_file}")
        
        # Read mapping file
        mapping = pd.read_csv(self.ensembl_mapping_file, sep='\t')
        
        # Create dictionary
        self.ensembl_to_symbol = dict(zip(
            mapping['ensembl_gene_id'], 
            mapping['hgnc_symbol']
        ))
        
        self.logger.info(f"Loaded {len(self.ensembl_to_symbol)} Ensembl → symbol mappings")
        
        return self.ensembl_to_symbol
    
    def parse_file(self, genes: List[str]) -> pd.DataFrame:
        """
        Parse HuRI file and extract PPIs involving input genes.
        
        Parameters:
        -----------
        genes : list
            List of gene symbols to filter for
        
        Returns:
        --------
        pd.DataFrame
            PPIs with columns: ['GeneA', 'GeneB', 'Source']
        """
        self.logger.info(f"Reading HuRI file: {self.file_path}")
        
        # Read HuRI file (Ensembl IDs, no header)
        huri = pd.read_csv(self.file_path, sep='\t', header=None, names=['GeneA', 'GeneB'])
        
        self.logger.info(f"Loaded {len(huri)} total interactions")
        
        # Load Ensembl mapping
        ensembl_to_symbol = self._load_ensembl_mapping()
        
        # Create symbol → Ensembl mapping for input genes
        symbol_to_ensembl = {v: k for k, v in ensembl_to_symbol.items() if v in genes}
        genes_ensembl = list(symbol_to_ensembl.values())
        
        self.logger.info(f"Mapped {len(genes_ensembl)} input genes to Ensembl IDs")
        
        # Filter for interactions involving input genes (in Ensembl ID space)
        mask = (
            huri['GeneA'].isin(genes_ensembl) | 
            huri['GeneB'].isin(genes_ensembl)
        )
        huri_filtered = huri[mask].copy()
        
        self.logger.info(f"Filtered to {len(huri_filtered)} interactions involving input genes")
        
        # Map Ensembl IDs back to gene symbols
        failed_count = 0
        mapped_rows = []
        
        for _, row in huri_filtered.iterrows():
            ensembl_a = row['GeneA']
            ensembl_b = row['GeneB']
            
            symbol_a = ensembl_to_symbol.get(ensembl_a)
            symbol_b = ensembl_to_symbol.get(ensembl_b)
            
            if symbol_a and symbol_b:
                mapped_rows.append({
                    'GeneA': symbol_a,
                    'GeneB': symbol_b
                })
            else:
                failed_count += 1
        
        if failed_count > 0:
            self.logger.warning(f"Failed to map {failed_count} Ensembl IDs to gene symbols")
        
        # Create DataFrame
        huri_mapped = pd.DataFrame(mapped_rows)
        
        if len(huri_mapped) == 0:
            self.logger.warning("No interactions after mapping!")
            return pd.DataFrame(columns=['GeneA', 'GeneB', 'Source'])
        
        # Add source column
        huri_mapped['Source'] = 'HuRI'
        
        # Remove duplicates
        huri_mapped = huri_mapped.drop_duplicates(subset=['GeneA', 'GeneB'])
        
        # Remove rows with NaN
        huri_mapped = huri_mapped.dropna()
        
        self.logger.info(f"Final: {len(huri_mapped)} interactions after mapping and deduplication")
        
        return huri_mapped
    
    def get_interactions(self, genes: List[str]) -> pd.DataFrame:
        """
        Get protein-protein interactions from HuRI.
        
        Parameters:
        -----------
        genes : list
            List of gene symbols
        
        Returns:
        --------
        pd.DataFrame
            PPIs with columns: ['GeneA', 'GeneB', 'Source']
        """
        self.logger.info(f"Querying HuRI for {len(genes)} genes...")
        
        # Parse file
        interactions = self.parse_file(genes)
        
        if len(interactions) == 0:
            self.logger.warning("No interactions found!")
            return pd.DataFrame(columns=['GeneA', 'GeneB', 'Source'])
        
        # Save to file
        output_file = config.DB_OUTPUT_FILES.get('HuRI',
                                                 os.path.join(config.OUTPUT_DIR, 'HuRI_interactions.csv'))
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        interactions.to_csv(output_file, index=False)
        self.logger.info(f"Saved to: {output_file}")
        
        return interactions


if __name__ == '__main__':
    # Test HuRI retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # Test with a small set of genes
    test_genes = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS']
    
    retriever = HuRIRetriever()
    interactions = retriever.get_interactions(test_genes)
    
    print(f"\nTest results:")
    print(f"  Total interactions: {len(interactions)}")
    print(f"  Unique genes: {len(set(interactions['GeneA']) | set(interactions['GeneB']))}")
    print(f"\nFirst 10 interactions:")
    print(interactions.head(10))
