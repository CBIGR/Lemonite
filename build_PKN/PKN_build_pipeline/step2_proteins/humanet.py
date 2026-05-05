"""
HumanNet retriever for protein-protein interactions.

Parses the local HumanNet HS-PI.symbol.tsv file and extracts interactions
involving genes from the metabolite-gene PKN.
"""

import pandas as pd
import logging
from typing import List
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config


class HumanNetRetriever:
    """
    Retrieve protein-protein interactions from HumanNet.
    """

    def __init__(self):
        self.logger = logging.getLogger('humanet')
        self.file_path = config.HUMANNET_LOCATION

        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"HumanNet file not found: {self.file_path}")

    def parse_file(self, genes: List[str]) -> pd.DataFrame:
        """
        Parse HumanNet file and extract PPIs involving input genes.
        """
        self.logger.info(f"Reading HumanNet file: {self.file_path}")

        humanet = pd.read_csv(self.file_path, sep='\t', header=None, usecols=[0, 1], names=['GeneA', 'GeneB'])
        self.logger.info(f"Loaded {len(humanet)} total interactions")

        # Filter for interactions involving any input gene
        mask = (
            humanet['GeneA'].isin(genes) |
            humanet['GeneB'].isin(genes)
        )
        filtered = humanet[mask].copy()
        self.logger.info(f"Filtered to {len(filtered)} interactions involving input genes")

        if len(filtered) == 0:
            return pd.DataFrame(columns=['GeneA', 'GeneB', 'Source'])

        filtered['Source'] = 'HumanNet'
        filtered = filtered.drop_duplicates(subset=['GeneA', 'GeneB'])
        filtered = filtered.dropna(subset=['GeneA', 'GeneB'])

        return filtered

    def get_interactions(self, genes: List[str]) -> pd.DataFrame:
        self.logger.info(f"Querying HumanNet for {len(genes)} genes...")

        interactions = self.parse_file(genes)

        if len(interactions) == 0:
            self.logger.warning("No HumanNet interactions found!")
            return pd.DataFrame(columns=['GeneA', 'GeneB', 'Source'])

        output_file = config.DB_OUTPUT_FILES.get('HumanNet', os.path.join(config.OUTPUT_DIR, 'HumanNet_interactions.csv'))
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        interactions.to_csv(output_file, index=False)
        self.logger.info(f"Saved to: {output_file}")

        return interactions


if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO)

    test_genes = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS']

    retriever = HumanNetRetriever()
    interactions = retriever.get_interactions(test_genes)

    print(f"\nTest results:")
    print(f"  Total interactions: {len(interactions)}")
    print(f"  Unique genes: {len(set(interactions['GeneA']) | set(interactions['GeneB']))}")
    print(f"\nFirst 10 interactions:")
    print(interactions.head(10))
