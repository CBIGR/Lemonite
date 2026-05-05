"""
STRING API retriever for protein-protein interactions.

Queries the STRING database API to retrieve PPI interactions for input genes.
Uses chunking to handle large gene lists efficiently.
"""

import pandas as pd
import requests
import time
import logging
from typing import List, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config
from utils.api_retry import retry_api_call


class STRINGRetriever:
    """
    Retrieve protein-protein interactions from STRING database via API.
    
    STRING provides known and predicted protein-protein interactions for various organisms.
    
    API: https://string-db.org/
    """
    
    def __init__(self, confidence_threshold=400):
        """
        Initialize STRING retriever.
        
        Parameters:
        -----------
        confidence_threshold : int
            Minimum combined confidence score (0-1000). Default 400 = medium confidence.
            STRING scores: low (150), medium (400), high (700), highest (900)
        """
        self.logger = logging.getLogger('string')
        self.api_url = "https://version-11-5.string-db.org/api"
        self.confidence_threshold = confidence_threshold
        self.species = 9606  # Homo sapiens
        self.caller_identity = "LemonIte"
    
    @retry_api_call('STRING')
    def _get_string_ids(self, gene_symbols: List[str]) -> Dict[str, str]:
        """
        Map gene symbols to STRING protein IDs.
        
        Parameters:
        -----------
        gene_symbols : list
            List of gene symbols to map
        
        Returns:
        --------
        dict
            Mapping of gene symbol → STRING ID
        """
        method = "get_string_ids"
        request_url = f"{self.api_url}/tsv/{method}"
        
        params = {
            "identifiers": "\r".join(gene_symbols),
            "species": self.species,
            "caller_identity": self.caller_identity
        }
        
        response = requests.post(request_url, data=params)
        response.raise_for_status()
        
        # Parse TSV response
        mapping = {}
        for line in response.text.strip().split("\n")[1:]:  # Skip header
            parts = line.split("\t")
            if len(parts) >= 3:
                query_symbol = parts[0]
                string_id = parts[2]
                mapping[query_symbol] = string_id
        
        return mapping
    
    @retry_api_call('STRING')
    def _get_interactions_for_string_ids(self, string_ids: List[str]) -> pd.DataFrame:
        """
        Get interaction partners for a list of STRING IDs.
        
        Parameters:
        -----------
        string_ids : list
            List of STRING protein IDs
        
        Returns:
        --------
        pd.DataFrame
            Interactions with columns: ['query_id', 'query_name', 'partner_name', 'score']
        """
        method = "interaction_partners"
        request_url = f"{self.api_url}/tsv/{method}"
        
        params = {
            "identifiers": "%0d".join(string_ids),
            "species": self.species,
            "caller_identity": self.caller_identity
        }
        
        response = requests.post(request_url, data=params)
        response.raise_for_status()
        
        # Parse response
        interactions = []
        for line in response.text.strip().split("\n")[1:]:  # Skip header
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                query_id = parts[0]
                query_name = parts[2]
                partner_name = parts[3]
                combined_score = int(parts[5])
                
                # Filter by confidence threshold
                if combined_score >= self.confidence_threshold:
                    interactions.append({
                        'query_id': query_id,
                        'query_name': query_name,
                        'partner_name': partner_name,
                        'score': combined_score
                    })
        
        return pd.DataFrame(interactions)
    
    def _process_gene_chunk(self, gene_chunk: List[str]) -> pd.DataFrame:
        """
        Process a chunk of genes to retrieve interactions.
        
        Parameters:
        -----------
        gene_chunk : list
            List of gene symbols to process
        
        Returns:
        --------
        pd.DataFrame
            Interactions for this chunk
        """
        try:
            # Get STRING IDs for genes in this chunk
            string_id_mapping = self._get_string_ids(gene_chunk)
            
            if not string_id_mapping:
                self.logger.warning(f"No STRING IDs found for chunk of {len(gene_chunk)} genes")
                return pd.DataFrame()
            
            # Get interactions
            string_ids = list(string_id_mapping.values())
            interactions = self._get_interactions_for_string_ids(string_ids)
            
            # Sleep to respect rate limits
            time.sleep(1)
            
            return interactions
        
        except Exception as e:
            self.logger.error(f"Failed to process gene chunk: {e}")
            return pd.DataFrame()
    
    def get_interactions(self, genes: List[str], chunk_size: int = 1000, 
                        max_workers: int = 10) -> pd.DataFrame:
        """
        Get protein-protein interactions from STRING for input genes.
        
        Parameters:
        -----------
        genes : list
            List of gene symbols
        chunk_size : int
            Number of genes per API request (default: 1000)
        max_workers : int
            Number of parallel workers (default: 10)
        
        Returns:
        --------
        pd.DataFrame
            Interactions with columns: ['GeneA', 'GeneB', 'combined_score', 'Source']
        """
        self.logger.info(f"Querying STRING for {len(genes)} genes...")
        self.logger.info(f"Using chunk size: {chunk_size}, max workers: {max_workers}")
        self.logger.info(f"Confidence threshold: {self.confidence_threshold}")
        
        # Split genes into chunks
        gene_chunks = [genes[i:i + chunk_size] for i in range(0, len(genes), chunk_size)]
        self.logger.info(f"Processing {len(gene_chunks)} chunks...")
        
        # Process chunks in parallel
        results = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(self._process_gene_chunk, chunk): idx 
                for idx, chunk in enumerate(gene_chunks)
            }
            
            for future in as_completed(futures):
                chunk_idx = futures[future]
                try:
                    result = future.result()
                    self.logger.info(f"✓ Chunk {chunk_idx + 1}/{len(gene_chunks)}: "
                                   f"{len(result)} interactions")
                    results.append(result)
                except Exception as e:
                    self.logger.error(f"✗ Chunk {chunk_idx + 1}/{len(gene_chunks)}: {e}")
        
        # Combine results
        if not results or all(len(df) == 0 for df in results):
            self.logger.warning("No interactions found!")
            return pd.DataFrame(columns=['GeneA', 'GeneB', 'combined_score', 'Source'])
        
        combined = pd.concat(results, ignore_index=True)
        
        # Standardize columns
        combined = combined.rename(columns={
            'query_name': 'GeneA',
            'partner_name': 'GeneB'
        })
        combined['Source'] = 'STRING'
        combined = combined[['GeneA', 'GeneB', 'combined_score', 'Source']]
        
        # Remove duplicates
        combined = combined.drop_duplicates(subset=['GeneA', 'GeneB'])
        
        self.logger.info(f"Retrieved {len(combined)} unique interactions from STRING")
        
        # Save to file
        output_file = config.DB_OUTPUT_FILES.get('STRING', 
                                                 os.path.join(config.OUTPUT_DIR, 'STRING_interactions.csv'))
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        combined.to_csv(output_file, index=False)
        self.logger.info(f"Saved to: {output_file}")
        
        return combined


if __name__ == '__main__':
    # Test STRING retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # Test with a small set of genes
    test_genes = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS']
    
    retriever = STRINGRetriever(confidence_threshold=400)
    interactions = retriever.get_interactions(test_genes, chunk_size=100, max_workers=2)
    
    print(f"\nTest results:")
    print(f"  Total interactions: {len(interactions)}")
    print(f"  Unique genes: {len(set(interactions['GeneA']) | set(interactions['GeneB']))}")
    print(f"\nFirst 10 interactions:")
    print(interactions.head(10))
