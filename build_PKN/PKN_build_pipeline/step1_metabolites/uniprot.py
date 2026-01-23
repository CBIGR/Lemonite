"""
UniProt metabolite-protein interaction retriever.

Uses UniProtKB REST API to search for proteins by InChIKey.
"""

import pandas as pd
import requests
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pipeline import APIRetriever
from utils.api_retry import retry_api_call
import config


class UniProtRetriever(APIRetriever):
    """
    Retrieve metabolite-gene interactions from UniProtKB.
    
    UniProt provides protein information with InChIKey identifiers.
    Uses REST API with retry logic for fault tolerance.
    """
    
    def __init__(self, reviewed_only=False, cache_file=None):
        """
        Initialize UniProt retriever.
        
        Parameters:
        -----------
        reviewed_only : bool
            If True, only return reviewed (SwissProt) entries
        """
        super().__init__(
            db_name='UniProtKB',
            cache_file=cache_file or config.DB_OUTPUT_FILES['UniProtKB'],
            max_workers=config.API_RETRY_CONFIG['UniProtKB']['max_workers']
        )
        self.reviewed_only = reviewed_only
    
    @retry_api_call(db_name='UniProtKB')
    def fetch_single(self, metabolite: Dict) -> List[str]:
        """
        Fetch UniProt interactions for a single metabolite.
        
        Parameters:
        -----------
        metabolite : dict
            Must contain 'InChIKey'
        
        Returns:
        --------
        list of str
            Gene symbols
        """
        inchikey = metabolite.get('InChIKey')
        
        if pd.isna(inchikey):
            return []
        
        # Get timeout from config
        timeout = config.API_RETRY_CONFIG['UniProtKB']['timeout']
        
        # Build query URL
        if self.reviewed_only:
            url = (
                f'https://rest.uniprot.org/uniprotkb/search?'
                f'query=(inchikey:{inchikey})%20AND%20(organism_id:9606)%20AND%20(reviewed:true)'
                f'&format=tsv&fields=accession,gene_names'
            )
        else:
            url = (
                f'https://rest.uniprot.org/uniprotkb/search?'
                f'query=(inchikey:{inchikey})%20AND%20(organism_id:9606)'
                f'&format=tsv&fields=accession,gene_names'
            )
        
        response = requests.get(url, timeout=timeout)
        
        if response.status_code != 200:
            raise requests.exceptions.RequestException(f"HTTP {response.status_code}")
        
        # Parse TSV response
        data = response.text.strip().split('\n')
        if len(data) < 2:  # Header + at least one row
            return []
        
        genes = []
        for line in data[1:]:  # Skip header
            cols = line.split('\t')
            
            if len(cols) < 2:
                continue
            
            # Gene names are space-separated
            gene_names = cols[1].split(' ') if cols[1] else []
            
            for gene in gene_names:
                gene = gene.strip()
                if gene:
                    genes.append(gene)
        
        return genes
    
    def generate_urls(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Generate UniProt URLs for metabolites.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'URL', 'UniProt_ID', 'InChIKey']
        """
        results = []
        
        for metabolite in metabolites:
            inchikey = metabolite.get('InChIKey')
            hmdb_id = metabolite.get('HMDB_ID')
            
            if pd.isna(inchikey):
                continue
            
            # Fetch data
            timeout = config.API_RETRY_CONFIG['UniProtKB']['timeout']
            
            if self.reviewed_only:
                url_query = (
                    f'https://rest.uniprot.org/uniprotkb/search?'
                    f'query=(inchikey:{inchikey})%20AND%20(organism_id:9606)%20AND%20(reviewed:true)'
                    f'&format=tsv&fields=accession,gene_names'
                )
            else:
                url_query = (
                    f'https://rest.uniprot.org/uniprotkb/search?'
                    f'query=(inchikey:{inchikey})%20AND%20(organism_id:9606)'
                    f'&format=tsv&fields=accession,gene_names'
                )
            
            try:
                response = requests.get(url_query, timeout=timeout)
                
                if response.status_code != 200:
                    continue
                
                data = response.text.strip().split('\n')
                if len(data) < 2:
                    continue
                
                for line in data[1:]:
                    cols = line.split('\t')
                    
                    if len(cols) < 2:
                        continue
                    
                    uniprot_id = cols[0]
                    gene_names = cols[1].split(' ') if cols[1] else []
                    
                    for gene in gene_names:
                        gene = gene.strip()
                        if gene:
                            results.append({
                                'HMDB_ID': hmdb_id,
                                'Gene': gene,
                                'URL': f'https://www.uniprot.org/uniprotkb/{uniprot_id}/entry',
                                'UniProt_ID': uniprot_id,
                                'InChIKey': inchikey
                            })
            
            except Exception as e:
                self.logger.warning(f"Failed to fetch URL data for {hmdb_id}: {e}")
                continue
        
        return pd.DataFrame(results)


if __name__ == '__main__':
    # Test the retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Test with first 10 metabolites
    test_metabolites = metabolites[:10]
    
    retriever = UniProtRetriever(reviewed_only=False)
    interactions = retriever.get_interactions(test_metabolites)
    print(f"\nUniProt interactions (test): {len(interactions)}")
    print(interactions.head())
