"""
IntAct metabolite-protein interaction retriever.

Uses IntAct REST API to search for molecular interactions by ChEBI ID.
"""

import pandas as pd
import requests
import json
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pipeline import APIRetriever
from utils.api_retry import retry_api_call
import config


class IntActRetriever(APIRetriever):
    """
    Retrieve metabolite-gene interactions from IntAct database.
    
    IntAct provides molecular interaction data with ChEBI identifiers.
    Uses REST API with MI score filtering.
    """
    
    def __init__(self, min_mi_score=0, cache_file=None):
        """
        Initialize IntAct retriever.
        
        Parameters:
        -----------
        min_mi_score : float
            Minimum MI (molecular interaction) score
        """
        super().__init__(
            db_name='IntAct',
            cache_file=cache_file or config.DB_OUTPUT_FILES['IntAct'],
            max_workers=config.API_RETRY_CONFIG['IntAct']['max_workers']
        )
        self.min_mi_score = min_mi_score
    
    @retry_api_call(db_name='IntAct')
    def fetch_single(self, metabolite: Dict) -> List[str]:
        """
        Fetch IntAct interactions for a single metabolite.
        
        Parameters:
        -----------
        metabolite : dict
            Must contain 'ChEBI' (numeric ChEBI ID without prefix)
        
        Returns:
        --------
        list of str
            Gene symbols
        """
        chebi_id = metabolite.get('ChEBI')
        
        if chebi_id == 'NA' or pd.isna(chebi_id):
            return []
        
        try:
            chebi_id = int(chebi_id)
        except (ValueError, TypeError):
            self.logger.warning(f"Invalid ChEBI_id value: {chebi_id}")
            return []
        
        # Get timeout from config
        timeout = config.API_RETRY_CONFIG['IntAct']['timeout']
        
        # Build query URL
        url = (
            f'https://www.ebi.ac.uk/intact/ws/interaction/list?'
            f'draw=1&maxMIScore=1&minMIScore={self.min_mi_score}'
            f'&negativeFilter=POSITIVE_ONLY&page=0&pageSize=10000'
            f'&query=CHEBI%3A{chebi_id}'
        )
        
        response = requests.post(url, timeout=timeout)
        
        if response.status_code != 200:
            raise requests.exceptions.RequestException(f"HTTP {response.status_code}")
        
        data = json.loads(response.text)
        
        if 'data' not in data or not data['data']:
            return []
        
        df = pd.DataFrame(data['data'])
        
        # Filter for human (taxId: 9606)
        df = df[
            df['taxIdAStyled'].str.contains('9606', na=False) | 
            df['taxIdBStyled'].str.contains('9606', na=False)
        ]
        
        if df.empty:
            return []
        
        genes = []
        for _, row in df.iterrows():
            # Check which molecule is the ChEBI compound
            if row['idA'] == f'CHEBI:{chebi_id} (chebi)':
                gene = row.get('moleculeB')
            elif row['idB'] == f'CHEBI:{chebi_id} (chebi)':
                gene = row.get('moleculeA')
            else:
                continue
            
            if gene and pd.notna(gene):
                genes.append(gene)
        
        return genes
    
    def generate_urls(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Generate IntAct URLs for metabolites.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'URL', 'Interaction_ID', 'ChEBI_ID']
        """
        results = []
        
        for metabolite in metabolites:
            chebi_id = metabolite.get('ChEBI')
            hmdb_id = metabolite.get('HMDB_ID')
            
            if chebi_id == 'NA' or pd.isna(chebi_id):
                continue
            
            try:
                chebi_id = int(chebi_id)
            except (ValueError, TypeError):
                continue
            
            # Fetch data
            timeout = config.API_RETRY_CONFIG['IntAct']['timeout']
            
            url_query = (
                f'https://www.ebi.ac.uk/intact/ws/interaction/list?'
                f'draw=1&maxMIScore=1&minMIScore={self.min_mi_score}'
                f'&negativeFilter=POSITIVE_ONLY&page=0&pageSize=10000'
                f'&query=CHEBI%3A{chebi_id}'
            )
            
            try:
                response = requests.post(url_query, timeout=timeout)
                
                if response.status_code != 200:
                    continue
                
                data = json.loads(response.text)
                
                if 'data' not in data or not data['data']:
                    continue
                
                df = pd.DataFrame(data['data'])
                
                # Filter for human
                df = df[
                    df['taxIdAStyled'].str.contains('9606', na=False) | 
                    df['taxIdBStyled'].str.contains('9606', na=False)
                ]
                
                for _, row in df.iterrows():
                    ebi_id = row.get('ac')
                    
                    if row['idA'] == f'CHEBI:{chebi_id} (chebi)':
                        gene = row.get('moleculeB')
                    elif row['idB'] == f'CHEBI:{chebi_id} (chebi)':
                        gene = row.get('moleculeA')
                    else:
                        continue
                    
                    if gene and pd.notna(gene):
                        results.append({
                            'HMDB_ID': hmdb_id,
                            'Gene': gene,
                            'URL': f'https://www.ebi.ac.uk/intact/details/interaction/{ebi_id}',
                            'Interaction_ID': ebi_id,
                            'ChEBI_ID': chebi_id
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
    
    retriever = IntActRetriever(min_mi_score=0)
    interactions = retriever.get_interactions(test_metabolites)
    print(f"\nIntAct interactions (test): {len(interactions)}")
    print(interactions.head())
