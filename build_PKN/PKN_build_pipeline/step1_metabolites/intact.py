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
    
    def get_url_for_metabolite(self, metabolite: Dict) -> str:
        """Return IntAct search URL for a metabolite by ChEBI ID (fallback)."""
        chebi_id = metabolite.get('ChEBI')
        if chebi_id and pd.notna(chebi_id):
            try:
                return f'https://www.ebi.ac.uk/intact/search?query=CHEBI:{int(float(chebi_id))}'
            except (ValueError, TypeError):
                pass
        return ''

    @retry_api_call(db_name='IntAct')
    def _fetch_with_urls(self, metabolite: Dict) -> List[Dict]:
        """
        Fetch IntAct interactions for a single metabolite, returning (gene, ebi_url) pairs.

        Returns:
        --------
        list of dict with 'gene' and 'url' keys
        """
        chebi_id = metabolite.get('ChEBI')

        if chebi_id == 'NA' or not chebi_id or pd.isna(chebi_id):
            return []

        try:
            chebi_id = int(chebi_id)
        except (ValueError, TypeError):
            self.logger.warning(f"Invalid ChEBI_id value: {chebi_id}")
            return []

        timeout = config.API_RETRY_CONFIG['IntAct']['timeout']
        url = (
            f'https://www.ebi.ac.uk/intact/ws/interaction/list?'
            f'draw=1&maxMIScore=1&minMIScore={self.min_mi_score}'
            f'&negativeFilter=POSITIVE_ONLY&page=0&pageSize=10000'
            f'&query=CHEBI%3A{chebi_id}'
        )

        response = requests.post(url, timeout=timeout)
        response.raise_for_status()

        data = json.loads(response.text)

        if 'data' not in data or not data['data']:
            return []

        df = pd.DataFrame(data['data'])
        df = df[
            df['taxIdAStyled'].str.contains('9606', na=False) |
            df['taxIdBStyled'].str.contains('9606', na=False)
        ]

        if df.empty:
            return []

        pairs = []
        for _, row in df.iterrows():
            ebi_id = row.get('ac')
            if row['idA'] == f'CHEBI:{chebi_id} (chebi)':
                gene = row.get('moleculeB')
            elif row['idB'] == f'CHEBI:{chebi_id} (chebi)':
                gene = row.get('moleculeA')
            else:
                continue
            if gene and pd.notna(gene):
                ebi_url = (
                    f'https://www.ebi.ac.uk/intact/details/interaction/{ebi_id}'
                    if ebi_id and pd.notna(ebi_id) else ''
                )
                pairs.append({'gene': gene, 'url': ebi_url})
        return pairs

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
        response.raise_for_status()
        
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

    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get IntAct interactions for metabolites, storing per-interaction EBI URLs.

        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'Source', 'url']
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from tqdm import tqdm
        from utils.api_retry import rate_limit_pause

        cached = self.load_cache()
        if cached is not None:
            return cached

        self.logger.info("Fetching interactions from IntAct API (with per-interaction EBI URLs)...")

        already_processed = self._load_checkpoint()
        results = []

        if already_processed and self.cache_file:
            partial_file = self.cache_file + '.partial'
            if os.path.exists(partial_file):
                try:
                    partial_df = pd.read_csv(partial_file)
                    results = partial_df.to_dict('records')
                    self.logger.info(
                        f"Resuming: {len(already_processed)} done, {len(results)} interactions loaded"
                    )
                except Exception:
                    already_processed = set()
                    results = []

        remaining = [m for m in metabolites if m.get('HMDB_ID') not in already_processed]
        self.logger.info(f"Processing {len(remaining)} metabolites ({len(already_processed)} already done)")

        processed_ids = set(already_processed)
        call_count = 0
        checkpoint_interval = 200

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self._fetch_with_urls, m): m for m in remaining}

            for future in tqdm(as_completed(futures), total=len(remaining), desc="IntAct API"):
                metabolite = futures[future]
                hmdb_id = metabolite.get('HMDB_ID')
                call_count += 1

                try:
                    pairs = future.result()
                    for pair in pairs:
                        results.append({
                            'HMDB_ID': hmdb_id,
                            'Gene': pair['gene'],
                            'Source': 'IntAct',
                            'url': pair['url'],
                        })
                except Exception as e:
                    self.logger.error(f"Failed to process {hmdb_id}: {e}")

                processed_ids.add(hmdb_id)

                if call_count % checkpoint_interval == 0:
                    self._save_checkpoint(processed_ids, results)
                    self.logger.info(f"Checkpoint: {len(processed_ids)}/{len(metabolites)} processed")

                rate_limit_pause(call_count, self.db_name)

        df = pd.DataFrame(results) if results else pd.DataFrame(
            columns=['HMDB_ID', 'Gene', 'Source', 'url']
        )
        # Keep first occurrence per (HMDB_ID, Gene) to avoid duplicates from multiple interactions
        if len(df) > 0:
            df = df.drop_duplicates(subset=['HMDB_ID', 'Gene'])
        self.logger.info(f"Retrieved {len(df)} interactions from IntAct")
        self.save_cache(df)
        self._clear_checkpoint()
        return df
    
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
