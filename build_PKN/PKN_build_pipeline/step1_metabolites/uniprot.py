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
    
    def get_url_for_metabolite(self, metabolite: Dict) -> str:
        """Return UniProtKB search URL by InChIKey (fallback; per-protein URLs stored in get_interactions)."""
        inchikey = metabolite.get('InChIKey')
        if inchikey and pd.notna(inchikey):
            return f'https://www.uniprot.org/uniprotkb/?query=inchikey:{inchikey}&organism_id=9606'
        return ''

    @retry_api_call(db_name='UniProtKB')
    def _fetch_with_urls(self, metabolite: Dict) -> List[Dict]:
        """
        Fetch UniProt interactions for a single metabolite, returning (gene, url) pairs.

        Returns:
        --------
        list of dict with 'gene' and 'url' keys
        """
        inchikey = metabolite.get('InChIKey')

        if not inchikey or pd.isna(inchikey):
            return []

        timeout = config.API_RETRY_CONFIG['UniProtKB']['timeout']

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
        response.raise_for_status()

        data = response.text.strip().split('\n')
        if len(data) < 2:
            return []

        pairs = []
        for line in data[1:]:
            cols = line.split('\t')
            if len(cols) < 2:
                continue
            uniprot_id = cols[0].strip()
            gene_names = cols[1].strip().split(' ') if cols[1].strip() else []
            entry_url = f'https://www.uniprot.org/uniprotkb/{uniprot_id}/entry'
            for gene in gene_names:
                gene = gene.strip()
                if gene:
                    pairs.append({'gene': gene, 'url': entry_url})
        return pairs

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
        response.raise_for_status()
        
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

    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get UniProtKB interactions for metabolites, storing per-protein entry URLs.

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

        self.logger.info("Fetching interactions from UniProtKB API (with per-protein URLs)...")

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

            for future in tqdm(as_completed(futures), total=len(remaining), desc="UniProtKB API"):
                metabolite = futures[future]
                hmdb_id = metabolite.get('HMDB_ID')
                call_count += 1

                try:
                    pairs = future.result()
                    for pair in pairs:
                        results.append({
                            'HMDB_ID': hmdb_id,
                            'Gene': pair['gene'],
                            'Source': 'UniProtKB',
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
        self.logger.info(f"Retrieved {len(df)} interactions from UniProtKB")
        self.save_cache(df)
        self._clear_checkpoint()
        return df
    
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
