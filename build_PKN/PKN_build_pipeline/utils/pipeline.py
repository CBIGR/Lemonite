"""
Base classes for database retrievers in the PKN pipeline.

Provides abstract base classes and utilities for implementing
consistent retriever interfaces across all data sources.
"""

from abc import ABC, abstractmethod
import pandas as pd
import logging
from typing import List, Dict, Optional
import os
import json
import time


class DatabaseRetriever(ABC):
    """
    Abstract base class for all database retrievers.
    
    All retrievers must implement get_interactions() which returns
    a standardized DataFrame with columns: ['HMDB_ID', 'Gene', 'Source']
    or ['Node1', 'Node2', 'Source'] for PPI networks.
    """
    
    def __init__(self, db_name: str, cache_file: Optional[str] = None):
        self.db_name = db_name
        self.cache_file = cache_file
        self.logger = logging.getLogger(self.db_name)
    
    @abstractmethod
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        pass
    
    def load_cache(self) -> Optional[pd.DataFrame]:
        """Load cached results if available. Logs cache age."""
        if self.cache_file and os.path.exists(self.cache_file):
            meta_file = self.cache_file + '.meta'
            age_info = ""
            if os.path.exists(meta_file):
                try:
                    with open(meta_file, 'r') as f:
                        meta = json.load(f)
                    created = meta.get('created_at', 'unknown')
                    age_info = f" (cached on {created})"
                except Exception:
                    pass
            self.logger.info(f"Loading cached data from {self.cache_file}{age_info}")
            return pd.read_csv(self.cache_file)
        return None
    
    def save_cache(self, df: pd.DataFrame):
        """Save results to cache file with metadata."""
        if self.cache_file:
            os.makedirs(os.path.dirname(self.cache_file), exist_ok=True)
            df.to_csv(self.cache_file, index=False)
            # Write metadata
            meta_file = self.cache_file + '.meta'
            meta = {
                'created_at': time.strftime('%Y-%m-%d %H:%M:%S'),
                'db_name': self.db_name,
                'row_count': len(df),
            }
            try:
                with open(meta_file, 'w') as f:
                    json.dump(meta, f, indent=2)
            except Exception:
                pass
            self.logger.info(f"Saved {len(df)} interactions to cache: {self.cache_file}")


class LocalFileRetriever(DatabaseRetriever):
    """
    Base class for retrievers that process local files (no API calls).
    """
    
    def __init__(self, db_name: str = '', file_path: str = '', cache_file: Optional[str] = None):
        super().__init__(db_name, cache_file)
        self.file_path = file_path
        
        if file_path and not os.path.exists(file_path):
            raise FileNotFoundError(f"Database file not found: {file_path}")
    
    @abstractmethod
    def parse_file(self) -> pd.DataFrame:
        pass
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        cached = self.load_cache()
        if cached is not None:
            return cached
        
        self.logger.info(f"Parsing {self.db_name} file: {self.file_path}")
        df = self.parse_file()
        
        metabolite_ids = {m.get('HMDB_ID') for m in metabolites if m.get('HMDB_ID')}
        df = df[df['HMDB_ID'].isin(metabolite_ids)]
        
        self.logger.info(f"Found {len(df)} interactions for {len(metabolite_ids)} metabolites")
        self.save_cache(df)
        return df


class APIRetriever(DatabaseRetriever):
    """
    Base class for retrievers that use external APIs.
    
    Implements retry logic, rate limiting, and checkpoint/resume support.
    """
    
    def __init__(self, db_name: str, cache_file: Optional[str] = None, 
                 max_workers: int = 4):
        super().__init__(db_name, cache_file)
        self.max_workers = max_workers
    
    @abstractmethod
    def fetch_single(self, metabolite: Dict) -> List[str]:
        pass
    
    def _get_checkpoint_file(self) -> str:
        """Return path for the checkpoint file used for resume."""
        if self.cache_file:
            return self.cache_file + '.checkpoint'
        return os.path.join('/tmp', f'{self.db_name}_checkpoint.json')
    
    def _load_checkpoint(self) -> set:
        """Load set of already-processed HMDB IDs from checkpoint."""
        cp_file = self._get_checkpoint_file()
        if os.path.exists(cp_file):
            try:
                with open(cp_file, 'r') as f:
                    data = json.load(f)
                return set(data.get('processed_ids', []))
            except Exception:
                return set()
        return set()
    
    def _save_checkpoint(self, processed_ids: set, partial_results: list):
        """Save checkpoint with processed IDs and partial results."""
        cp_file = self._get_checkpoint_file()
        os.makedirs(os.path.dirname(cp_file) if os.path.dirname(cp_file) else '.', exist_ok=True)
        with open(cp_file, 'w') as f:
            json.dump({'processed_ids': list(processed_ids)}, f)
        # Also save partial results CSV
        if partial_results and self.cache_file:
            partial_file = self.cache_file + '.partial'
            pd.DataFrame(partial_results).to_csv(partial_file, index=False)
    
    def _clear_checkpoint(self):
        """Remove checkpoint files after successful completion."""
        cp_file = self._get_checkpoint_file()
        for f in [cp_file, self.cache_file + '.partial' if self.cache_file else '']:
            if f and os.path.exists(f):
                os.remove(f)
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from tqdm import tqdm
        from utils.api_retry import rate_limit_pause
        
        # Try full cache first
        cached = self.load_cache()
        if cached is not None:
            return cached
        
        self.logger.info(f"Fetching interactions from {self.db_name} API...")
        
        # Load checkpoint for resume
        already_processed = self._load_checkpoint()
        results = []
        
        # Load partial results if resuming
        if already_processed and self.cache_file:
            partial_file = self.cache_file + '.partial'
            if os.path.exists(partial_file):
                try:
                    partial_df = pd.read_csv(partial_file)
                    results = partial_df.to_dict('records')
                    self.logger.info(
                        f"Resuming from checkpoint: {len(already_processed)} metabolites "
                        f"already processed, {len(results)} interactions loaded"
                    )
                except Exception:
                    already_processed = set()
                    results = []
        
        # Filter to unprocessed metabolites
        remaining = [m for m in metabolites if m.get('HMDB_ID') not in already_processed]
        self.logger.info(f"Processing {len(remaining)} metabolites ({len(already_processed)} already done)")
        
        processed_ids = set(already_processed)
        call_count = 0
        checkpoint_interval = 200
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self.fetch_single, m): m for m in remaining}
            
            for future in tqdm(as_completed(futures), total=len(remaining), 
                              desc=f"{self.db_name} API"):
                metabolite = futures[future]
                hmdb_id = metabolite.get('HMDB_ID')
                call_count += 1
                
                try:
                    genes = future.result()
                    if genes is not None:
                        for gene in genes:
                            results.append({
                                'HMDB_ID': hmdb_id,
                                'Gene': gene,
                                'Source': self.db_name
                            })
                except Exception as e:
                    self.logger.error(f"Failed to process {hmdb_id}: {e}")
                
                processed_ids.add(hmdb_id)
                
                # Periodic checkpoint save
                if call_count % checkpoint_interval == 0:
                    self._save_checkpoint(processed_ids, results)
                    self.logger.info(f"Checkpoint saved: {len(processed_ids)}/{len(metabolites)} processed")
                
                # Rate limiting
                rate_limit_pause(call_count, self.db_name)
        
        df = pd.DataFrame(results)
        self.logger.info(f"Retrieved {len(df)} interactions from {self.db_name}")
        
        self.save_cache(df)
        self._clear_checkpoint()
        
        return df
