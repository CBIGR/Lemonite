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


class DatabaseRetriever(ABC):
    """
    Abstract base class for all database retrievers.
    
    All retrievers must implement get_interactions() which returns
    a standardized DataFrame with columns: ['HMDB_ID', 'Gene', 'Source']
    or ['Node1', 'Node2', 'Source'] for PPI networks.
    """
    
    def __init__(self, db_name: str, cache_file: Optional[str] = None):
        """
        Initialize the retriever.
        
        Parameters:
        -----------
        db_name : str
            Name of the database (used for logging)
        cache_file : str, optional
            Path to cache file for storing intermediate results
        """
        self.db_name = db_name
        self.cache_file = cache_file
        self.logger = logging.getLogger(self.db_name)
    
    @abstractmethod
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Retrieve interactions from the database.
        
        Parameters:
        -----------
        metabolites : list of dict
            List of metabolite records with identifiers (HMDB_ID, InChIKey, etc.)
        
        Returns:
        --------
        pd.DataFrame
            Standardized interaction data with columns:
            - HMDB_ID: Metabolite identifier
            - Gene: Gene symbol
            - Source: Database name
        """
        pass
    
    def load_cache(self) -> Optional[pd.DataFrame]:
        """Load cached results if available."""
        if self.cache_file and os.path.exists(self.cache_file):
            self.logger.info(f"Loading cached data from {self.cache_file}")
            return pd.read_csv(self.cache_file)
        return None
    
    def save_cache(self, df: pd.DataFrame):
        """Save results to cache file."""
        if self.cache_file:
            os.makedirs(os.path.dirname(self.cache_file), exist_ok=True)
            df.to_csv(self.cache_file, index=False)
            self.logger.info(f"Saved {len(df)} interactions to cache: {self.cache_file}")


class LocalFileRetriever(DatabaseRetriever):
    """
    Base class for retrievers that process local files (no API calls).
    
    Examples: BioGRID, STITCH, UniProt XML, Human-GEM files
    """
    
    def __init__(self, db_name: str, file_path: str, cache_file: Optional[str] = None):
        """
        Initialize local file retriever.
        
        Parameters:
        -----------
        db_name : str
            Database name
        file_path : str
            Path to local database file
        cache_file : str, optional
            Path to cache processed results
        """
        super().__init__(db_name, cache_file)
        self.file_path = file_path
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Database file not found: {file_path}")
    
    @abstractmethod
    def parse_file(self) -> pd.DataFrame:
        """Parse the local file into a standardized DataFrame."""
        pass
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get interactions by parsing local file and filtering to metabolites.
        
        Default implementation:
        1. Check cache
        2. Parse file
        3. Filter to metabolites of interest
        4. Save to cache
        """
        # Try cache first
        cached = self.load_cache()
        if cached is not None:
            return cached
        
        # Parse file
        self.logger.info(f"Parsing {self.db_name} file: {self.file_path}")
        df = self.parse_file()
        
        # Filter to metabolites of interest
        metabolite_ids = {m.get('HMDB_ID') for m in metabolites if m.get('HMDB_ID')}
        df = df[df['HMDB_ID'].isin(metabolite_ids)]
        
        self.logger.info(f"Found {len(df)} interactions for {len(metabolite_ids)} metabolites")
        
        # Save cache
        self.save_cache(df)
        
        return df


class APIRetriever(DatabaseRetriever):
    """
    Base class for retrievers that use external APIs.
    
    Examples: UniProtKB, IntAct, chEMBL, STRING
    Implements retry logic and rate limiting.
    """
    
    def __init__(self, db_name: str, cache_file: Optional[str] = None, 
                 max_workers: int = 4):
        """
        Initialize API retriever.
        
        Parameters:
        -----------
        db_name : str
            Database name
        cache_file : str, optional
            Path to cache results
        max_workers : int
            Number of concurrent API workers
        """
        super().__init__(db_name, cache_file)
        self.max_workers = max_workers
    
    @abstractmethod
    def fetch_single(self, metabolite: Dict) -> List[str]:
        """
        Fetch interactions for a single metabolite.
        
        Should be decorated with @retry_api_call for fault tolerance.
        
        Returns:
        --------
        list of str
            List of gene symbols that interact with the metabolite
        """
        pass
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get interactions using concurrent API calls with retry logic.
        
        Default implementation:
        1. Check cache
        2. Use ThreadPoolExecutor for concurrent calls
        3. Collect results
        4. Save to cache
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from tqdm import tqdm
        
        # Try cache first
        cached = self.load_cache()
        if cached is not None:
            return cached
        
        self.logger.info(f"Fetching interactions from {self.db_name} API...")
        
        results = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {executor.submit(self.fetch_single, m): m for m in metabolites}
            
            for future in tqdm(as_completed(futures), total=len(metabolites), 
                              desc=f"{self.db_name} API"):
                metabolite = futures[future]
                try:
                    genes = future.result()
                    for gene in genes:
                        results.append({
                            'HMDB_ID': metabolite.get('HMDB_ID'),
                            'Gene': gene,
                            'Source': self.db_name
                        })
                except Exception as e:
                    self.logger.error(f"Failed to process {metabolite.get('HMDB_ID')}: {e}")
        
        df = pd.DataFrame(results)
        self.logger.info(f"Retrieved {len(df)} interactions from {self.db_name}")
        
        # Save cache
        self.save_cache(df)
        
        return df
