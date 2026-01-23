"""
chEMBL bioactivity data retriever.

Uses chEMBL web resource client to search for compound-target interactions.
"""

import pandas as pd
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pipeline import APIRetriever
from utils.api_retry import retry_api_call
import config


class ChEMBLRetriever(APIRetriever):
    """
    Retrieve metabolite-gene interactions from chEMBL database.
    
    chEMBL provides bioactivity data for compounds against targets.
    Uses chembl_webresource_client library with retry logic.
    """
    
    def __init__(self, cache_file=None):
        """Initialize chEMBL retriever."""
        super().__init__(
            db_name='chEMBL',
            cache_file=cache_file or config.DB_OUTPUT_FILES['chEMBL'],
            max_workers=config.API_RETRY_CONFIG['chEMBL']['max_workers']
        )
        
        # Suppress chEMBL client logs
        logging.getLogger("chembl_webresource_client").setLevel(logging.WARNING)
    
    @retry_api_call(db_name='chEMBL')
    def fetch_single(self, metabolite: Dict) -> List[str]:
        """
        Fetch chEMBL interactions for a single metabolite.
        
        Parameters:
        -----------
        metabolite : dict
            Must contain 'ChEMBL_id'
        
        Returns:
        --------
        list of str
            Gene symbols
        """
        from chembl_webresource_client.new_client import new_client
        
        chembl_id = metabolite.get('ChEMBL_id')
        
        if chembl_id == 'none' or pd.isna(chembl_id) or chembl_id == '' or chembl_id == 'NA':
            return []
        
        # Query activities for this compound in human (tax_id 9606)
        activity = new_client.activity
        act = activity.filter(molecule_chembl_id=chembl_id, target_tax_id='9606')
        
        dat = pd.DataFrame(act)
        
        if dat.empty:
            return []
        
        # Filter by interaction confidence
        if 'activity_comment' in dat.columns:
            dat = dat[dat['activity_comment'].str.lower().isin(['active', 'substrate'])]
        else:
            self.logger.warning(
                f"No 'activity_comment' column found in chEMBL data for ChEMBL_ID: {chembl_id}"
            )
            return []
        
        if dat.empty:
            return []
        
        # Get unique target IDs
        target_ids = dat['target_chembl_id'].unique()
        
        if target_ids.size == 0:
            return []
        
        genes = []
        
        # Query each target to get gene symbols
        for target_id in target_ids:
            try:
                target_client = new_client.target
                target_data = target_client.filter(chembl_id=target_id)
                
                if not target_data or len(target_data) == 0:
                    continue
                
                target_info = target_data[0]
                target_components = target_info.get('target_components', [])
                
                if not target_components:
                    continue
                
                # Extract gene symbols from target components
                for component in target_components:
                    synonyms = component.get('target_component_synonyms', [])
                    
                    if not synonyms:
                        continue
                    
                    info_df = pd.DataFrame(synonyms)
                    
                    if info_df.empty:
                        continue
                    
                    if 'syn_type' not in info_df.columns or 'component_synonym' not in info_df.columns:
                        continue
                    
                    gene_symbols = info_df[
                        info_df['syn_type'] == 'GENE_SYMBOL'
                    ]['component_synonym'].values
                    
                    for gene in gene_symbols:
                        if pd.notna(gene):
                            genes.append(gene)
            
            except Exception as e:
                self.logger.warning(f"Error fetching target {target_id}: {e}")
                continue
        
        return genes
    
    def generate_urls(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Generate chEMBL URLs for metabolites.
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'URL', 'ChEMBL_ID', 'Target_ID']
        """
        from chembl_webresource_client.new_client import new_client
        
        results = []
        
        for metabolite in metabolites:
            chembl_id = metabolite.get('ChEMBL_id')
            hmdb_id = metabolite.get('HMDB_ID')
            
            if chembl_id == 'none' or pd.isna(chembl_id) or chembl_id == '':
                continue
            
            try:
                activity = new_client.activity
                act = activity.filter(molecule_chembl_id=chembl_id, target_tax_id='9606')
                
                dat = pd.DataFrame(act)
                
                if dat.empty:
                    continue
                
                if 'activity_comment' in dat.columns:
                    dat = dat[dat['activity_comment'].str.lower().isin(['active', 'substrate'])]
                
                if dat.empty:
                    continue
                
                target_ids = dat['target_chembl_id'].unique()
                
                for target_id in target_ids:
                    target_client = new_client.target
                    target_data = target_client.filter(chembl_id=target_id)
                    
                    if not target_data or len(target_data) == 0:
                        continue
                    
                    target_info = target_data[0]
                    target_components = target_info.get('target_components', [])
                    
                    for component in target_components:
                        synonyms = component.get('target_component_synonyms', [])
                        
                        if not synonyms:
                            continue
                        
                        info_df = pd.DataFrame(synonyms)
                        
                        if info_df.empty:
                            continue
                        
                        if 'syn_type' not in info_df.columns or 'component_synonym' not in info_df.columns:
                            continue
                        
                        gene_symbols = info_df[
                            info_df['syn_type'] == 'GENE_SYMBOL'
                        ]['component_synonym'].values
                        
                        for gene in gene_symbols:
                            if pd.notna(gene):
                                results.append({
                                    'HMDB_ID': hmdb_id,
                                    'Gene': gene,
                                    'URL': f'https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}',
                                    'ChEMBL_ID': chembl_id,
                                    'Target_ID': target_id
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
    
    # Note: chEMBL requires ChEMBL IDs which need to be pre-computed
    print("\nchEMBL retriever requires ChEMBL IDs to be added to metabolites first")
    print("Run preprocessing.py to add ChEMBL mappings before using chEMBL")
