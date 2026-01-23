"""
Integration module for combining protein-protein interactions from multiple databases.

Handles:
- Combining results from STRING, BioGRID, and HuRI
- Creating final PPI network
- Generating Venn diagram visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from typing import Dict, Set, Tuple
import os
import sys

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config


def load_metabolite_gene_network() -> pd.DataFrame:
    """
    Load metabolite-gene PKN from Step 1.
    
    Returns:
    --------
    pd.DataFrame
        Metabolite-gene network with gene list
    """
    logger = logging.getLogger('ppi_integration')
    logger.info(f"Loading metabolite-gene PKN from: {config.METABOLITE_GENE_PKN}")
    
    pkn = pd.read_csv(config.METABOLITE_GENE_PKN, sep='\t')
    genes = sorted(set(pkn['Gene'].unique()))
    
    logger.info(f"Extracted {len(genes)} unique genes from metabolite-gene PKN")
    
    return genes


def integrate_ppi_databases(results: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Integrate PPI results from all databases.
    
    Parameters:
    -----------
    results : dict
        Mapping of database name → interactions DataFrame
        Each DataFrame must have columns: ['GeneA', 'GeneB', 'Source']
    
    Returns:
    --------
    pd.DataFrame
        Combined PPI network with columns: ['GeneA', 'GeneB', 'Source']
        (may include 'combined_score' if STRING is present)
    """
    logger = logging.getLogger('ppi_integration')
    logger.info("Integrating PPI database results...")
    
    # Filter out None/empty results
    valid_results = {}
    for db_name, df in results.items():
        if df is not None and len(df) > 0:
            valid_results[db_name] = df
            logger.info(f"  {db_name}: {len(df)} interactions")
        else:
            logger.warning(f"  {db_name}: No interactions (skipped)")
    
    if not valid_results:
        logger.error("No valid PPI results to integrate!")
        return pd.DataFrame(columns=['GeneA', 'GeneB', 'Source'])
    
    # Combine all dataframes
    combined = pd.concat(valid_results.values(), ignore_index=True)
    
    # Ensure GeneA and GeneB are strings
    combined['GeneA'] = combined['GeneA'].astype(str)
    combined['GeneB'] = combined['GeneB'].astype(str)
    
    # Remove rows with NaN values
    combined = combined[~combined['GeneA'].str.contains('nan', na=False)]
    combined = combined[~combined['GeneB'].str.contains('nan', na=False)]
    combined = combined.dropna(subset=['GeneA', 'GeneB'])
    
    # Get statistics
    all_genes = set(combined['GeneA']) | set(combined['GeneB'])
    
    logger.info(f"\nIntegration summary:")
    logger.info(f"  Total interactions: {len(combined)}")
    logger.info(f"  Unique genes: {len(all_genes)}")
    
    # Count overlaps by grouping
    combined_grouped = combined.groupby(['GeneA', 'GeneB'])['Source'].apply(
        lambda x: '|'.join(sorted(set(x)))
    ).reset_index()
    
    logger.info(f"  Unique gene pairs: {len(combined_grouped)}")
    
    # For STRING, keep the combined_score if available
    if 'combined_score' in combined.columns:
        # For each gene pair, keep the max score from STRING
        score_map = combined[combined['Source'] == 'STRING'].groupby(
            ['GeneA', 'GeneB']
        )['combined_score'].max().to_dict()
        
        combined_grouped['combined_score'] = combined_grouped.apply(
            lambda row: score_map.get((row['GeneA'], row['GeneB']), np.nan),
            axis=1
        )
    
    # Save final network
    output_file = config.PPI_NETWORK
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    combined_grouped.to_csv(output_file, index=False, sep='\t')
    logger.info(f"  Saved to: {output_file}")
    
    return combined_grouped


def create_venn_diagram(results: Dict[str, pd.DataFrame], output_dir: str = None):
    """
    Create Venn diagram showing overlap between PPI databases.
    
    Parameters:
    -----------
    results : dict
        Mapping of database name → interactions DataFrame
    output_dir : str, optional
        Output directory (defaults to config.OUTPUT_DIR)
    """
    if output_dir is None:
        output_dir = config.OUTPUT_DIR
    
    logger = logging.getLogger('ppi_integration')
    logger.info("Creating Venn diagram...")
    
    try:
        from matplotlib_venn import venn3, venn3_circles
        
        # Create sets of unordered (GeneA, GeneB) pairs for each database
        sets = {}
        for db_name, df in results.items():
            if df is not None and len(df) > 0:
                # Create unordered pairs (sort to make A-B same as B-A)
                pairs = set(
                    tuple(sorted([row['GeneA'], row['GeneB']]))
                    for _, row in df.iterrows()
                )
                sets[db_name] = pairs
        
        if len(sets) < 2:
            logger.warning("Need at least 2 databases for Venn diagram, skipping")
            return
        
        # Get the three main sets (STRING, BioGRID_PPI, HuRI)
        set_string = sets.get('STRING', set())
        set_biogrid = sets.get('BioGRID_PPI', set())
        set_huri = sets.get('HuRI', set())
        
        # Calculate overlap counts
        only_string = len(set_string - set_biogrid - set_huri)
        only_biogrid = len(set_biogrid - set_string - set_huri)
        only_huri = len(set_huri - set_string - set_biogrid)
        string_biogrid = len(set_string & set_biogrid - set_huri)
        string_huri = len(set_string & set_huri - set_biogrid)
        biogrid_huri = len(set_biogrid & set_huri - set_string)
        all_three = len(set_string & set_biogrid & set_huri)
        
        # Create Venn diagram
        fig, ax = plt.subplots(figsize=(10, 8))
        
        venn = venn3(
            subsets=(
                only_string, only_biogrid, string_biogrid,
                only_huri, string_huri, biogrid_huri, all_three
            ),
            set_labels=('STRING', 'BioGRID', 'HuRI'),
            ax=ax
        )
        
        # Customize colors
        venn3_circles(subsets=(
            only_string, only_biogrid, string_biogrid,
            only_huri, string_huri, biogrid_huri, all_three
        ), ax=ax, linewidth=1.5)
        
        ax.set_title('PPI Database Overlap', fontsize=14, fontweight='bold')
        
        # Save
        os.makedirs(os.path.join(output_dir, 'figures'), exist_ok=True)
        output_file = os.path.join(output_dir, 'figures', 'ppi_venn_diagram.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Saved Venn diagram to: {output_file}")
        
        # Print statistics
        logger.info(f"\nOverlap statistics:")
        logger.info(f"  STRING only: {only_string}")
        logger.info(f"  BioGRID only: {only_biogrid}")
        logger.info(f"  HuRI only: {only_huri}")
        logger.info(f"  STRING + BioGRID: {string_biogrid}")
        logger.info(f"  STRING + HuRI: {string_huri}")
        logger.info(f"  BioGRID + HuRI: {biogrid_huri}")
        logger.info(f"  All three: {all_three}")
    
    except ImportError:
        logger.warning("matplotlib-venn package not installed, skipping Venn diagram")
    except Exception as e:
        logger.error(f"Failed to create Venn diagram: {e}")


def build_ppi_network() -> pd.DataFrame:
    """
    Build complete PPI network by querying all databases.
    
    Returns:
    --------
    pd.DataFrame
        Combined PPI network
    """
    from step2_proteins.string_api import STRINGRetriever
    from step2_proteins.biogrid_ppi import BioGRIDPPIRetriever
    from step2_proteins.huri import HuRIRetriever
    
    logger = logging.getLogger('ppi_integration')
    logger.info("="*80)
    logger.info("BUILDING PPI NETWORK")
    logger.info("="*80)
    
    # Load genes from metabolite-gene PKN
    genes = load_metabolite_gene_network()
    
    # Query databases
    results = {}
    
    # STRING
    try:
        logger.info("\n" + "="*80)
        logger.info("Querying STRING")
        logger.info("="*80)
        string_retriever = STRINGRetriever(confidence_threshold=400)
        results['STRING'] = string_retriever.get_interactions(genes, chunk_size=1000, max_workers=10)
    except Exception as e:
        logger.error(f"STRING failed: {e}", exc_info=True)
        results['STRING'] = None
    
    # BioGRID
    try:
        logger.info("\n" + "="*80)
        logger.info("Querying BioGRID PPI")
        logger.info("="*80)
        biogrid_retriever = BioGRIDPPIRetriever()
        results['BioGRID_PPI'] = biogrid_retriever.get_interactions(genes)
    except Exception as e:
        logger.error(f"BioGRID failed: {e}", exc_info=True)
        results['BioGRID_PPI'] = None
    
    # HuRI
    try:
        logger.info("\n" + "="*80)
        logger.info("Querying HuRI")
        logger.info("="*80)
        huri_retriever = HuRIRetriever()
        results['HuRI'] = huri_retriever.get_interactions(genes)
    except Exception as e:
        logger.error(f"HuRI failed: {e}", exc_info=True)
        results['HuRI'] = None
    
    # Integrate results
    logger.info("\n" + "="*80)
    logger.info("INTEGRATING PPI RESULTS")
    logger.info("="*80)
    
    ppi_network = integrate_ppi_databases(results)
    create_venn_diagram(results)
    
    logger.info("\n" + "="*80)
    logger.info("PPI NETWORK COMPLETE")
    logger.info("="*80)
    logger.info(f"Final network: {len(ppi_network)} interactions")
    logger.info(f"Saved to: {config.PPI_NETWORK}")
    
    return ppi_network


if __name__ == '__main__':
    # Test PPI integration
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # Build full network
    ppi_network = build_ppi_network()
    
    print(f"\nFinal PPI network:")
    print(f"  Total interactions: {len(ppi_network)}")
    print(f"  Unique genes: {len(set(ppi_network['GeneA']) | set(ppi_network['GeneB']))}")
    print(f"\nFirst 10 interactions:")
    print(ppi_network.head(10))
