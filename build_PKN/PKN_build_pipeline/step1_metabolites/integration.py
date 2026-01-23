"""
Integration module for combining metabolite-gene interactions from all databases.

Handles:
- Combining results from all retrievers
- Creating final metabolite-gene PKN
- Generating visualizations (UpSet plot, barplot)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from typing import Dict
import os
import sys

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config


def integrate_databases(results: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Integrate results from all database retrievers.
    
    Parameters:
    -----------
    results : dict
        Mapping of database name → interactions DataFrame
        Each DataFrame must have columns: ['HMDB_ID', 'Gene', 'Source']
    
    Returns:
    --------
    pd.DataFrame
        Combined metabolite-gene network with columns:
        ['HMDB_ID', 'Gene', 'Source']
    """
    logger = logging.getLogger('integration')
    logger.info("Integrating database results...")
    
    # Filter out None/empty results
    valid_results = {}
    for db_name, df in results.items():
        if df is not None and len(df) > 0:
            valid_results[db_name] = df
            logger.info(f"  {db_name}: {len(df)} interactions")
        else:
            logger.warning(f"  {db_name}: No interactions (skipped)")
    
    if not valid_results:
        logger.error("No valid results to integrate!")
        return pd.DataFrame(columns=['HMDB_ID', 'Gene', 'Source'])
    
    # Combine all dataframes
    combined = pd.concat(valid_results.values(), ignore_index=True)
    
    # Remove duplicates (same metabolite-gene pair from multiple sources)
    # Keep all source information by grouping
    combined_grouped = combined.groupby(['HMDB_ID', 'Gene'])['Source'].apply(
        lambda x: '|'.join(sorted(set(x)))
    ).reset_index()
    
    logger.info(f"\nIntegration summary:")
    logger.info(f"  Total interactions: {len(combined)}")
    logger.info(f"  Unique metabolite-gene pairs: {len(combined_grouped)}")
    logger.info(f"  Unique metabolites: {combined_grouped['HMDB_ID'].nunique()}")
    logger.info(f"  Unique genes: {combined_grouped['Gene'].nunique()}")
    
    # Save final network
    output_file = config.METABOLITE_GENE_PKN
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    combined_grouped.to_csv(output_file, index=False, sep='\t')
    logger.info(f"  Saved to: {output_file}")
    
    return combined_grouped


def create_database_matrix(results: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Create binary matrix showing which databases have interactions for each metabolite.
    
    Parameters:
    -----------
    results : dict
        Mapping of database name → interactions DataFrame
    
    Returns:
    --------
    pd.DataFrame
        Binary matrix with metabolites as rows, databases as columns
    """
    logger = logging.getLogger('integration')
    logger.info("Creating database matrix...")
    
    # Get all unique metabolites
    all_metabolites = set()
    for df in results.values():
        if df is not None and len(df) > 0:
            all_metabolites.update(df['HMDB_ID'].unique())
    
    all_metabolites = sorted(all_metabolites)
    
    # Create binary matrix
    matrix_data = {}
    
    for db_name, df in results.items():
        if df is None or len(df) == 0:
            matrix_data[db_name] = [0] * len(all_metabolites)
        else:
            metabolites_with_interactions = set(df['HMDB_ID'].unique())
            matrix_data[db_name] = [
                1 if met in metabolites_with_interactions else 0
                for met in all_metabolites
            ]
    
    matrix = pd.DataFrame(matrix_data, index=all_metabolites)
    
    logger.info(f"Created matrix: {len(matrix)} metabolites × {len(matrix.columns)} databases")
    
    return matrix


def create_upset_plot(matrix: pd.DataFrame, output_dir: str):
    """
    Create UpSet plot showing database overlap.
    
    Parameters:
    -----------
    matrix : pd.DataFrame
        Binary matrix from create_database_matrix()
    output_dir : str
        Directory to save plot
    """
    logger = logging.getLogger('integration')
    logger.info("Creating UpSet plot...")
    
    try:
        from upsetplot import plot as upset_plot, from_indicators
        
        # Convert to UpSet format
        upset_data = from_indicators(matrix.columns, data=matrix)
        
        # Create plot
        fig = plt.figure(figsize=(12, 6))
        upset_plot(upset_data, fig=fig, show_counts=True)
        plt.suptitle('Database Coverage Overlap', fontsize=14, y=0.98)
        
        # Save
        os.makedirs(os.path.join(output_dir, 'figures'), exist_ok=True)
        output_file = os.path.join(output_dir, 'figures', 'database_upset_plot.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Saved UpSet plot to: {output_file}")
    
    except ImportError:
        logger.warning("upsetplot package not installed, skipping UpSet plot")
    except Exception as e:
        logger.error(f"Failed to create UpSet plot: {e}")


def create_coverage_barplot(matrix: pd.DataFrame, output_dir: str):
    """
    Create barplot showing number of metabolites per database.
    
    Parameters:
    -----------
    matrix : pd.DataFrame
        Binary matrix from create_database_matrix()
    output_dir : str
        Directory to save plot
    """
    logger = logging.getLogger('integration')
    logger.info("Creating coverage barplot...")
    
    try:
        # Calculate counts per database
        counts = matrix.sum(axis=0).sort_values(ascending=False)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        counts.plot(kind='bar', ax=ax, color='steelblue')
        ax.set_xlabel('Database', fontsize=12)
        ax.set_ylabel('Number of Metabolites', fontsize=12)
        ax.set_title('Metabolite Coverage by Database', fontsize=14)
        ax.tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for i, v in enumerate(counts):
            ax.text(i, v + max(counts)*0.01, str(int(v)), 
                   ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        
        # Save
        os.makedirs(os.path.join(output_dir, 'figures'), exist_ok=True)
        output_file = os.path.join(output_dir, 'figures', 'database_coverage_barplot.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Saved barplot to: {output_file}")
    
    except Exception as e:
        logger.error(f"Failed to create coverage barplot: {e}")


def create_visualizations(results: Dict[str, pd.DataFrame], output_dir: str = None):
    """
    Create all visualizations for database integration.
    
    Parameters:
    -----------
    results : dict
        Mapping of database name → interactions DataFrame
    output_dir : str, optional
        Output directory (defaults to config.OUTPUT_DIR)
    """
    if output_dir is None:
        output_dir = config.OUTPUT_DIR
    
    logger = logging.getLogger('integration')
    logger.info("Creating visualizations...")
    
    # Create database matrix
    matrix = create_database_matrix(results)
    
    # Create plots
    create_upset_plot(matrix, output_dir)
    create_coverage_barplot(matrix, output_dir)
    
    logger.info("Visualizations complete!")


if __name__ == '__main__':
    # Test integration
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # Create dummy results for testing
    results = {
        'BioGRID': pd.DataFrame({
            'HMDB_ID': ['HMDB0000001', 'HMDB0000002'],
            'Gene': ['TP53', 'BRCA1'],
            'Source': ['BioGRID', 'BioGRID']
        }),
        'STITCH': pd.DataFrame({
            'HMDB_ID': ['HMDB0000001', 'HMDB0000003'],
            'Gene': ['TP53', 'EGFR'],
            'Source': ['STITCH', 'STITCH']
        }),
        'MetalinksDB': pd.DataFrame({
            'HMDB_ID': ['HMDB0000002', 'HMDB0000003'],
            'Gene': ['BRCA1', 'EGFR'],
            'Source': ['MetalinksDB', 'MetalinksDB']
        })
    }
    
    # Test integration
    final_network = integrate_databases(results)
    print(f"\nFinal network: {len(final_network)} interactions")
    print(final_network.head())
    
    # Test visualizations
    create_visualizations(results, config.OUTPUT_DIR)
