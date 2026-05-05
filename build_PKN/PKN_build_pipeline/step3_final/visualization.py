"""
PKN visualization module for creating network statistics and plots.

Handles:
- Network statistics plots
- Database comparison charts
- Degree distribution plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config


def plot_database_comparison(output_dir: str = None):
    """
    Create comparison plot showing metabolites, genes, and interactions per database.
    
    Parameters:
    -----------
    output_dir : str, optional
        Output directory (defaults to config.OUTPUT_DIR)
    """
    if output_dir is None:
        output_dir = config.OUTPUT_DIR
    
    logger = logging.getLogger('visualization')
    logger.info("Creating database comparison plot...")
    
    # Load metabolite-gene PKN
    pkn = pd.read_csv(config.METABOLITE_GENE_PKN, sep='\t')
    
    # Calculate statistics per database
    stats = []
    
    for source in pkn['Source'].unique():
        source_df = pkn[pkn['Source'].str.contains(source, na=False)]
        stats.append({
            'Database': source,
            'Metabolites': source_df['HMDB_ID'].nunique(),
            'Genes': source_df['Gene'].nunique(),
            'Interactions': len(source_df)
        })
    
    stats_df = pd.DataFrame(stats)
    stats_df = stats_df.sort_values('Interactions', ascending=False)
    
    # Create plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: Number of metabolites
    axes[0].barh(stats_df['Database'], stats_df['Metabolites'], color='steelblue')
    axes[0].set_xlabel('Number of Metabolites', fontsize=12)
    axes[0].set_title('Metabolite Coverage', fontsize=14, fontweight='bold')
    axes[0].invert_yaxis()
    
    # Plot 2: Number of genes
    axes[1].barh(stats_df['Database'], stats_df['Genes'], color='coral')
    axes[1].set_xlabel('Number of Genes', fontsize=12)
    axes[1].set_title('Gene Coverage', fontsize=14, fontweight='bold')
    axes[1].invert_yaxis()
    
    # Plot 3: Number of interactions
    axes[2].barh(stats_df['Database'], stats_df['Interactions'], color='mediumseagreen')
    axes[2].set_xlabel('Number of Interactions', fontsize=12)
    axes[2].set_title('Total Interactions', fontsize=14, fontweight='bold')
    axes[2].invert_yaxis()
    
    plt.tight_layout()
    
    # Save
    os.makedirs(os.path.join(output_dir, 'figures'), exist_ok=True)
    output_file = os.path.join(output_dir, 'figures', 'database_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  Saved database comparison to: {output_file}")


def plot_network_statistics(output_dir: str = None):
    """
    Create network statistics plots (degree distribution, etc.).
    
    Parameters:
    -----------
    output_dir : str, optional
        Output directory (defaults to config.OUTPUT_DIR)
    """
    if output_dir is None:
        output_dir = config.OUTPUT_DIR
    
    logger = logging.getLogger('visualization')
    logger.info("Creating network statistics plots...")
    
    # Load final PKN
    pkn = pd.read_csv(config.FINAL_PKN_FILE, sep='\t')
    
    # Calculate degree distribution
    metabolite_gene = pkn[pkn['Type'] == 'metabolite-gene']
    
    # Metabolite degrees (how many genes per metabolite)
    metabolite_degrees = metabolite_gene.groupby('Node1').size()
    
    # Gene degrees (how many metabolites per gene)
    gene_degrees = metabolite_gene.groupby('Node2').size()
    
    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Metabolite degree distribution
    axes[0].hist(metabolite_degrees, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Number of Gene Partners', fontsize=12)
    axes[0].set_ylabel('Number of Metabolites', fontsize=12)
    axes[0].set_title('Metabolite Degree Distribution', fontsize=14, fontweight='bold')
    axes[0].set_yscale('log')
    
    # Plot 2: Gene degree distribution
    axes[1].hist(gene_degrees, bins=50, color='coral', edgecolor='black', alpha=0.7)
    axes[1].set_xlabel('Number of Metabolite Partners', fontsize=12)
    axes[1].set_ylabel('Number of Genes', fontsize=12)
    axes[1].set_title('Gene Degree Distribution', fontsize=14, fontweight='bold')
    axes[1].set_yscale('log')
    
    plt.tight_layout()
    
    # Save
    os.makedirs(os.path.join(output_dir, 'figures'), exist_ok=True)
    output_file = os.path.join(output_dir, 'figures', 'network_statistics.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  Saved network statistics to: {output_file}")


def create_all_visualizations():
    """Create all visualization plots."""
    logger = logging.getLogger('visualization')
    logger.info("="*80)
    logger.info("CREATING VISUALIZATIONS")
    logger.info("="*80)
    
    plot_database_comparison()
    plot_network_statistics()
    
    logger.info("\nAll visualizations complete!")


if __name__ == '__main__':
    # Test visualization
    import logging
    logging.basicConfig(level=logging.INFO)
    
    create_all_visualizations()
