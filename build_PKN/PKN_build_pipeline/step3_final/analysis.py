"""
PKN analysis module for comparing against reference databases.

Handles:
- Comparison with MetalinksDB
- Comparison with MEBOCOST
- Coverage statistics
- Overlap analysis
"""

import pandas as pd
import logging
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config


def analyze_coverage() -> dict:
    """
    Analyze PKN coverage compared to reference databases.
    
    Returns:
    --------
    dict
        Statistics about coverage
    """
    logger = logging.getLogger('analysis')
    logger.info("="*80)
    logger.info("ANALYZING PKN COVERAGE")
    logger.info("="*80)
    
    # Load final PKN
    logger.info(f"\nLoading final PKN from: {config.FINAL_PKN_FILE}")
    pkn = pd.read_csv(config.FINAL_PKN_FILE, sep='\t')
    
    # Get basic statistics
    metabolite_gene = pkn[pkn['Type'] == 'metabolite-gene']
    ppi = pkn[pkn['Type'] == 'PPI']
    
    metabolites = set(metabolite_gene['Node1'])
    genes = set(metabolite_gene['Node2']) | set(ppi['Node1']) | set(ppi['Node2'])
    
    stats = {
        'total_interactions': len(pkn),
        'metabolite_gene_interactions': len(metabolite_gene),
        'ppi_interactions': len(ppi),
        'unique_metabolites': len(metabolites),
        'unique_genes': len(genes),
        'unique_nodes': len(metabolites) + len(genes)
    }
    
    logger.info(f"\nFinal PKN Statistics:")
    logger.info(f"  Total interactions: {stats['total_interactions']}")
    logger.info(f"  Metabolite-gene: {stats['metabolite_gene_interactions']}")
    logger.info(f"  PPI: {stats['ppi_interactions']}")
    logger.info(f"  Unique metabolites: {stats['unique_metabolites']}")
    logger.info(f"  Unique genes: {stats['unique_genes']}")
    logger.info(f"  Total nodes: {stats['unique_nodes']}")
    
    # Compare with individual databases
    logger.info(f"\nCoverage by database:")
    
    for source in metabolite_gene['Source'].unique():
        source_df = metabolite_gene[metabolite_gene['Source'].str.contains(source, na=False)]
        logger.info(f"  {source}: {len(source_df)} interactions, "
                   f"{source_df['Node1'].nunique()} metabolites, "
                   f"{source_df['Node2'].nunique()} genes")
    
    return stats


if __name__ == '__main__':
    # Test analysis
    import logging
    logging.basicConfig(level=logging.INFO)
    
    stats = analyze_coverage()
    
    print(f"\nAnalysis complete:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
