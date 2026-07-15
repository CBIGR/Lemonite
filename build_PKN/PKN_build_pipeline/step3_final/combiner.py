"""
PKN combiner module for integrating metabolite-gene and PPI networks.

Handles:
- Loading metabolite-gene PKN from Step 1
- Loading PPI network from Step 2
- Combining into unified network format
- Saving final combined PKN
"""

import pandas as pd
import logging
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config


def combine_networks() -> pd.DataFrame:
    """
    Combine metabolite-gene PKN and PPI network into final unified network.
    
    Returns:
    --------
    pd.DataFrame
        Combined network with columns: ['Node1', 'Node2', 'Type', 'Source']
        Type is 'metabolite-gene' or 'PPI'
    """
    logger = logging.getLogger('combiner')
    logger.info("="*80)
    logger.info("COMBINING METABOLITE-GENE PKN AND PPI NETWORK")
    logger.info("="*80)
    
    # Load metabolite-gene PKN
    logger.info(f"\nLoading metabolite-gene PKN from: {config.METABOLITE_GENE_PKN}")
    metabolite_gene = pd.read_csv(config.METABOLITE_GENE_PKN, sep='\t')
    
    logger.info(f"  Metabolite-gene PKN:")
    logger.info(f"    Interactions: {len(metabolite_gene)}")
    # Support both 'Metabolite' (new format) and 'HMDB_ID' (legacy format)
    metabolite_col = 'Metabolite' if 'Metabolite' in metabolite_gene.columns else 'HMDB_ID'
    logger.info(f"    Unique metabolites: {metabolite_gene[metabolite_col].nunique()}")
    logger.info(f"    Unique genes: {metabolite_gene['Gene'].nunique()}")

    # Rename columns to standardized format
    metabolite_gene = metabolite_gene.rename(columns={
        metabolite_col: 'Node1',
        'Gene': 'Node2'
    })
    metabolite_gene['Type'] = 'metabolite-gene'

    # Load PPI network
    logger.info(f"\nLoading PPI network from: {config.PPI_NETWORK}")
    ppi = pd.read_csv(config.PPI_NETWORK, sep='\t')

    # Rename columns to standardized format (handle both GeneA/GeneB and Node1/Node2)
    if 'GeneA' in ppi.columns and 'GeneB' in ppi.columns:
        ppi = ppi.rename(columns={'GeneA': 'Node1', 'GeneB': 'Node2'})

    ppi['Type'] = 'PPI'

    logger.info(f"  PPI network:")
    logger.info(f"    Interactions: {len(ppi)}")
    logger.info(f"    Unique genes: {len(set(ppi['Node1']) | set(ppi['Node2']))}")

    # Concatenate, preserving optional columns (url, reaction_id, combined_score) with NaN fill
    combined = pd.concat([metabolite_gene, ppi], ignore_index=True, sort=False)

    # Keep only the standardized output columns
    output_cols = ['Node1', 'Node2', 'Type', 'Source']
    for opt_col in ['url', 'reaction_id', 'combined_score']:
        if opt_col in combined.columns:
            output_cols.append(opt_col)
    combined = combined[output_cols]
    
    # Get statistics
    metabolite_nodes = set(combined[combined['Type'] == 'metabolite-gene']['Node1'])
    gene_nodes = set(combined['Node2']) | set(combined[combined['Type'] == 'PPI']['Node1'])
    
    logger.info(f"\nCombined network statistics:")
    logger.info(f"  Total interactions: {len(combined)}")
    logger.info(f"  Unique metabolites: {len(metabolite_nodes)}")
    logger.info(f"  Unique genes: {len(gene_nodes)}")
    logger.info(f"  Total unique nodes: {len(metabolite_nodes) + len(gene_nodes)}")
    
    # Save combined network
    output_file = config.FINAL_PKN_FILE
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    combined.to_csv(output_file, sep='\t', index=False)
    logger.info(f"  Saved to: {output_file}")
    
    return combined


if __name__ == '__main__':
    # Test combiner
    import logging
    logging.basicConfig(level=logging.INFO)
    
    combined = combine_networks()
    
    print(f"\nFinal combined PKN:")
    print(f"  Total interactions: {len(combined)}")
    print(combined.head(10))
