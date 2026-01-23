#!/usr/bin/env python3
"""
Parse GSEA results from MixOmics and MOFA analyses and summarize top enriched pathways.
Extracts top 3 pathways per component/factor for Biological Process (BP) terms,
ranked by absolute normalized enrichment score (NES).
"""

import pandas as pd
import os
from pathlib import Path

def parse_combined_gsea_file(filepath, component_or_factor='Component'):
    """
    Parse combined GSEA file that contains all components/factors.
    
    Args:
        filepath: Path to the combined GSEA results file
        component_or_factor: Either 'Component' or 'Factor'
    
    Returns:
        Dictionary with component/factor as key and list of top pathways as value
    """
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return {}
    
    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return {}
    
    # Filter for Biological Process only
    if 'Database' in df.columns:
        df_bp = df[df['Database'] == 'BP'].copy()
    else:
        print(f"Warning: No 'Database' column in {filepath}")
        return {}
    
    if df_bp.empty:
        print(f"No BP (Biological Process) entries found in {filepath}")
        return {}
    
    # Add absolute NES column
    df_bp['abs_NES'] = df_bp['NES'].abs()
    
    # Group by component/factor
    results = {}
    component_col = component_or_factor
    
    if component_col not in df_bp.columns:
        print(f"Warning: Column '{component_col}' not found in {filepath}")
        return {}
    
    for comp_factor in df_bp[component_col].unique():
        comp_df = df_bp[df_bp[component_col] == comp_factor].copy()
        
        # Sort by absolute NES (descending)
        comp_df_sorted = comp_df.sort_values('abs_NES', ascending=False)
        
        # Get top 3
        top3 = comp_df_sorted.head(3)
        
        pathways = []
        for _, row in top3.iterrows():
            # Add + or - prefix based on NES sign
            sign = '+' if row['NES'] > 0 else '-'
            pathway_name = f"{sign} {row['Description']}"
            pathways.append({
                'pathway': pathway_name,
                'NES': row['NES'],
                'abs_NES': row['abs_NES'],
                'pvalue': row['pvalue'],
                'p.adjust': row['p.adjust']
            })
        
        results[comp_factor] = pathways
    
    return results


def parse_individual_factor_files(gsea_dir, n_factors):
    """
    Parse individual factor GSEA files from MOFA analysis.
    
    Args:
        gsea_dir: Directory containing Factor_X subdirectories
        n_factors: Number of factors to parse
    
    Returns:
        Dictionary with factor as key and list of top pathways as value
    """
    if not os.path.exists(gsea_dir):
        print(f"Directory not found: {gsea_dir}")
        return {}
    
    results = {}
    
    for factor_num in range(1, n_factors + 1):
        factor_dir = os.path.join(gsea_dir, f"Factor_{factor_num}")
        bp_file = os.path.join(factor_dir, f"Factor{factor_num}_gseGO_BP_results.txt")
        
        if not os.path.exists(bp_file):
            print(f"File not found: {bp_file}")
            continue
        
        try:
            df = pd.read_csv(bp_file, sep='\t')
        except Exception as e:
            print(f"Error reading {bp_file}: {e}")
            continue
        
        if df.empty:
            continue
        
        # Add absolute NES column
        df['abs_NES'] = df['NES'].abs()
        
        # Sort by absolute NES (descending)
        df_sorted = df.sort_values('abs_NES', ascending=False)
        
        # Get top 3
        top3 = df_sorted.head(3)
        
        pathways = []
        for _, row in top3.iterrows():
            # Add + or - prefix based on NES sign
            sign = '+' if row['NES'] > 0 else '-'
            pathway_name = f"{sign} {row['Description']}"
            pathways.append({
                'pathway': pathway_name,
                'NES': row['NES'],
                'abs_NES': row['abs_NES'],
                'pvalue': row['pvalue'],
                'p.adjust': row['p.adjust']
            })
        
        results[factor_num] = pathways
    
    return results


def write_summary(output_file, all_results):
    """
    Write summary of top pathways to a text file.
    
    Args:
        output_file: Path to output file
        all_results: Dictionary with analysis names and their results
    """
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("TOP 3 ENRICHED BIOLOGICAL PROCESS PATHWAYS PER FACTOR/COMPONENT\n")
        f.write("Ranked by Absolute Normalized Enrichment Score (NES)\n")
        f.write("=" * 80 + "\n\n")
        
        for analysis_name, results in all_results.items():
            f.write(f"\n{'=' * 80}\n")
            f.write(f"{analysis_name}\n")
            f.write(f"{'=' * 80}\n\n")
            
            if not results:
                f.write("No results found.\n")
                continue
            
            for comp_factor, pathways in sorted(results.items()):
                if isinstance(comp_factor, int):
                    f.write(f"Factor {comp_factor}:\n")
                else:
                    f.write(f"Component {comp_factor}:\n")
                f.write("-" * 80 + "\n")
                
                for i, pathway_info in enumerate(pathways, 1):
                    f.write(f"  {i}. {pathway_info['pathway']}\n")
                    f.write(f"     NES: {pathway_info['NES']:.3f} (|NES|: {pathway_info['abs_NES']:.3f})\n")
                    f.write(f"     p-value: {pathway_info['pvalue']:.2e}, p.adjust: {pathway_info['p.adjust']:.2e}\n")
                    if i < len(pathways):
                        f.write("\n")
                
                f.write("\n")


def main():
    # Define paths to GSEA result files
    lloyd_price_base = "/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results"
    wang_gbm_base = "/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results"
    
    all_results = {}
    
    # 1. Lloyd-Price MixOmics
    print("Parsing Lloyd-Price MixOmics GSEA results...")
    lloyd_mixomics_file = os.path.join(lloyd_price_base, "MixOmics/GSEA_all_components/Combined_GSEA_all_components.txt")
    lloyd_mixomics_results = parse_combined_gsea_file(lloyd_mixomics_file, component_or_factor='Component')
    if lloyd_mixomics_results:
        all_results["Lloyd-Price - MixOmics Analysis"] = lloyd_mixomics_results
    
    # 2. Lloyd-Price MOFA (individual files)
    print("Parsing Lloyd-Price MOFA GSEA results...")
    lloyd_mofa_dir = os.path.join(lloyd_price_base, "MOFA/GSEA_all_factors")
    lloyd_mofa_results = parse_individual_factor_files(lloyd_mofa_dir, n_factors=15)
    if lloyd_mofa_results:
        all_results["Lloyd-Price - MOFA Analysis"] = lloyd_mofa_results
    
    # 3. Wang_GBM MixOmics
    print("Parsing Wang_GBM MixOmics GSEA results...")
    wang_mixomics_file = os.path.join(wang_gbm_base, "MixOmics/GSEA_all_components/Combined_GSEA_all_components.txt")
    wang_mixomics_results = parse_combined_gsea_file(wang_mixomics_file, component_or_factor='Component')
    if wang_mixomics_results:
        all_results["Wang_GBM - MixOmics Analysis"] = wang_mixomics_results
    
    # 4. Wang_GBM MOFA
    print("Parsing Wang_GBM MOFA GSEA results...")
    wang_mofa_file = os.path.join(wang_gbm_base, "MOFA_with_lipidomics/GSEA_all_factors/Combined_GSEA_GO_BP_MF_CC_all_factors.txt")
    wang_mofa_results = parse_combined_gsea_file(wang_mofa_file, component_or_factor='Factor')
    if wang_mofa_results:
        all_results["Wang_GBM - MOFA Analysis"] = wang_mofa_results
    
    # Write summary to file
    output_file = "/home/borisvdm/repo/LemonIte/GSEA_top_pathways_summary.txt"
    print(f"\nWriting summary to {output_file}...")
    write_summary(output_file, all_results)
    
    print(f"Done! Summary written to {output_file}")


if __name__ == "__main__":
    main()
