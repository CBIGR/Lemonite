#!/usr/bin/env python3

"""
Lemonite Pipeline Summary Report Generator

This script generates a comprehensive HTML summary report of the Lemonite pipeline run,
including input data statistics, preprocessing details, parameter documentation,
and results summaries with interactive tables and visualizations.

Author: Boris Vandemoortele
"""

import os
import sys
import argparse
import json
import glob
import re
import gc
from datetime import datetime
from collections import defaultdict, OrderedDict
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np

# Try to import optional dependencies
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: plotly not available. Some visualizations will be disabled.")


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Generate Lemonite Pipeline Summary Report')
    parser.add_argument('--input_dir', required=True, help='Pipeline input directory')
    parser.add_argument('--output_dir', required=True, help='Pipeline output/results directory')
    parser.add_argument('--run_id', required=True, help='Pipeline run identifier')
    parser.add_argument('--regulator_types', default='TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt',
                       help='Regulator types configuration')
    parser.add_argument('--parameters_file', default=None, help='Path to pipeline_parameters_log.txt')
    parser.add_argument('--organism', default='human', help='Organism (human/mouse)')
    args = parser.parse_args()
    return args


def safe_read_file(filepath, sep='\t', nrows=None, **kwargs):
    """Safely read a file, returning None if it doesn't exist"""
    if filepath and os.path.exists(filepath):
        try:
            # Default to limiting rows unless explicitly specified
            if nrows is None and 'nrows' not in kwargs:
                nrows = 10000  # Limit to 10k rows by default for memory safety
            return pd.read_csv(filepath, sep=sep, nrows=nrows, **kwargs)
        except Exception as e:
            print(f"Warning: Could not read {filepath}: {e}")
    return None


def count_lines(filepath):
    """Count non-empty lines in a file"""
    if not filepath or not os.path.exists(filepath):
        return 0
    try:
        with open(filepath, 'r') as f:
            return sum(1 for line in f if line.strip())
    except:
        return 0


def parse_regulator_types(regulator_types_str):
    """Parse regulator types string into list of dicts
    
    Format: 'Prefix1:DataFile1[:DataType1],Prefix2:DataFile2[:DataType2]'
    DataType is optional: 'c' for continuous (default), 'd' for discrete/binary
    """
    configs = []
    for config_str in regulator_types_str.split(','):
        config_str = config_str.strip()
        if ':' in config_str:
            parts = config_str.split(':')
            prefix = parts[0].strip()
            data_file = parts[1].strip() if len(parts) > 1 else ''
            data_type = parts[2].strip().lower() if len(parts) > 2 else 'c'
            configs.append({'prefix': prefix, 'data_file': data_file, 'data_type': data_type})
    return configs


def collect_input_statistics(input_dir, regulator_configs, output_dir='.'):
    """Collect statistics about input data
    
    Args:
        input_dir: Original input directory (may not be accessible in container)
        regulator_configs: List of regulator configuration dictionaries
        output_dir: Output directory where LemonTree results are staged (default: current dir)
    """
    stats = {
        'data_types': [],
        'samples': 0,
        'features_per_type': {},
        'metadata_columns': [],
        'files_found': []
    }
    
    # Try multiple possible data directory locations
    # IMPORTANT: Include local LemonTree/Preprocessing directory which is always staged by Nextflow
    search_dirs = [
        os.path.join(output_dir, 'LemonTree', 'Preprocessing'),  # Staged preprocessing data (always accessible)
        os.path.join('.', 'LemonTree', 'Preprocessing'),  # Current directory staging
        os.path.join(input_dir, 'data'),  # Original input location
        input_dir,
        os.path.join(input_dir, '..', 'data')  # One level up
    ]
    
    # DEBUG: Print what we're searching for
    print(f"[DEBUG] Input dir: {input_dir}", file=sys.stderr)
    print(f"[DEBUG] Output dir: {output_dir}", file=sys.stderr)
    print(f"[DEBUG] Search dirs: {search_dirs}", file=sys.stderr)
    
    data_dir = None
    for d in search_dirs:
        if os.path.exists(d):
            data_dir = d
            break
    
    if data_dir is None:
        data_dir = input_dir
    
    # Find expression file - look for raw count data in input directories first,
    # only fall back to preprocessed files if originals are not accessible
    # Search raw input files only in input data directories (not preprocessing dirs)
    input_data_dirs = [
        os.path.join(input_dir, 'data'),
        input_dir,
        os.path.join(input_dir, '..', 'data')
    ]
    raw_expr_patterns = [
        '*counts*.tsv', '*Counts*.tsv', '*expression*.tsv', '*host_tx_counts.tsv'
    ]
    # Preprocessed files as fallback (only in preprocessing dirs)
    preproc_dirs = [
        os.path.join(output_dir, 'LemonTree', 'Preprocessing'),
        os.path.join('.', 'LemonTree', 'Preprocessing')
    ]
    preproc_patterns = [
        'LemonPreprocessed_expression.txt',
        'LemonPreprocessed_complete.txt'
    ]
    expr_file = None
    # First: search for raw count files in input data directories
    for pattern in raw_expr_patterns:
        print(f"[DEBUG] Searching for raw pattern: {pattern}", file=sys.stderr)
        for search_dir in input_data_dirs:
            if not os.path.exists(search_dir):
                print(f"[DEBUG]   Dir doesn't exist: {search_dir}", file=sys.stderr)
                continue
            print(f"[DEBUG]   Searching in: {search_dir}", file=sys.stderr)
            matches = glob.glob(os.path.join(search_dir, pattern))
            print(f"[DEBUG]   Found {len(matches)} matches", file=sys.stderr)
            if matches:
                expr_file = matches[0]
                print(f"[DEBUG]   Selected: {expr_file}", file=sys.stderr)
                break
        if expr_file:
            break
    # Fallback: search for preprocessed files
    if not expr_file:
        for pattern in preproc_patterns:
            print(f"[DEBUG] Searching for preprocessed pattern: {pattern}", file=sys.stderr)
            for search_dir in preproc_dirs:
                if not os.path.exists(search_dir):
                    continue
                matches = glob.glob(os.path.join(search_dir, pattern))
                if matches:
                    expr_file = matches[0]
                    print(f"[DEBUG]   Fallback selected: {expr_file}", file=sys.stderr)
                    break
            if expr_file:
                break
    
    if expr_file:
        stats['files_found'].append(('Transcriptomics', os.path.basename(expr_file)))
        stats['data_types'].append('Transcriptomics')
        
        df = safe_read_file(expr_file)
        if df is not None:
            # First column is usually gene names
            stats['features_per_type']['Transcriptomics (genes)'] = {
                'input': len(df),
                'file': os.path.basename(expr_file)
            }
            # Count samples (columns minus gene column)
            stats['samples'] = len(df.columns) - 1
            print(f"[DEBUG] Expression file found: {expr_file}", file=sys.stderr)
            print(f"[DEBUG] DataFrame shape: {df.shape}, Columns: {len(df.columns)}, Samples: {stats['samples']}", file=sys.stderr)
    
    # Find metadata file - include DESeq_groups.txt from preprocessing
    meta_patterns = [
        'DESeq_groups.txt',  # Staged from preprocessing (highest priority)
        '*metadata*.txt', '*Metadata*.txt', 'metadata.txt', '*metadata*.tsv'
    ]
    meta_file = None
    for pattern in meta_patterns:
        for search_dir in search_dirs:
            if not os.path.exists(search_dir):
                continue
            matches = glob.glob(os.path.join(search_dir, pattern))
            if matches:
                meta_file = matches[0]
                break
        if meta_file:
            break
    
    if meta_file:
        stats['files_found'].append(('Metadata', os.path.basename(meta_file)))
        
        df = safe_read_file(meta_file)
        if df is not None:
            stats['metadata_columns'] = list(df.columns)
            if stats['samples'] == 0:
                stats['samples'] = len(df)
    
    # Find regulator data files
    for config in regulator_configs:
        prefix = config['prefix']
        data_file = config['data_file']
        
        # Handle TF lists - these are reference files, not user input data files
        # Do not count them in the "Input Files" statistic
        if 'TF' in prefix.upper() and 'list' in data_file.lower():
            # Count TFs in the list for feature statistics only
            tf_file = os.path.join(data_dir, data_file)
            if not os.path.exists(tf_file):
                # Check PKN directory
                tf_file = os.path.join(input_dir, 'PKN', data_file)
            
            if os.path.exists(tf_file):
                n_tfs = count_lines(tf_file)
                stats['features_per_type'][f'{prefix} (regulators)'] = {
                    'input': n_tfs,
                    'file': data_file
                }
            continue
        
        # Look for omics data file
        omics_file = None
        # Search in all directories for the omics file
        for search_dir in search_dirs:
            potential_file = os.path.join(search_dir, data_file)
            if os.path.exists(potential_file):
                omics_file = potential_file
                break
        
        if omics_file:
            stats['data_types'].append(prefix)
            stats['files_found'].append((prefix, data_file))
            
            df = safe_read_file(omics_file)
            if df is not None:
                # First column is usually feature names
                stats['features_per_type'][f'{prefix} (features)'] = {
                    'input': len(df),
                    'file': data_file
                }
        else:
            # Even if the physical file isn't found (e.g. running in container),
            # the regulator_types config tells us this omics data was provided as input
            stats['data_types'].append(prefix)
            stats['files_found'].append((prefix, data_file))
            
            # Fallback: try to count features from preprocessed file in staged Preprocessing dir
            # Preprocessing creates LemonPreprocessed_{prefix_lower}.txt for each omics type
            preproc_fallback_name = f'LemonPreprocessed_{prefix.lower()}.txt'
            preproc_search_dirs = [
                os.path.join(output_dir, 'LemonTree', 'Preprocessing'),
                os.path.join('.', 'LemonTree', 'Preprocessing')
            ]
            for pdir in preproc_search_dirs:
                preproc_path = os.path.join(pdir, preproc_fallback_name)
                if os.path.exists(preproc_path):
                    df = safe_read_file(preproc_path)
                    if df is not None:
                        stats['features_per_type'][f'{prefix} (features)'] = {
                            'input': len(df),
                            'file': data_file
                        }
                        print(f"[DEBUG] Fallback: counted {len(df)} features for {prefix} from {preproc_path}", file=sys.stderr)
                    break
    
    return stats


def collect_preprocessing_statistics(output_dir, run_id):
    """Collect statistics from preprocessing outputs"""
    stats = {
        'genes_after_hvg': 0,
        'genes_retained': 0,
        'scaling_methods': {},
        'tfa_performed': False,
        'deseq_groups': [],
        'pca_variance': None,
        'input_features': {}  # Track input features per omics type
    }
    
    # Find preprocessing directory - when run_id=".", look directly in output_dir
    if run_id == '.':
        preproc_dir = os.path.join(output_dir, 'LemonTree', 'Preprocessing')
    else:
        preproc_dir = os.path.join(output_dir, run_id, 'LemonTree', 'Preprocessing')
    if not os.path.exists(preproc_dir):
        preproc_dir = os.path.join(output_dir, 'LemonTree', 'Preprocessing')
    
    # Read preprocessed expression file
    expr_file = os.path.join(preproc_dir, 'LemonPreprocessed_expression.txt')
    if os.path.exists(expr_file):
        df = safe_read_file(expr_file)
        if df is not None:
            # First column is 'symbol', rest are samples
            # Subtract 1 for header line
            stats['genes_retained'] = len(df)
            stats['genes_after_hvg'] = len(df)
            stats['input_features']['Genes'] = len(df)
    
    # Read TF activity file
    tfa_file = os.path.join(preproc_dir, 'tfs.txt')
    if os.path.exists(tfa_file):
        n_tfs = count_lines(tfa_file) - 1  # Subtract header
        stats['input_features']['TFs'] = n_tfs
        stats['tfa_performed'] = True
    
    # Read metabolites file
    metab_file = os.path.join(preproc_dir, 'metabolites.txt')
    if os.path.exists(metab_file):
        n_metab = count_lines(metab_file) - 1  # Subtract header
        stats['input_features']['Metabolites'] = n_metab
    
    # Read proteins file
    prot_file = os.path.join(preproc_dir, 'proteins.txt')
    if os.path.exists(prot_file):
        n_prot = count_lines(prot_file) - 1  # Subtract header
        stats['input_features']['Proteins'] = n_prot
    
    # Read complete preprocessed file (includes all omics)
    complete_file = os.path.join(preproc_dir, 'LemonPreprocessed_complete.txt')
    if os.path.exists(complete_file):
        df = safe_read_file(complete_file)
        if df is not None:
            stats['total_features_integrated'] = len(df)
    
    # Read DESeq groups
    deseq_file = os.path.join(preproc_dir, 'DESeq_groups.txt')
    if os.path.exists(deseq_file):
        df = safe_read_file(deseq_file)
        if df is not None:
            stats['deseq_groups'] = list(df.columns)
            stats['n_samples'] = len(df)
            # Get group counts
            stats['group_counts'] = {}
            for col in df.columns:
                stats['group_counts'][col] = df[col].value_counts().to_dict()
    
    # Check for TFA results
    tfa_dir = os.path.join(output_dir, run_id, 'TFA')
    if not os.path.exists(tfa_dir):
        tfa_dir = os.path.join(output_dir, 'TFA')
    
    if os.path.exists(tfa_dir) and len(os.listdir(tfa_dir)) > 0:
        stats['tfa_performed'] = True
        tfa_files = glob.glob(os.path.join(tfa_dir, '*TFA*.txt'))
        stats['tfa_files'] = len(tfa_files)
    
    return stats


def collect_clustering_statistics(output_dir, run_id):
    """Collect statistics from LemonTree clustering"""
    stats = {
        'n_clustering_runs': 0,
        'n_tight_clusters': 0,
        'cluster_sizes': []
    }
    
    # Find Lemon_out directory with multiple search paths - when run_id=".", look directly in output_dir
    if run_id == '.':
        lemon_out_paths = [
            os.path.join(output_dir, 'LemonTree', 'Lemon_out'),
            os.path.join(output_dir, 'Lemon_out')
        ]
    else:
        lemon_out_paths = [
            os.path.join(output_dir, run_id, 'LemonTree', 'Lemon_out'),
            os.path.join(output_dir, 'LemonTree', 'Lemon_out'),
            os.path.join(output_dir, run_id, 'Lemon_out')
        ]
    
    lemon_out = None
    for path in lemon_out_paths:
        if os.path.exists(path):
            lemon_out = path
            break
    
    if lemon_out:
        # Count clustering runs - check for both cluster_* and Cluster_* patterns
        if os.path.exists(lemon_out):
            # Try Lemon_results subdirectory first
            lemon_results = os.path.join(lemon_out, 'Lemon_results')
            if os.path.exists(lemon_results):
                cluster_dirs = glob.glob(os.path.join(lemon_results, '[Cc]luster_*'))
                stats['n_clustering_runs'] = len(cluster_dirs)
            else:
                # Check directly in lemon_out
                cluster_dirs = glob.glob(os.path.join(lemon_out, '[Cc]luster_*'))
                stats['n_clustering_runs'] = len(cluster_dirs)
        
        # Read tight clusters - count unique module IDs directly to avoid safe_read_file nrows limit
        tight_clusters_file = os.path.join(lemon_out, 'tight_clusters.txt')
        if os.path.exists(tight_clusters_file):
            try:
                module_gene_map = defaultdict(list)
                with open(tight_clusters_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            gene, mod_id = parts[0], parts[1]
                            module_gene_map[mod_id].append(gene)
                if module_gene_map:
                    module_sizes_series = pd.Series({k: len(v) for k, v in module_gene_map.items()})
                    stats['n_tight_clusters'] = len(module_sizes_series)
                    stats['cluster_sizes'] = module_sizes_series.tolist()
                    stats['mean_cluster_size'] = module_sizes_series.mean()
                    stats['median_cluster_size'] = module_sizes_series.median()
                    stats['min_cluster_size'] = module_sizes_series.min()
                    stats['max_cluster_size'] = module_sizes_series.max()
            except Exception as e:
                print(f"Warning: Could not read tight_clusters.txt: {e}")
    
    # Get initial module count from Module_coherence_scores.txt (before filtering)
    networks_paths = [
        os.path.join(output_dir, run_id, 'LemonTree', 'Networks'),
        os.path.join(output_dir, 'LemonTree', 'Networks')
    ]
    for networks_dir in networks_paths:
        coherence_file = os.path.join(networks_dir, 'Module_coherence_scores.txt')
        if os.path.exists(coherence_file):
            # Subtract 1 for header
            stats['n_initial_modules'] = count_lines(coherence_file) - 1
            break
    
    # Fallback: Try to get stats from Module_Overview.csv
    if stats['n_tight_clusters'] == 0:
        overview_paths = [
            os.path.join(output_dir, run_id, 'LemonTree', 'Module_Overview', 'Module_Overview.csv'),
            os.path.join(output_dir, 'LemonTree', 'Module_Overview', 'Module_Overview.csv')
        ]
        for overview_file in overview_paths:
            if os.path.exists(overview_file):
                df = safe_read_file(overview_file, sep='\t')  # Tab-separated
                if df is not None and 'Module' in df.columns:
                    stats['n_tight_clusters'] = len(df)
                    # Get module sizes from Module_genes column
                    if 'Module_genes' in df.columns:
                        sizes = df['Module_genes'].apply(lambda x: len(str(x).split('|')) if pd.notna(x) else 0)
                        stats['cluster_sizes'] = sizes.tolist()
                        stats['mean_cluster_size'] = sizes.mean()
                        stats['median_cluster_size'] = sizes.median()
                        stats['min_cluster_size'] = sizes.min()
                        stats['max_cluster_size'] = sizes.max()
                    break
    
    return stats


def collect_network_statistics(output_dir, run_id, regulator_configs):
    """Collect statistics from network generation"""
    stats = {
        'modules_total': 0,
        'modules_after_coherence': 0,
        'coherence_scores': {},
        'regulators_per_type': {},
        'network_edges': 0,
        'genes_in_network': 0,
        'de_modules_count': 0,
        'ppi_enriched_modules_count': 0
    }
    
    # Find Module_Overview directory (primary source) - when run_id=".", look directly in output_dir
    if run_id == '.':
        overview_paths = [
            os.path.join(output_dir, 'LemonTree', 'Module_Overview'),
            os.path.join(output_dir, 'Module_Overview')
        ]
    else:
        overview_paths = [
            os.path.join(output_dir, run_id, 'LemonTree', 'Module_Overview'),
            os.path.join(output_dir, 'LemonTree', 'Module_Overview')
        ]
    
    overview_dir = None
    for path in overview_paths:
        if os.path.exists(path):
            overview_dir = path
            break
    
    # Read Module_Overview.csv (primary source for all stats)
    if overview_dir:
        module_overview_file = os.path.join(overview_dir, 'Module_Overview.csv')
        if os.path.exists(module_overview_file):
            df = safe_read_file(module_overview_file, sep='\t')  # Tab-separated
            if df is not None:
                # Total modules
                if 'Module' in df.columns:
                    stats['modules_after_coherence'] = len(df)
                
                # Coherence scores
                if 'Coherence' in df.columns:
                    coherence_values = pd.to_numeric(df['Coherence'], errors='coerce').dropna()
                    if len(coherence_values) > 0:
                        stats['coherence_scores'] = dict(zip(df['Module'].astype(str), df['Coherence']))
                        stats['mean_coherence'] = coherence_values.mean()
                        stats['median_coherence'] = coherence_values.median()
                
                # Count genes in network from Module_genes column
                if 'Module_genes' in df.columns:
                    all_genes = set()
                    for genes_str in df['Module_genes']:
                        if pd.notna(genes_str):
                            genes = str(genes_str).split('|')
                            all_genes.update([g.strip() for g in genes if g.strip()])
                    stats['genes_in_network'] = len(all_genes)
                
                # Count PPI enriched modules
                if 'PPI_FDR' in df.columns:
                    ppi_fdr = pd.to_numeric(df['PPI_FDR'], errors='coerce')
                    stats['ppi_enriched_modules_count'] = int((ppi_fdr < 0.05).sum())
        
        # Read DE analysis results - but only count modules that are in Module_Overview
        de_analysis_file = os.path.join(overview_dir, 'module_expression_analysis.csv')
        if os.path.exists(de_analysis_file) and 'Module' in df.columns:
            de_df = safe_read_file(de_analysis_file, sep=',')  # Actually comma-separated
            if de_df is not None and 'Module' in de_df.columns:
                # Get modules that are in the final filtered set
                final_modules = set(df['Module'].astype(str))
                
                # Look for adjusted p-value column
                pval_cols = [c for c in de_df.columns if 'adjusted' in c.lower() and 'pval' in c.lower()]
                if not pval_cols:
                    pval_cols = [c for c in de_df.columns if 'p_adjusted' in c.lower()]
                if pval_cols:
                    pval_col = pval_cols[0]
                    # Filter to only final modules
                    de_df_filtered = de_df[de_df['Module'].astype(str).isin(final_modules)]
                    pvals = pd.to_numeric(de_df_filtered[pval_col], errors='coerce')
                    stats['de_modules_count'] = int((pvals < 0.05).sum())
    
    # Find Networks directory (for network file) - when run_id=".", look directly in output_dir
    if run_id == '.':
        networks_paths = [
            os.path.join(output_dir, 'LemonTree', 'Networks'),
            os.path.join(output_dir, 'Networks')
        ]
    else:
        networks_paths = [
            os.path.join(output_dir, run_id, 'LemonTree', 'Networks'),
            os.path.join(output_dir, 'LemonTree', 'Networks')
        ]
    
    networks_dir = None
    for path in networks_paths:
        if os.path.exists(path):
            networks_dir = path
            break
    
    # Read initial module count from Module_coherence_scores.txt (before coherence filtering)
    if networks_dir and stats['modules_total'] == 0:
        coherence_file = os.path.join(networks_dir, 'Module_coherence_scores.txt')
        if os.path.exists(coherence_file):
            stats['modules_total'] = count_lines(coherence_file) - 1  # Subtract header
    
    # Read network edges from Cytoscape network file or LemonNetwork file
    if networks_dir:
        # Try Cytoscape network file first
        network_files = glob.glob(os.path.join(networks_dir, 'Cytoscape_network_*.txt'))
        if not network_files:
            network_files = glob.glob(os.path.join(networks_dir, 'LemonNetwork_*.txt'))
        
        if network_files:
            df = safe_read_file(network_files[0])
            if df is not None:
                stats['network_edges'] = len(df)
        
        # Fallback: Read specific modules file
        if stats['modules_after_coherence'] == 0:
            specific_modules_file = os.path.join(networks_dir, 'specific_modules.txt')
            if os.path.exists(specific_modules_file):
                with open(specific_modules_file, 'r') as f:
                    modules = [line.strip() for line in f if line.strip()]
                    stats['modules_after_coherence'] = len(modules)
                    stats['filtered_modules'] = modules
    
    # Find ModuleViewer_files directory (for genes/regulators if needed)
    viewer_paths = [
        os.path.join(output_dir, run_id, 'LemonTree', 'ModuleViewer_files'),
        os.path.join(output_dir, 'LemonTree', 'ModuleViewer_files')
    ]
    
    viewer_dir = None
    for path in viewer_paths:
        if os.path.exists(path):
            viewer_dir = path
            break
    
    # Read clusters list for gene counts (fallback if not in Module_Overview)
    if stats['genes_in_network'] == 0 and viewer_dir:
        clusters_file = os.path.join(viewer_dir, 'clusters_list.txt')
        if os.path.exists(clusters_file):
            df = safe_read_file(clusters_file, header=None)
            if df is not None and len(df.columns) >= 2:
                all_genes = set()
                for genes_str in df[1]:
                    if pd.notna(genes_str):
                        genes = str(genes_str).split('|')
                        all_genes.update([g.strip() for g in genes if g.strip()])
                stats['genes_in_network'] = len(all_genes)
    
    # Count regulators per type from Module_Overview.csv if available
    if overview_dir:
        module_overview_file = os.path.join(overview_dir, 'Module_Overview.csv')
        if os.path.exists(module_overview_file):
            df = safe_read_file(module_overview_file, sep='\t')  # Tab-separated
            if df is not None:
                # Count TFs, Metabolites, Proteins regulators
                for reg_type in ['TFs', 'Metabolites', 'Proteins']:
                    col_name = f'{reg_type}_regulators'
                    if col_name in df.columns:
                        all_regs = set()
                        modules_with_regs = 0
                        for regs_str in df[col_name]:
                            if pd.notna(regs_str) and str(regs_str).strip() and str(regs_str) not in ['nan', 'NA', '']:
                                regs = str(regs_str).split('|')
                                all_regs.update([r.strip() for r in regs if r.strip()])
                                modules_with_regs += 1
                        if len(all_regs) > 0:
                            stats['regulators_per_type'][reg_type] = {
                                'selected': len(all_regs),
                                'modules_with_regulators': modules_with_regs
                            }
    
    # Fallback: Count regulators from viewer_dir files
    if not stats['regulators_per_type'] and viewer_dir:
        for config in regulator_configs:
            prefix = config['prefix']
            reg_file = os.path.join(viewer_dir, f'{prefix}.selected_regs_list.txt')
            if os.path.exists(reg_file):
                df = safe_read_file(reg_file, header=None)
                if df is not None and len(df.columns) >= 2:
                    all_regs = set()
                    for regs_str in df[1]:
                        if pd.notna(regs_str):
                            regs = str(regs_str).split('|')
                            all_regs.update([r.strip() for r in regs if r.strip()])
                    stats['regulators_per_type'][prefix] = {
                        'selected': len(all_regs),
                        'modules_with_regulators': len(df)
                    }
    
    return stats


def collect_pkn_evaluation_statistics(output_dir, run_id):
    """Collect statistics from PKN evaluation"""
    stats = {
        'pkn_evaluation_performed': False,
        'metrics': {},
        'ppi_enrichment': {}
    }
    
    # Find PKN_Evaluation directory - check multiple locations, when run_id=".", look directly in output_dir
    if run_id == '.':
        pkn_paths = [
            os.path.join(output_dir, 'LemonTree', 'PKN_Evaluation'),
            os.path.join(output_dir, 'PKN_Evaluation')
        ]
    else:
        pkn_paths = [
            os.path.join(output_dir, run_id, 'LemonTree', 'PKN_Evaluation'),
            os.path.join(output_dir, 'LemonTree', 'PKN_Evaluation'),
            os.path.join(output_dir, run_id, 'PKN_Evaluation'),
            os.path.join(output_dir, 'PKN_Evaluation')
        ]
    
    pkn_dir = None
    for path in pkn_paths:
        if os.path.exists(path):
            pkn_dir = path
            break
    
    if pkn_dir and os.path.exists(pkn_dir):
        eval_files = glob.glob(os.path.join(pkn_dir, '*evaluation*.txt'))
        if eval_files:
            stats['pkn_evaluation_performed'] = True
            
            # Try to parse evaluation metrics
            for eval_file in eval_files:
                df = safe_read_file(eval_file)
                if df is not None:
                    stats['evaluation_file'] = os.path.basename(eval_file)
                    # Store basic stats
                    if 'precision' in df.columns.str.lower():
                        stats['metrics']['precision'] = df['precision'].mean() if 'precision' in df.columns else None
                    if 'recall' in df.columns.str.lower():
                        stats['metrics']['recall'] = df['recall'].mean() if 'recall' in df.columns else None
    
    # Also check for enrichment results in ModuleViewer_files as indicator of PKN evaluation
    if not stats['pkn_evaluation_performed']:
        if run_id == '.':
            viewer_paths = [
                os.path.join(output_dir, 'LemonTree', 'ModuleViewer_files'),
                os.path.join(output_dir, 'ModuleViewer_files')
            ]
        else:
            viewer_paths = [
                os.path.join(output_dir, run_id, 'LemonTree', 'ModuleViewer_files'),
                os.path.join(output_dir, 'LemonTree', 'ModuleViewer_files'),
                os.path.join(output_dir, run_id, 'ModuleViewer_files'),
                os.path.join(output_dir, 'ModuleViewer_files')
            ]
        
        for viewer_dir in viewer_paths:
            if os.path.exists(viewer_dir):
                viewer_files = os.listdir(viewer_dir)
                # Check for PPI or metabolite enrichment files
                if any(f in viewer_files for f in ['PPI_enrichment_results.csv', 'Metabolite_Gene_enrichment_results.csv']):
                    stats['pkn_evaluation_performed'] = True
                    break
    
    return stats


def collect_enrichment_statistics(output_dir, run_id):
    """Collect statistics from enrichment analysis - memory optimized"""
    stats = {
        'enrichment_performed': False,
        'method': None,
        'modules_with_enrichment': 0,
        'top_pathways': [],
        'databases_used': []
    }
    
    # Find Enrichment directory - search only in top level, when run_id=".", look directly in output_dir
    if run_id == '.':
        enrich_paths = [
            os.path.join(output_dir, 'LemonTree', 'Enrichment'),
            os.path.join(output_dir, 'Enrichment'),
        ]
    else:
        enrich_paths = [
            os.path.join(output_dir, run_id, 'LemonTree', 'Enrichment'),
            os.path.join(output_dir, 'LemonTree', 'Enrichment'),
            os.path.join(output_dir, 'Enrichment'),
        ]
    
    enrich_dir = None
    for path in enrich_paths:
        if os.path.exists(path):
            enrich_dir = path
            break
    
    if enrich_dir:
        stats['enrichment_performed'] = True
        
        # Collect enrichment data across ALL method directories (EnrichR + GSEA)
        csv_count = 0
        all_modules = set()
        all_terms = {}  # term -> count
        for root, dirs, files in os.walk(enrich_dir):
            # Allow up to 3 levels deep to find CSVs in Enrichment/Modules_enrichr/*.csv
            depth = root[len(enrich_dir):].count(os.sep)
            if depth > 2:
                dirs[:] = []  # Don't recurse deeper
                continue
            
            csv_files = [f for f in files if f.endswith('.csv')]
            csv_count += len(csv_files)
            
            if csv_files:
                # Read all summary/top files to aggregate across EnrichR and GSEA
                priority_files = [f for f in csv_files if 'summary' in f.lower() or 'top_' in f.lower()]
                files_to_read = priority_files if priority_files else csv_files[:1]
                
                for target_file in files_to_read:
                    filepath = os.path.join(root, target_file)
                    try:
                        # Read only specific columns to save memory
                        df = pd.read_csv(filepath, sep=',', usecols=lambda x: x in ['Module', 'Term', 'module', 'term'], nrows=10000)
                        if df is not None and len(df) > 0:
                            # Collect unique modules
                            mod_col = 'Module' if 'Module' in df.columns else ('module' if 'module' in df.columns else None)
                            if mod_col:
                                all_modules.update(df[mod_col].dropna().unique())
                            
                            # Collect term frequencies
                            term_col = 'Term' if 'Term' in df.columns else ('term' if 'term' in df.columns else None)
                            if term_col:
                                for term, count in df[term_col].value_counts().items():
                                    all_terms[term] = all_terms.get(term, 0) + count
                            
                            del df
                            gc.collect()
                    except Exception as e:
                        print(f"Warning: Could not read enrichment file {filepath}: {e}")
        
        if all_modules:
            stats['modules_with_enrichment'] = len(all_modules)
        if all_terms:
            sorted_terms = sorted(all_terms.items(), key=lambda x: x[1], reverse=True)[:20]
            stats['top_pathways'] = dict(sorted_terms)
        
        # Detect enrichment method from directory structure
        all_items = os.listdir(enrich_dir)
        has_enrichr = any('enrichr' in str(d).lower() for d in all_items if os.path.isdir(os.path.join(enrich_dir, d)))
        has_gsea = any('gsea' in str(d).lower() for d in all_items if os.path.isdir(os.path.join(enrich_dir, d)))
        if has_enrichr and has_gsea:
            stats['method'] = 'EnrichR + GSEA'
        elif has_enrichr or 'EnrichR' in str(enrich_dir):
            stats['method'] = 'EnrichR'
        elif has_gsea or 'GSEA' in str(enrich_dir):
            stats['method'] = 'GSEA'
        else:
            stats['method'] = 'EnrichR/GSEA' if csv_count > 0 else None
    
    # Fallback: if no enrichment CSVs were found or modules_with_enrichment is 0,
    # try to count modules with pathway data from Module_Overview.csv
    if stats['modules_with_enrichment'] == 0:
        if run_id == '.':
            overview_paths = [
                os.path.join(output_dir, 'LemonTree', 'Module_Overview', 'Module_Overview.csv'),
                os.path.join(output_dir, 'Module_Overview', 'Module_Overview.csv'),
            ]
        else:
            overview_paths = [
                os.path.join(output_dir, run_id, 'LemonTree', 'Module_Overview', 'Module_Overview.csv'),
                os.path.join(output_dir, 'LemonTree', 'Module_Overview', 'Module_Overview.csv'),
            ]
        for overview_file in overview_paths:
            if os.path.exists(overview_file):
                try:
                    df = safe_read_file(overview_file, sep='\t', nrows=500)
                    if df is not None:
                        pathway_cols = [c for c in df.columns if 'Top_3_pathways' in c]
                        if pathway_cols:
                            has_enrichment = df[pathway_cols].apply(
                                lambda row: any(pd.notna(v) and str(v).strip() != '' for v in row), axis=1
                            )
                            stats['modules_with_enrichment'] = int(has_enrichment.sum())
                            if stats['modules_with_enrichment'] > 0:
                                stats['enrichment_performed'] = True
                    del df
                    gc.collect()
                except Exception as e:
                    print(f"Warning: Could not read Module_Overview for enrichment fallback: {e}")
                break
    
    return stats


def collect_module_rankings(output_dir, run_id):
    """Collect module rankings from various sources"""
    rankings = {
        'by_coherence': [],
        'by_de': [],
        'by_size': [],
        'by_regulator_count': [],
        'by_ppi_enrichment': [],
        'by_metgene_enrichment': [],
        'de_modules_count': 0,
        'ppi_enriched_modules_count': 0,
        'metgene_enriched_modules_count': 0,
        'modules_with_enrichment_count': 0
    }
    
    # Find Networks directory, when run_id=".", look directly in output_dir
    if run_id == '.':
        networks_dir = os.path.join(output_dir, 'LemonTree', 'Networks')
        if not os.path.exists(networks_dir):
            networks_dir = os.path.join(output_dir, 'Networks')
    else:
        networks_dir = os.path.join(output_dir, run_id, 'LemonTree', 'Networks')
        if not os.path.exists(networks_dir):
            networks_dir = os.path.join(output_dir, 'LemonTree', 'Networks')
    
    # Find Module_Overview directory, when run_id=".", look directly in output_dir
    if run_id == '.':
        overview_dir = os.path.join(output_dir, 'LemonTree', 'Module_Overview')
        if not os.path.exists(overview_dir):
            overview_dir = os.path.join(output_dir, 'Module_Overview')
    else:
        overview_dir = os.path.join(output_dir, run_id, 'LemonTree', 'Module_Overview')
        if not os.path.exists(overview_dir):
            overview_dir = os.path.join(output_dir, 'LemonTree', 'Module_Overview')
    
    # Find ModuleViewer_files directory, when run_id=".", look directly in output_dir
    if run_id == '.':
        viewer_dir = os.path.join(output_dir, 'LemonTree', 'ModuleViewer_files')
        if not os.path.exists(viewer_dir):
            viewer_dir = os.path.join(output_dir, 'ModuleViewer_files')
    else:
        viewer_dir = os.path.join(output_dir, run_id, 'LemonTree', 'ModuleViewer_files')
        if not os.path.exists(viewer_dir):
            viewer_dir = os.path.join(output_dir, 'LemonTree', 'ModuleViewer_files')
    
    # Read coherence scores from Module_Overview.csv first
    if overview_dir:
        overview_file = os.path.join(overview_dir, 'Module_Overview.csv')
        if os.path.exists(overview_file):
            df = safe_read_file(overview_file, sep='\t', nrows=500)  # Limit rows to save memory
            if df is not None and 'Module' in df.columns and 'Coherence' in df.columns:
                df_sorted = df.sort_values('Coherence', ascending=False)
                # Convert only top 20 to dict to save memory
                rankings['by_coherence'] = df_sorted.head(20)[['Module', 'Coherence']].rename(columns={'Coherence': 'Coherence_Score'}).to_dict('records')
                del df_sorted
                gc.collect()
            del df
            gc.collect()
    
    # Fallback: Try networks_dir
    if not rankings['by_coherence'] and networks_dir:
        coherence_file = os.path.join(networks_dir, 'Module_coherence_scores.txt')
        if os.path.exists(coherence_file):
            df = safe_read_file(coherence_file)
            if df is not None and 'Module' in df.columns and 'Coherence_Score' in df.columns:
                df_sorted = df.sort_values('Coherence_Score', ascending=False)
                rankings['by_coherence'] = df_sorted.head(20).to_dict('records')
    
    # Read Module_Overview.csv - single read, multiple uses
    summary_file = os.path.join(overview_dir, 'Module_Overview.csv')
    if os.path.exists(summary_file):
        df = safe_read_file(summary_file, sep='\t', nrows=500)  # Limit rows to save memory
        if df is not None:
            # Size rankings - use Module_genes column
            if 'Module_genes' in df.columns:
                df['n_genes'] = df['Module_genes'].apply(lambda x: len(str(x).split('|')) if pd.notna(x) else 0)
                df_sorted = df.sort_values('n_genes', ascending=False)
                rankings['by_size'] = df_sorted.head(20)[['Module', 'n_genes', 'Coherence']].to_dict('records')
                del df_sorted
                gc.collect()
            
            # PPI enrichment rankings - use PPI_fold_enrichment as the display score
            if 'PPI_fold_enrichment' in df.columns:
                df_numeric = df.copy()
                df_numeric['PPI_fold_enrichment'] = pd.to_numeric(df_numeric['PPI_fold_enrichment'], errors='coerce')
                df_numeric['PPI_FDR'] = pd.to_numeric(df_numeric.get('PPI_FDR', pd.NA), errors='coerce')
                # Sort by fold enrichment (higher is better)
                df_sorted = df_numeric.dropna(subset=['PPI_fold_enrichment']).sort_values('PPI_fold_enrichment', ascending=False)
                # Use PPI_fold_enrichment as the display column, rename to ppi_pvalue for consistency
                display_df = df_sorted.head(20)[['Module', 'PPI_fold_enrichment', 'PPI_FDR']].copy()
                display_df = display_df.rename(columns={'PPI_fold_enrichment': 'ppi_pvalue'})
                rankings['by_ppi_enrichment'] = display_df.to_dict('records')
                
                # Count modules with significant PPI enrichment (FDR < 0.05)
                rankings['ppi_enriched_modules_count'] = int((df_numeric['PPI_FDR'] < 0.05).sum())
                del df_numeric, df_sorted, display_df
                gc.collect()
            
            # Metabolite-gene interaction enrichment rankings - show adjusted p-value (FDR)
            if 'MetGene_FDR' in df.columns:
                df_numeric = df.copy()
                df_numeric['MetGene_FDR'] = pd.to_numeric(df_numeric['MetGene_FDR'], errors='coerce')
                # Sort by FDR ascending (best adjusted p-values first)
                df_sorted = df_numeric.dropna(subset=['MetGene_FDR']).sort_values('MetGene_FDR', ascending=True)
                display_df = df_sorted.head(20)[['Module', 'MetGene_FDR']].copy()
                display_df = display_df.rename(columns={'MetGene_FDR': 'metgene_adjusted_pvalue'})
                rankings['by_metgene_enrichment'] = display_df.to_dict('records')
                
                # Count modules with significant metabolite-gene enrichment (FDR < 0.05)
                rankings['metgene_enriched_modules_count'] = int((df_numeric['MetGene_FDR'] < 0.05).sum())
                del df_numeric, df_sorted, display_df
                gc.collect()
            
            # DE rankings from Module_Overview.csv
            if 'Expression_adjusted_pval' in df.columns:
                df_numeric = df[[col for col in df.columns if col in ['Module', 'Expression_adjusted_pval']]].copy()
                df_numeric['Expression_adjusted_pval'] = pd.to_numeric(df_numeric['Expression_adjusted_pval'], errors='coerce')
                df_sorted = df_numeric.dropna(subset=['Expression_adjusted_pval']).sort_values('Expression_adjusted_pval', ascending=True)
                
                # Rename Expression_adjusted_pval to de_pvalue for display consistency
                display_df = df_sorted.head(20)[['Module', 'Expression_adjusted_pval']].copy()
                display_df = display_df.rename(columns={'Expression_adjusted_pval': 'de_pvalue'})
                rankings['by_de'] = display_df.to_dict('records')
                
                # Count significant DE modules (p < 0.05)
                rankings['de_modules_count'] = int((df_numeric['Expression_adjusted_pval'] < 0.05).sum())
                del df_numeric, df_sorted, display_df
                gc.collect()
            
            del df
            gc.collect()
    
    # Fallback: Try module_expression_analysis.csv if available
    if not rankings.get('by_de'):
        de_analysis_file = os.path.join(overview_dir, 'module_expression_analysis.csv')
        if os.path.exists(de_analysis_file):
            df = safe_read_file(de_analysis_file, sep=',')
            if df is not None:
                # DE rankings - look for adjusted p-value columns
                pval_cols = [c for c in df.columns if 'adjusted' in c.lower() and 'pval' in c.lower()]
                if pval_cols:
                    pval_col = pval_cols[0]
                    df_numeric = df.copy()
                    df_numeric[pval_col] = pd.to_numeric(df_numeric[pval_col], errors='coerce')
                    df_sorted = df_numeric.dropna(subset=[pval_col]).sort_values(pval_col, ascending=True)
                    
                    # Get relevant columns for display
                    display_cols = ['Module', pval_col]
                    # Add rank column if available
                    rank_cols = [c for c in df.columns if 'rank' in c.lower()]
                    if rank_cols:
                        display_cols.append(rank_cols[0])
                    
                    rankings['by_de'] = df_sorted.head(20)[display_cols].to_dict('records')
                    
                    # Count significant DE modules (p < 0.05)
                    rankings['de_modules_count'] = int((df_numeric[pval_col] < 0.05).sum())
    
    # Also try to find prioritization files (limit search scope)
    for parent_dir in [os.path.join(output_dir, run_id), os.path.join(output_dir, run_id, 'LemonTree'), output_dir]:
        if os.path.exists(parent_dir):
            try:
                for fname in os.listdir(parent_dir):
                    if 'prioritization' in fname.lower() and fname.endswith('.txt'):
                        pf = os.path.join(parent_dir, fname)
                        df = safe_read_file(pf)
                        if df is not None:
                            rankings['prioritization_data'] = df.head(20).to_dict('records')
                            break
            except:
                pass
    
    return rankings


def parse_parameters_file(params_file):
    """Parse the pipeline_parameters_log.txt file"""
    params = {}
    
    if not params_file or not os.path.exists(params_file):
        return params
    
    try:
        with open(params_file, 'r') as f:
            content = f.read()
        
        # Parse key-value pairs
        for line in content.split('\n'):
            if '=' in line and not line.strip().startswith('#'):
                parts = line.split('=', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    params[key] = value
    except Exception as e:
        print(f"Warning: Could not parse parameters file: {e}")
    
    return params


def collect_regulator_rankings(output_dir, run_id, regulator_configs):
    """Collect regulator-module pair rankings from selected_regulators_scores files.
    
    Reads {Type}.selected_regulators_scores.txt files from ModuleViewer_files/
    and returns ranked regulator-module pairs and per-regulator summaries,
    filtered to only include coherence-filtered modules.
    """
    regulator_rankings = {}
    
    # Find ModuleViewer_files directory
    if run_id == '.':
        viewer_candidates = [
            os.path.join(output_dir, 'LemonTree', 'ModuleViewer_files'),
            os.path.join(output_dir, 'ModuleViewer_files'),
        ]
    else:
        viewer_candidates = [
            os.path.join(output_dir, run_id, 'LemonTree', 'ModuleViewer_files'),
            os.path.join(output_dir, 'LemonTree', 'ModuleViewer_files'),
        ]
    
    viewer_dir = None
    for path in viewer_candidates:
        if os.path.exists(path):
            viewer_dir = path
            break
    
    if not viewer_dir:
        print("ModuleViewer_files directory not found — skipping regulator rankings")
        return regulator_rankings
    
    # Find coherence-filtered modules from Module_Overview.csv
    filtered_modules = None
    if run_id == '.':
        overview_candidates = [
            os.path.join(output_dir, 'LemonTree', 'Module_Overview', 'Module_Overview.csv'),
            os.path.join(output_dir, 'Module_Overview', 'Module_Overview.csv'),
        ]
    else:
        overview_candidates = [
            os.path.join(output_dir, run_id, 'LemonTree', 'Module_Overview', 'Module_Overview.csv'),
            os.path.join(output_dir, 'LemonTree', 'Module_Overview', 'Module_Overview.csv'),
        ]
    
    for ov_path in overview_candidates:
        if os.path.exists(ov_path):
            ov_df = safe_read_file(ov_path, sep='\t')
            if ov_df is not None and 'Module' in ov_df.columns:
                filtered_modules = set(ov_df['Module'].astype(str))
                print(f"Filtering regulator rankings to {len(filtered_modules)} coherence-filtered modules")
            break
    
    # Discover score files: {Type}.selected_regulators_scores.txt
    score_files = glob.glob(os.path.join(viewer_dir, '*.selected_regulators_scores.txt'))
    
    for score_file in score_files:
        basename = os.path.basename(score_file)
        # Extract type name: e.g. "TFs.selected_regulators_scores.txt" -> "TFs"
        reg_type = basename.replace('.selected_regulators_scores.txt', '')
        
        try:
            df = pd.read_csv(score_file, sep='\t')
            if df.empty:
                continue
            
            # Standardize column names
            if 'Target' in df.columns:
                df = df.rename(columns={'Target': 'Module'})
            df['Module'] = df['Module'].astype(str)
            
            # Filter to coherence-filtered modules only
            if filtered_modules is not None:
                df = df[df['Module'].isin(filtered_modules)]
            
            if df.empty:
                print(f"  {reg_type}: no regulator-module pairs after filtering to coherence modules")
                continue
            
            # Sort by score descending
            pairs_df = df.sort_values('Score', ascending=False).reset_index(drop=True)
            
            # Build per-regulator summary
            summary_rows = []
            for reg_name, grp in df.groupby('Regulator'):
                total_score = grp['Score'].sum()
                n_modules = grp['Module'].nunique()
                module_list = ', '.join(sorted(grp['Module'].unique(), key=lambda x: int(x) if x.isdigit() else x))
                summary_rows.append({
                    'Regulator': reg_name,
                    'Total_Score': total_score,
                    'N_Modules': n_modules,
                    'Target_Modules': module_list,
                })
            summary_df = pd.DataFrame(summary_rows).sort_values('Total_Score', ascending=False).reset_index(drop=True)
            
            regulator_rankings[reg_type] = {
                'pairs': pairs_df.to_dict('records'),
                'summary': summary_df.to_dict('records'),
            }
            print(f"  {reg_type}: {len(pairs_df)} regulator-module pairs, {len(summary_df)} unique regulators")
            
        except Exception as e:
            print(f"Warning: Could not load regulator scores from {score_file}: {e}")
    
    return regulator_rankings


def generate_regulator_ranking_section(regulator_rankings):
    """Generate HTML for regulator ranking tables section."""
    if not regulator_rankings:
        return '<p>No regulator score data available.</p>'
    
    parts = []
    
    for reg_type, data in regulator_rankings.items():
        pairs = data.get('pairs', [])
        summary = data.get('summary', [])
        safe_id = reg_type.replace(' ', '_')
        
        # --- Pairs table ---
        parts.append(f'<div class="subsection">')
        parts.append(f'<h3>🔬 {reg_type} — Regulator–Module Pairs (ranked by score)</h3>')
        parts.append(f'<input type="text" class="search-box" id="search_{safe_id}_pairs" '
                     f'onkeyup="filterRegTable(\'{safe_id}_pairs_table\', this.value)" '
                     f'placeholder="Search regulators or modules...">')
        parts.append(f'<div class="ranking-table" style="max-height:500px;overflow-y:auto;">')
        parts.append(f'<table id="{safe_id}_pairs_table">')
        parts.append('<thead><tr>'
                     '<th>Rank</th><th>Regulator</th><th>Module</th><th>Score</th><th>Overall Rank</th>'
                     '</tr></thead><tbody>')
        
        for rank, item in enumerate(pairs, 1):
            score = item.get('Score', 0)
            score_style = 'color:#27ae60;font-weight:600;' if score > 0 else 'color:#e74c3c;'
            overall_rank = item.get('Overall_rank', 'N/A')
            parts.append(
                f'<tr><td>{rank}</td>'
                f'<td>{item.get("Regulator", "N/A")}</td>'
                f'<td>Module {item.get("Module", "N/A")}</td>'
                f'<td style="{score_style}">{score:.0f}</td>'
                f'<td>{overall_rank}</td></tr>'
            )
        
        parts.append('</tbody></table></div>')
        
        # --- Summary table ---
        parts.append(f'<h3 style="margin-top:1.5rem;">📊 {reg_type} — Regulator Summary (aggregated across modules)</h3>')
        parts.append(f'<input type="text" class="search-box" id="search_{safe_id}_summary" '
                     f'onkeyup="filterRegTable(\'{safe_id}_summary_table\', this.value)" '
                     f'placeholder="Search regulators...">')
        parts.append(f'<div class="ranking-table" style="max-height:500px;overflow-y:auto;">')
        parts.append(f'<table id="{safe_id}_summary_table">')
        parts.append('<thead><tr>'
                     '<th>Rank</th><th>Regulator</th><th>Sum of Scores</th><th>N Target Modules</th><th>Target Modules</th>'
                     '</tr></thead><tbody>')
        
        for rank, item in enumerate(summary, 1):
            total_score = item.get('Total_Score', 0)
            score_style = 'color:#27ae60;font-weight:600;' if total_score > 0 else 'color:#e74c3c;'
            parts.append(
                f'<tr><td>{rank}</td>'
                f'<td>{item.get("Regulator", "N/A")}</td>'
                f'<td style="{score_style}">{total_score:.0f}</td>'
                f'<td>{item.get("N_Modules", 0)}</td>'
                f'<td style="font-size:0.9em;max-width:400px;word-wrap:break-word;">{item.get("Target_Modules", "")}</td></tr>'
            )
        
        parts.append('</tbody></table></div>')
        parts.append('</div>')  # close subsection
    
    return '\n'.join(parts)


def generate_html_report(all_stats, output_path):
    """Generate the comprehensive HTML report"""
    
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Lemonite Pipeline Summary Report</title>
    <style>
        :root {
            --primary-color: #2E7D32;
            --secondary-color: #81C784;
            --accent-color: #FDD835;
            --bg-color: #FAFAFA;
            --card-bg: #FFFFFF;
            --text-color: #333333;
            --text-light: #666666;
            --border-color: #E0E0E0;
            --success-color: #4CAF50;
            --warning-color: #FF9800;
            --error-color: #F44336;
            --info-color: #2196F3;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-color);
            line-height: 1.6;
        }
        
        /* Header */
        .header {
            background: linear-gradient(135deg, var(--primary-color) 0%, #1B5E20 100%);
            color: white;
            padding: 2rem 0;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        
        .header h1 {
            font-size: 2.5rem;
            margin-bottom: 0.5rem;
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 1rem;
        }
        
        .header .subtitle {
            font-size: 1.1rem;
            opacity: 0.9;
        }
        
        .header .run-info {
            margin-top: 1rem;
            font-size: 0.95rem;
            opacity: 0.85;
        }
        
        /* Navigation */
        .nav {
            background: var(--card-bg);
            padding: 1rem;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            display: flex;
            justify-content: center;
            flex-wrap: wrap;
            gap: 0.5rem;
        }
        
        .nav a {
            color: var(--primary-color);
            text-decoration: none;
            padding: 0.5rem 1rem;
            border-radius: 20px;
            transition: all 0.3s ease;
            font-weight: 500;
        }
        
        .nav a:hover {
            background: var(--secondary-color);
            color: white;
        }
        
        /* Main container */
        .container {
            max-width: 1400px;
            margin: 0 auto;
            padding: 2rem;
        }
        
        /* Section */
        .section {
            background: var(--card-bg);
            border-radius: 12px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }
        
        .section-header {
            display: flex;
            align-items: center;
            gap: 1rem;
            margin-bottom: 1.5rem;
            padding-bottom: 1rem;
            border-bottom: 2px solid var(--secondary-color);
        }
        
        .section-header h2 {
            color: var(--primary-color);
            font-size: 1.5rem;
        }
        
        .section-icon {
            font-size: 1.8rem;
        }
        
        /* Stats grid */
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1.5rem;
            margin: 1.5rem 0;
        }
        
        .stat-card {
            background: linear-gradient(135deg, var(--bg-color) 0%, #fff 100%);
            border-radius: 10px;
            padding: 1.5rem;
            text-align: center;
            border: 1px solid var(--border-color);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }
        
        .stat-card:hover {
            transform: translateY(-4px);
            box-shadow: 0 6px 12px rgba(0,0,0,0.1);
        }
        
        .stat-value {
            font-size: 2.5rem;
            font-weight: 700;
            color: var(--primary-color);
        }
        
        .stat-label {
            font-size: 0.9rem;
            color: var(--text-light);
            margin-top: 0.5rem;
        }
        
        .stat-detail {
            font-size: 0.8rem;
            color: var(--text-light);
            margin-top: 0.25rem;
        }
        
        /* Tables */
        .table-container {
            overflow-x: auto;
            margin: 1rem 0;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.95rem;
        }
        
        th, td {
            padding: 0.75rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }
        
        th {
            background: var(--primary-color);
            color: white;
            font-weight: 600;
            position: sticky;
            top: 0;
        }
        
        tr:nth-child(even) {
            background: var(--bg-color);
        }
        
        tr:hover {
            background: #E8F5E9;
        }
        
        /* Badges */
        .badge {
            display: inline-block;
            padding: 0.25rem 0.75rem;
            border-radius: 20px;
            font-size: 0.8rem;
            font-weight: 500;
        }
        
        .badge-success { background: #C8E6C9; color: #2E7D32; }
        .badge-warning { background: #FFE082; color: #F57F17; }
        .badge-info { background: #B3E5FC; color: #0277BD; }
        .badge-error { background: #FFCDD2; color: #C62828; }
        
        /* Progress bars */
        .progress-container {
            margin: 1rem 0;
        }
        
        .progress-label {
            display: flex;
            justify-content: space-between;
            margin-bottom: 0.5rem;
            font-size: 0.9rem;
        }
        
        .progress-bar {
            height: 12px;
            background: var(--border-color);
            border-radius: 6px;
            overflow: hidden;
        }
        
        .progress-fill {
            height: 100%;
            background: linear-gradient(90deg, var(--primary-color), var(--secondary-color));
            border-radius: 6px;
            transition: width 0.5s ease;
        }
        
        /* Info boxes */
        .info-box {
            padding: 1rem 1.5rem;
            border-radius: 8px;
            margin: 1rem 0;
            display: flex;
            align-items: flex-start;
            gap: 1rem;
        }
        
        .info-box-success { background: #E8F5E9; border-left: 4px solid var(--success-color); }
        .info-box-warning { background: #FFF8E1; border-left: 4px solid var(--warning-color); }
        .info-box-info { background: #E3F2FD; border-left: 4px solid var(--info-color); }
        .info-box-error { background: #FFEBEE; border-left: 4px solid var(--error-color); }
        
        /* Two column layout */
        .two-col {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 2rem;
        }
        
        @media (max-width: 900px) {
            .two-col {
                grid-template-columns: 1fr;
            }
        }
        
        /* Subsection */
        .subsection {
            margin: 1.5rem 0;
        }
        
        .subsection h3 {
            color: var(--primary-color);
            font-size: 1.1rem;
            margin-bottom: 1rem;
            padding-bottom: 0.5rem;
            border-bottom: 1px solid var(--border-color);
        }
        
        /* Parameter list */
        .param-list {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 0.5rem;
        }
        
        .param-item {
            display: flex;
            padding: 0.5rem;
            background: var(--bg-color);
            border-radius: 4px;
        }
        
        .param-name {
            font-weight: 600;
            min-width: 200px;
            color: var(--text-light);
        }
        
        .param-value {
            color: var(--primary-color);
            font-family: 'Consolas', monospace;
        }
        
        /* Collapsible */
        .collapsible {
            cursor: pointer;
            padding: 1rem;
            background: var(--bg-color);
            border: none;
            text-align: left;
            width: 100%;
            font-size: 1rem;
            border-radius: 8px;
            margin: 0.5rem 0;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        
        .collapsible:hover {
            background: #E8F5E9;
        }
        
        .collapsible-content {
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.3s ease;
        }
        
        .collapsible-content.active {
            max-height: 2000px;
        }
        
        /* Footer */
        .footer {
            text-align: center;
            padding: 2rem;
            color: var(--text-light);
            font-size: 0.9rem;
            border-top: 1px solid var(--border-color);
            margin-top: 2rem;
        }
        
        /* Charts placeholder */
        .chart-container {
            width: 100%;
            height: 300px;
            background: var(--bg-color);
            border-radius: 8px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin: 1rem 0;
        }
        
        /* Module ranking table */
        .ranking-table {
            max-height: 400px;
            overflow-y: auto;
        }
        
        .rank-badge {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 28px;
            height: 28px;
            border-radius: 50%;
            background: var(--primary-color);
            color: white;
            font-weight: 700;
            font-size: 0.85rem;
        }
        
        .rank-badge.gold { background: #FFD700; color: #333; }
        .rank-badge.silver { background: #C0C0C0; color: #333; }
        .rank-badge.bronze { background: #CD7F32; color: white; }
        
        /* Search box for regulator tables */
        .search-box {
            padding: 8px 12px;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            width: 300px;
            margin: 10px 0;
            font-size: 14px;
        }
        .search-box:focus {
            outline: none;
            border-color: var(--primary-color);
            box-shadow: 0 0 3px rgba(46,125,50,0.3);
        }
        
        /* Print styles */
        @media print {
            .nav { display: none; }
            .section { break-inside: avoid; }
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>🍋🌳 Lemonite Pipeline Summary Report</h1>
        <div class="subtitle">Comprehensive Multi-Omics Integration & Network Analysis</div>
        <div class="run-info">
            <strong>Run ID:</strong> {run_id} | 
            <strong>Generated:</strong> {timestamp}
        </div>
    </div>
    
    <nav class="nav">
        <a href="#summary-section">📋 Summary</a>
        <a href="#input-section">📊 Input Data</a>
        <a href="#preprocessing-section">⚙️ Preprocessing</a>
        <a href="#parameters-section">🎛️ Parameters</a>
        <a href="#clustering-section">🧩 Clustering</a>
        <a href="#network-section">🔗 Network</a>
        <a href="#evaluation-section">✅ Evaluation</a>
        <a href="#enrichment-section">🔬 Enrichment</a>
        <a href="#files-section">📂 Output Files</a>
        <a href="#rankings-section">🏆 Module Rankings</a>
        <a href="#regulator-rankings-section">🧬 Regulator Rankings</a>
    </nav>
    
    <div class="container">
        {content}
    </div>
    
    <div class="footer">
        <p>🍋🌳 <strong>Lemonite Pipeline</strong> | Developed by Boris Vandemoortele, CBIGR Lab @ Ghent University</p>
        <p>Report generated on {timestamp}</p>
    </div>
    
    <script>
        // Collapsible functionality
        document.querySelectorAll('.collapsible').forEach(button => {
            button.addEventListener('click', () => {
                const content = button.nextElementSibling;
                content.classList.toggle('active');
                button.querySelector('.toggle-icon').textContent = 
                    content.classList.contains('active') ? '▼' : '▶';
            });
        });
        
        // Smooth scroll for navigation
        document.querySelectorAll('.nav a').forEach(anchor => {
            anchor.addEventListener('click', function(e) {
                e.preventDefault();
                const target = document.querySelector(this.getAttribute('href'));
                target.scrollIntoView({ behavior: 'smooth', block: 'start' });
            });
        });
        
        // Filter function for regulator ranking tables
        function filterRegTable(tableId, query) {
            var table = document.getElementById(tableId);
            if (!table) return;
            var rows = table.querySelectorAll('tbody tr');
            var q = query.toLowerCase();
            rows.forEach(function(row) {
                row.style.display = row.textContent.toLowerCase().indexOf(q) > -1 ? '' : 'none';
            });
        }
    </script>
</body>
</html>"""
    
    # Build content sections
    content_parts = []
    
    # === EXECUTIVE SUMMARY SECTION ===
    input_stats = all_stats.get('input', {})
    preproc_stats = all_stats.get('preprocessing', {})
    network_stats = all_stats.get('network', {})
    rankings = all_stats.get('rankings', {})
    enrich_stats = all_stats.get('enrichment', {})
    
    # Count all omics data types including regulators (TFs, Metabolites, Proteins)
    omics_types = list(input_stats.get('data_types', []))  # Start with transcriptomics/metadata
    input_features = preproc_stats.get('input_features', {})
    if 'TFs' in input_features:
        omics_types.append('TFs')
    if 'Metabolites' in input_features:
        omics_types.append('Metabolites')
    if 'Proteins' in input_features:
        omics_types.append('Proteins')
    # Remove duplicates while preserving order
    omics_types = list(dict.fromkeys(omics_types))
    # Filter out 'Metadata' for omics count
    omics_only = [t for t in omics_types if t != 'Metadata']
    
    # Get initial and final module counts
    n_initial_modules = all_stats.get('clustering', {}).get('n_initial_modules', 'N/A')
    n_final_modules = network_stats.get('modules_after_coherence', 'N/A')
    
    content_parts.append(f"""
    <section class="section" id="summary-section" style="background: linear-gradient(135deg, #E8F5E9 0%, #fff 100%);">
        <div class="section-header">
            <span class="section-icon">📋</span>
            <h2>Executive Summary</h2>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card" style="background: linear-gradient(135deg, #C8E6C9 0%, #fff 100%);">
                <div class="stat-value">{input_stats.get('samples', 'N/A')}</div>
                <div class="stat-label">Samples Analyzed</div>
            </div>
            <div class="stat-card" style="background: linear-gradient(135deg, #B3E5FC 0%, #fff 100%);">
                <div class="stat-value">{preproc_stats.get('genes_retained', 'N/A')}</div>
                <div class="stat-label">Genes Analyzed</div>
                <div class="stat-detail">After HVG selection + re-introduction of TFs</div>
            </div>
            <div class="stat-card" style="background: linear-gradient(135deg, #FFE082 0%, #fff 100%);">
                <div class="stat-value">{network_stats.get('modules_after_coherence', 'N/A')}</div>
                <div class="stat-label">Final Modules</div>
                <div class="stat-detail">Passed coherence filter</div>
            </div>
            <div class="stat-card" style="background: linear-gradient(135deg, #CE93D8 0%, #fff 100%);">
                <div class="stat-value">{network_stats.get('de_modules_count', 'N/A')}</div>
                <div class="stat-label">DE Modules</div>
                <div class="stat-detail">Differentially expressed</div>
            </div>
        </div>
        
        <div class="info-box info-box-success">
            <span>✅</span>
            <div>
                <strong>Analysis completed successfully!</strong> This report summarizes the Lemonite multi-omics integration analysis.
                <br><br>
                <strong>Key Results:</strong>
                <ul style="margin-top: 0.5rem; margin-left: 1rem;">
                    <li>{len(omics_only)} omics data types integrated ({', '.join(omics_only)})</li>
                    <li>{n_initial_modules} consensus modules identified, {n_final_modules} passed quality filtering</li>
                    <li>{sum(d.get('selected', 0) for d in network_stats.get('regulators_per_type', {}).values())} regulators assigned across {len(network_stats.get('regulators_per_type', {}))} regulator type(s)</li>
                    <li>{network_stats.get('genes_in_network', 'N/A')} unique genes in the regulatory network</li>
                    <li>{'✓ Enrichment analysis completed' if enrich_stats.get('enrichment_performed') else '✗ Enrichment analysis not available'}</li>
                </ul>
            </div>
        </div>
    </section>
    """)
    
    # === INPUT DATA SECTION ===
    preproc_stats = all_stats.get('preprocessing', {})
    input_features = preproc_stats.get('input_features', {})
    
    content_parts.append(f"""
    <section class="section" id="input-section">
        <div class="section-header">
            <span class="section-icon">📊</span>
            <h2>Input Data Overview</h2>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{input_stats.get('samples', 'N/A')}</div>
                <div class="stat-label">Samples</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{len(input_stats.get('data_types', []))}</div>
                <div class="stat-label">Omics Data Types</div>
                <div class="stat-detail">{', '.join(input_stats.get('data_types', []))}</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{len(input_stats.get('files_found', []))}</div>
                <div class="stat-label">Input Files</div>
            </div>
        </div>
        
        {'<div class="subsection"><h3>📊 Input Features per Omics Type</h3><div class="stats-grid">' + ''.join([f'<div class="stat-card"><div class="stat-value">{count}</div><div class="stat-label">{omics_type}</div></div>' for omics_type, count in input_features.items()]) + '</div></div>' if input_features else ''}
        
        <div class="subsection">
            <h3>📁 Input Files</h3>
            <div class="table-container">
                <table>
                    <tr><th>Data Type</th><th>File</th><th>Features</th></tr>
                    {generate_input_files_table(input_stats)}
                </table>
            </div>
        </div>
        
        <div class="subsection">
            <h3>📋 Metadata Columns</h3>
            <p>{', '.join(input_stats.get('metadata_columns', ['N/A'])) if input_stats.get('metadata_columns') else 'No metadata columns found'}</p>
        </div>
    </section>
    """)
    
    # === PREPROCESSING SECTION ===
    content_parts.append(f"""
    <section class="section" id="preprocessing-section">
        <div class="section-header">
            <span class="section-icon">⚙️</span>
            <h2>Preprocessing Summary</h2>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{preproc_stats.get('genes_retained', 'N/A')}</div>
                <div class="stat-label">Genes After HVG Selection</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{preproc_stats.get('total_features_integrated', 'N/A')}</div>
                <div class="stat-label">Total Features Integrated</div>
                <div class="stat-detail">All omics combined</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{preproc_stats.get('n_samples', 'N/A')}</div>
                <div class="stat-label">Samples</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{'✓' if preproc_stats.get('tfa_performed', False) else '✗'}</div>
                <div class="stat-label">TFA Analysis</div>
                <div class="stat-detail">{f"Completed" if preproc_stats.get('tfa_performed') else 'Not performed'}</div>
            </div>
        </div>
        
        <div class="info-box info-box-info">
            <span>ℹ️</span>
            <div>
                <strong>Preprocessing Steps Applied:</strong>
                <ul style="margin-top: 0.5rem; margin-left: 1rem;">
                    <li>DESeq2 normalization (transcriptomics)</li>
                    <li>Log transformation (auto-detected)</li>
                    <li>Omics-specific scaling (Pareto for metabolomics/lipidomics, Z-score for transcriptomics)</li>
                    <li>Highly variable gene selection</li>
                    <li>Vertical integration of all omics</li>
                </ul>
            </div>
        </div>
        
        {generate_group_counts_section(preproc_stats)}
    </section>
    """)
    
    # === PARAMETERS SECTION ===
    params = all_stats.get('parameters', {})
    content_parts.append(f"""
    <section class="section" id="parameters-section">
        <div class="section-header">
            <span class="section-icon">🎛️</span>
            <h2>Pipeline Parameters</h2>
        </div>
        
        <button class="collapsible">
            <span><strong>🔧 Preprocessing Parameters</strong></span>
            <span class="toggle-icon">▶</span>
        </button>
        <div class="collapsible-content">
            <div class="param-list" style="padding: 1rem;">
                {generate_param_items(params, ['top_n_genes', 'perform_tfa', 'use_omics_specific_scaling', 'expression_col', 'sample_id_col', 'metadata_columns', 'design_formula', 'deseq_contrast1', 'organism'])}
            </div>
        </div>
        
        <button class="collapsible">
            <span><strong>🧩 Clustering Parameters</strong></span>
            <span class="toggle-icon">▶</span>
        </button>
        <div class="collapsible-content">
            <div class="param-list" style="padding: 1rem;">
                {generate_param_items(params, ['n_clusters', 'coherence_threshold', 'use_deseq_priors', 'min_cluster_size', 'tight_clusters_only', 'max_n_iterations'])}
            </div>
        </div>
        
        <button class="collapsible">
            <span><strong>🔗 Regulator & Network Parameters</strong></span>
            <span class="toggle-icon">▶</span>
        </button>
        <div class="collapsible-content">
            <div class="param-list" style="padding: 1rem;">
                {generate_param_items(params, ['regulator_types', 'regulator_selection_method', 'top_n_percent_regulators', 'regulator_fold_cutoff', 'min_regulator_size', 'max_regulator_size', 'min_module_size', 'min_targets'])}
            </div>
        </div>
        
        <button class="collapsible">
            <span><strong>🔬 Analysis Parameters</strong></span>
            <span class="toggle-icon">▶</span>
        </button>
        <div class="collapsible-content">
            <div class="param-list" style="padding: 1rem;">
                {generate_param_items(params, ['enrichment_method', 'enrichr_libraries', 'pkn_network', 'use_megago', 'prioritize_by_expression', 'overview_n_clusters'])}
            </div>
        </div>
    </section>
    """)
    
    # === CLUSTERING SECTION ===
    cluster_stats = all_stats.get('clustering', {})
    content_parts.append(f"""
    <section class="section" id="clustering-section">
        <div class="section-header">
            <span class="section-icon">🧩</span>
            <h2>LemonTree Clustering Results</h2>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{cluster_stats.get('n_clustering_runs', 'N/A')}</div>
                <div class="stat-label">Clustering Runs</div>
                <div class="stat-detail">Parallel Gibbs samplers</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{cluster_stats.get('n_tight_clusters', 'N/A')}</div>
                <div class="stat-label">Consensus Modules</div>
                <div class="stat-detail">Tight clusters</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{cluster_stats.get('mean_cluster_size', 'N/A') if not isinstance(cluster_stats.get('mean_cluster_size'), (int, float)) else f"{cluster_stats.get('mean_cluster_size'):.1f}"}</div>
                <div class="stat-label">Mean Module Size</div>
                <div class="stat-detail">genes per module</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{cluster_stats.get('median_cluster_size', 'N/A') if not isinstance(cluster_stats.get('median_cluster_size'), (int, float)) else f"{cluster_stats.get('median_cluster_size'):.0f}"}</div>
                <div class="stat-label">Median Module Size</div>
                <div class="stat-detail">Range: {cluster_stats.get('min_cluster_size', '?')}-{cluster_stats.get('max_cluster_size', '?')}</div>
            </div>
        </div>
    </section>
    """)
    
    # === NETWORK SECTION ===
    network_stats = all_stats.get('network', {})
    rankings = all_stats.get('rankings', {})
    
    # Format coherence values
    mean_coh = network_stats.get('mean_coherence', 'N/A')
    mean_coh_str = f"{mean_coh:.3f}" if isinstance(mean_coh, (int, float)) else 'N/A'
    median_coh = network_stats.get('median_coherence', 'N/A')
    median_coh_str = f"{median_coh:.3f}" if isinstance(median_coh, (int, float)) else 'N/A'
    
    content_parts.append(f"""
    <section class="section" id="network-section">
        <div class="section-header">
            <span class="section-icon">🔗</span>
            <h2>Network Generation Results</h2>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{cluster_stats.get('n_initial_modules', 'N/A')}</div>
                <div class="stat-label">Initial Modules</div>
                <div class="stat-detail">Before coherence filtering</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{network_stats.get('modules_after_coherence', 'N/A')}</div>
                <div class="stat-label">Final Modules</div>
                <div class="stat-detail">After coherence filtering</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{network_stats.get('de_modules_count', 'N/A')}</div>
                <div class="stat-label">DE Modules</div>
                <div class="stat-detail">p-value &lt; 0.05</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{rankings.get('ppi_enriched_modules_count', network_stats.get('ppi_enriched_modules_count', 'N/A'))}</div>
                <div class="stat-label">PPI Enriched Modules</div>
                <div class="stat-detail">FDR &lt; 0.05</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{rankings.get('metgene_enriched_modules_count', 'N/A')}</div>
                <div class="stat-label">Met-Gene Enriched Modules</div>
                <div class="stat-detail">FDR &lt; 0.05</div>
            </div>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{network_stats.get('genes_in_network', 'N/A')}</div>
                <div class="stat-label">Genes in Network</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{network_stats.get('regulators_per_type', {}).get('TFs', {}).get('selected', 'N/A')}</div>
                <div class="stat-label">TFs in Network</div>
                <div class="stat-detail">Out of {preproc_stats.get('input_features', {}).get('TFs', '?')} input</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{network_stats.get('regulators_per_type', {}).get('Metabolites', {}).get('selected', 'N/A')}</div>
                <div class="stat-label">Metabolites in Network</div>
                <div class="stat-detail">Out of {preproc_stats.get('input_features', {}).get('Metabolites', '?')} input</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{network_stats.get('network_edges', 'N/A')}</div>
                <div class="stat-label">Network Edges</div>
            </div>
        </div>
        
        {generate_coherence_progress(network_stats)}
        
        <div class="subsection">
            <h3>🎯 Regulators per Type</h3>
            {generate_regulators_table(network_stats.get('regulators_per_type', {}))}
        </div>
        
        <div class="info-box info-box-success">
            <span>✓</span>
            <div>
                <strong>Module Coherence Statistics:</strong>
                Mean coherence: {mean_coh_str} | 
                Median coherence: {median_coh_str}
            </div>
        </div>
        
        <div class="subsection">
            <h3>🌐 Interactive Module Network</h3>
            <p>Visualize module-regulator interactions, functional clusters, and enrichment annotations.</p>
            <div style="display:flex; gap:1rem; flex-wrap:wrap; margin-top:0.75rem;">
                <a href="LemonTree/Module_Overview/Module_Overview_Comprehensive.html" target="_blank"
                   style="display:inline-block; padding:0.6rem 1.2rem; background:#2E7D32; color:white;
                          border-radius:6px; text-decoration:none; font-weight:600;">
                    📊 Open Comprehensive Report
                </a>
                <a href="LemonTree/Module_Overview/interactive_module_network.html" target="_blank"
                   style="display:inline-block; padding:0.6rem 1.2rem; background:#1565C0; color:white;
                          border-radius:6px; text-decoration:none; font-weight:600;">
                    🌐 Open Network Visualization
                </a>
                <a href="LemonTree/Module_Overview/interactive_module_network_movable.html" target="_blank"
                   style="display:inline-block; padding:0.6rem 1.2rem; background:#6A1B9A; color:white;
                          border-radius:6px; text-decoration:none; font-weight:600;">
                    🖱️ Open Draggable Network
                </a>
            </div>
        </div>
    </section>
    """)
    
    # === EVALUATION SECTION ===
    pkn_stats = all_stats.get('pkn_evaluation', {})
    content_parts.append(f"""
    <section class="section" id="evaluation-section">
        <div class="section-header">
            <span class="section-icon">✅</span>
            <h2>PKN Evaluation</h2>
        </div>
        
        {'<div class="info-box info-box-success"><span>✓</span><div>PKN evaluation was performed successfully.</div></div>' if pkn_stats.get('pkn_evaluation_performed') else '<div class="info-box info-box-warning"><span>⚠️</span><div>PKN evaluation was not performed or no results found.</div></div>'}
        
        {generate_pkn_metrics(pkn_stats)}
    </section>
    """)
    
    # === ENRICHMENT SECTION ===
    enrich_stats = all_stats.get('enrichment', {})
    content_parts.append(f"""
    <section class="section" id="enrichment-section">
        <div class="section-header">
            <span class="section-icon">🔬</span>
            <h2>Pathway Enrichment Analysis</h2>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{enrich_stats.get('method', 'N/A')}</div>
                <div class="stat-label">Method Used</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{enrich_stats.get('modules_with_enrichment', 'N/A')}</div>
                <div class="stat-label">Modules with Enrichment</div>
            </div>
        </div>
        
        {'<div class="info-box info-box-success"><span>✓</span><div>Enrichment analysis completed successfully.</div></div>' if enrich_stats.get('enrichment_performed') else '<div class="info-box info-box-warning"><span>⚠️</span><div>Enrichment analysis was not performed or no results found.</div></div>'}
        
        {generate_top_pathways(enrich_stats.get('top_pathways', {}))}
    </section>
    """)
    
    # === KEY OUTPUT FILES SECTION ===
    content_parts.append(f"""
    <section class="section" id="files-section">
        <div class="section-header">
            <span class="section-icon">📂</span>
            <h2>Key Output Files</h2>
        </div>
        
        <div class="info-box info-box-info">
            <span>📁</span>
            <div>
                <strong>Your results are organized in the following structure:</strong>
            </div>
        </div>
        
        <div class="table-container">
            <table>
                <tr><th>Category</th><th>File/Directory</th><th>Description</th></tr>
                <tr><td><span class="badge badge-info">Preprocessing</span></td><td>LemonTree/Preprocessing/</td><td>Normalized expression data, DESeq groups, PCA plots</td></tr>
                <tr><td><span class="badge badge-info">Clustering</span></td><td>LemonTree/Lemon_out/</td><td>Consensus modules, regulator scores</td></tr>
                <tr><td><span class="badge badge-success">Networks</span></td><td>LemonTree/Networks/LemonNetwork_*.txt</td><td>Main regulatory network (Cytoscape-compatible)</td></tr>
                <tr><td><span class="badge badge-success">Networks</span></td><td>LemonTree/Networks/Module_coherence_scores.txt</td><td>Quality scores for each module</td></tr>
                <tr><td><span class="badge badge-warning">Viewer</span></td><td>LemonTree/ModuleViewer_files/</td><td>Module-gene and module-regulator mappings</td></tr>
                <tr><td><span class="badge badge-warning">Heatmaps</span></td><td>LemonTree/module_heatmaps/</td><td>Per-module expression heatmaps</td></tr>
                <tr><td><span class="badge badge-info">Enrichment</span></td><td>LemonTree/Enrichment/</td><td>Pathway enrichment results (GO, KEGG, Reactome)</td></tr>
                <tr><td><span class="badge badge-info">PKN</span></td><td>LemonTree/PKN_Evaluation/</td><td>Validation against prior knowledge network</td></tr>
                <tr><td><span class="badge badge-success">Overview</span></td><td>LemonTree/Module_Overview/</td><td>Interactive module overview and summary tables</td></tr>
                <tr><td><span class="badge badge-success">Overview</span></td><td><a href="LemonTree/Module_Overview/Module_Overview_Comprehensive.html" target="_blank">LemonTree/Module_Overview/Module_Overview_Comprehensive.html</a></td><td>📊 Combined report: interactive network + regulator rankings</td></tr>
                <tr><td><span class="badge badge-success">Overview</span></td><td><a href="LemonTree/Module_Overview/interactive_module_network.html" target="_blank">LemonTree/Module_Overview/interactive_module_network.html</a></td><td>🌐 Interactive module-regulator network (standalone)</td></tr>
                <tr><td><span class="badge badge-success">Overview</span></td><td><a href="LemonTree/Module_Overview/interactive_module_network_movable.html" target="_blank">LemonTree/Module_Overview/interactive_module_network_movable.html</a></td><td>🖱️ Draggable network visualization</td></tr>
                <tr><td><span class="badge badge-success">Report</span></td><td>Lemonite_Summary_Report.html</td><td>This summary report</td></tr>
            </table>
        </div>
    </section>
    """)
    
    # === MODULE RANKINGS SECTION ===
    rankings = all_stats.get('rankings', {})
    content_parts.append(f"""
    <section class="section" id="rankings-section">
        <div class="section-header">
            <span class="section-icon">🏆</span>
            <h2>Module Rankings</h2>
        </div>
        
        <div class="two-col">
            <div class="subsection">
                <h3>🎯 Top Modules by Coherence</h3>
                <div class="ranking-table">
                    {generate_ranking_table(rankings.get('by_coherence', []), 'Coherence_Score')}
                </div>
            </div>
            <div class="subsection">
                <h3>📊 Top Modules by Size</h3>
                <div class="ranking-table">
                    {generate_ranking_table(rankings.get('by_size', []), 'n_genes')}
                </div>
            </div>
        </div>
        
        {f'''<div class="two-col">
            <div class="subsection">
                <h3>📈 Top Differentially Expressed Modules</h3>
                <div class="ranking-table">
                    {generate_ranking_table(rankings.get('by_de', []), 'de_pvalue')}
                </div>
            </div>
            <div class="subsection">
                <h3>🔗 Top Modules by PPI Enrichment</h3>
                <div class="ranking-table">
                    {generate_ranking_table(rankings.get('by_ppi_enrichment', []), 'ppi_pvalue')}
                </div>
            </div>
        </div>''' if rankings.get('by_de') or rankings.get('by_ppi_enrichment') else ''}
        
        {f'''<div class="two-col">
            <div class="subsection">
                <h3>🧬 Top Modules by Metabolite-Gene Interaction Enrichment</h3>
                <div class="ranking-table">
                    {generate_ranking_table(rankings.get('by_metgene_enrichment', []), 'metgene_adjusted_pvalue')}
                </div>
            </div>
            <div class="subsection">
                <p></p>
            </div>
        </div>''' if rankings.get('by_metgene_enrichment') else ''}
    </section>
    """)
    
    # === REGULATOR RANKINGS SECTION ===
    regulator_rankings = all_stats.get('regulator_rankings', {})
    content_parts.append(f"""
    <section class="section" id="regulator-rankings-section">
        <div class="section-header">
            <span class="section-icon">🧬</span>
            <h2>Regulator Rankings</h2>
        </div>
        <p style="margin-bottom:1rem;">Regulator-module pairs ranked by LemonTree association score. Higher scores indicate stronger predicted regulatory relationships. Only coherence-filtered modules are included.</p>
        {generate_regulator_ranking_section(regulator_rankings)}
    </section>
    """)
    
    # Combine all sections
    content = '\n'.join(content_parts)
    
    # Fill template
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    run_id = all_stats.get('run_id', 'unknown')
    
    # Build final HTML using string replacement instead of .format() to avoid CSS variable conflicts
    html = html_template.replace('{run_id}', str(run_id))
    html = html.replace('{timestamp}', str(timestamp))
    html = html.replace('{content}', str(content))
    
    # Clean up memory before writing
    del content_parts
    del content
    gc.collect()
    
    # Write output to file
    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8', buffering=65536) as f:
            f.write(html)
        print(f"Summary report generated: {output_path}")
    except Exception as e:
        print(f"Error writing report: {e}")
    finally:
        del html
        gc.collect()
    return output_path


def generate_input_files_table(input_stats):
    """Generate table rows for input files"""
    rows = []
    features = input_stats.get('features_per_type', {})
    
    for data_type, file_name in input_stats.get('files_found', []):
        # Try multiple key formats: (features), (genes), (regulators)
        feature_info = (features.get(f'{data_type} (features)') or 
                       features.get(f'{data_type} (genes)') or 
                       features.get(f'{data_type} (regulators)') or {})
        count = feature_info.get('input', 'N/A') if isinstance(feature_info, dict) else 'N/A'
        rows.append(f'<tr><td><span class="badge badge-info">{data_type}</span></td><td>{file_name}</td><td>{count}</td></tr>')
    
    return '\n'.join(rows) if rows else '<tr><td colspan="3">No input files found</td></tr>'


def generate_group_counts_section(preproc_stats):
    """Generate sample group counts section"""
    group_counts = preproc_stats.get('group_counts', {})
    if not group_counts:
        return ''
    
    html = '<div class="subsection"><h3>📊 Sample Distribution by Groups</h3>'
    
    for col, counts in group_counts.items():
        html += f'<p><strong>{col}:</strong> '
        items = [f'{k}: {v}' for k, v in counts.items()]
        html += ', '.join(items)
        html += '</p>'
    
    html += '</div>'
    return html


def generate_param_items(params, keys):
    """Generate parameter items HTML"""
    items = []
    for key in keys:
        value = params.get(key, 'N/A')
        items.append(f'<div class="param-item"><span class="param-name">{key}</span><span class="param-value">{value}</span></div>')
    return '\n'.join(items)


def generate_coherence_progress(network_stats):
    """Generate coherence filtering progress bar"""
    total = network_stats.get('modules_total', 0)
    # If modules_total equals modules_after_coherence, it wasn't set correctly - use n_initial_modules
    filtered = network_stats.get('modules_after_coherence', 0)
    if total == filtered and total > 0:
        # modules_total wasn't differentiated; this is already the filtered count
        pass
    
    if total == 0:
        return ''
    
    percentage = (filtered / total) * 100
    
    return f'''
    <div class="progress-container">
        <div class="progress-label">
            <span>Module Retention After Coherence Filtering</span>
            <span>{filtered} / {total} ({percentage:.1f}%)</span>
        </div>
        <div class="progress-bar">
            <div class="progress-fill" style="width: {percentage}%"></div>
        </div>
    </div>
    '''


def generate_regulators_table(regulators_per_type):
    """Generate regulators per type table"""
    if not regulators_per_type:
        return '<p>No regulator data available</p>'
    
    rows = []
    for reg_type, data in regulators_per_type.items():
        rows.append(f'''
        <tr>
            <td><span class="badge badge-success">{reg_type}</span></td>
            <td>{data.get('selected', 'N/A')}</td>
            <td>{data.get('modules_with_regulators', 'N/A')}</td>
        </tr>
        ''')
    
    return f'''
    <div class="table-container">
        <table>
            <tr><th>Regulator Type</th><th>Selected Regulators</th><th>Modules with Regulators</th></tr>
            {''.join(rows)}
        </table>
    </div>
    '''


def generate_pkn_metrics(pkn_stats):
    """Generate PKN evaluation metrics"""
    if not pkn_stats.get('pkn_evaluation_performed'):
        return ''
    
    metrics = pkn_stats.get('metrics', {})
    if not metrics:
        return '<p>No detailed metrics available</p>'
    
    return f'''
    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{metrics.get('precision', 'N/A'):.3f if isinstance(metrics.get('precision'), (int, float)) else 'N/A'}</div>
            <div class="stat-label">Precision</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{metrics.get('recall', 'N/A'):.3f if isinstance(metrics.get('recall'), (int, float)) else 'N/A'}</div>
            <div class="stat-label">Recall</div>
        </div>
    </div>
    '''


def generate_top_pathways(top_pathways):
    """Generate top pathways list"""
    if not top_pathways:
        return ''
    
    items = []
    for i, (pathway, count) in enumerate(list(top_pathways.items())[:10], 1):
        items.append(f'<tr><td>{i}</td><td>{pathway}</td><td>{count}</td></tr>')
    
    return f'''
    <div class="subsection">
        <h3>🔝 Most Frequently Enriched Pathways</h3>
        <div class="table-container">
            <table>
                <tr><th>#</th><th>Pathway</th><th>Modules Enriched</th></tr>
                {''.join(items)}
            </table>
        </div>
    </div>
    '''


def generate_ranking_table(rankings, score_col):
    """Generate module ranking table"""
    if not rankings:
        return '<p>No ranking data available</p>'
    
    rows = []
    for i, item in enumerate(rankings, 1):  # Show all modules, not just top 15
        # Determine rank badge class
        badge_class = ''
        if i == 1:
            badge_class = 'gold'
        elif i == 2:
            badge_class = 'silver'
        elif i == 3:
            badge_class = 'bronze'
        
        module = item.get('Module', item.get('module', 'N/A'))
        # Try multiple column name variations
        score = item.get(score_col, item.get(score_col.lower(), item.get(score_col.replace('_', ''), 'N/A')))
        
        # For coherence, also try 'Coherence' column
        if score == 'N/A' and 'coherence' in score_col.lower():
            score = item.get('Coherence', item.get('coherence', 'N/A'))
        
        # Detect if this is a p-value column for scientific notation formatting
        is_pvalue = any(kw in score_col.lower() for kw in ['pvalue', 'pval', 'fdr', 'adjusted'])
        if isinstance(score, float):
            score = f'{score:.4e}' if is_pvalue else f'{score:.4f}'
        elif isinstance(score, str) and score != 'N/A':
            try:
                score = f'{float(score):.4e}' if is_pvalue else f'{float(score):.4f}'
            except:
                pass
        
        rows.append(f'''
        <tr>
            <td><span class="rank-badge {badge_class}">{i}</span></td>
            <td>Module {module}</td>
            <td>{score}</td>
        </tr>
        ''')
    
    return f'''
    <table>
        <tr><th>Rank</th><th>Module</th><th>{score_col.replace('_', ' ').title()}</th></tr>
        {''.join(rows)}
    </table>
    '''


def main():
    args = parse_arguments()
    
    print("="*60)
    print("LEMONITE SUMMARY REPORT GENERATOR")
    print("="*60)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Run ID: {args.run_id}")
    print("="*60)
    
    # Parse regulator types
    regulator_configs = parse_regulator_types(args.regulator_types)
    
    # Collect all statistics with garbage collection
    print("Collecting input statistics...")
    input_stats = collect_input_statistics(args.input_dir, regulator_configs, args.output_dir)
    print(f"[DEBUG MAIN] Input stats collected: samples={input_stats.get('samples')}, data_types={input_stats.get('data_types')}", file=sys.stderr)
    all_stats = {
        'run_id': args.run_id,
        'organism': args.organism,
        'input': input_stats,
    }
    gc.collect()
    
    print("Collecting preprocessing statistics...")
    all_stats['preprocessing'] = collect_preprocessing_statistics(args.output_dir, args.run_id)
    gc.collect()
    
    print("Collecting clustering statistics...")
    all_stats['clustering'] = collect_clustering_statistics(args.output_dir, args.run_id)
    gc.collect()
    
    print("Collecting network statistics...")
    all_stats['network'] = collect_network_statistics(args.output_dir, args.run_id, regulator_configs)
    gc.collect()
    
    print("Collecting PKN evaluation statistics...")
    all_stats['pkn_evaluation'] = collect_pkn_evaluation_statistics(args.output_dir, args.run_id)
    gc.collect()
    
    print("Collecting enrichment statistics...")
    all_stats['enrichment'] = collect_enrichment_statistics(args.output_dir, args.run_id)
    gc.collect()
    
    print("Collecting module rankings...")
    all_stats['rankings'] = collect_module_rankings(args.output_dir, args.run_id)
    gc.collect()
    
    print("Collecting regulator rankings...")
    all_stats['regulator_rankings'] = collect_regulator_rankings(args.output_dir, args.run_id, regulator_configs)
    gc.collect()
    
    print("Parsing parameters...")
    all_stats['parameters'] = {}
    
    # Parse parameters file if available
    params_file = args.parameters_file
    if not params_file:
        params_file = os.path.join(args.output_dir, args.run_id, 'pipeline_parameters_log.txt')
        if not os.path.exists(params_file):
            params_file = os.path.join(args.output_dir, 'pipeline_parameters_log.txt')
    
    all_stats['parameters'] = parse_parameters_file(params_file)
    gc.collect()
    
    # Generate HTML report
    print("Generating HTML report...")
    output_path = os.path.join(args.output_dir, args.run_id, 'Lemonite_Summary_Report.html')
    
    # Create directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    generate_html_report(all_stats, output_path)
    gc.collect()
    
    print("\n" + "="*60)
    print("REPORT GENERATION COMPLETE")
    print("="*60)
    print(f"Report saved to: {output_path}")
    

if __name__ == '__main__':
    main()
