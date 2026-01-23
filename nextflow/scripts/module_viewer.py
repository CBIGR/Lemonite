#!/usr/bin/env python3
"""
ModuleViewer: Create comprehensive heatmap visualizations for LemonTree modules

This script generates detailed heatmap visualizations for each module, including:
- Gene expression data with eigengene-based sample ordering
- Regulator information (transcription factors, metabolites) with optional scores
- Sample annotations and legends
- Metabolite-gene interaction overlays
- Publication-ready figures in multiple formats

Usage:
    python module_viewer.py --input_dir <path> --output_dir <path> [options]
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import gridspec
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

def load_sample_mapping(folder, filename):
    """Load sample mapping and legend from .mvf file
    
    Supports both old single-metadata format and new multi-metadata format.
    
    Returns:
        dict: Dictionary of metadata types, each containing 'legend', 'samples', and 'label'
              e.g., {'diagnosis': {'legend': {...}, 'samples': {...}, 'label': 'Clinical Status'}}
              
              For backward compatibility with old format, returns single 'diagnosis' entry.
    """
    file_path = os.path.join(folder, filename)
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found")
        print("This file is required for the module viewer to show sample annotations")
        sys.exit(1)
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Check if this is new multi-metadata format (contains '---' separator)
    if '\n---\n' in content:
        # New format: multiple metadata sections separated by '---'
        sections = content.split('\n---\n')
        metadata_dict = {}
        
        for section in sections:
            lines = section.strip().split('\n')
            if len(lines) < 6:  # Need at least TYPE, GLOBAL, VALUES, OBJECT, LEGEND, mapping
                continue
                
            metadata_type = None
            legend_dict = {}
            sample_color_dict = {}
            label = None
            
            for line in lines:
                if line.startswith('::TYPE='):
                    metadata_type = line.split('=', 1)[1].strip()
                elif line.startswith('::LEGEND='):
                    legend_part = line.split('=', 1)[1]
                    parts = legend_part.split('\t')
                    legend_items = parts[0].split('|')
                    label = parts[1] if len(parts) > 1 else metadata_type.replace('_', ' ').title()
                    
                    for item in legend_items:
                        if ':' in item:
                            lab, color = item.split(':', 1)
                            legend_dict[lab.strip()] = color.strip().lower()
                            
                elif line.startswith('|'):
                    for item in line.strip('| \n').split('|'):
                        if ':' in item:
                            sample, color = item.split(':', 1)
                            sample_color_dict[sample.strip()] = color.strip().lower()
            
            if metadata_type and legend_dict and sample_color_dict:
                metadata_dict[metadata_type] = {
                    'legend': legend_dict,
                    'samples': sample_color_dict,
                    'label': label
                }
        
        if metadata_dict:
            return metadata_dict
        else:
            print(f"Error: Could not parse multi-metadata format in {filename}")
            sys.exit(1)
    
    else:
        # Old format: single metadata type (backward compatibility)
        lines = content.split('\n')
        legend_line = None
        mapping_line = None
        
        for line in lines:
            if line.startswith('::LEGEND='):
                legend_line = line
            elif line.startswith('|'):
                mapping_line = line
        
        if not legend_line or not mapping_line:
            print(f"Error: Invalid format in {filename}")
            print("The file exists but does not contain the expected LEGEND and mapping lines")
            sys.exit(1)
        
        # Parse legend
        legend_dict = {}
        legend_items = legend_line.split('=', 1)[1].split('\t')[0].split('|')
        for item in legend_items:
            if ':' in item:
                lab, color = item.split(':', 1)
                legend_dict[lab.strip()] = color.strip().lower()
        
        # Parse sample mapping
        sample_color_dict = {}
        for item in mapping_line.strip('| \n').split('|'):
            if ':' in item:
                sample, color = item.split(':', 1)
                sample_color_dict[sample.strip()] = color.strip().lower()
        
        # Get label from legend line if available
        legend_parts = legend_line.split('\t')
        label = legend_parts[1].strip() if len(legend_parts) > 1 else 'Clinical Status'
        
        # Return in new format for consistency
        return {
            'diagnosis': {
                'legend': legend_dict,
                'samples': sample_color_dict,
                'label': label
            }
        }

def create_subset(main_df, module, cluster_data, gene_column='genes'):
    """Create ordered subset for a module"""
    mask = cluster_data['module'] == module
    if not mask.any():
        mask = cluster_data['module'].astype(str) == str(module)
    
    if not mask.any():
        return pd.DataFrame()
    
    genes = cluster_data.loc[mask, gene_column].values[0]
    if isinstance(genes, str):
        genes = [g.strip() for g in genes.split('|') if g.strip()]

    # Ensure main_df has a 'symbol' column (handle common alternatives)
    if 'symbol' not in main_df.columns:
        for alt in ['Gene', 'gene', 'Symbol', 'GENE']:
            if alt in main_df.columns:
                main_df = main_df.rename(columns={alt: 'symbol'})
                break

    # Create a fast lookup of available symbols
    available_symbols = set(main_df['symbol'].astype(str).str.strip().values)

    # Preserve requested order but only keep genes present in main_df
    genes_present = [g for g in genes if g in available_symbols]

    # Try case-insensitive matching for missing genes and map back to actual symbols
    if len(genes_present) < len(genes):
        missing = [g for g in genes if g not in genes_present]
        upper_map = {s.upper(): s for s in available_symbols}
        for gm in missing:
            gm_up = gm.upper()
            if gm_up in upper_map:
                genes_present.append(upper_map[gm_up])

    if not genes_present:
        # Nothing to return
        return pd.DataFrame()

    # Subset and preserve the requested order
    subset = main_df[main_df['symbol'].astype(str).str.strip().isin(genes_present)].copy()
    order_map = {g: i for i, g in enumerate(genes_present)}
    subset['__order'] = subset['symbol'].astype(str).map(lambda s: order_map.get(s, 999999))
    subset = subset.sort_values('__order').drop(columns=['__order'])
    subset = subset.reset_index(drop=True)
    
    # Diagnostic print
    requested = len(genes)
    found = len(genes_present)
    if found < requested:
        missing_final = [g for g in genes if g not in genes_present]
        print(f"Warning: Module {module} requested {requested} genes, found {found} in expression data. Missing examples: {missing_final[:10]}")

    return subset

def create_subset_with_scores(main_df, module, cluster_data, score_data=None, gene_column='genes', show_scores=False):
    """Create ordered subset for a module with optional scores"""
    subset = create_subset(main_df, module, cluster_data, gene_column)
    
    if subset.empty or not show_scores or score_data is None:
        return subset
    
    # Add scores if available
    module_scores = score_data[score_data['Target'].astype(str) == str(module)]
    score_mapping = dict(zip(module_scores['Regulator'], module_scores['Score']))
    
    # Update symbol column to include scores
    new_symbols = []
    for symbol in subset['symbol']:
        if symbol in score_mapping:
            score = score_mapping[symbol]
            new_symbols.append(f"{symbol} ({score:.2f})")
        else:
            new_symbols.append(symbol)
    
    subset['symbol'] = new_symbols
    return subset

def load_metabolite_interactions(folder, filename):
    """Load metabolite-gene interactions from .mvf file"""
    file_path = os.path.join(folder, filename)
    
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found")
        print("This file is required for the module viewer to show metabolite-gene interactions")
        sys.exit(1)
    
    metadata = {}
    data_rows = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('::'):
                key_value = line[2:].split('=', 1) if '=' in line else line[2:].split(':', 1)
                if len(key_value) == 2:
                    metadata[key_value[0].strip()] = key_value[1].strip()
            elif line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    data_rows.append(parts[:3])
    
    if data_rows:
        return pd.DataFrame(data_rows, columns=['Module', 'Genes', 'Metabolite']), metadata
    else:
        print(f"Warning: {file_path} exists but contains no data")
        return pd.DataFrame(columns=['Module', 'Genes', 'Metabolite']), {}

def load_ppi_interactions(folder, filename):
    """Load PPI interactions from .mvf file
    
    Args:
        folder (str): Directory containing the PPI file
        filename (str): Name of the PPI file
    
    Returns:
        tuple: (DataFrame with columns ['Module', 'GenePair'], metadata dict)
    """
    metadata = {}
    data_rows = []
    file_path = os.path.join(folder, filename)
    
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('::'):
                    key_value = line[2:].split('=', 1) if '=' in line else line[2:].split(':', 1)
                    if len(key_value) == 2:
                        metadata[key_value[0].strip()] = key_value[1].strip()
                elif line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        data_rows.append(parts[:2])
        
        if data_rows:
            return pd.DataFrame(data_rows, columns=['Module', 'GenePair']), metadata
        else:
            print(f"Info: {file_path} exists but contains no PPI data")
            return pd.DataFrame(columns=['Module', 'GenePair']), metadata
    else:
        print(f"Info: {file_path} not found - PPI visualization will be skipped")
        return pd.DataFrame(columns=['Module', 'GenePair']), {}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Generate ModuleViewer heatmaps')
    parser.add_argument('--input_dir', required=True, help='Input directory with LemonTree results')
    parser.add_argument('--output_dir', required=True, help='Output directory for figures')
    parser.add_argument('--regulator_files', required=True,
                       help='Comma-separated list of regulator files (format: Type:Path,Type:Path)')
    # Note: viewer files are always expected in input_dir/ModuleViewer_files
    # The --viewer_files_dir option is intentionally ignored to avoid duplicate locations
    parser.add_argument('--viewer_files_dir', help='(deprecated) Directory containing ModuleViewer files (ignored)')
    parser.add_argument('--expression_file', default='LemonPreprocessed_expression.txt', help='Expression data file for main heatmap')
    parser.add_argument('--complete_file', default='LemonPreprocessed_complete.txt', help='Complete omics data file for regulator blocks')
    parser.add_argument('--show_regulator_scores', action='store_true', help='Show regulator scores in labels')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for output figures')
    parser.add_argument('--modules', help='Comma-separated list of specific modules to process')
    parser.add_argument('--annotation_types', default='diagnosis',
                       help='Comma-separated list of metadata columns to display as annotation bars (default: diagnosis)')
    parser.add_argument('--annotation_labels', default=None,
                       help='Optional: Comma-separated custom labels for annotation bars (must match --annotation_types count)')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Setup paths
    input_dir = args.input_dir
    output_dir = args.output_dir
    # Canonical ModuleViewer files directory (do not allow alternate dirs)
    viewer_files_dir = os.path.join(input_dir, 'ModuleViewer_files')
    if not os.path.exists(viewer_files_dir):
        print(f"Error: ModuleViewer files directory not found: {viewer_files_dir}")
        print("Ensure that the pipeline produced ModuleViewer_files with sample_mapping.mvf and metabolite_LemoniteKG_interactions.mvf")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Processing ModuleViewer data from: {viewer_files_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Regulator files: {args.regulator_files}")
    
    # Configuration
    SHOW_REGULATOR_SCORES = args.show_regulator_scores
    
    # Parse regulator files dynamically
    regulator_configs = {}
    if args.regulator_files.strip():  # Only parse if not empty
        for reg_config in args.regulator_files.split(','):
            if reg_config.strip():  # Skip empty entries
                reg_type, reg_path = reg_config.split(':', 1)
                regulator_configs[reg_type] = reg_path.strip()
    
    print(f"Found regulator types: {list(regulator_configs.keys())}")
    
    # Define regulators dynamically
    REGULATORS = []
    for reg_type, reg_path in regulator_configs.items():
        # Extract filename from path
        reg_filename = os.path.basename(reg_path)
        
        # Determine score file based on regulator type prefix
        # The score files from LemonTree are: <prefix>.topreg.txt, <prefix>.allreg.txt, <prefix>.randomreg.txt
        # Extract prefix from the list file (e.g., TFs.selected_regs_list.txt -> TFs)
        # Or Metabolites.selected_regs_list.txt -> Metabolites
        prefix = reg_filename.split('.')[0]  # Get first part before dot
        score_file = f"{prefix}.topreg.txt"  # LemonTree generates .topreg.txt for selected regulators
        
        # Determine column name based on regulator type
        if reg_type.lower() in ['metabolite', 'lipid']:
            column_name = 'metabolites'
        else:
            column_name = 'genes'
        
        REGULATORS.append({
            "name": reg_type,
            "file": reg_filename,
            "score_file": score_file,
            "column": column_name
        })
    
    print(f"Configured {len(REGULATORS)} regulator types: {[r['name'] for r in REGULATORS]}")
    
    # Helper function to find file in multiple locations
    def find_data_file(filename, file_description):
        locations = [
            os.path.join(viewer_files_dir, '..', 'Preprocessing', filename),
            os.path.join(input_dir, filename),
            os.path.join(input_dir, 'Preprocessing', filename),
            os.path.join(input_dir, 'results', 'test_fixed', 'LemonTree', 'Preprocessing', filename)
        ]
        
        for location in locations:
            if os.path.exists(location):
                return location
        
        # Try searching for any results directory
        import glob
        potential_paths = glob.glob(os.path.join(input_dir, 'results', '*', 'LemonTree', 'Preprocessing', filename))
        if potential_paths:
            return potential_paths[0]
        
        print(f"Error: {file_description} file not found: {filename}")
        print(f"Searched in: {viewer_files_dir}/../Preprocessing/, {input_dir}/, {input_dir}/Preprocessing/, {input_dir}/results/*/LemonTree/Preprocessing/")
        sys.exit(1)
    
    # Load expression data for main heatmap (gene expression only, no other omics)
    expression_file = find_data_file(args.expression_file, 'Expression data')
    print(f"Loading expression data from: {expression_file}")
    expression_data = pd.read_csv(expression_file, sep='\t')
    
    # Load complete omics data for regulator blocks (includes metabolomics, proteomics, etc.)
    complete_file = find_data_file(args.complete_file, 'Complete omics data')
    print(f"Loading complete omics data from: {complete_file}")
    complete_data = pd.read_csv(complete_file, sep='\t')
    
    # Load cluster data
    cluster_file = os.path.join(viewer_files_dir, 'clusters_list.txt')
    if not os.path.exists(cluster_file):
        print(f"Error: Clusters file not found: {cluster_file}")
        sys.exit(1)
    
    cluster = pd.read_csv(cluster_file, sep='\t', header=None, names=['module', 'genes'])
    cluster['genes'] = cluster['genes'].str.split('|')
    
    # Filter modules if specific_modules.txt exists
    specific_modules_file = os.path.join(input_dir, 'Networks', 'specific_modules.txt')
    if os.path.exists(specific_modules_file):
        with open(specific_modules_file, 'r') as f:
            specific_modules = [line.strip() for line in f.readlines()]
        
        print(f"Found {len(specific_modules)} filtered modules from LemonTree_to_network")
        
        # Convert to same type as cluster module column
        if cluster['module'].dtype == 'object':
            specific_modules = [str(m) for m in specific_modules]
        else:
            specific_modules = [int(m) if m.isdigit() else m for m in specific_modules]
        
        original_count = len(cluster['module'].unique())
        cluster = cluster[cluster['module'].isin(specific_modules)]
        filtered_count = len(cluster['module'].unique())
        
        print(f"Modules before filtering: {original_count}")
        print(f"Modules after filtering: {filtered_count}")
    else:
        print("Warning: specific_modules.txt not found - processing all modules")
    
    # Filter for specific modules if requested
    if args.modules:
        requested_modules = [m.strip() for m in args.modules.split(',')]
        if cluster['module'].dtype != 'object':
            requested_modules = [int(m) if m.isdigit() else m for m in requested_modules]
        cluster = cluster[cluster['module'].isin(requested_modules)]
        print(f"Processing only requested modules: {requested_modules}")
    
    # Load regulators
    loaded_regulators = []
    for reg in REGULATORS:
        file_path = os.path.join(viewer_files_dir, reg['file'])
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep='\t', header=None, names=['module', reg['column']])
            df[reg['column']] = df[reg['column']].str.split('|')
            reg['data'] = df
            
            # Load score data if available
            # Try multiple locations for score files
            score_file_locations = [
                os.path.join(input_dir, 'Lemon_out', reg['score_file']),
                os.path.join(viewer_files_dir, '..', 'Lemon_out', reg['score_file']),
                os.path.join(input_dir, reg['score_file']),
                os.path.join(viewer_files_dir, reg['score_file'])
            ]
            
            score_file_path = None
            for location in score_file_locations:
                if os.path.exists(location):
                    score_file_path = location
                    break
                    
            if score_file_path:
                score_df = pd.read_csv(score_file_path, sep='\t')
                reg['score_data'] = score_df
                if SHOW_REGULATOR_SCORES:
                    print(f"Loaded scores for {reg['name']} from {score_file_path}")
            else:
                reg['score_data'] = None
                if SHOW_REGULATOR_SCORES:
                    print(f"Warning: Score file not found for {reg['name']} in any location")
            
            loaded_regulators.append(reg)
        else:
            print(f"Warning: Regulator file not found: {file_path}")
    
    # Load sample mappings (now returns dict of metadata types)
    all_metadata = load_sample_mapping(viewer_files_dir, 'sample_mapping.mvf')
    
    # Parse requested annotation types
    requested_annotations = [a.strip() for a in args.annotation_types.split(',')]
    
    # Filter to only requested annotations that are available
    selected_metadata = {}
    for ann_type in requested_annotations:
        if ann_type in all_metadata:
            selected_metadata[ann_type] = all_metadata[ann_type]
        else:
            print(f"Warning: Requested annotation type '{ann_type}' not found in sample_mapping.mvf")
            print(f"Available types: {', '.join(all_metadata.keys())}")
    
    if not selected_metadata:
        print("Error: No valid annotation types found. Using default 'diagnosis' if available.")
        if 'diagnosis' in all_metadata:
            selected_metadata['diagnosis'] = all_metadata['diagnosis']
        else:
            # Use first available
            first_key = list(all_metadata.keys())[0]
            selected_metadata[first_key] = all_metadata[first_key]
    
    # Parse custom labels if provided
    custom_labels = {}
    if args.annotation_labels:
        labels_list = [l.strip() for l in args.annotation_labels.split(',')]
        if len(labels_list) == len(selected_metadata):
            for ann_type, custom_label in zip(selected_metadata.keys(), labels_list):
                custom_labels[ann_type] = custom_label
        else:
            print(f"Warning: Number of custom labels ({len(labels_list)}) doesn't match annotations ({len(selected_metadata)}), using defaults")
    
    print(f"Displaying {len(selected_metadata)} annotation track(s): {', '.join(selected_metadata.keys())}")
    
    # Load metabolite interactions
    metabo, meta_metadata = load_metabolite_interactions(viewer_files_dir, 'metabolite_LemoniteKG_interactions.mvf')
    
    # Load PPI interactions
    ppis, ppi_metadata = load_ppi_interactions(viewer_files_dir, 'PPI_interactions.mvf')
    if not ppis.empty:
        print(f"Loaded {len(ppis)} PPI interactions")
    else:
        print("No PPI interactions loaded - visualization will skip PPI connections")
    
    # Process each module
    modules_processed = 0
    for module_number in cluster['module'].unique():
        try:
            print(f"Processing module {module_number}...")
            
            # Create eigengene for sample ordering using expression data
            subset_cluster = create_subset(expression_data, module_number, cluster)
            if subset_cluster.empty:
                print(f"Warning: No genes found for module {module_number}, skipping")
                continue
            
            subset_cluster = subset_cluster.sort_values('symbol')
            module_expr = subset_cluster.iloc[:, 2:]  # Skip 'symbol' and first data column
            
            if module_expr.shape[0] < 2:
                print(f"Warning: Module {module_number} has too few genes for PCA, using mean")
                eigengene_series = module_expr.mean(axis=0)
            else:
                pca = PCA(n_components=1)
                eigengene = pca.fit_transform(module_expr.T)
                eigengene_series = pd.Series(eigengene.flatten(), index=module_expr.columns)
            
            sorted_samples = eigengene_series.sort_values().index.tolist()
            
            # Prepare subsets for all regulators (use complete_data for metabolomics, proteomics, etc.)
            subsets = []
            titles = []
            for reg in loaded_regulators:
                try:
                    subset = create_subset_with_scores(complete_data, module_number, reg['data'], 
                                                     score_data=reg.get('score_data'), 
                                                     gene_column=reg['column'],
                                                     show_scores=SHOW_REGULATOR_SCORES)
                    if not subset.empty:
                        subsets.append(subset)
                        titles.append(reg['name'])
                except Exception as e:
                    print(f"  Error loading {reg['name']}: {str(e)}")
                    continue
            
            # Add main expression data
            subsets.append(subset_cluster)
            titles.append("Expression Data")
            
            if not subsets:
                print(f"Warning: No data found for module {module_number}, skipping")
                continue
            
            # Create figure
            total_rows = sum(len(s) for s in subsets)
            fig = plt.figure(figsize=(15, max(total_rows/2 + 8, 12)))
            fig.suptitle(f'Module {module_number}', fontsize=18, fontweight='bold', y=0.92)
            
            # Grid layout - dynamic based on number of annotation tracks
            num_annotations = len(selected_metadata)
            total_extra_rows = num_annotations + 2  # annotations + legend row + colorbar row
            annotation_heights = [1] * num_annotations  # Each annotation bar gets height of 1
            
            gs = gridspec.GridSpec(len(subsets) + total_extra_rows, 2,
                                 height_ratios=[len(s) for s in subsets] + annotation_heights + [2.0, 1.5],
                                 width_ratios=[1, 1],
                                 hspace=0.25, wspace=0.05)
            
            # Plot each heatmap
            for idx, (subset_df, title) in enumerate(zip(subsets, titles)):
                ax = fig.add_subplot(gs[idx, :])
                numeric_df = subset_df.iloc[:, 2:].reindex(columns=sorted_samples)
                labels = subset_df['symbol'].values
                
                # Scale TF data (z-score normalization per TF across samples)
                # TFA scores are not pre-scaled in LemonPreprocessed_complete.txt
                if title.upper() in ['TFS', 'TRANSCRIPTION FACTORS']:
                    # Apply row-wise z-score normalization
                    numeric_df = numeric_df.apply(lambda row: (row - row.mean()) / row.std() if row.std() > 0 else row, axis=1)
                    print(f"  Scaled {title} data using z-score normalization")
                
                sns.heatmap(numeric_df,
                           cmap=LinearSegmentedColormap.from_list('custom', ['blue', 'black', 'yellow']),
                           vmin=-2.0, vmax=2.0,
                           yticklabels=labels,
                           xticklabels=False,
                           cbar=False,
                           ax=ax,
                           linewidths=0.5,
                           linecolor='black')
                
                # Optimize font size
                ax_height_inches = ax.get_position().height * fig.get_figheight()
                row_height_points = (ax_height_inches * 72) / len(labels)
                optimal_fontsize = max(6, min(row_height_points * 0.9, 14))
                
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position('right')
                ax.set_yticklabels(labels, rotation=0, va='center', 
                                 fontsize=optimal_fontsize, fontweight='bold')
                ax.set_title(title, fontsize=14, fontweight='bold', pad=8)
                
                # Add PPI interactions for expression data (on the left side)
                if title == "Expression Data" and not ppis.empty:
                    module_ppis = ppis[ppis['Module'] == str(module_number)]
                    
                    if len(module_ppis) > 0:
                        # Get the position of the main heatmap
                        pos = ax.get_position()
                        
                        # Create a new axis on the left side for PPI connections
                        ppi_width = 0.03  # Width of the PPI connection area
                        ax_ppi = fig.add_axes([
                            pos.x0 - ppi_width - 0.01,  # Position to the left of heatmap
                            pos.y0,
                            ppi_width,
                            pos.height
                        ])
                        
                        ax_ppi.set_xlim(0, 1)
                        ax_ppi.set_ylim(0, len(labels))
                        ax_ppi.axis('off')
                        
                        # Get PPI color from metadata (default is blue)
                        ppi_color = ppi_metadata.get('COLOR', 'blue').lower()
                        
                        # Draw lines connecting genes that have PPIs
                        ppi_count = 0
                        for _, ppi_row in module_ppis.iterrows():
                            gene_pair = ppi_row['GenePair'].split('|')
                            if len(gene_pair) == 2:
                                gene1, gene2 = gene_pair[0].strip(), gene_pair[1].strip()
                                
                                # Find positions of genes in the labels
                                try:
                                    idx1 = list(labels).index(gene1)
                                    idx2 = list(labels).index(gene2)
                                    
                                    # Draw a curved line connecting the two genes
                                    # Y coordinates are centered on each row (0.5 offset)
                                    y1 = len(labels) - idx1 - 0.5
                                    y2 = len(labels) - idx2 - 0.5
                                    
                                    # Draw the connection line
                                    ax_ppi.plot([0.2, 0.8], [y1, y2], 
                                              color=ppi_color, 
                                              linewidth=1.5, 
                                              alpha=0.7,
                                              solid_capstyle='round')
                                    
                                    # Add small dots at connection points
                                    ax_ppi.plot([0.2], [y1], 'o', 
                                              color=ppi_color, 
                                              markersize=3, 
                                              alpha=0.8)
                                    ax_ppi.plot([0.8], [y2], 'o', 
                                              color=ppi_color, 
                                              markersize=3, 
                                              alpha=0.8)
                                    
                                    ppi_count += 1
                                    
                                except ValueError:
                                    # Gene not found in labels (might have been filtered)
                                    continue
                        
                        if ppi_count > 0:
                            print(f"  Added {ppi_count} PPI connections to visualization")
                
                # Add metabolite interactions for expression data
                if title == "Expression Data" and not metabo.empty:
                    module_data = metabo[metabo['Module'] == str(module_number)]
                    metabolites = module_data['Metabolite'].unique()
                    
                    if len(metabolites) > 0:
                        match_matrix = np.zeros((len(labels), len(metabolites)))
                        for i, gene in enumerate(labels):
                            for j, metabolite in enumerate(metabolites):
                                gene_lists = module_data[module_data['Metabolite'] == metabolite]['Genes'].values
                                if any(gene in gl.split('|') for gl in gene_lists):
                                    match_matrix[i, j] = 1
                        
                        # Add metabolite panel
                        meta_color = meta_metadata.get('COLOR', 'yellow').lower()
                        pos = ax.get_position()
                        num_cols_main = numeric_df.shape[1]
                        cell_width_main = pos.width / num_cols_main
                        metabolite_panel_width = cell_width_main * len(metabolites)
                        
                        ax_meta = fig.add_axes([
                            pos.x1 + 0.10,
                            pos.y0,
                            metabolite_panel_width,
                            pos.height
                        ])
                        
                        sns.heatmap(
                            match_matrix,
                            cmap=LinearSegmentedColormap.from_list('custom_meta', ['white', meta_color]),
                            cbar=False,
                            xticklabels=metabolites,
                            yticklabels=False,
                            linewidths=0.5,
                            linecolor='black',
                            ax=ax_meta,
                            square=False,
                            clip_on=False
                        )
                        ax_meta.set_xticklabels(metabolites, rotation=90, fontsize=10, fontweight='bold')
            
            # Annotation bars - one per metadata type
            for anno_idx, (ann_type, ann_data) in enumerate(selected_metadata.items()):
                # Get color mapping for this annotation type
                sample_colors = ann_data['samples']
                legend_colors = ann_data['legend']
                color_map = {color.upper(): color.lower() for color in legend_colors.values()}
                
                # Create color data for samples
                color_data = []
                for s in sorted_samples:
                    sample_color = sample_colors.get(s, 'grey')
                    final_color = color_map.get(sample_color.upper(), sample_color.lower())
                    color_data.append(final_color)
                
                # Create annotation bar
                ax_anno = fig.add_subplot(gs[len(subsets) + anno_idx, :])
                for j, color in enumerate(color_data):
                    rect = plt.Rectangle((j, 0), 1, 1, facecolor=color, edgecolor='black', linewidth=0.5)
                    ax_anno.add_patch(rect)
                ax_anno.set_xlim(0, len(color_data))
                ax_anno.set_ylim(0, 1)
                
                # Add sample labels to the bottom annotation bar
                is_last_annotation = (anno_idx == len(selected_metadata) - 1)
                if is_last_annotation:
                    ax_anno.set_xticks([j + 0.5 for j in range(len(sorted_samples))])
                    ax_anno.set_xticklabels(sorted_samples, rotation=45, ha='right', fontsize=8)
                    ax_anno.tick_params(axis='x', which='both', length=0)  # Hide tick marks
                    # Keep only bottom spine visible
                    ax_anno.spines['top'].set_visible(False)
                    ax_anno.spines['left'].set_visible(False)
                    ax_anno.spines['right'].set_visible(False)
                    ax_anno.spines['bottom'].set_visible(False)
                    ax_anno.set_yticks([])
                else:
                    ax_anno.axis('off')
                
                # Use custom label if provided, otherwise use default from MVF or formatted type name
                label = custom_labels.get(ann_type, ann_data.get('label', ann_type.replace('_', ' ').title()))
                ax_anno.text(-0.01, 0.5, label, 
                            transform=ax_anno.transAxes,
                            rotation=0, va='center', ha='right',
                            fontsize=10, fontweight='bold')
            
            # Legend row - spans both columns, positioned below annotation bars
            ax_leg = fig.add_subplot(gs[len(subsets) + num_annotations, :])
            ax_leg.axis('off')
            
            # Create separate legend for each metadata type and place them side by side
            legend_objects = []
            for ann_idx, (ann_type, ann_data) in enumerate(selected_metadata.items()):
                legend_dict = ann_data['legend']
                color_map = {color.upper(): color.lower() for color in legend_dict.values()}
                label_text = custom_labels.get(ann_type, ann_data.get('label', ann_type.replace('_', ' ').title()))
                
                # Create legend elements for this metadata type
                legend_elements = []
                for label, color in legend_dict.items():
                    final_color = color_map.get(color.upper(), color.lower())
                    legend_elements.append(
                        Patch(facecolor=final_color, label=label, edgecolor='black', linewidth=0.5)
                    )
                
                # Calculate horizontal position for this legend
                num_legends = len(selected_metadata)
                x_position = ann_idx / num_legends
                
                # Add legend for this metadata type
                leg = ax_leg.legend(handles=legend_elements, 
                                   loc='upper left',
                                   frameon=True,
                                   fontsize=9,
                                   title=label_text,
                                   title_fontsize=10,
                                   ncol=min(len(legend_elements), 4),  # Max 4 columns per legend
                                   bbox_to_anchor=(x_position, 1.0),
                                   columnspacing=1.0,
                                   handletextpad=0.5)
                leg.get_frame().set_linewidth(0.5)
                legend_objects.append(leg)
                
                # Add previous legends back (matplotlib only keeps the last one by default)
                for prev_leg in legend_objects[:-1]:
                    ax_leg.add_artist(prev_leg)
            
            # Color bar row
            ax_colorbar = fig.add_subplot(gs[len(subsets) + num_annotations + 1, :])
            ax_colorbar.axis('off')
            cax = fig.add_axes([0.35, 0.05, 0.3, 0.02])
            cbar = plt.colorbar(plt.cm.ScalarMappable(
                norm=plt.Normalize(-2.0, 2.0),
                cmap=LinearSegmentedColormap.from_list('custom', ['blue', 'black', 'yellow'])
            ), cax=cax, orientation='horizontal')
            cbar.set_ticks([-2.0, 0, 2.0])
            cbar.set_ticklabels(['<= -2.0', '0', '>= 2.0'])
            cbar.outline.set_linewidth(0.5)
            cbar.ax.set_title('Expression ratios', fontsize=11, pad=10)
            
            # Save figure
            plt.savefig(f'{output_dir}/Module_{module_number}_heatmap.png', 
                       dpi=args.dpi, bbox_inches='tight')
            plt.savefig(f'{output_dir}/Module_{module_number}_heatmap.pdf', 
                       dpi=args.dpi, bbox_inches='tight')
            plt.close()
            
            modules_processed += 1
            print(f"  Saved module {module_number}")
            
        except Exception as e:
            print(f"Error processing module {module_number}: {str(e)}")
            import traceback
            traceback.print_exc()
            plt.close('all')
    
    print(f"\nProcessing complete! Generated heatmaps for {modules_processed} modules")
    print(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()
