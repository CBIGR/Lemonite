#!/usr/bin/env python3

"""
LemonTree to Network Conversion Script
Adapted from LemonTree_to_network.ipynb for Nextflow pipeline
Supports flexible regulator types configuration
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse
from scipy.stats import mannwhitneyu, pearsonr
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx


class RegulatorConfig:
    """Class to hold regulator type configuration"""
    def __init__(self, prefix, data_file, data_type='c'):
        """
        Args:
            prefix: Prefix used for LemonTree output files (e.g., 'TFs', 'Metabolites', 'Proteins')
            data_file: Original input data filename (stored for reference, not used in network generation)
            data_type: Data type - 'c' for continuous (default), 'd' for discrete/binary
        """
        self.type = prefix  # Used as edge type in network (e.g., 'TFs-gene')
        self.prefix = prefix  # Used for finding LemonTree output files
        self.data_file = data_file  # Store for reference/logging
        self.data_type = data_type  # 'c' for continuous, 'd' for discrete/binary
        self.allreg = None
        self.randomreg = None
        self.selected_file = None
        
    def set_files(self, input_dir, lemon_out_dir='Lemon_out'):
        """Set file paths for this regulator type"""
        # Try multiple possible file patterns:
        # 1. Files with .txt extension (expected format)
        # 2. Files without .txt extension (LemonTree output format)
        # 3. Files with different extensions
        
        possible_allreg = [
            os.path.join(input_dir, lemon_out_dir, f'{self.prefix}.allreg.txt'),
            os.path.join(input_dir, lemon_out_dir, f'{self.prefix}.allreg'),
            os.path.join(input_dir, lemon_out_dir, f'{self.prefix}.allreg.tsv'),
        ]
        
        possible_randomreg = [
            os.path.join(input_dir, lemon_out_dir, f'{self.prefix}.randomreg.txt'),
            os.path.join(input_dir, lemon_out_dir, f'{self.prefix}.randomreg'),
            os.path.join(input_dir, lemon_out_dir, f'{self.prefix}.randomreg.tsv'),
        ]
        
        # Set paths to the first existing file, or default to .txt version if none exist
        self.allreg = next((f for f in possible_allreg if os.path.exists(f)), possible_allreg[0])
        self.randomreg = next((f for f in possible_randomreg if os.path.exists(f)), possible_randomreg[0])
        
    def files_exist(self):
        """Check if required files exist"""
        return (self.allreg is not None and os.path.exists(self.allreg) and 
                self.randomreg is not None and os.path.exists(self.randomreg))
    
    def __repr__(self):
        return f"RegulatorConfig(prefix='{self.prefix}', data_file='{self.data_file}', data_type='{self.data_type}')"


def parse_regulator_types(regulator_types_str):
    """
    Parse regulator types string into RegulatorConfig objects
    
    Format: 'Prefix1:DataFile1[:DataType1],Prefix2:DataFile2[:DataType2]'
    Example: 'TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,ClinicalParams:clinical.txt:d'
    
    DataType is optional: 'c' for continuous (default), 'd' for discrete/binary.
    
    Note: DataFile is stored for reference but not used in network generation.
          The prefix is used to find LemonTree output files (e.g., TFs.allreg.txt)
    """
    configs = []
    for config_str in regulator_types_str.split(','):
        config_str = config_str.strip()
        if ':' not in config_str:
            print(f"Warning: Invalid regulator config format: '{config_str}' (expected Prefix:DataFile[:DataType])")
            continue
        
        parts = config_str.split(':')
        prefix = parts[0].strip()
        data_file = parts[1].strip() if len(parts) > 1 else ''
        data_type = parts[2].strip().lower() if len(parts) > 2 else 'c'
        if data_type not in ('c', 'd'):
            print(f"Warning: Invalid data type '{data_type}' for {prefix} - using 'c' (continuous). Valid: 'c' or 'd'")
            data_type = 'c'
        configs.append(RegulatorConfig(prefix.strip(), data_file, data_type))
    
    return configs

def main():
    parser = argparse.ArgumentParser(description='Convert LemonTree results to networks')
    parser.add_argument('--input_dir', required=True, help='Input directory containing LemonTree results')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--run_id', required=True, help='Unique run identifier for tracking')
    parser.add_argument('--coherence_threshold', type=float, default=0.6, help='Coherence threshold for module filtering')
    
    # Regulator configuration
    parser.add_argument('--regulator_types', type=str, 
                       default='TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt',
                       help='Comma-separated list of Prefix:DataFile pairs (e.g., TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Lipids:Lipidomics.txt)')
    
    # Regulator selection method and thresholds
    parser.add_argument('--regulator_selection_method', type=str, default='percentage', 
                       choices=['percentage', 'fold_per_module'],
                       help='Method for selecting regulators')
    parser.add_argument('--top_n_percent_regulators', type=float, default=2.0, 
                       help='Top N percent of regulators (for percentage method)')
    parser.add_argument('--regulator_fold_cutoff', type=float, default=2.0,
                       help='Fold cutoff above random max (for fold_per_module method)')
    parser.add_argument('--skip_coherence_filtering', action='store_true', default=False, 
                       help='Skip coherence-based module filtering and use all modules')
    
    args = parser.parse_args()
    
    # Parse regulator types
    regulator_configs = parse_regulator_types(args.regulator_types)
    
    # Determine method suffix and threshold based on selection method
    if args.regulator_selection_method == 'percentage':
        threshold_value = args.top_n_percent_regulators
        method_suffix = f"top{threshold_value}pct"
    elif args.regulator_selection_method == 'fold_per_module':
        threshold_value = args.regulator_fold_cutoff
        method_suffix = f"fold{threshold_value}x_permodule"
    
    # Set working directory
    os.chdir(args.output_dir)
    
    # Create output directories
    os.makedirs('ModuleViewer_files', exist_ok=True)
    os.makedirs('Networks', exist_ok=True)
    
    # Define file paths
    # file used for module coherence (must not contain duplicate symbols)
    expression_dataset = os.path.join(args.input_dir, 'Preprocessing', 'LemonPreprocessed_expression.txt')
    # full omics dataset used elsewhere (not for coherence)
    expression_complete = os.path.join(args.input_dir, 'Preprocessing', 'LemonPreprocessed_complete.txt')
    DESeq_groups = os.path.join(args.input_dir, 'Preprocessing', 'DESeq_groups.txt')
    clusterfile = os.path.join(args.input_dir, 'Lemon_out', 'tight_clusters.txt')
    
    coherence_threshold = args.coherence_threshold
    skip_coherence_filtering = args.skip_coherence_filtering
    
    # Set file paths for each regulator type and filter to available ones
    for config in regulator_configs:
        config.set_files(args.input_dir)
    
    available_configs = [c for c in regulator_configs if c.files_exist()]
    missing_configs = [c for c in regulator_configs if not c.files_exist()]
    
    # Print configuration
    print(f"\n{'='*80}")
    print(f"LEMONTREE TO NETWORK CONVERSION")
    print(f"{'='*80}")
    print(f"Run ID: {args.run_id}")
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Regulator selection method: {args.regulator_selection_method}")
    print(f"Method suffix: {method_suffix}")
    print(f"Threshold value: {threshold_value}")
    print(f"Coherence threshold: {coherence_threshold}")
    print(f"Skip coherence filtering: {skip_coherence_filtering}")
    
    print(f"\n{'='*80}")
    print(f"REGULATOR CONFIGURATION")
    print(f"{'='*80}")
    print(f"Available regulator types: {len(available_configs)}")
    for config in available_configs:
        print(f"  [OK] {config.type}: {config.prefix}")
    
    if missing_configs:
        print(f"\nMissing regulator types: {len(missing_configs)}")
        for config in missing_configs:
            print(f"  [MISSING] {config.type}: {config.prefix}")
            print(f"    Expected: {config.allreg}")
    
    if not available_configs:
        print("\nERROR: No valid regulator files found!")
        sys.exit(1)
    
    # Step 1: Create sample mapping file
    print(f"\n{'='*80}")
    print("STEP 1: Creating sample mapping file")
    print(f"{'='*80}")
    create_sample_mapping(DESeq_groups)
    
    # Step 2: Select regulators using unified method for all regulator types (on ORIGINAL scores)
    print(f"\n{'='*80}")
    print(f"STEP 2: REGULATOR SELECTION ({args.regulator_selection_method.upper()})")
    print("Selection performed on ORIGINAL LemonTree scores before normalization")
    print(f"{'='*80}")
    
    reg_files = []
    for config in available_configs:
        print(f"\nProcessing {config.type} regulators ({config.prefix})...")
        config.selected_file = f'./Lemon_out/{config.prefix}.{method_suffix}.txt'
        
        if args.regulator_selection_method == 'percentage':
            selectregulators_percentage(
                config.allreg,
                config.randomreg,
                config.selected_file,
                threshold_value
            )
        elif args.regulator_selection_method == 'fold_per_module':
            selectregulators_scorecutoff_permodule(
                config.allreg,
                config.randomreg,
                threshold_value,
                config.selected_file
            )
        
        # Create list file
        list_output = f'./ModuleViewer_files/{config.prefix}.selected_regs_list.txt'
        createListFile(config.selected_file, list_output)
        
        # Create distribution plot
        create_reg_distribution(list_output)
        
        reg_files.append(config.selected_file)
    
    # Step 3: Normalize selected regulator scores (global sum normalization)
    print(f"\n{'='*80}")
    print(f"STEP 3: NORMALIZING SELECTED REGULATOR SCORES (SUM NORMALIZATION)")
    print("Global sum normalization divides scores by the TOTAL sum of selected regulator scores,")
    print("then scales by 100 to make numbers more interpretable")
    print("Normalized scores OVERWRITE original selection files for downstream use")
    print(f"{'='*80}")
    
    for config in available_configs:
        print(f"\nNormalizing selected {config.type} regulators ({config.prefix})...")
        normalize_selected_regulators_sum(config.selected_file)
    
    # Step 4: Create list files and distribution plots for selected regulators
    print(f"\n{'='*80}")
    print("STEP 4: Creating list files and distribution plots")
    print(f"{'='*80}")
    # Note: List files and plots were already created during selection (Step 2)
    # This step label is for documentation purposes only
    
    # Step 5: Create list file for clusters
    print(f"\n{'='*80}")
    print("STEP 5: Creating cluster list files")
    print(f"{'='*80}")
    createListFile(clusterfile, './ModuleViewer_files/clusters_list.txt')
    
    # Copy selected regulator score files to ModuleViewer_files for downstream use
    for config in available_configs:
        if config.selected_file and os.path.exists(config.selected_file):
            dest = f'./ModuleViewer_files/{config.prefix}.selected_regulators_scores.txt'
            shutil.copy2(config.selected_file, dest)
            print(f"  Copied selected regulator scores: {dest}")
    
    # Step 6: Read clusters and create cluster2gene mapping
    print(f"\n{'='*80}")
    print("STEP 6: Loading module-gene mappings")
    print(f"{'='*80}")
    clusters_list = pd.read_csv('./ModuleViewer_files/clusters_list.txt', sep='\t', header=None)
    clusters_list[0] = clusters_list[0].astype(str)
    cluster2gene = dict(zip(clusters_list[0], clusters_list[1]))
    print(f"Found {len(cluster2gene)} modules")
    
    # Step 7: Calculate module coherence
    print(f"\n{'='*80}")
    print("STEP 7: Calculating module coherence")
    print(f"{'='*80}")
    if not skip_coherence_filtering:
        # use expression-only file for coherence calculation
        low_coherence_modules, coherence_scores, module_eigengenes = calculate_module_coherence(
            expression_dataset, cluster2gene, coherence_threshold
        )
        
        # Filter cluster2gene dictionary
        print(f"Original number of modules: {len(cluster2gene)}")
        filtered_cluster2gene = {k: v for k, v in cluster2gene.items() if str(k) not in low_coherence_modules}
        cluster2gene = filtered_cluster2gene
        print(f"Modules after coherence filtering: {len(cluster2gene)}")
    else:
        print(f"Skipping coherence filtering - using all {len(cluster2gene)} modules")
        # Still calculate coherence scores for informational purposes, but don't filter
        low_coherence_modules, coherence_scores, module_eigengenes = calculate_module_coherence(
            expression_dataset, cluster2gene, coherence_threshold
        )
        # Keep all modules (no filtering)
        low_coherence_modules = []
    
    # Save coherence scores
    coherence_df = pd.DataFrame(list(coherence_scores.items()), columns=['Module', 'Coherence_Score'])
    coherence_df = coherence_df.sort_values('Coherence_Score', ascending=False)
    coherence_df.to_csv(f'./Networks/Module_coherence_scores.txt', sep='\t', index=False)
    print(f"Saved coherence scores to: ./Networks/Module_coherence_scores.txt")
    
    # Step 8: Prioritize modules
    filtered_module_names = list(cluster2gene.keys())
    order = Prioritize_modules_coherence(coherence_scores, filtered_module_names)

    
    # Step 9: Define specific modules to use
    print(f"\n{'='*80}")
    print("STEP 9: Defining modules for network construction")
    print(f"{'='*80}")
    specific_modules = list(order.keys())
    specific_modules = [str(module) for module in specific_modules]
    print(f'{len(specific_modules)} modules in specific_modules')
    
    # Save specific modules
    with open(f'./Networks/specific_modules.txt', 'w') as handle:
        for module in specific_modules:
            handle.write(module + '\n')
    
    # Step 10: Print summary
    print_summary(coherence_scores, low_coherence_modules, cluster2gene, specific_modules, 
                  coherence_threshold, method_suffix, args.regulator_selection_method)
    
    # Step 11: Build network
    print(f"\n{'='*80}")
    print("STEP 11: Building regulatory network")
    print(f"{'='*80}")
    # Build list of selected regulator files
    reg_files = [config.selected_file for config in available_configs]
    network_df = build_network(reg_files, cluster2gene, available_configs,
                               specific_modules_list=specific_modules, 
                               method_suffix=method_suffix)
    
    # Step 12: Create regulator-target mappings
    print(f"\n{'='*80}")
    print("STEP 12: Creating regulator-target mappings")
    print(f"{'='*80}")
    regulator_mappings, regulator2module = create_regulator_mappings(network_df)
    
    save_regulator_targets(regulator_mappings, method_suffix, 
                          len(specific_modules), regulator_configs)
    
    # Step 13: Create ranked regulator files
    print(f"\n{'='*80}")
    print("STEP 13: Creating ranked regulator files")
    print(f"{'='*80}")
    for reg_type in network_df['Type'].unique():
        type_network = network_df[network_df['Type'] == reg_type]
        create_ranked_reg_file(type_network, len(specific_modules), regulator_mappings, 
                              regulator2module, method_suffix)
    
    # Step 14: Create Cytoscape files
    print(f"\n{'='*80}")
    print("STEP 14: Creating Cytoscape network files")
    print(f"{'='*80}")
    createCytoscapeFiles(
        modules_to_use=specific_modules,
        reg_files={config.prefix: config.selected_file for config in available_configs},
        method_suffix=method_suffix,
        regulator_mappings=regulator_mappings
    )
    
    # Step 15: Create filtered distribution plots
    print(f"\n{'='*80}")
    print("STEP 15: Creating filtered distribution plots")
    print(f"{'='*80}")
    for config in available_configs:
        list_file = f'{config.prefix}.selected_regs_list.txt'
        create_reg_distribution_top_filtered(list_file, specific_modules, method_suffix)
    
    print(f"\n{'='*80}")
    print(f"Network generation completed successfully for run: {args.run_id}")
    print(f"{'='*80}\n")

def create_sample_mapping(DESeq_groups):
    """Create sample mapping file for ModuleViewer - supports multiple metadata annotations
    
    Generates an MVF file with multiple metadata types (diagnosis, sex, age, batch, etc.)
    that can be displayed as separate annotation tracks in module heatmaps.
    """
    annotations = pd.read_csv(DESeq_groups, sep='\t', index_col=0)
    
    # Get sample IDs
    sample_ids = annotations.index.tolist()
    print(f"{len(sample_ids)} samples included in the analysis.")
    
    # Define color palettes for different metadata types
    colors_categorical = ['RED', 'BLUE', 'GREEN', 'ORANGE', 'PURPLE', 'YELLOW', 'PINK', 'CYAN', 
                         'BROWN', 'MAGENTA', 'TEAL', 'LIME', 'NAVY', 'MAROON', 'OLIVE', 'CORAL']
    
    # Color palettes for specific known metadata types
    colors_sex = {'M': 'BLUE', 'F': 'PINK', 'Male': 'BLUE', 'Female': 'PINK', 
                  'male': 'BLUE', 'female': 'PINK', 'Unknown': 'GREY'}
    
    # Nice labels for common metadata columns
    default_labels = {
        'diagnosis': 'Clinical Status',
        'sex': 'Biological Sex',
        'age': 'Age Group',
        'batch': 'Batch',
        'biopsy_location': 'Biopsy Location',
        'tissue': 'Tissue Type',
        'condition': 'Condition'
    }
    
    # Process all available metadata columns
    mvf_sections = []
    available_columns = []
    
    for col in annotations.columns:
        # Skip if all values are NA or empty
        if annotations[col].isna().all():
            print(f"Skipping column '{col}': all values are NA")
            continue
            
        # Get unique values (excluding NAs)
        unique_values = sorted([str(v) for v in annotations[col].dropna().unique()])
        
        if len(unique_values) == 0:
            print(f"Skipping column '{col}': no valid values")
            continue
        
        if len(unique_values) > 20:
            print(f"Skipping column '{col}': too many unique values ({len(unique_values)} > 20)")
            continue
            
        available_columns.append(col)
        print(f"Found metadata column '{col}' with {len(unique_values)} unique values: {', '.join(unique_values)}")
        
        # Determine color mapping based on column type
        if col.lower() in ['sex', 'gender']:
            color_map = {}
            for val in unique_values:
                color_map[val] = colors_sex.get(val, 'GREY')
        else:
            color_map = {}
            for i, val in enumerate(unique_values):
                color_map[val] = colors_categorical[i % len(colors_categorical)]
        
        # Create MVF section for this metadata type
        section = []
        section.append(f'::TYPE={col}')
        section.append('::GLOBAL')
        section.append('::VALUES=color')
        section.append('::OBJECT=CONDITIONS')
        
        # Create legend
        legend_items = [f"{val}:{color_map[val]}" for val in unique_values]
        label = default_labels.get(col, col.replace('_', ' ').title())
        legend_str = '|'.join(legend_items) + f'\t{label}'
        section.append(f'::LEGEND={legend_str}')
        
        # Create sample-to-color mapping
        sample_to_map = ""
        for sample_id in sample_ids:
            value = str(annotations.loc[sample_id, col]) if pd.notna(annotations.loc[sample_id, col]) else 'Unknown'
            color = color_map.get(value, 'GREY')
            sample_to_map += '|' + str(sample_id) + ':' + color
        section.append(sample_to_map)
        
        mvf_sections.append('\n'.join(section))
    
    # Write complete MVF file with all metadata types
    with open('./ModuleViewer_files/sample_mapping.mvf', 'w') as handle:
        handle.write('\n---\n'.join(mvf_sections))
    
    print(f"Created sample mapping with {len(available_columns)} metadata types: {', '.join(available_columns)}")

def selectregulators_scorecutoff_permodule(allregfile, randomregfile, fold_cutoff, outfile, ensemble_to_symbol=False):
    
    # Create a dictionary with the highest random score per module
    with open(randomregfile, 'r') as handle:
        random_scores = {}
        for line in handle:
            module = line.rstrip().split('\t')[1]
            if module in random_scores:
                if float(line.rstrip().split('\t')[2]) > random_scores[module]:
                    random_scores[module] = float(line.rstrip().split('\t')[2])
            else:
                random_scores[module] = float(line.rstrip().split('\t')[2])
        
    # Print max random score across all modules
    print('Max random score across all modules:')
    max_random_overall = max(random_scores.values())
    print(max_random_overall)

    # print(random_scores)
    
    with open(allregfile, 'r') as handle2:
        module2regulators = {} # This will be a dictionary with as key the module number and value the regulator name
        
        all_scores = []
        for line in handle2:
            
            line = line.strip().split('\t')
            module = line[1]
            max_random = random_scores[module]
            all_scores.append(float(line[2])) # Save score in list
            if float(line[2]) >= fold_cutoff*float(max_random):
                if module in module2regulators: # Add regulator
                    to_append = line[0] + '|' + line[2]
                    module2regulators[module].append(to_append)
                else:
                    to_append = line[0] + '|' + line[2]
                    module2regulators[module] = []
                    module2regulators[module].append(to_append)


    all_scores.sort(reverse=True)
    print(all_scores)
    
    # Convert regulator ensembl id's to gene symbols, store mapping in dictionary
    ensemble_mapping = {}
    try:
        with open('./Preprocessing/LemonPreprocessed_complete.txt', 'r') as handle5:
            next(handle5)
            for line in handle5:
                line = line.rstrip().split('\t')
                ensembl = line[0]
                symbol = line[1]
                ensemble_mapping[ensembl] = symbol
    except:
        pass
    
    # Create output regulator file and print some diagnosticsModules
    with open(outfile, 'w') as handle3:
        out = 'Regulator' + '\t' + 'Target' + '\t' + 'Score' + '\t' + 'Overall_rank' + '\n'
        for element in module2regulators:
            regs_for_module = module2regulators[element]
            print(f'{len(regs_for_module)} remain for module {element} with a cutoff of {fold_cutoff} > random')
            for reg in regs_for_module:
                regulator = reg.split('|')[0]
                score = reg.split('|')[1]
                rank = str(all_scores.index(float(score)) + 1)
                if ensemble_to_symbol == False:
                    out = out + regulator + '\t' + element + '\t' + score + '\t' + rank +'\n'
                else:
                    out = out + ensemble_mapping[regulator] + '\t' + element + '\t' + score + '\t' + rank + '\n'
        handle3.write(out)


def normalize_selected_regulators_sum(selected_file):
    """
    Normalize already-selected regulator scores by dividing by the total sum of scores
    across all modules, then scale by 100.
    This is applied AFTER selection to make scores comparable across regulator types.
    
    Global sum normalization converts scores to percentages that sum to ~100 across all modules:
    - Each score represents its proportion of the total regulatory activity across all modules
    - Scores across all modules sum to ~100
    - Preserves relative differences between regulators across modules
    - Values are in [0, 100] range
    
    Parameters:
    - selected_file: path to selected regulator file (4 columns: Regulator, Target, Score, Overall_rank)
    
    Returns:
    - normalized_df
    
    The normalized scores OVERWRITE the original file so all downstream steps use normalized scores.
    """
    # Read the selected regulator file
    selected = pd.read_csv(selected_file, sep='\t')
    
    # Store original score range for reporting
    orig_min = selected['Score'].min()
    orig_max = selected['Score'].max()
    orig_sum = selected['Score'].sum()
    
    # Apply GLOBAL sum normalization: divide by total sum across all modules and scale by 100
    total_sum = selected['Score'].sum()
    if total_sum > 0:
        selected['Score'] = selected['Score'] / total_sum * 1000.0
        # round to nearest integer
        selected['Score'] = selected['Score'].round()
    else:
        # Edge case: if total sum is zero, distribute 100 equally across entries
        selected['Score'] = 1000.0 / len(selected)
    
    # Save normalized scores back to the SAME file
    selected.to_csv(selected_file, sep='\t', index=False)
    
    print(f"  Normalized {len(selected)} selected regulator-module pairs")
    print(f"  Original range: [{orig_min:.4f}, {orig_max:.4f}]")
    print(f"  Original sum (all modules): {orig_sum:.4f}")
    print(f"  Normalized sum (all modules): {selected['Score'].sum():.4f} (scaled to sum=1000 across all modules)")
    print(f"  Normalized range: [{selected['Score'].min():.6f}, {selected['Score'].max():.6f}]")
    print(f"  Note: Scores now sum to ~100 across all modules (global normalization and scaling)")
    print(f"  Overwritten: {selected_file}")
    
    return selected


def selectregulators_percentage(allregfile, randomregfile, outfile, percentage, ensemble_to_symbol=False):
    
    # Create a dictionary with the highest random score per module
    with open(randomregfile, 'r') as handle:
        random_scores = {}
        for line in handle:
            module = line.rstrip().split('\t')[1]
            if module in random_scores:
                if float(line.rstrip().split('\t')[2]) > random_scores[module]:
                    random_scores[module] = float(line.rstrip().split('\t')[2])
            else:
                random_scores[module] = float(line.rstrip().split('\t')[2])
        
    # Print max random score across all modules
    print('Max random score across all modules:')
    max_random_overall = max(random_scores.values())
    print(max_random_overall)


    allreg = pd.read_csv(allregfile, sep = '\t', header = None)
    allreg.columns = ['Regulator', 'Module', 'Score']
    allreg_filtered = allreg[allreg['Score'] >= allreg['Score'].quantile(1 - (percentage/100))]
    allreg_filtered = allreg_filtered[allreg_filtered['Score'] > max_random_overall]

    module2regulators = {}
    all_scores = []
    # create dict from Module to list of regulators
    for index, row in allreg_filtered.iterrows():
        module = str(row['Module'])
        regulator = row['Regulator']
        score = row['Score']
        all_scores.append(score)
        to_append = f"{regulator}|{score}"
        if module in module2regulators:
            module2regulators[module].append(to_append)
        else:
            module2regulators[module] = [to_append]

    
    all_scores.sort(reverse=True)
    print(all_scores)
    
    # Convert regulator ensembl id's to gene symbols, store mapping in dictionary
    ensemble_mapping = {}
    try:
        with open('./Preprocessing/LemonPreprocessed_complete.txt', 'r') as handle5:
            next(handle5)
            for line in handle5:
                line = line.rstrip().split('\t')
                ensembl = line[0]
                symbol = line[1]
                ensemble_mapping[ensembl] = symbol
    except:
        pass
    
    # Create output regulator file and print some diagnosticsModules
    with open(outfile, 'w') as handle3:
        out = 'Regulator' + '\t' + 'Target' + '\t' + 'Score' + '\t' + 'Overall_rank' + '\n'
        for element in module2regulators:
            regs_for_module = module2regulators[element]
            print(f'{len(regs_for_module)} remain for module {element} with a cutoff of {percentage} > random')
            for reg in regs_for_module:
                regulator = reg.split('|')[0]
                score = reg.split('|')[1]
                rank = str(all_scores.index(float(score)) + 1)
                if ensemble_to_symbol == False:
                    out = out + regulator + '\t' + element + '\t' + str(score) + '\t' + rank +'\n'
                else:
                    out = out + ensemble_mapping[regulator] + '\t' + element + '\t' + score + '\t' + rank + '\n'
        handle3.write(out)

        
def create_regulator_heatmap(method_suffix):
    """Create regulator heatmap"""
    try:
        # Load regulator-module score data
        regulator_file = f'./Lemon_out/Metabolite.{method_suffix}.txt'
        df = pd.read_csv(regulator_file, sep='\t')

        # Select only relevant columns
        df = df[['Regulator', 'Target', 'Score']]

        # Get top 15 regulators based on highest scores
        top_regulators = df.groupby('Regulator')['Score'].max().nlargest(15).index

        # Filter data to keep only top regulators
        df = df[df['Regulator'].isin(top_regulators)]

        # Pivot to create a matrix
        heatmap_df = df.pivot(index='Target', columns='Regulator', values='Score')

        # Fill NaN values with 0
        heatmap_df = heatmap_df.fillna(0)

        # Create a custom color map
        cmap = LinearSegmentedColormap.from_list('custom_cmap', ['white', 'red'])

        # Set annotations: show value only if it's not zero
        annot = heatmap_df.applymap(lambda x: '' if x == 0 else f'{x:.2f}')

        # Plot heatmap
        plt.figure(figsize=(9, 12))
        sns.heatmap(heatmap_df, cmap=cmap, annot=annot, fmt='', linewidths=0.5, linecolor='black', cbar=False, 
                    annot_kws={'color': 'black', 'size': 8})

        plt.title('Top 15 Metabolite Regulators Across Gene Modules')
        plt.xlabel('Regulators')
        plt.ylabel('Gene Modules')
        plt.xticks(rotation=45, ha='right')

        plt.savefig(f'./Metabolite.{method_suffix}_top15_module_scores.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Warning: Could not create regulator heatmap: {e}")

def createListFile(input_file, output_file):
    """Create list files for visualization in moduleviewer - Cell 6 from notebook"""
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Warning: {input_file} not found, skipping createListFile")
        return
    
    # Reading the input file into a Pandas DataFrame
    df = pd.read_csv(input_file, delimiter='\t', header=None)
    print(f"Processing {input_file}: {df.shape[1]} columns")
    
    # Check if df has 2 columns, then it is the cluster file
    if df.shape[1] == 2:
        df.columns = ['gene', 'cluster']
    elif df.shape[1] == 3: # Then it's a topreg file
        df.columns = ['regulator', 'target', 'score']
        df['cluster'] = df['target']
        df['gene'] = df['regulator']
        df = df[['gene', 'cluster']]
    elif df.shape[1] == 4: # Then it's a fold reg file
        df.columns = ['regulator', 'target', 'score', 'overall_rank']
        df['cluster'] = df['target']
        df['gene'] = df['regulator']
        df = df[['gene', 'cluster']]

    # Grouping by cluster and aggregating genes into a '|' separated string
    clusters_and_genes = df.groupby('cluster')['gene'].agg(lambda x: '|'.join(x)).to_dict()
    if 'Target' in clusters_and_genes:
        del clusters_and_genes['Target']
    
    # Formatting the output file
    out = [f"{cluster}\t{genes}\n" for cluster, genes in clusters_and_genes.items()]

    # Writing the output to a new file
    with open(output_file, 'w') as handle:
        handle.write(''.join(out))
    
    print(f"Written {len(clusters_and_genes)} clusters to {output_file}")

def create_reg_distribution(reg_file_list):
    """Create regulator distribution plot - Cell 7 from notebook"""
    try:
        with open('./ModuleViewer_files/' + reg_file_list, 'r') as handle:
            regs = {}
            for line in handle:
                line = line.rstrip().split('\t')
                if len(line) >= 2:
                    regulators = line[1].split('|')
                    for regulator in regulators:
                        if regulator in regs:
                            regs[regulator] += 1
                        else:
                            regs[regulator] = 1

        # Sort keys based on counts in decreasing order
        sorted_keys = sorted(regs, key=regs.get, reverse=True)

        # Only keep the top 15
        sorted_keys = sorted_keys[:15]

        # Extract sorted keys and counts
        sorted_counts = [regs[key] for key in sorted_keys]

        # Plotting
        bars = plt.bar(sorted_keys, sorted_counts, width=0.5)

        plt.xlabel('Keys')
        plt.ylabel('Counts')
        plt.title('Distribution of Counts for Each Key')

        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.tight_layout()

        # Save the plot as a PNG in Networks directory
        plot_filename = os.path.basename(reg_file_list)[:-4] + '_distribution_plot.png'
        plot_path = os.path.join('Networks', plot_filename)
        plt.savefig(plot_path)
        plt.close()
        print(f"Saved distribution plot to: {plot_path}")
        
    except Exception as e:
        print(f"Warning: Could not create distribution plot for {reg_file_list}: {e}")

def calculate_module_coherence(expression_file, clusters_dict, coherence_threshold=0.6):
    """Calculate eigengenes, kME, and module coherence scores

    The input should be the *expression-only* dataset (LemonPreprocessed_expression.txt)
    rather than the complete file, because the complete version may include duplicate
    gene symbols (one row per omics data type).  Using the expression file avoids
    pandas returning multiple rows for a single symbol (which led to mismatched
    vector lengths during Pearson correlation).
    """
    
    # Load expression data
    expression_data = pd.read_csv(expression_file, sep='\t')
    print(f"Loaded expression data: {expression_data.shape}")
    
    # Set gene symbols as index for easier subsetting
    expression_data_indexed = expression_data.set_index('symbol')
    # drop any duplicated symbols (should not happen for expression-only file)
    if expression_data_indexed.index.duplicated().any():
        dupes = expression_data_indexed.index[expression_data_indexed.index.duplicated()]
        print(f"Warning: {len(dupes)} duplicated gene symbols in expression file - keeping first occurrence")
        expression_data_indexed = expression_data_indexed[~expression_data_indexed.index.duplicated()]
    
    # Remove non-numeric columns (keep only sample columns)
    numeric_cols = expression_data_indexed.select_dtypes(include=[np.number]).columns
    expression_matrix = expression_data_indexed[numeric_cols]
    
    # collapse duplicates in expression matrix by averaging, since pearsonr fails otherwise
    if expression_matrix.index.duplicated().any():
        n_dups = expression_matrix.index.duplicated().sum()
        print(f"Warning: {n_dups} duplicated gene symbols found in expression data; collapsing by mean")
        expression_matrix = expression_matrix.groupby(expression_matrix.index).mean()
    
    coherence_scores = {}
    eigengenes = {}
    modules_to_remove = []
    
    print("Calculating module coherence scores...")
    
    for module, genes in clusters_dict.items():
        genes = genes.split('|')  # Ensure genes are split correctly
        try:
            # Get expression data for genes in this module
            module_genes_in_data = [gene for gene in genes if gene in expression_matrix.index]
            
            if len(module_genes_in_data) < 3:  # Need at least 3 genes for meaningful analysis
                print(f"Module {module}: Too few genes ({len(module_genes_in_data)}), marking for removal")
                modules_to_remove.append(str(module))
                coherence_scores[module] = 0.0
                continue
            
            # Extract expression data for module genes
            module_expression = expression_matrix.loc[module_genes_in_data]
            
            # Calculate eigengene using PCA (first principal component)
            pca = PCA(n_components=1)
            eigengene = pca.fit_transform(module_expression.T)
            eigengene_series = pd.Series(eigengene.flatten(), index=module_expression.columns)
            eigengenes[module] = eigengene_series
            
            # Calculate kME (correlation of each gene with eigengene)
            kme_values = []
            for gene in module_genes_in_data:
                gene_expression = module_expression.loc[gene]
                correlation, _ = pearsonr(gene_expression, eigengene_series)
                kme_values.append(abs(correlation))  # Use absolute correlation
            
            # Module coherence = mean absolute kME
            coherence_score = np.mean(kme_values)
            coherence_scores[module] = coherence_score
            
            # Check if module should be removed
            if coherence_score < coherence_threshold:
                modules_to_remove.append(str(module))
                print(f"Module {module}: Low coherence ({coherence_score:.3f} < {coherence_threshold}), marking for removal")

        except Exception as e:
            print(f"Error processing module {module}: {str(e)}")
            modules_to_remove.append(str(module))
            coherence_scores[module] = 0.0
    
    return modules_to_remove, coherence_scores, eigengenes

def Prioritize_modules_coherence(coherence_scores, filtered_modules):
    """Prioritize modules based on coherence scores"""
    
    # Filter coherence scores to only include modules that passed filtering
    filtered_coherence = {k: v for k, v in coherence_scores.items() if str(k) in filtered_modules}
    
    # Sort by coherence score (descending - higher coherence is better)
    modules_ordered = {k: v for k, v in sorted(filtered_coherence.items(), key=lambda item: item[1], reverse=True)}
    
    # Convert to df and sort
    df = pd.DataFrame(list(modules_ordered.items()), columns=['Module', 'coherence_score'])
    print(df)
    df.to_csv(f'Module_prioritization_coherence.txt', sep='\t', index=False)
    return modules_ordered

def print_summary(coherence_scores, low_coherence_modules, cluster2gene, specific_modules, coherence_threshold, method_suffix, regulator_selection_method):
    """Print summary of filtering and prioritization results"""
    print("="*80)
    print("SUMMARY OF MODULE FILTERING AND PRIORITIZATION")
    print("="*80)

    print(f"Original number of modules: {len(coherence_scores)}")
    print(f"Modules removed due to low coherence (< {coherence_threshold}): {len(low_coherence_modules)}")
    print(f"Modules remaining after coherence filtering: {len(cluster2gene)}")

    print(f"\nPrioritization method: Module coherence scores")
    print(f"All coherence-filtered modules used for network construction: {len(specific_modules)}")

    print(f"\nFinal network will be built using {len(specific_modules)} modules")
    print(f"Coherence threshold used: {coherence_threshold}")
    print(f"Regulator selection method: {regulator_selection_method}")
    print(f"Method suffix: {method_suffix}")
    print("\nIMPORTANT: No modules are removed based on expression significance.")
    print("Modules are only filtered by coherence, then ordered by coherence scores.")

    print("="*80)

def build_network(reg_files, cluster2gene, regulator_configs, specific_modules_list=None, method_suffix='top2pct'):
    """Build the network with selected regulators and modules"""
    # Initialize an empty list to store DataFrames
    dfs = []

    if specific_modules_list is not None:
        modules_to_keep = specific_modules_list
        modules_name = f'{len(modules_to_keep)}'
        print(f'Building network for {len(modules_to_keep)} selected modules')
    else:
        # Default: use all filtered modules
        modules_to_keep = list(cluster2gene.keys())
        modules_to_keep = [str(module) for module in modules_to_keep]
        modules_name = f'{len(modules_to_keep)}'
        print(f'Building network for all {len(modules_to_keep)} filtered modules')

    print(f'Modules to keep: {modules_to_keep[:10]}...')  # Show first 10
    
    # Create a mapping from file prefix to regulator type for quick lookup
    prefix_to_type = {config.prefix: config.type for config in regulator_configs}

    for reg_file in reg_files: # Loop over regulator files
        data = pd.read_csv(reg_file, sep='\t')
        if len(data.columns) == 3:
            data.columns = ['Regulator', 'Target', 'Score']
        elif len(data.columns) == 4:
            data.columns = ['Regulator', 'Target', 'Score', 'Overall_rank']

        data['Target'] = data['Target'].astype(str)
        
        # Detect regulator type from filename
        reg_type = 'Unknown'
        for prefix, rtype in prefix_to_type.items():
            if prefix in reg_file:
                reg_type = rtype
                break
        
        # Loop over this dataframe
        for index, row in data.iterrows():
            # Get Regulator
            regulator = row['Regulator']
            # Get module
            module = row['Target']
            
            if str(module) in modules_to_keep:
                genes = cluster2gene.get(module, [])
                if isinstance(genes, str):
                    genes = genes.split('|')  # Ensure genes are split correctly
                
                for gene in genes:
                    # Create a DataFrame for each gene with detected type
                    gene_df = pd.DataFrame({
                        'Regulator': [regulator], 
                        'Target': [gene], 
                        'Score': [row['Score']], 
                        'Lemon_module': str(module), 
                        'Type': [f'{reg_type}-gene']
                    })
                    dfs.append(gene_df)  # Append the DataFrame to the list
    
    # Concatenate all DataFrames in the list
    if dfs:
        network_df = pd.concat(dfs, ignore_index=True)
    else:
        print("Warning: No network edges found!")
        network_df = pd.DataFrame(columns=['Regulator', 'Target', 'Score', 'Lemon_module', 'Type'])

    # Write to file with method_suffix
    network_df.to_csv(f'./Networks/LemonNetwork_{method_suffix}_{modules_name}modules.txt', sep='\t', index=False)
    print(f"Network saved to: ./Networks/LemonNetwork_{method_suffix}_{modules_name}modules.txt")
    
    # Print network statistics
    print(f'The network contains {network_df["Target"].nunique()} unique genes.')
    for reg_type in network_df['Type'].unique():
        print(f'The network contains {network_df[network_df["Type"] == reg_type]["Regulator"].nunique()} unique {reg_type.split("-")[0]} regulators.')
    
    return network_df

def create_regulator_mappings(network_df):
    """Create dictionaries that map regulators to targets for all regulator types"""
    # Create a mapping for each regulator type dynamically
    regulator_mappings = {}
    
    for reg_type in network_df['Type'].unique():
        type_name = reg_type.split('-')[0]  # e.g., 'TF', 'Metabolite', 'Lipid'
        subset = network_df[network_df['Type'] == reg_type]
        
        reg2targets = {}
        for index, row in subset.iterrows():
            regulator = row['Regulator']
            target = row['Target']
            if regulator in reg2targets:
                reg2targets[regulator].append(target)
            else:
                reg2targets[regulator] = [target]
        
        regulator_mappings[type_name] = reg2targets
        print(f"  {type_name}: {len(reg2targets)} regulators")
    
    # Create regulator-to-module mapping (used by all types)
    regulator2module = {}
    for index, row in network_df.iterrows():
        regulator = row['Regulator']
        module = row['Lemon_module']
        if regulator in regulator2module:
            if module not in regulator2module[regulator]:
                regulator2module[regulator].append(module)
        else:
            regulator2module[regulator] = [module]

    return regulator_mappings, regulator2module

def save_regulator_targets(regulator_mappings, method_suffix, actual_n_modules, regulator_configs):
    """Save regulator-target mappings for all regulator types using the correct file prefix"""
    # Create a mapping from type to prefix
    type_to_prefix = {config.type: config.prefix for config in regulator_configs}
    
    for reg_type, reg2targets in regulator_mappings.items():
        df = pd.DataFrame(list(reg2targets.items()), columns=['Regulator', 'Targets'])
        df['Targets'] = df['Targets'].apply(lambda x: '|'.join(x))
        
        # Use the file prefix from configuration instead of the type name
        file_prefix = type_to_prefix.get(reg_type, reg_type)
        output_file = f'./Networks/{file_prefix}2targets_{method_suffix}_{actual_n_modules}_modules.txt'
        df.to_csv(output_file, sep='\t', index=False)
        print(f"  {reg_type}: {len(df)} regulators -> {output_file}")

def create_ranked_reg_file(network, n_modules_actual, regulator_mappings, regulator2module, method_suffix):
    """Create ranked regulator file for a specific regulator type"""
    type_name = network['Type'].unique()[0].split('-')[0]
    # remove Target and Type columns
    network = network.drop(['Target', 'Type'], axis=1)
    
    # Group by Regulator and calculate the mean score
    grouped_df = network.groupby(['Regulator', 'Lemon_module'])[['Score']].first().reset_index()
    
    # Group by 'Regulator' and sum the 'Score' column
    sum_scores_df = grouped_df.groupby('Regulator')[['Score']].sum().reset_index()

    # Sort the dataframe based on the sum of scores
    sum_scores_df = sum_scores_df.sort_values(by='Score', ascending=False)

    # Add columns with number of target modules, number of target genes, target modules and target genes
    reg2targets = regulator_mappings.get(type_name, {})
    
    for row in sum_scores_df.iterrows():
        regulator = row[1]['Regulator']
        modules = regulator2module.get(regulator, [])
        sum_scores_df.at[row[0], 'N_modules'] = int(len(modules))
        sum_scores_df.at[row[0], 'N_targets'] = int(len(reg2targets.get(regulator, [])))
        sum_scores_df.at[row[0], 'Modules'] = '|'.join(modules)
        sum_scores_df.at[row[0], 'Targets'] = '|'.join(reg2targets.get(regulator, []))
    
    print(f"  Top {type_name} regulators:")
    print(sum_scores_df.head())
    
    # Write to file
    output_file = f'./Networks/Network_{method_suffix}_{n_modules_actual}_modules_{type_name}_ranked_regulators.txt'
    sum_scores_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved: {output_file}")

def createCytoscapeFiles(modules_to_use=None, reg_files=None, method_suffix='top2pct', regulator_mappings=None):
    """Create Cytoscape files
    
    Args:
        modules_to_use: List of module IDs to include
        reg_files: Dict mapping prefix to file path for each regulator type
        method_suffix: Method suffix for output files
        regulator_mappings: Dict mapping regulator type to {regulator: [targets]} dict
        run_id: Run identifier
    """
    if modules_to_use is None or reg_files is None:
        return
    
    # Ensure modules are integers for consistency
    modules = [int(i) if str(i).isdigit() else str(i) for i in modules_to_use]
    print(f"Creating Cytoscape files for {len(modules)} modules: {modules[:10]}...")
    
    n_modules = len(modules)
    
    with open(f'./Networks/Cytoscape_network_{n_modules}_modules_{method_suffix}_filtered.txt', 'w') as handle:
        data_frames = []
        
        # Read all regulator files
        for prefix, file_path in reg_files.items():
            data = pd.read_csv(file_path, sep='\t')
            # Select rows in which the value of column 'Target' is present in the list of modules
            data = data.loc[data['Target'].isin(modules)]
            print(f"{prefix} regulators: {len(data)} interactions for filtered modules")
            data_frames.append(data)
        
        if data_frames:
            # concatenate the dataframes
            combined_data = pd.concat(data_frames, ignore_index=True)
            # Write the dataframe to a file
            combined_data.to_csv(handle, sep='\t', index=False)
            
            # Now we need to create an attributes file
            regulators = combined_data['Regulator'].unique().tolist()
            targets = combined_data['Target'].unique().tolist()
            
            with open(f'./Networks/Cytoscape_attributes_{n_modules}_modules_filtered.txt', 'w') as handle2:
                handle2.write('Name' + '\t' + 'Type' + '\n')
                
                # Write regulator types based on regulator_mappings
                for regulator in regulators:
                    # Find which type this regulator belongs to
                    reg_type_found = False
                    if regulator_mappings:
                        for reg_type, reg_dict in regulator_mappings.items():
                            if regulator in reg_dict.keys():
                                # Extract type name (e.g., 'TF-gene' -> 'TF')
                                type_name = reg_type.replace('-gene', '')
                                handle2.write(str(regulator) + '\t' + type_name + '\n')
                                reg_type_found = True
                                break
                    if not reg_type_found:
                        handle2.write(str(regulator) + '\t' + 'Unknown' + '\n')
                
                # Write module types
                for target in targets:
                    handle2.write(str(target) + '\t' + 'Module' + '\n')
            
            print(f"Created Cytoscape files for {len(regulators)} regulators and {len(targets)} modules")
        else:
            print("Warning: No data to write to Cytoscape files")

def create_reg_distribution_top_filtered(reg_file_list, modules_to_keep, method_suffix):
    """Create distribution plot for top regulators in filtered modules"""
    try:
        with open('./ModuleViewer_files/' + reg_file_list, 'r') as handle:
            regs = {}
            for line in handle:
                line = line.rstrip().split('\t')
                if len(line) < 2:
                    continue
                module = line[0]
                if str(module) not in [str(m) for m in modules_to_keep]:
                    continue
                regulators = line[1].split('|')
                for regulator in regulators:
                    if regulator in regs:
                        regs[regulator] += 1
                    else:
                        regs[regulator] = 1

        if not regs:
            print(f"Warning: No regulators found for file {reg_file_list}")
            return
        
        # Sort keys based on counts in decreasing order
        sorted_keys = sorted(regs, key=regs.get, reverse=True)

        # Extract sorted keys and counts
        sorted_counts = [regs[key] for key in sorted_keys]

        # Plotting
        fig, ax = plt.subplots()
        bars = ax.bar(sorted_keys[:15], sorted_counts[:15], width=0.5)

        ax.set_xlabel('Regulators')
        ax.set_ylabel('Number of Modules')
        ax.set_title(f'Distribution of Top 15 Regulators Across {len(modules_to_keep)} Filtered Modules')

        ax.set_xticks(range(len(sorted_keys[:15])))
        ax.set_xticklabels(sorted_keys[:15], rotation=45, ha='right', fontsize=6)
        plt.tight_layout()

        # Explicitly set the face color
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        # Save the plot with a white background
        plt.savefig(f'./Networks/{reg_file_list[:-4]}_distribution_plot_top30_filtered_{len(modules_to_keep)}modules.png', 
                    facecolor=fig.get_facecolor(), edgecolor='none', dpi=300)

        plt.close()
        print(f"Saved plot for {len(sorted_keys)} unique regulators")
        
    except Exception as e:
        print(f"Warning: Could not create filtered distribution plot for {reg_file_list}: {e}")

if __name__ == "__main__":
    main()
