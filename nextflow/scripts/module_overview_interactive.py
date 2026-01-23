#!/usr/bin/env python3

"""
Interactive Module Overview Generator Script
Enhanced version with megago clustering and interactive visualization

This script generates a comprehensive overview of LemonTree modules by integrating:
- Gene assignments per module
- Regulator assignments (TF, metabolites) per module  
- Pathway enrichment analysis results
- Expression-based module prioritization (differential analysis)
- Megago-based module clustering using functional similarity
- Interactive network visualization with plotly
- Module meta-clustering visualization
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse
import json
import warnings
import re
import subprocess
import glob
from scipy.stats import mannwhitneyu, kruskal, rankdata
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import plotly.graph_objects as go
import plotly.express as px
import plotly.offline as pyo
import plotly.subplots as sp
import networkx as nx
from collections import defaultdict, Counter
import traceback
import itertools
from concurrent.futures import ThreadPoolExecutor, as_completed

# Try to import matplotlib and seaborn for heatmap generation
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib/seaborn not available. Module expression heatmap will not be generated.")

# Try to import statsmodels for multiple testing correction
try:
    from statsmodels.stats.multitest import multipletests
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False
    print("Warning: statsmodels not available. Using fallback for multiple testing correction.")

# Try to import megago for functional clustering
try:
    import megago
    MEGAGO_AVAILABLE = True
    print("Megago available for functional clustering")
except ImportError:
    MEGAGO_AVAILABLE = False
    print("Warning: Megago not available. Functional clustering will use basic similarity metrics.")
    print("To enable megago clustering, install with: pip install megago")

warnings.filterwarnings('ignore')

def get_regulators(regfile):
    """
    Parse regulator/gene list files into dictionaries.
    
    Parameters:
    regfile (str): Path to the regulator/gene list file
    
    Returns:
    dict: Dictionary mapping module IDs to lists of regulators/genes
    """
    regs = {}
    try:
        with open(regfile) as f:
            for line in f:
                if line.strip():  # Skip empty lines
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        module = parts[0].strip()
                        regulators = parts[1].rstrip().split('|')
                        regs[module] = regulators
    except FileNotFoundError:
        print(f"Warning: File {regfile} not found. Creating empty dictionary.")
        regs = {}
    except Exception as e:
        print(f"Error reading {regfile}: {e}")
        regs = {}
    
    return regs

def load_enrichment_data(enrichment_dir, method='EnrichR'):
    """
    Load pathway enrichment results from either EnrichR or GSEA methods
    
    Parameters:
    enrichment_dir (str): Directory containing enrichment results
    method (str): 'EnrichR' or 'GSEA'
    
    Returns:
    dict: Dictionary with pathway type as keys and DataFrames as values
    """
    enrichment_data = {}

    try:
        # quick exit if enrichment dir doesn't exist
        if not os.path.isdir(enrichment_dir):
            return enrichment_data

        # collect CSV files under enrichment_dir
        csv_files = []
        for root, dirs, files in os.walk(enrichment_dir):
            for f in files:
                if f.lower().endswith('.csv'):
                    csv_files.append(os.path.join(root, f))

        # Debug: list discovered enrichment csv files
        if csv_files:
            print(f"Found enrichment CSV files ({len(csv_files)}):")
            for p in csv_files:
                print(f"  - {p}")
        else:
            print(f"No enrichment CSV files found under: {enrichment_dir}")

        if not csv_files:
            # no enrichment files found
            return enrichment_data

        filenames_lower = [os.path.basename(p).lower() for p in csv_files]

        # detect method if possible
        enrichr_candidates = [p for p,fn in zip(csv_files, filenames_lower) if 'enrichr' in fn or 'top_10_enriched_pathways' in fn]
        gsea_candidates = [p for p,fn in zip(csv_files, filenames_lower) if 'gsea' in fn or fn.startswith('gsea_')]

        if enrichr_candidates and not gsea_candidates:
            method = 'EnrichR'
            print('Auto-detected EnrichR method')
        elif gsea_candidates and not enrichr_candidates:
            method = 'GSEA'
            print('Auto-detected GSEA method')
        else:
            print(f'Both or ambiguous enrichment outputs found; using requested method: {method}')

        def _standardize_df(df):
            df = df.copy()
            # Module column
            mod_col = next((c for c in df.columns if c.lower() in ('module','cluster','moduleid','mod')), None)
            if mod_col and mod_col != 'Module':
                df = df.rename(columns={mod_col: 'Module'})
            # Term column
            term_col = next((c for c in df.columns if c.lower() in ('term','description','pathway','name')), None)
            if term_col and term_col != 'Term':
                df = df.rename(columns={term_col: 'Term'})
            # p.adjust alternatives
            padj_col = next((c for c in df.columns if c.lower() in ('padj','p.adjust','fdr','qvalue','p_adj','p.value')), None)
            if padj_col and padj_col != 'p.adjust':
                df = df.rename(columns={padj_col: 'p.adjust'})
            if 'Module' in df.columns:
                df['Module'] = df['Module'].astype(str).str.strip()
            return df

        # If user-provided EnrichR combined files exist, load them (preferred)
        up_name = 'enrichr_top_10_enriched_pathways_up_per_module.csv'
        down_name = 'enrichr_top_10_enriched_pathways_down_per_module.csv'
        # Also check for capitalized versions
        up_name_cap = 'Enrichr_top_10_enriched_pathways_up_per_module.csv'
        down_name_cap = 'Enrichr_top_10_enriched_pathways_down_per_module.csv'
        
        # Look for files that match the expected pattern
        up_candidates = [p for p in csv_files if 'enrichr' in os.path.basename(p).lower() and 'up' in os.path.basename(p).lower() and 'per_module' in os.path.basename(p).lower()]
        down_candidates = [p for p in csv_files if 'enrichr' in os.path.basename(p).lower() and 'down' in os.path.basename(p).lower() and 'per_module' in os.path.basename(p).lower()]
        
        # If we have multiple candidates, try to find the one that matches current module count
        up_path = None
        down_path = None
        
        if up_candidates:
            if len(up_candidates) == 1:
                up_path = up_candidates[0]
            else:
                print(f"Multiple up enrichment files found ({len(up_candidates)}), trying to match module count...")
                # Try to determine current module count from other files if possible
                current_module_count = None
                # Look for module assignment files or other indicators
                for root, dirs, files in os.walk(os.path.dirname(enrichment_dir)):
                    for f in files:
                        if 'module' in f.lower() and f.endswith('.txt'):
                            try:
                                with open(os.path.join(root, f), 'r') as file:
                                    lines = file.readlines()
                                    if lines:
                                        # Count unique module IDs
                                        modules = set()
                                        for line in lines[1:]:  # Skip header
                                            parts = line.strip().split('\t')
                                            if len(parts) >= 2:
                                                modules.add(parts[1])  # Assuming module ID is in second column
                                        if modules:
                                            current_module_count = len(modules)
                                            print(f"Detected {current_module_count} modules from {f}")
                                            break
                            except:
                                continue
                
                if current_module_count:
                    # Find file that contains the module count in filename
                    for candidate in up_candidates:
                        basename = os.path.basename(candidate)
                        if str(current_module_count) in basename:
                            up_path = candidate
                            print(f"Selected up file matching {current_module_count} modules: {basename}")
                            break
                
                if not up_path:
                    up_path = up_candidates[0]  # Fallback to first file
                    print(f"Using first up file as fallback: {os.path.basename(up_path)}")
        
        if down_candidates:
            if len(down_candidates) == 1:
                down_path = down_candidates[0]
            else:
                print(f"Multiple down enrichment files found ({len(down_candidates)}), trying to match module count...")
                if 'current_module_count' in locals() and current_module_count:
                    # Find file that contains the module count in filename
                    for candidate in down_candidates:
                        basename = os.path.basename(candidate)
                        if str(current_module_count) in basename:
                            down_path = candidate
                            print(f"Selected down file matching {current_module_count} modules: {basename}")
                            break
                
                if not down_path:
                    down_path = down_candidates[0]  # Fallback to first file
                    print(f"Using first down file as fallback: {os.path.basename(down_path)}")
        
        # Fallback to original method if no candidates found
        if not up_path:
            up_path = next((p for p in csv_files if os.path.basename(p).lower() == up_name or os.path.basename(p).lower() == up_name_cap.lower()), None)
        if not down_path:
            down_path = next((p for p in csv_files if os.path.basename(p).lower() == down_name or os.path.basename(p).lower() == down_name_cap.lower()), None)

        if up_path or down_path:
            parts = []
            for pth in (up_path, down_path):
                if pth:
                    try:
                        df = pd.read_csv(pth)
                        # Set direction based on filename
                        filename = os.path.basename(pth).lower()
                        if 'up' in filename:
                            df['__direction__'] = 'Up'
                        elif 'down' in filename:
                            df['__direction__'] = 'Down'
                        else:
                            df['__direction__'] = 'Up'  # Default to up if unclear
                        parts.append(df)
                    except Exception as e:
                        print(f"Warning: could not read enrichment file {pth}: {e}")

            if parts:
                combined = pd.concat(parts, ignore_index=True)
                combined = _standardize_df(combined)

                # If Database column present, split by database
                db_col = next((c for c in combined.columns if c.lower() in ('database','db','source')), None)
                if db_col:
                    for db_val, sub in combined.groupby(db_col):
                        v = str(db_val).lower()
                        if 'bp' in v or 'go' in v or 'biol' in v:
                            key = 'bp'
                        elif 'mf' in v:
                            key = 'mf'
                        elif 'cc' in v:
                            key = 'cc'
                        elif 'kegg' in v:
                            key = 'kegg'
                        elif 'react' in v:
                            key = 'reactome'
                        else:
                            key = 'other'

                        enrichment_data.setdefault(key, pd.DataFrame())
                        enrichment_data[key] = pd.concat([enrichment_data[key], _standardize_df(sub)], ignore_index=True)
                else:
                    # fallback: put into 'bp' to ensure downstream has something
                    enrichment_data['bp'] = combined

                # ensure keys exist
                for k in ('bp','mf','cc','kegg','reactome'):
                    enrichment_data.setdefault(k, pd.DataFrame(columns=['Module','Term','p.adjust']))

                return enrichment_data

        # General fallback: try to load any EnrichR/GSEA-like files and bucket by filename
        for pth in csv_files:
            fn = os.path.basename(pth).lower()
            try:
                df = pd.read_csv(pth)
            except Exception:
                continue

            std = _standardize_df(df)
            key = 'other'
            if 'gsea' in fn:
                # For GSEA files, check if they contain biological process terms
                # GSEA files have BP terms under 'BP' or 'ALL' database category
                if 'Database' in std.columns:
                    # Check if this file contains biological process terms
                    db_values = std['Database'].str.lower().unique()
                    if any('bp' in db_val or 'all' in db_val or 'go:' in db_val for db_val in db_values if isinstance(db_val, str)):
                        # Extract BP terms from GSEA file - prefer 'BP' over 'ALL'
                        bp_terms = pd.DataFrame()
                        if 'bp' in [str(x).lower() for x in db_values]:
                            bp_terms = std[std['Database'].str.lower() == 'bp'].copy()
                        elif 'all' in [str(x).lower() for x in db_values]:
                            bp_terms = std[std['Database'].str.lower() == 'all'].copy()

                        if not bp_terms.empty:
                            enrichment_data.setdefault('bp', pd.DataFrame())
                            enrichment_data['bp'] = pd.concat([enrichment_data['bp'], bp_terms], ignore_index=True)

                        # Extract other database categories
                        for db_val in db_values:
                            if isinstance(db_val, str):
                                db_lower = db_val.lower()
                                if db_lower == 'cc':
                                    key = 'cc'
                                elif db_lower == 'mf':
                                    key = 'mf'
                                elif 'kegg' in db_lower:
                                    key = 'kegg'
                                elif 'react' in db_lower:
                                    key = 'reactome'
                                else:
                                    continue

                                db_terms = std[std['Database'].str.lower() == db_lower]
                                if not db_terms.empty:
                                    enrichment_data.setdefault(key, pd.DataFrame())
                                    enrichment_data[key] = pd.concat([enrichment_data[key], db_terms], ignore_index=True)
                        continue  # Skip the general key assignment below

                # Fallback for GSEA files without Database column or other patterns
                if 'bp' in fn or 'biol' in fn:
                    key = 'bp'
                elif 'mf' in fn:
                    key = 'mf'
                elif 'cc' in fn:
                    key = 'cc'
                elif 'kegg' in fn:
                    key = 'kegg'
                elif 'react' in fn:
                    key = 'reactome'
            elif 'enrichr' in fn or 'top_10_enriched_pathways' in fn:
                # heuristic
                if 'bp' in fn or 'biol' in fn:
                    key = 'bp'
                elif 'mf' in fn:
                    key = 'mf'
                elif 'cc' in fn:
                    key = 'cc'
                elif 'kegg' in fn:
                    key = 'kegg'
                elif 'react' in fn:
                    key = 'reactome'

            enrichment_data.setdefault(key, pd.DataFrame())
            enrichment_data[key] = pd.concat([enrichment_data[key], std], ignore_index=True)

        # ensure expected keys exist
        for k in ('bp','mf','cc','kegg','reactome'):
            enrichment_data.setdefault(k, pd.DataFrame(columns=['Module','Term','p.adjust']))

    except Exception as e:
        print(f"Error loading enrichment data: {e}")

    return enrichment_data

def calculate_pathway_similarity_matrix(module_pathways):
    """
    Calculate similarity matrix between modules based on shared pathways
    
    Parameters:
    module_pathways (dict): Dictionary mapping module IDs to sets of pathway terms
    
    Returns:
    pandas.DataFrame: Similarity matrix
    """
    modules = list(module_pathways.keys())
    n_modules = len(modules)
    similarity_matrix = np.zeros((n_modules, n_modules))
    
    for i, mod1 in enumerate(modules):
        for j, mod2 in enumerate(modules):
            if i == j:
                similarity_matrix[i, j] = 1.0
            else:
                pathways1 = set(module_pathways[mod1])
                pathways2 = set(module_pathways[mod2])
                
                if len(pathways1) == 0 and len(pathways2) == 0:
                    similarity = 0.0
                elif len(pathways1) == 0 or len(pathways2) == 0:
                    similarity = 0.0
                else:
                    # Jaccard similarity
                    intersection = len(pathways1.intersection(pathways2))
                    union = len(pathways1.union(pathways2))
                    similarity = intersection / union if union > 0 else 0.0
                
                similarity_matrix[i, j] = similarity
    
    return pd.DataFrame(similarity_matrix, index=modules, columns=modules)

def create_megago_files(enrichment_data, output_dir):
    """
    Create MegaGO files with GO BP terms for each module
    
    Parameters:
    enrichment_data (dict): Dictionary with pathway type as keys and DataFrames as values
    output_dir (str): Output directory for megaGO files
    
    Returns:
    str: Directory path where megaGO files are created
    """
    megago_dir = os.path.join(output_dir, 'megaGO_files')
    os.makedirs(megago_dir, exist_ok=True)
    
    files_created = 0
    
    # Use biological process enrichment data for megaGO
    if 'bp' in enrichment_data and not enrichment_data['bp'].empty:
        print("Creating MegaGO files with GO BP terms for each module...")
        
        bp_data = enrichment_data['bp']
        
        for module in bp_data['Module'].unique():
            module_data = bp_data[bp_data['Module'] == module]
            bp_terms = module_data['Term'].tolist()
            
            if bp_terms:  # Only create file if module has BP terms
                # Clean GO terms by extracting GO IDs robustly or falling back to a cleaned term label
                import re
                cleaned_terms = []
                for term in bp_terms:
                    if not isinstance(term, str):
                        continue
                    # Try to find a GO ID anywhere in the string (e.g. "(...GO:0006954)" or "GO:0006954 - term")
                    m = re.search(r'GO:(\d+)', term)
                    if m:
                        go_id = m.group(1).zfill(7)
                        cleaned_terms.append(f"GO:{go_id}")
                        continue

                    # Fallback: take text before common separators like ' - ' or '(' or ';' or ','
                    fallback = re.split(r' - |\(|;|,', term)[0].strip()
                    if fallback:
                        cleaned_terms.append(fallback)
                
                if cleaned_terms:
                    filepath = os.path.join(megago_dir, f"{module}_BP_terms.txt")
                    with open(filepath, "w") as f:
                        f.write('GO_TERM\n')
                        f.write("\n".join(cleaned_terms))
                    
                    files_created += 1
        
        print(f"Created {files_created} GO BP term files for MegaGO clustering")
    else:
        print("Warning: No biological process enrichment data available - cannot create MegaGO files")
    
    return megago_dir if files_created > 0 else None

def run_single_megago_pair(file1, file2, megago_dir):
    """
    Run MegaGO on a single pair of files
    
    Parameters:
    file1 (str): Path to first BP terms file
    file2 (str): Path to second BP terms file
    megago_dir (str): Directory containing the files
    
    Returns:
    tuple: (file1_name, file2_name, output_text, error_text, returncode)
    """
    file1_name = os.path.basename(file1)
    file2_name = os.path.basename(file2)
    
    megago_command = f"megago {file1_name} {file2_name}"
    
    try:
        result = subprocess.run(
            megago_command,
            shell=True,
            capture_output=True,
            text=True,
            cwd=megago_dir,
            timeout=300  # 5 minute timeout per pair
        )
        return file1_name, file2_name, result.stdout, result.stderr, result.returncode
    except subprocess.TimeoutExpired:
        return file1_name, file2_name, "", "Timeout after 5 minutes", -1
    except Exception as e:
        return file1_name, file2_name, "", str(e), -1

def parse_single_pair_output(output_text, file1_name, file2_name):
    """
    Parse MegaGO output for a single pair to extract biological_process similarity score
    
    Parameters:
    output_text (str): Raw MegaGO output for a single pair
    file1_name (str): Name of first file
    file2_name (str): Name of second file
    
    Returns:
    float: Similarity score or None if parsing failed
    """
    if not output_text:
        return None
    
    # Pattern to extract biological_process score
    pattern = r'biological_process,([0-9.]+)'
    match = re.search(pattern, output_text)
    
    if match:
        return float(match.group(1))
    else:
        return None

def run_megago_clustering(megago_dir):
    """
    Run MegaGO command-line tool in parallel on all pairwise combinations
    
    Parameters:
    megago_dir (str): Directory containing megaGO files
    
    Returns:
    tuple: (similarity_matrix, module_ids) or (None, None) if failed
    """
    import glob
    
    # Find all BP_terms.txt files
    bp_files = glob.glob(os.path.join(megago_dir, "*_BP_terms.txt"))
    bp_files.sort()  # Sort to ensure consistent order
    
    if not bp_files:
        print("No BP_terms.txt files found for MegaGO analysis")
        return None, None
    
    print(f"Found {len(bp_files)} BP term files for MegaGO analysis")
    
    # Extract module IDs from file names
    module_ids = []
    for filepath in bp_files:
        filename = os.path.basename(filepath)
        module_id = filename.replace('_BP_terms.txt', '')
        module_ids.append(module_id)
    
    n_modules = len(module_ids)
    print(f"Processing {n_modules} modules: {module_ids}")
    
    # Generate all pairwise combinations
    file_pairs = list(itertools.combinations(bp_files, 2))
    print(f"Will run {len(file_pairs)} pairwise MegaGO comparisons")
    
    # Initialize similarity matrix
    similarity_matrix = np.eye(n_modules)  # Identity matrix (diagonal = 1.0)
    
    # Create mapping from filename to index
    file_to_index = {os.path.basename(f): i for i, f in enumerate(bp_files)}
    
    # Run MegaGO comparisons in parallel
    successful_comparisons = 0
    failed_comparisons = 0
    
    with ThreadPoolExecutor(max_workers=min(8, len(file_pairs))) as executor:
        # Submit all pairwise comparisons
        future_to_pair = {
            executor.submit(run_single_megago_pair, file1, file2, megago_dir): (file1, file2)
            for file1, file2 in file_pairs
        }
        
        # Process completed comparisons
        for future in as_completed(future_to_pair):
            file1, file2 = future_to_pair[future]
            try:
                file1_name, file2_name, stdout, stderr, returncode = future.result()
                
                if returncode == 0:
                    # Parse the output for this pair
                    score = parse_single_pair_output(stdout, file1_name, file2_name)
                    if score is not None:
                        idx1 = file_to_index[file1_name]
                        idx2 = file_to_index[file2_name]
                        similarity_matrix[idx1, idx2] = score
                        similarity_matrix[idx2, idx1] = score  # Make symmetric
                        print(f"  [OK] {module_ids[idx1]} vs {module_ids[idx2]}: {score:.4f}")
                        successful_comparisons += 1
                    else:
                        print(f"  [WARNING] Failed to parse output for {file1_name} vs {file2_name}")
                        failed_comparisons += 1
                else:
                    print(f"  [ERROR] MegaGO failed for {file1_name} vs {file2_name}: {stderr}")
                    failed_comparisons += 1
                    
            except Exception as e:
                print(f"  [ERROR] Error processing {os.path.basename(file1)} vs {os.path.basename(file2)}: {e}")
                failed_comparisons += 1
    
    print(f"MegaGO parallel execution completed:")
    print(f"  Successful comparisons: {successful_comparisons}")
    print(f"  Failed comparisons: {failed_comparisons}")
    
    if successful_comparisons == 0:
        print("No successful MegaGO comparisons - cannot create similarity matrix")
        return None, None
    
    # Save similarity matrix to file for inspection
    try:
        output_file = os.path.join(megago_dir, '../megago_similarity_matrix.csv')
        similarity_df = pd.DataFrame(similarity_matrix, index=module_ids, columns=module_ids)
        similarity_df.to_csv(output_file)
        print(f"\nSaved MegaGO similarity matrix to: {output_file}")
    except Exception as e:
        print(f"Warning: Could not save similarity matrix: {e}")
    
    return similarity_matrix, module_ids

def parse_megago_output(output_text, bp_files):
    """
    Parse megago output to extract biological_process similarity scores
    
    Parameters:
    output_text (str): Raw megago output
    bp_files (list): List of BP term files used
    
    Returns:
    tuple: (similarity_matrix, module_ids) or (None, None) if failed
    """
    if not output_text:
        print("No megago output available to parse.")
        return None, None
    
    # Pattern to extract sample comparisons and biological_process scores
    pattern = r'Results for sample (\d+) and (\d+).*?biological_process,([0-9.]+)'
    
    # Find all matches
    matches = re.findall(pattern, output_text, re.DOTALL)
    
    if not matches:
        print("No similarity scores found in the output.")
        print("Output preview:")
        print(output_text[:500] + "..." if len(output_text) > 500 else output_text)
        return None, None
    
    print(f"Found {len(matches)} similarity score pairs")
    
    # Extract module IDs from file names
    module_ids = []
    for filepath in bp_files:
        filename = os.path.basename(filepath)
        module_id = filename.replace('_BP_terms.txt', '')
        module_ids.append(module_id)
    
    n_modules = len(module_ids)
    print(f"Processing {n_modules} modules: {module_ids}")
    
    # Initialize similarity matrix
    similarity_matrix = np.eye(n_modules)  # Identity matrix (diagonal = 1.0)
    
    # Fill the matrix with extracted scores
    for match in matches:
        sample1, sample2, score = int(match[0]), int(match[1]), float(match[2])
        
        # Convert sample indices to module indices
        if sample1 < n_modules and sample2 < n_modules:
            similarity_matrix[sample1, sample2] = score
            similarity_matrix[sample2, sample1] = score  # Make symmetric
            print(f"  Module {module_ids[sample1]} vs Module {module_ids[sample2]}: {score:.4f}")
    
    return similarity_matrix, module_ids

def megago_cluster_modules(module_pathways, n_clusters=5, use_megago=True, enrichment_data=None, output_dir='.'):
    """
    Cluster modules using megaGO command-line tool if available, otherwise use pathway similarity
    
    Parameters:
    module_pathways (dict): Dictionary mapping module IDs to lists of pathway terms
    n_clusters (int): Number of clusters to create
    use_megago (bool): Whether to use megago for clustering
    enrichment_data (dict): Enrichment data for creating megaGO files
    output_dir (str): Output directory for temporary files
    
    Returns:
    tuple: (module_clusters dict, similarity_matrix) where module_clusters maps module IDs to cluster labels
    """
    if not module_pathways:
        print("No pathway data available for clustering")
        return {}, None
    
    modules = list(module_pathways.keys())
    
    # Skip clustering if requested or not enough clusters
    if n_clusters <= 1:
        print("Functional clustering disabled - assigning all modules to single cluster")
        return {mod: 1 for mod in modules}, None
    
    if not use_megago:
        print("MegaGO clustering disabled by user - using pathway similarity")
    else:
        print("\nAttempting MegaGO clustering...")
        
        # Try to use actual megaGO command-line tool
        if enrichment_data and 'bp' in enrichment_data:
            print(f"   Found biological process enrichment data with {len(enrichment_data['bp'])} entries")
            # Create megaGO files
            megago_dir = create_megago_files(enrichment_data, output_dir)
            
            if megago_dir:
                print(f"   Created megaGO files in: {megago_dir}")
                # Run megaGO clustering
                similarity_matrix, megago_module_ids = run_megago_clustering(megago_dir)
                
                if similarity_matrix is not None:
                    print("Successfully obtained MegaGO similarity matrix")
                    
                    # Convert similarity to distance
                    distance_matrix = 1 - similarity_matrix
                    
                    # Perform hierarchical clustering
                    if len(megago_module_ids) >= 2:
                        # Convert to condensed distance matrix for scipy
                        condensed_dist = squareform(distance_matrix, checks=False)
                        
                        # Hierarchical clustering
                        linkage_matrix = linkage(condensed_dist, method='ward')
                        
                        # Determine optimal number of clusters
                        if n_clusters == 'auto':
                            n_clusters = min(5, max(2, len(megago_module_ids) // 3))
                        
                        cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
                        
                        # Create module-to-cluster mapping
                        module_clusters = {}
                        for i, module_id in enumerate(megago_module_ids):
                            module_clusters[module_id] = cluster_labels[i]
                        
                        # Assign cluster 0 to modules without MegaGO data
                        for module in modules:
                            if module not in module_clusters:
                                module_clusters[module] = 0
                        
                        print(f"Created {n_clusters} module clusters using MegaGO semantic similarity")
                        
                        # Print cluster summary
                        cluster_counts = Counter(module_clusters.values())
                        for cluster_id, count in sorted(cluster_counts.items()):
                            print(f"  Cluster {cluster_id}: {count} modules")
                        
                        return module_clusters, similarity_matrix
                    
                else:
                    print("MegaGO clustering failed - falling back to pathway similarity")
            else:
                print("Could not create MegaGO files - falling back to pathway similarity")
        else:
            if not enrichment_data:
                print("No enrichment data available - falling back to pathway similarity")
            elif 'bp' not in enrichment_data:
                print(f"No biological process enrichment data (available keys: {list(enrichment_data.keys())}) - falling back to pathway similarity")
            else:
                print("No biological process enrichment data for MegaGO - falling back to pathway similarity")
    
    # Fall back to enhanced pathway similarity clustering
    print("Using pathway similarity for clustering...")
    
    # Method 1: Jaccard similarity on pathway terms
    jaccard_sim = calculate_pathway_similarity_matrix(module_pathways)

    
    # Combine similarities (average of Jaccard and weighted)
    similarity_matrix = (jaccard_sim.values)
    valid_modules = list(module_pathways.keys())

    print("Using enhanced pathway similarity clustering (Jaccard)")

    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix
    
    # Perform hierarchical clustering
    if len(valid_modules) < 2:
        return {mod: 0 for mod in modules}, None
    
    # Convert to condensed distance matrix for scipy
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method='ward')
    
    # Determine optimal number of clusters if not specified
    if n_clusters == 'auto':
        n_clusters = min(5, max(2, len(valid_modules) // 3))
    
    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    
    # Create module-to-cluster mapping
    module_clusters = {}
    for i, module in enumerate(valid_modules):
        module_clusters[module] = cluster_labels[i]
    
    # Assign cluster 0 to modules without pathways
    for module in modules:
        if module not in module_clusters:
            module_clusters[module] = 0
    
    print(f"Created {n_clusters} module clusters based on functional similarity")
    
    # Print cluster summary
    cluster_counts = Counter(module_clusters.values())
    for cluster_id, count in sorted(cluster_counts.items()):
        print(f"  Cluster {cluster_id}: {count} modules")
    
    return module_clusters, similarity_matrix

def create_interactive_network_visualization(module_data, module_clusters, output_dir, enrichment_data=None, go_similarity_matrix=None):
    """
    Create interactive network visualization using plotly with MegaGO clustering layout
    
    Parameters:
    module_data (list): List of module data dictionaries
    module_clusters (dict): Dictionary mapping module IDs to cluster labels
    output_dir (str): Output directory for saving visualization
    enrichment_data (dict): Enrichment data for enhanced hover information
    go_similarity_matrix (numpy.ndarray): GO similarity matrix for spatial clustering
    """
    print("Creating interactive network visualization with MegaGO clustering...")
    
    # Create network graph
    G = nx.Graph()
    
    # Add nodes (modules and regulators)
    nodes = []
    edges = []
    
    # Track regulators to avoid duplicates (dynamically populated)
    regulator_modules = {}
    
    # Process each module
    for module_info in module_data:
        module_id = str(module_info['Module'])
        cluster_id = module_clusters.get(module_id, 0)
        
        # Calculate node attributes
        n_genes = len(module_info.get('Module_genes', '').split('|')) if module_info.get('Module_genes', '') != 'NA' else 0
        
        # Count regulators dynamically for all regulator types
        regulator_counts = {}
        regulator_columns = [col for col in module_info.keys() if col.endswith('_regulators')]
        for reg_col in regulator_columns:
            reg_type = reg_col.replace('_regulators', '')
            reg_count = len(module_info.get(reg_col, '').split('|')) if module_info.get(reg_col, '') != 'NA' else 0
            regulator_counts[reg_type] = reg_count
        
        # Keep legacy variables for backwards compatibility (use first two regulator types if available)
        all_reg_types = sorted(regulator_counts.keys())
        n_tf_regs = regulator_counts.get(all_reg_types[0], 0) if len(all_reg_types) > 0 else 0
        n_met_regs = regulator_counts.get(all_reg_types[1], 0) if len(all_reg_types) > 1 else 0
        
        # Count pathways
        pathway_counts = 0
        for pathway_type in ['Top_3_pathways_bio_process', 'Top_3_pathways_molecular_function', 
                           'Top_3_pathways_cellular_component', 'Top_3_pathways_KEGG', 'Top_3_pathways_Reactome']:
            if module_info.get(pathway_type, 'NA') != 'NA':
                pathway_counts += len(module_info[pathway_type].split('|'))
        
        # Create module node
        module_node = {
            'id': f"Module_{module_id}",
            'label': f"Module {module_id}",
            'type': 'module',
            'cluster': cluster_id,
            'n_genes': n_genes,
            'n_tf_regs': n_tf_regs,
            'n_met_regs': n_met_regs,
            'n_pathways': pathway_counts,
            'expression_rank': module_info.get('Expression_rank', 'NA'),
            'hover_info': create_module_hover_info(module_info, enrichment_data)
        }
        nodes.append(module_node)
        G.add_node(module_node['id'], **module_node)
        
        # Collect regulators dynamically for all regulator types
        regulator_columns = [col for col in module_info.keys() if col.endswith('_regulators')]
        for reg_col in regulator_columns:
            if module_info.get(reg_col, 'NA') != 'NA':
                # Extract regulator type name
                reg_type = reg_col.replace('_regulators', '')
                regs = module_info[reg_col].split('|')
                
                # Initialize regulator type dict if needed
                if reg_type not in regulator_modules:
                    regulator_modules[reg_type] = {}
                
                for reg in regs:
                    reg = reg.strip()
                    if reg:
                        if reg not in regulator_modules[reg_type]:
                            regulator_modules[reg_type][reg] = []
                        regulator_modules[reg_type][reg].append(module_id)
    
    # Create regulator nodes
    for reg_type, regulators in regulator_modules.items():
        for regulator, target_modules in regulators.items():
            modules_str = ', '.join(sorted(target_modules))
            regulator_node = {
                'id': f"{reg_type}_{regulator}",
                'label': regulator,
                'type': reg_type,
                'hover_info': f"<b>{regulator}</b><br>Type: {reg_type.upper()}<br>Regulates: Modules {modules_str}<br>Module count: {len(target_modules)}"
            }
            nodes.append(regulator_node)
            G.add_node(regulator_node['id'], **regulator_node)
            
            # Create edges
            for module_id in target_modules:
                edge = {
                    'source': f"{reg_type}_{regulator}",
                    'target': f"Module_{module_id}",
                    'type': f"{reg_type}_regulation"
                }
                edges.append(edge)
                G.add_edge(edge['source'], edge['target'])
    
    # Create MegaGO clustering layout
    module_list = list(module_clusters.keys())
    pos = create_megago_clustering_layout(nodes, edges, go_similarity_matrix, module_list)
    
    # Prepare data for plotly
    edge_x = []
    edge_y = []
    
    for edge in edges:
        if edge['source'] in pos and edge['target'] in pos:
            x0, y0 = pos[edge['source']]
            x1, y1 = pos[edge['target']]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
    
    # Create edge trace
    edge_trace = go.Scatter(x=edge_x, y=edge_y,
                           line=dict(width=1, color='rgba(150,150,150,0.5)'),
                           hoverinfo='none',
                           mode='lines',
                           name='Connections',
                           showlegend=False)
    
    # Create single node trace for all modules (no cluster-based coloring)
    all_module_nodes = [n for n in nodes if n['type'] == 'module' and n['id'] in pos]
    
    if all_module_nodes:
        node_x = [pos[n['id']][0] for n in all_module_nodes]
        node_y = [pos[n['id']][1] for n in all_module_nodes]
        
        # Get hover texts and sizes
        hover_texts = []
        node_sizes = []
        
        for node in all_module_nodes:
            hover_texts.append(node.get('hover_info', node['id']))
            # Size based on number of genes with increased base size for better visibility
            n_genes = node.get('n_genes', 10)
            size = max(15, min(80, n_genes * 3))  # Increased base size and scaling factor
            node_sizes.append(size)
        
        # Single trace for all modules with uniform color
        node_trace = go.Scatter(x=node_x, y=node_y,
                               mode='markers+text',
                               hoverinfo='text',
                               text=[n['label'] for n in all_module_nodes],
                               textposition='middle center',
                               textfont=dict(size=10, color='black'),  # Increased text size
                               hovertext=hover_texts,
                               name='Modules',
                               marker=dict(size=node_sizes,
                                         color='lightblue',  # Single uniform color
                                         line=dict(width=2, color='darkblue'),
                                         opacity=0.8),
                               showlegend=True)
        node_traces = [node_trace]
    else:
        node_traces = []
    
    # Create regulator traces dynamically for all regulator types
    regulator_traces = []
    
    # Define default colors and symbols for known regulator types
    default_styles = {
        'TFs': {'color': 'green', 'symbol': 'triangle-up'},
        'Metabolites': {'color': 'red', 'symbol': 'circle'},
        'Lipids': {'color': 'gold', 'symbol': 'square'},
        'Proteins': {'color': 'purple', 'symbol': 'diamond'}
    }
    
    # Define colors and symbols to use for unknown regulator types dynamically
    fallback_colors = ['orange', 'cyan', 'magenta', 'brown', 'pink', 
                       'lime', 'navy', 'teal', 'coral']
    fallback_symbols = ['circle', 'square', 'diamond', 'cross', 'x', 'triangle-up', 'triangle-down', 
                        'star', 'hexagon']
    
    # Build regulator_styles dynamically based on found regulator types
    regulator_styles = {}
    for idx, reg_type in enumerate(sorted(regulator_modules.keys())):
        if reg_type in default_styles:
            # Use predefined style for known types
            regulator_styles[reg_type] = {
                'color': default_styles[reg_type]['color'],
                'symbol': default_styles[reg_type]['symbol'],
                'name': f'{reg_type} Regulators'
            }
        else:
            # Dynamically assign color and symbol for unknown types
            color_idx = idx % len(fallback_colors)
            symbol_idx = idx % len(fallback_symbols)
            regulator_styles[reg_type] = {
                'color': fallback_colors[color_idx],
                'symbol': fallback_symbols[symbol_idx],
                'name': f'{reg_type} Regulators'
            }
    
    # Create traces for all regulator types found in the data
    for reg_type in sorted(regulator_modules.keys()):
        reg_nodes = [n for n in nodes if n['type'] == reg_type and n['id'] in pos]
        
        if reg_nodes:
            reg_x = [pos[n['id']][0] for n in reg_nodes]
            reg_y = [pos[n['id']][1] for n in reg_nodes]
            reg_hover = [n['hover_info'] for n in reg_nodes]
            reg_labels = [n['label'] for n in reg_nodes]
            
            # Get style for this regulator type (guaranteed to exist now)
            style = regulator_styles[reg_type]
            
            reg_trace = go.Scatter(x=reg_x, y=reg_y,
                                  mode='markers+text',
                                  hoverinfo='text',
                                  text=reg_labels,
                                  textposition='middle center',
                                  textfont=dict(size=7, color='black'),
                                  hovertext=reg_hover,
                                  name=style['name'],
                                  marker=dict(size=18,  # Increased regulator size
                                            color=style['color'],
                                            symbol=style['symbol'],
                                            line=dict(width=1.5, color='darkgray'),
                                            opacity=0.8),
                                  showlegend=True)
            regulator_traces.append(reg_trace)
    
    # Create the figure with larger dimensions and better scaling
    fig = go.Figure(data=[edge_trace] + node_traces + regulator_traces,
                   layout=go.Layout(
                        title=dict(
                            text='Interactive Module-Regulator Network',
                            font=dict(size=18),
                            x=0.5
                        ),
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=60,l=40,r=40,t=80),
                        width=1200,  # Explicit width for larger network display
                        height=800,  # Explicit height for larger network display
                        annotations=[
                            dict(
                                text="Modules positioned by GO biological process similarity.<br>Node size = number of genes, all modules shown in uniform color.<br>Hover over nodes for details, zoom and pan to explore the network.",
                                showarrow=False,
                                xref="paper", yref="paper",
                                x=0.005, y=-0.002,
                                xanchor='left', yanchor='bottom',
                                font=dict(size=12)
                            )
                        ],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, scaleanchor="y", scaleratio=1),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        plot_bgcolor='white'
                        ))
    
    # Save interactive plot
    output_file = os.path.join(output_dir, 'interactive_module_network.html')
    pyo.plot(fig, filename=output_file, auto_open=False)
    print(f"Interactive network visualization saved to: {output_file}")
    
    return fig

def create_cluster_heatmap(module_data, module_clusters, output_dir):
    """
    Create heatmap showing module characteristics by cluster
    
    Parameters:
    module_data (list): List of module data dictionaries  
    module_clusters (dict): Dictionary mapping module IDs to cluster labels
    output_dir (str): Output directory for saving visualization
    """
    print("Creating cluster characteristics heatmap...")
    
    # Prepare data matrix
    clusters = sorted(set(module_clusters.values()))
    characteristics = ['n_genes', 'n_tf_regs', 'n_met_regs', 'n_pathways_bp', 'n_pathways_mf', 
                      'n_pathways_cc', 'n_pathways_kegg', 'n_pathways_reactome']
    
    cluster_stats = []
    
    for cluster_id in clusters:
        cluster_modules = [mod for mod, clust in module_clusters.items() if clust == cluster_id]
        
        # Calculate statistics for this cluster
        stats = {'cluster': cluster_id, 'n_modules': len(cluster_modules)}
        
        # Initialize counters
        genes_list = []
        # Dynamic regulator lists - initialize based on available regulator types
        regulator_lists = {}
        
        # Get all regulator types from the first module
        if cluster_modules:
            first_mod = next((m for m in module_data if str(m['Module']) == cluster_modules[0]), {})
            regulator_columns = [col for col in first_mod.keys() if col.endswith('_regulators')]
            for reg_col in regulator_columns:
                reg_type = reg_col.replace('_regulators', '')
                regulator_lists[reg_type] = []
        
        pathway_counts = {pt: [] for pt in ['bp', 'mf', 'cc', 'kegg', 'reactome']}
        
        for mod_id in cluster_modules:
            mod_data = next((m for m in module_data if str(m['Module']) == mod_id), {})
            
            # Count genes
            genes = mod_data.get('Module_genes', 'NA')
            n_genes = len(genes.split('|')) if genes != 'NA' else 0
            genes_list.append(n_genes)
            
            # Count regulators dynamically for all types
            for reg_type in regulator_lists.keys():
                reg_col = f'{reg_type}_regulators'
                regs = mod_data.get(reg_col, 'NA')
                n_regs = len(regs.split('|')) if regs != 'NA' else 0
                regulator_lists[reg_type].append(n_regs)
            
            # Count pathways by type
            pathway_mapping = {
                'bp': 'Top_3_pathways_bio_process',
                'mf': 'Top_3_pathways_molecular_function',
                'cc': 'Top_3_pathways_cellular_component',
                'kegg': 'Top_3_pathways_KEGG',
                'reactome': 'Top_3_pathways_Reactome'
            }
            
            for pt_short, pt_long in pathway_mapping.items():
                pathways = mod_data.get(pt_long, 'NA')
                n_pathways = len(pathways.split('|')) if pathways != 'NA' else 0
                pathway_counts[pt_short].append(n_pathways)
        
        # Calculate means
        stats['mean_genes'] = np.mean(genes_list) if genes_list else 0
        
        # Calculate mean regulators dynamically for all types
        for reg_type, reg_list in regulator_lists.items():
            stats[f'mean_{reg_type}_regs'] = np.mean(reg_list) if reg_list else 0
        
        for pt in ['bp', 'mf', 'cc', 'kegg', 'reactome']:
            stats[f'mean_pathways_{pt}'] = np.mean(pathway_counts[pt]) if pathway_counts[pt] else 0
        
        cluster_stats.append(stats)
    
    # Create DataFrame
    stats_df = pd.DataFrame(cluster_stats)
    
    # Prepare data for heatmap
    heatmap_data = []
    heatmap_labels = []
    
    for _, row in stats_df.iterrows():
        heatmap_data.append([
            row['mean_genes'],
            row['mean_tf_regs'], 
            row['mean_met_regs'],
            row['mean_pathways_bp'],
            row['mean_pathways_mf'],
            row['mean_pathways_cc'],
            row['mean_pathways_kegg'],
            row['mean_pathways_reactome']
        ])
        heatmap_labels.append(f"Cluster {row['cluster']}<br>({row['n_modules']} modules)")
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=['Mean Genes', 'Mean TF Regs', 'Mean Met Regs', 'Mean BP Pathways',
           'Mean MF Pathways', 'Mean CC Pathways', 'Mean KEGG Pathways', 'Mean Reactome Pathways'],
        y=heatmap_labels,
        colorscale='Viridis',
        showscale=True
    ))
    
    fig.update_layout(
        title='Module Cluster Characteristics',
        xaxis_title='Characteristics',
        yaxis_title='Functional Clusters',
        height=max(400, len(clusters) * 50),
        margin=dict(l=150)
    )
    
    # COMMENTED OUT: Save heatmap - not needed in output
    # output_file = os.path.join(output_dir, 'cluster_characteristics_heatmap.html')
    # pyo.plot(fig, filename=output_file, auto_open=False)
    # print(f"Cluster characteristics heatmap saved to: {output_file}")
    
    # COMMENTED OUT: Save cluster statistics - not needed in output
    # stats_file = os.path.join(output_dir, 'cluster_statistics.csv')
    # stats_df.to_csv(stats_file, index=False)
    # print(f"Cluster statistics saved to: {stats_file}")
    
    print("Cluster characteristics heatmap creation skipped (not needed in output)")
    
    return fig, stats_df

def auto_prioritize_modules_expression(group_column='diagnosis', modules_dict=None, 
                                      expression_file=None, deseq_groups_file=None):
    """
    Prioritize modules based on differential expression using auto-detected groups.
    Auto-detects whether to use Mann-Whitney U (2 groups) or Kruskal-Wallis (3+ groups).
    """
    # Basic argument checks
    if modules_dict is None or expression_file is None or deseq_groups_file is None:
        print("Warning: Missing required files for expression prioritization")
        return {}, pd.DataFrame()

    try:
        # Load expression data (flexible format)
        print(f"Loading expression data from: {expression_file}")
        expr_df_raw = pd.read_csv(expression_file, sep='\t', header=0)

        # Determine whether genes are in a 'symbol' column or index
        if 'symbol' in expr_df_raw.columns.str.lower():
            # normalize column name
            sym_col = next(c for c in expr_df_raw.columns if c.lower() == 'symbol')
            expr_df = expr_df_raw.set_index(sym_col)
        else:
            # assume first column is gene names if not labelled
            # if the file already has genes as index, try to set index
            if expr_df_raw.columns[0].lower() in ('gene','symbol'):
                expr_df = expr_df_raw.set_index(expr_df_raw.columns[0])
            else:
                # if the file already uses index, try reading again with index_col=0
                expr_df = pd.read_csv(expression_file, sep='\t', index_col=0)

        # Load metadata
        print(f"Loading metadata from: {deseq_groups_file}")
        metadata = pd.read_csv(deseq_groups_file, sep='\t', header=0, index_col=0)

        if group_column not in metadata.columns:
            print(f"Error: Column '{group_column}' not found in metadata")
            return {}, pd.DataFrame()

        # Align samples between expression and metadata
        sample_ids = [c for c in expr_df.columns if c in metadata.index]
        if len(sample_ids) < 2:
            print(f"Error: Not enough common samples between expression and metadata: {len(sample_ids)}")
            return {}, pd.DataFrame()

        expr_df = expr_df[sample_ids]
        metadata = metadata.loc[sample_ids]

        unique_groups = metadata[group_column].dropna().unique()
        n_groups = len(unique_groups)
        print(f"Auto-detected {n_groups} groups in '{group_column}': {list(unique_groups)}")

        # Build module mean expression per sample
        module_scores = {}
        for module_id, genes in (modules_dict or {}).items():
            if isinstance(genes, str):
                genes = [g for g in genes.split('|') if g]
            elif not isinstance(genes, (list, tuple)):
                continue

            present_genes = [g for g in genes if g in expr_df.index]
            if len(present_genes) == 0:
                continue

            module_mean = expr_df.loc[present_genes].mean(axis=0)
            module_scores[module_id] = module_mean

        if not module_scores:
            print("No modules with expression data found")
            return {}, pd.DataFrame()

        # Perform tests
        results = []
        for module_id, series in module_scores.items():
            # series indexed by sample ids
            groups_data = []
            group_means = {}
            for grp in unique_groups:
                samples = metadata.index[metadata[group_column] == grp].tolist()
                vals = series.reindex(samples).dropna().values
                group_means[str(grp)] = float(np.mean(vals)) if len(vals) > 0 else np.nan
                if len(vals) > 0:
                    groups_data.append(vals)

            if len(groups_data) < 2:
                # not enough data
                continue

            try:
                if len(unique_groups) == 2:
                    # Mann-Whitney U
                    g1, g2 = groups_data[0], groups_data[1]
                    stat, p = mannwhitneyu(g1, g2, alternative='two-sided')
                    test_name = 'Mann-Whitney U'
                    # effect size: rank-biserial approximation
                    n1, n2 = len(g1), len(g2)
                    effect = 1 - (2 * stat) / (n1 * n2) if (n1 * n2) > 0 else 0.0
                else:
                    # Kruskal-Wallis
                    stat, p = kruskal(*groups_data)
                    test_name = 'Kruskal-Wallis'
                    n_total = sum(len(g) for g in groups_data)
                    effect = (stat - len(groups_data) + 1) / (n_total - len(groups_data)) if (n_total - len(groups_data)) > 0 else 0.0
                    effect = max(0.0, effect)

                results.append({
                    'Module': module_id,
                    'Test': test_name,
                    'Statistic': float(stat),
                    'P_value': float(p),
                    'Effect_size': float(effect),
                    'N_samples': int((series.notna()).sum()),
                    **{f'Mean_{g}': group_means.get(str(g), np.nan) for g in unique_groups}
                })
            except Exception as e:
                print(f"Warning: test failed for module {module_id}: {e}")
                continue

        if not results:
            print("No modules could be analyzed after grouping")
            return {}, pd.DataFrame()

        results_df = pd.DataFrame(results)

        # Multiple testing correction using statsmodels if available
        try:
            from statsmodels.stats.multitest import multipletests
            pvals = results_df['P_value'].values
            adj = multipletests(pvals, method='fdr_bh')[1]
            results_df['P_adjusted'] = adj
        except Exception:
            # fallback BH
            from scipy.stats import rankdata
            p_values = results_df['P_value'].values
            ranked = rankdata(p_values)
            n = len(p_values)
            results_df['P_adjusted'] = np.minimum(p_values * n / ranked, 1.0)

        # Rank modules by adjusted p-value and effect size
        results_df['Combined_score'] = -np.log10(results_df['P_adjusted'] + 1e-12) * np.abs(results_df['Effect_size'])
        results_df = results_df.sort_values(['P_adjusted', 'Combined_score'], ascending=[True, False]).reset_index(drop=True)
        results_df['Expression_rank'] = range(1, len(results_df) + 1)

        priority_dict = dict(zip(results_df['Module'].astype(str), results_df['Expression_rank']))

        print(f"Expression analysis completed: {len(results_df)} modules tested; {sum(results_df['P_adjusted'] < 0.05)} significant (FDR<0.05)")

        return priority_dict, results_df

    except Exception as e:
        print(f"Error in expression prioritization: {e}")
        traceback.print_exc()
        return {}, pd.DataFrame()
        return {}, pd.DataFrame()

def load_coherence_filtered_modules(coherence_file, threshold=0.6):
    """
    Load module coherence scores and filter modules based on threshold
    """
    # Implementation from original script
    # [Previous implementation would go here - keeping it the same]
    try:
        coherence_df = pd.read_csv(coherence_file, sep='\t')
        print(f"Loaded coherence scores for {len(coherence_df)} modules")
        
        # Filter modules that meet coherence threshold
        filtered_modules = coherence_df[coherence_df['Coherence_Score'] >= threshold]
        removed_modules = coherence_df[coherence_df['Coherence_Score'] < threshold]
        
        print(f"Coherence filtering (threshold = {threshold}):")
        print(f"- Modules passing filter: {len(filtered_modules)}")
        print(f"- Modules removed: {len(removed_modules)}")
        
        if len(removed_modules) > 0:
            print(f"- Removed modules: {sorted(removed_modules['Module'].tolist())}")
            
        # Return list of module IDs that pass the filter
        return filtered_modules['Module'].tolist(), coherence_df
        
    except FileNotFoundError:
        print(f"Warning: Coherence scores file not found: {coherence_file}")
        print("Proceeding without coherence filtering...")
        return None, None
    except Exception as e:
        print(f"Error loading coherence scores: {e}")
        return None, None

def load_filtered_modules_from_networks(input_dir):
    """
    Load the list of coherence-filtered modules from Networks/specific_modules.txt
    This file is created by lemontree_to_network.py and contains only modules that passed filtering
    
    Parameters:
    input_dir (str): Input directory containing Networks subdirectory
    
    Returns:
    set: Set of module IDs that passed coherence filtering, or None if file not found
    """
    # Try Networks/specific_modules.txt first (standard location)
    networks_dir = os.path.join(input_dir, 'Networks')
    specific_modules_file = os.path.join(networks_dir, 'specific_modules.txt')
    
    # If not found, try specific_modules.txt directly in input_dir (Nextflow symlink location)
    if not os.path.exists(specific_modules_file):
        specific_modules_file = os.path.join(input_dir, 'specific_modules.txt')
        
    if not os.path.exists(specific_modules_file):
        print(f"Warning: specific_modules.txt not found at: {networks_dir}/specific_modules.txt or {input_dir}/specific_modules.txt")
        return None
    
    try:
        with open(specific_modules_file, 'r') as f:
            # Read module IDs from file (one per line, remove whitespace)
            filtered_modules = set()
            for line in f:
                module_id = line.strip()
                if module_id:  # Skip empty lines
                    filtered_modules.add(module_id)
        
        print(f"Loaded {len(filtered_modules)} coherence-filtered modules from Networks/specific_modules.txt")
        print(f"   Filtered modules: {sorted(list(filtered_modules))}")
        return filtered_modules
        
    except Exception as e:
        print(f"Error reading Networks/specific_modules.txt: {e}")
        return None

def create_megago_clustering_layout(nodes, edges, go_similarity_matrix=None, module_list=None):
    """
    Create network layout using megaGO-inspired GO semantic clustering
    Uses hierarchical clustering on GO semantic similarity for BIOLOGICAL PROCESS (BP) terms only
    """
    print("Creating megaGO-inspired GO clustering layout (BP terms only)...")

    # Get module nodes
    module_nodes = [n for n in nodes if n['type'] == 'module']

    if go_similarity_matrix is None or len(module_nodes) < 2:
        print("No GO similarity data or insufficient modules, using spring layout")
        G_temp = nx.Graph()
        for node in nodes:
            G_temp.add_node(node['id'])
        for edge in edges:
            G_temp.add_edge(edge['source'], edge['target'])
        return nx.spring_layout(G_temp, k=8, iterations=100, seed=42, scale=10)

    # Step 1: Apply hierarchical clustering to GO similarity matrix
    distance_matrix = 1 - go_similarity_matrix

    try:
        # Use hierarchical clustering (similar to megaGO's approach)
        n_clusters = min(5, len(module_nodes) // 3 + 1)  # More clusters for finer GO grouping

        # Perform hierarchical clustering
        linkage_matrix = linkage(distance_matrix, method='ward')
        cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
        cluster_labels = cluster_labels - 1  # Convert to 0-based indexing

    except Exception as e:
        print(f"Hierarchical clustering failed: {e}, using single cluster")
        cluster_labels = [0] * len(module_nodes)
        n_clusters = 1

    # Step 2: Create spatial layout based on GO clusters
    print("Step 2: Creating BP GO-informed spatial layout...")
    pos = {}

    # Position modules based on their GO clusters using concentric circles with larger scale
    if n_clusters > 1:
        # Create concentric arrangement for better visualization with larger radius
        if n_clusters <= 4:
            # Single ring for small number of clusters - increased radius
            cluster_centers = [(8 * np.cos(2*np.pi*i/n_clusters), 8 * np.sin(2*np.pi*i/n_clusters))
                              for i in range(n_clusters)]
        else:
            # Multiple concentric rings for more clusters - increased radii
            inner_clusters = n_clusters // 2
            outer_clusters = n_clusters - inner_clusters

            cluster_centers = []
            # Inner ring - increased radius
            for i in range(inner_clusters):
                angle = 2*np.pi*i/inner_clusters
                cluster_centers.append((6 * np.cos(angle), 6 * np.sin(angle)))

            # Outer ring - increased radius
            for i in range(outer_clusters):
                angle = 2*np.pi*i/outer_clusters + np.pi/outer_clusters  # Offset for better spacing
                cluster_centers.append((12 * np.cos(angle), 12 * np.sin(angle)))
    else:
        cluster_centers = [(0, 0)]

    np.random.seed(123)  # Different seed for GO clustering

    # Map module IDs to positions with increased jitter for better spread
    module_id_to_pos = {}
    for i, (module_id, cluster) in enumerate(zip(module_list, cluster_labels)):
        center = cluster_centers[cluster % len(cluster_centers)]
        # Add GO-specific jitter pattern with larger spread
        jitter_x = np.random.normal(0, 1.2)  # Increased jitter for better visibility
        jitter_y = np.random.normal(0, 1.2)  # Increased jitter for better visibility
        module_id_to_pos[f"Module_{module_id}"] = (center[0] + jitter_x, center[1] + jitter_y)

    # Position all module nodes
    for node in module_nodes:
        if node['id'] in module_id_to_pos:
            pos[node['id']] = module_id_to_pos[node['id']]
        else:
            # Fallback position with larger spread
            pos[node['id']] = (np.random.normal(0, 5), np.random.normal(0, 5))

    # Step 3: Position regulators with GO-aware placement
    for node in nodes:
        if node['type'] != 'module' and node['id'] not in pos:
            # Find connected modules
            connected_modules = []
            for edge in edges:
                if edge['source'] == node['id'] and edge['target'].startswith('Module_'):
                    connected_modules.append(edge['target'])
                elif edge['target'] == node['id'] and edge['source'].startswith('Module_'):
                    connected_modules.append(edge['source'])

            if connected_modules:
                # Position based on connected modules' GO clusters
                connected_positions = [pos[mod] for mod in connected_modules if mod in pos]

                if connected_positions:
                    avg_x = np.mean([p[0] for p in connected_positions])
                    avg_y = np.mean([p[1] for p in connected_positions])

                    # Add offset based on regulator type with GO-specific spacing
                    offset_distance = 2.5  # Increased offset for better visibility
                    if node['type'] == 'TF':
                        angle_offset = 0
                    elif node['type'] == 'metabolite':
                        angle_offset = 2*np.pi/3
                    else:  # hPTM
                        angle_offset = 4*np.pi/3

                    # Add more variation for GO-based layout
                    angle = angle_offset + np.random.normal(0, 0.6)  # Increased variation
                    offset_x = offset_distance * np.cos(angle)
                    offset_y = offset_distance * np.sin(angle)

                    pos[node['id']] = (avg_x + offset_x, avg_y + offset_y)
                else:
                    pos[node['id']] = (np.random.normal(0, 4), np.random.normal(0, 4))
            else:
                # Random position if no connections with larger spread
                pos[node['id']] = (np.random.normal(0, 6), np.random.normal(0, 6))
                pos[node['id']] = (np.random.normal(0, 1.5), np.random.normal(0, 1.5))

    print(f"Successfully created megaGO-inspired BP GO clustering layout")
    return pos

def create_module_hover_info(row, enrichment_data=None):
    """Create detailed hover information for module nodes with enhanced pathway details"""
    module_id = row['Module']

    # Start with basic info
    hover_text = f"<b>Module {module_id}</b><br>"

    # Add expression prioritization information if available
    if 'Expression_p_value' in row and row['Expression_p_value'] != 'NA':
        hover_text += f"<b>Expression Analysis:</b><br>"
        hover_text += f"  • p-value: {row['Expression_p_value']:.2e}<br>"
        hover_text += f"  • Rank: {row['Expression_rank']}<br>"
        hover_text += f"  • Significant: {row['Expression_significant']}<br>"
        hover_text += "<br>"

    # Add gene count
    if row['Module_genes'] != 'NA' and pd.notna(row['Module_genes']):
        gene_count = len(row['Module_genes'].split('|'))
        hover_text += f"<b>Genes:</b> {gene_count} genes<br>"

    # Add detailed regulator information dynamically for all regulator types
    # Find all columns ending with '_regulators'
    # Handle both dict and pandas Series types
    if isinstance(row, dict):
        regulator_columns = [col for col in row.keys() if col.endswith('_regulators')]
    else:
        regulator_columns = [col for col in row.index if col.endswith('_regulators')]
    
    for reg_col in sorted(regulator_columns):  # Sort for consistent order
        if row[reg_col] != 'NA' and pd.notna(row[reg_col]):
            # Extract regulator type name (e.g., 'TFs_regulators' -> 'TFs')
            reg_type = reg_col.replace('_regulators', '')
            regs = row[reg_col].split('|')
            hover_text += f"<br><b>{reg_type} Regulators ({len(regs)}):</b><br>"
            # Show up to 5 regulators
            shown_regs = regs[:5]
            for reg in shown_regs:
                hover_text += f"  • {reg}<br>"
            if len(regs) > 5:
                hover_text += f"  ... and {len(regs) - 5} more<br>"

    # Add detailed pathway information grouped by database with up/down regulation
    if enrichment_data:
        hover_text += f"<br><b>Top Enriched Pathways:</b><br>"

        # Define databases to show
        databases = ['bp', 'mf', 'cc', 'kegg', 'reactome']

        for database in databases:
            # Get up-regulated pathways for this module and database
            if database in enrichment_data and not enrichment_data[database].empty:
                module_up_db = enrichment_data[database][
                    (enrichment_data[database]['Module'] == str(module_id)) &
                    (enrichment_data[database]['__direction__'] == 'Up')
                ]

                # Get down-regulated pathways for this module and database
                module_down_db = enrichment_data[database][
                    (enrichment_data[database]['Module'] == str(module_id)) &
                    (enrichment_data[database]['__direction__'] == 'Down')
                ]

                # Show top 2 pathways per database per direction
                if not module_up_db.empty or not module_down_db.empty:
                    hover_text += f"<br><b>{database.upper()}:</b><br>"

                    # Up-regulated pathways
                    if not module_up_db.empty:
                        top_up = module_up_db.nsmallest(2, 'p.adjust')
                        for _, pathway in top_up.iterrows():
                            term = pathway['Term']
                            # Truncate long pathway names
                            if len(term) > 45:
                                term = term[:42] + "..."
                            hover_text += f"  ↑ {term}<br>"

                    # Down-regulated pathways
                    if not module_down_db.empty:
                        top_down = module_down_db.nsmallest(2, 'p.adjust')
                        for _, pathway in top_down.iterrows():
                            term = pathway['Term']
                            # Truncate long pathway names
                            if len(term) > 45:
                                term = term[:42] + "..."
                            hover_text += f"  ↓ {term}<br>"
    else:
        # Fallback to the old format if enrichment data is not available
        hover_text += f"<br><b>Top Enriched Pathways:</b><br>"

        pathways = [
            ('BP', row.get('Top_3_pathways_bio_process', 'NA')),
            ('MF', row.get('Top_3_pathways_molecular_function', 'NA')),
            ('CC', row.get('Top_3_pathways_cellular_component', 'NA')),
            ('KEGG', row.get('Top_3_pathways_KEGG', 'NA')),
            ('Reactome', row.get('Top_3_pathways_Reactome', 'NA'))
        ]

        for pathway_type, pathway_data in pathways:
            if pathway_data != 'NA' and pd.notna(pathway_data) and pathway_data != '':
                pathways_list = str(pathway_data).split('|')
                if len(pathways_list) > 0:
                    hover_text += f"<br><b>{pathway_type}:</b><br>"
                    # Show up to 2 pathways per database
                    shown_pathways = pathways_list[:2]
                    for pathway in shown_pathways:
                        # Truncate long pathway names
                        if len(pathway) > 45:
                            pathway = pathway[:42] + "..."
                        hover_text += f"  • {pathway}<br>"
                    if len(pathways_list) > 2:
                        hover_text += f"  ... and {len(pathways_list) - 2} more<br>"

    return hover_text

def create_module_expression_heatmap(module_genes, modules_to_process, expression_file, metadata_file, 
                                     module_pvalues, output_dir, group_column='diagnosis'):
    """
    Create a heatmap showing average module expression across sample subtypes.
    Only creates the heatmap if differential expression analysis was performed.
    
    Parameters:
    -----------
    module_genes : dict
        Dictionary mapping module IDs to lists of gene symbols
    modules_to_process : set or list
        List of module IDs to include in the heatmap
    expression_file : str
        Path to expression data file (LemonPreprocessed_expression.txt)
    metadata_file : str
        Path to metadata file (DESeq_groups.txt)
    module_pvalues : dict
        Dictionary mapping module IDs to p-values from differential expression
    output_dir : str
        Directory to save the heatmap
    group_column : str
        Column name for sample groups in metadata file
        
    Returns:
    --------
    bool : True if heatmap was created successfully, False otherwise
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Matplotlib/seaborn not available. Skipping Module_Expression_Heatmap.png generation.")
        return False
    
    # Skip if no expression file or metadata file
    if not expression_file or not os.path.exists(expression_file):
        print("Expression file not found. Skipping Module_Expression_Heatmap.png generation.")
        return False
        
    if not metadata_file or not os.path.exists(metadata_file):
        print("Metadata file not found. Skipping Module_Expression_Heatmap.png generation.")
        return False
    
    # Skip if no p-values (meaning differential expression was not performed)
    if not module_pvalues or len(module_pvalues) == 0:
        print("Differential expression not performed. Skipping Module_Expression_Heatmap.png generation.")
        return False
    
    try:
        print("\nCreating Module Expression Heatmap...")
        
        # Load expression data
        print("  Loading expression data...")
        expression_data = pd.read_csv(expression_file, sep='\t')
        print(f"     Expression data loaded: {expression_data.shape}")
        
        # Load sample annotations
        print("  Loading sample annotations...")
        annotations = pd.read_csv(metadata_file, sep='\t')
        
        # Auto-detect the sample group column
        # Try common column names
        possible_group_columns = [group_column, 'multiomic', 'diagnosis', 'condition', 'group', 'subtype']
        group_col = None
        for col in possible_group_columns:
            if col in annotations.columns:
                group_col = col
                break
        
        if group_col is None:
            print(f"  Could not find group column. Tried: {possible_group_columns}")
            return False
        
        print(f"     Using group column: '{group_col}'")
        
        # Create sample mapping
        if 'Sample_ID' not in annotations.columns:
            # Use index as Sample_ID if not present
            annotations['Sample_ID'] = annotations.index.astype(str)
        
        sample_mapping = dict(zip(annotations['Sample_ID'].astype(str), annotations[group_col]))
        
        # Get unique groups
        unique_groups = sorted(annotations[group_col].dropna().unique())
        print(f"     Detected {len(unique_groups)} sample groups: {unique_groups}")
        
        # Prepare data for heatmap
        heatmap_data = []
        heatmap_module_labels = []
        
        print(f"  Processing {len(modules_to_process)} modules...")
        
        for module in modules_to_process:
            module_str = str(module)
            if module_str in module_genes:
                genes = module_genes[module_str]
                
                # Filter expression data for genes in this module
                module_expression = expression_data[expression_data['symbol'].isin(genes)]
                
                if len(module_expression) == 0:
                    continue
                
                # Remove non-numeric columns (keep only sample columns)
                numeric_cols = module_expression.select_dtypes(include=[np.number]).columns
                module_expression = module_expression[numeric_cols]
                
                # Calculate average expression per sample
                sample_means = module_expression.mean(axis=0)
                
                # Calculate average expression per group
                group_averages = {}
                for group in unique_groups:
                    group_samples = [sample for sample, group_name in sample_mapping.items() if group_name == group]
                    group_samples_in_data = [s for s in group_samples if s in sample_means.index]
                    
                    if group_samples_in_data:
                        group_avg = sample_means[group_samples_in_data].mean()
                        group_averages[group] = group_avg
                    else:
                        group_averages[group] = np.nan
                
                # Add to heatmap data
                row_data = [group_averages[group] for group in unique_groups]
                heatmap_data.append(row_data)
                
                # Create module label (add * if significant)
                p_value = module_pvalues.get(module_str, 1.0)
                if isinstance(p_value, (int, float)) and p_value < 0.05:
                    module_label = f"Module {module} *"
                else:
                    module_label = f"Module {module}"
                heatmap_module_labels.append(module_label)
        
        if not heatmap_data:
            print("  No data available for heatmap")
            return False
        
        # Create DataFrame for heatmap
        heatmap_df = pd.DataFrame(heatmap_data,
                                  index=heatmap_module_labels,
                                  columns=unique_groups)
        
        # Remove rows with all NaN values
        heatmap_df = heatmap_df.dropna(how='all')
        
        if heatmap_df.empty:
            print("  No valid data for heatmap after removing NaN values")
            return False
        
        print(f"     Heatmap data prepared: {heatmap_df.shape}")
        
        # Create custom colormap (blue to white to red)
        cmap = LinearSegmentedColormap.from_list('custom_cmap', ['blue', 'white', 'red'])
        
        # Create the heatmap
        fig_height = max(8, len(heatmap_df) * 0.3)
        plt.figure(figsize=(8, fig_height))
        
        # Create heatmap
        ax = sns.heatmap(heatmap_df,
                        cmap=cmap,
                        annot=True,
                        fmt='.2f',
                        linewidths=0.5,
                        linecolor='gray',
                        cbar_kws={'label': 'Average Expression'},
                        square=False)
        
        # Customize the plot
        plt.title('Module Expression Across Sample Subtypes', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Sample Subtypes', fontsize=14, labelpad=10)
        plt.ylabel('Modules', fontsize=14, labelpad=10)
        
        # Rotate x-axis labels
        plt.xticks(rotation=45, ha='right', fontsize=12)
        plt.yticks(fontsize=10)
        
        # Add colorbar label
        cbar = ax.collections[0].colorbar
        cbar.set_label('Average Expression', fontsize=12, labelpad=10)
        
        # Tight layout
        plt.tight_layout()
        
        # Save the plot
        output_file = os.path.join(output_dir, 'Module_Expression_Heatmap.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Heatmap created successfully!")
        print(f"     - Modules: {len(heatmap_df)}")
        print(f"     - Sample groups: {len(unique_groups)}")
        print(f"     - Significant modules marked with *")
        print(f"     - Saved as: {output_file}")
        print(f"     - Figure size: 8 x {fig_height:.1f} inches")
        
        return True
        
    except Exception as e:
        print(f"  Error creating heatmap: {e}")
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(description='Generate interactive module overview with functional clustering')
    parser.add_argument('--input_dir', type=str, required=True,
                       help='Input directory containing module files')
    parser.add_argument('--output_dir', type=str, default='.',
                       help='Output directory for results')
    parser.add_argument('--regulator_files', type=str, required=False, default='',
                       help='Comma-separated list of regulator files (format: Type:Path,Type:Path)')
    parser.add_argument('--enrichment_method', type=str, default='auto',
                       choices=['EnrichR', 'GSEA', 'auto'],
                       help='Enrichment analysis method to use')
    parser.add_argument('--n_clusters', type=int, default=5,
                       help='Number of functional clusters to create')
    parser.add_argument('--use_megago', action='store_true', default=True,
                       help='Use megago for advanced functional clustering (default: enabled)')
    parser.add_argument('--no_megago', dest='use_megago', action='store_false',
                       help='Disable megago clustering and use simple pathway similarity')
    parser.add_argument('--prioritize_by_expression', action='store_true', default=True,
                       help='Enable expression-based module prioritization (default: enabled)')
    parser.add_argument('--no_prioritize_by_expression', dest='prioritize_by_expression', action='store_false',
                       help='Disable expression-based module prioritization')
    parser.add_argument('--group_column', type=str, default='diagnosis',
                       help='Metadata column name for grouping samples (auto-detects statistical test)')
    parser.add_argument('--coherence_threshold', type=float, default=0.6,
                       help='Coherence threshold for module filtering')
    parser.add_argument('--expression_file', type=str, default=None,
                       help='Path to expression data file (LemonPreprocessed_expression.txt)')
    parser.add_argument('--metadata_file', type=str, default=None,
                       help='Path to metadata file (DESeq_groups.txt)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = os.path.join(args.output_dir, 'Module_Overview')
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*60)
    print("INTERACTIVE MODULE OVERVIEW GENERATOR")
    print("="*60)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Regulator files: {args.regulator_files}")
    print(f"Enrichment method: {args.enrichment_method}")
    print(f"Number of clusters: {args.n_clusters}")
    print(f"Coherence threshold: {args.coherence_threshold}")
    print(f"Expression prioritization: {args.prioritize_by_expression}")
    print(f"Use megago clustering: {args.use_megago}")
    if args.prioritize_by_expression:
        print(f"Grouping column: {args.group_column}")
    
    # Load module data
    print("\nLoading module data...")
    
    # Parse and load regulator files dynamically
    regulator_configs = {}
    if args.regulator_files:
        for reg_config in args.regulator_files.split(','):
            if ':' in reg_config:
                reg_type, reg_path = reg_config.split(':', 1)
                regulator_configs[reg_type] = reg_path.strip()
    
    print(f"Found regulator types: {list(regulator_configs.keys())}")
    
    # Load regulator files dynamically
    regulators_dict = {}
    for reg_type, reg_path in regulator_configs.items():
        print(f"Reading {reg_type} regulators from {reg_path}")
        regulators_dict[reg_type] = get_regulators(reg_path)
    
    module_genes = get_regulators(os.path.join(args.input_dir, 'clusters_list.txt'))
    
    # Count regulators by type
    regulator_counts = {reg_type: len(regs) for reg_type, regs in regulators_dict.items()}
    total_regulator_modules = sum(regulator_counts.values())
    
    print(f"Loaded data for {len(module_genes)} modules")
    for reg_type, count in regulator_counts.items():
        print(f"- {reg_type} regulators: {count} modules")
    
    # Load enrichment data
    enrichment_dir = os.path.join(args.input_dir, 'Enrichment')
    # Also check for nested enrichment directories
    if not os.path.exists(enrichment_dir):
        # Look for enrichment directories recursively
        for root, dirs, files in os.walk(args.input_dir):
            if 'Enrichment' in dirs or 'enrichment' in dirs:
                enrichment_dir = os.path.join(root, [d for d in dirs if 'Enrichment' in d or 'enrichment' in d][0])
                break
    
    enrichment_data = load_enrichment_data(enrichment_dir, args.enrichment_method)
    
    # Add direction column to enrichment data for hover information
    for key in enrichment_data:
        if not enrichment_data[key].empty and '__direction__' not in enrichment_data[key].columns:
            # If no direction column, assume all are up-regulated for compatibility
            enrichment_data[key]['__direction__'] = 'Up'
    
    # Load coherence filtering - prioritize Networks/specific_modules.txt over coherence scores file
    
    print("\nLoading coherence-filtered modules...")    # First, try to load from Networks/specific_modules.txt (created by lemontree_to_network.py)
    filtered_modules_from_networks = load_filtered_modules_from_networks(args.input_dir)
    
    if filtered_modules_from_networks is not None:
        # Use the pre-filtered modules from Networks/specific_modules.txt
        modules_to_process = filtered_modules_from_networks.intersection(set(module_genes.keys()))
        print(f"Using coherence-filtered modules from Networks/specific_modules.txt")
        print(f"   Processing {len(modules_to_process)} modules that passed coherence filtering")
        
        # Load coherence scores file for informational purposes only
        coherence_file = os.path.join(args.input_dir, 'Module_coherence_scores.txt')
        _, coherence_df = load_coherence_filtered_modules(coherence_file, args.coherence_threshold)
        
    else:
        # Fallback to coherence scores file if Networks/specific_modules.txt is not available
        print("Networks/specific_modules.txt not found, falling back to coherence scores file")
        coherence_file = os.path.join(args.input_dir, 'Module_coherence_scores.txt')
        filtered_module_ids, coherence_df = load_coherence_filtered_modules(coherence_file, args.coherence_threshold)
        
        if filtered_module_ids is not None:
            print(f"Applying coherence filtering from coherence scores file...")
            all_modules = set(module_genes.keys())
            filtered_modules = set(str(m) for m in filtered_module_ids)
            modules_to_process = all_modules.intersection(filtered_modules)
            print(f"Processing {len(modules_to_process)} modules after coherence filtering")
        else:
            modules_to_process = set(module_genes.keys())
            print(f"No coherence filtering available - processing all {len(modules_to_process)} modules")
            print("   Note: Module_Overview.csv will contain all modules, not just coherence-filtered ones")
    
    # Create module data structure
    module_data = []
    module_pathways = {}  # For clustering
    
    # Create coherence lookup dictionary if available
    coherence_lookup = {}
    if coherence_df is not None and not coherence_df.empty:
        coherence_lookup = dict(zip(
            coherence_df['Module'].astype(str), 
            coherence_df['Coherence_Score']
        ))
    
    for module_id in sorted(modules_to_process, key=lambda x: int(x) if x.isdigit() else float('inf')):
        module_entry = {
            'Module': module_id,
            'Module_genes': '|'.join(module_genes.get(module_id, [])) if module_genes.get(module_id) else 'NA',
            'Coherence': coherence_lookup.get(str(module_id), 'NA')
        }
        
        # Add regulator data dynamically
        for reg_type, regulators in regulators_dict.items():
            column_name = f'{reg_type}_regulators'
            module_entry[column_name] = '|'.join(regulators.get(module_id, [])) if regulators.get(module_id) else 'NA'
        
        # Add pathway enrichment data
        pathway_types = {
            'mf': 'Top_3_pathways_molecular_function', 
            'cc': 'Top_3_pathways_cellular_component',
            'bp': 'Top_3_pathways_biological_process',
            'reactome': 'Top_3_pathways_Reactome',
            'kegg': 'Top_3_pathways_KEGG'
        }
        
        all_pathways = set()
        
        for pathway_type, column_name in pathway_types.items():
            # Default NA
            module_entry[column_name] = 'NA'
            if pathway_type not in enrichment_data:
                continue

            df = enrichment_data[pathway_type]
            if df is None or df.empty:
                continue

            # Accept several possible term column names
            term_col = next((c for c in df.columns if c.lower() in ('term','description','pathway','name')), None)
            mod_col = next((c for c in df.columns if c.lower() in ('module','cluster','moduleid','mod')), 'Module')

            # Ensure Module column exists and compare as string
            if mod_col not in df.columns:
                # cannot match modules
                continue

            try:
                mod_series = df[mod_col].astype(str).str.strip()
            except Exception:
                mod_series = df[mod_col].astype(str)

            # Use term_col if present, otherwise try common alternatives
            if term_col and term_col in df.columns:
                terms_series = df[term_col].astype(str)
            elif 'Term' in df.columns:
                terms_series = df['Term'].astype(str)
            else:
                # fallback: take first non-numeric column besides Module
                other_cols = [c for c in df.columns if c != mod_col]
                if other_cols:
                    terms_series = df[other_cols[0]].astype(str)
                else:
                    continue

            # Select rows matching module_id (string compare)
            mask = (mod_series == str(module_id))
            top_terms = terms_series[mask].head(3).tolist()

            if top_terms:
                module_entry[column_name] = '|'.join([t for t in top_terms if t and t.lower() != 'nan'])
                all_pathways.update([t for t in top_terms if t and t.lower() != 'nan'])
        
        # Store pathways for clustering
        module_pathways[module_id] = list(all_pathways)
        
        module_data.append(module_entry)
    
    print(f"Prepared data for {len(module_data)} modules")
    
    # Perform expression-based prioritization if requested
    expression_priority = {}
    expression_results_df = pd.DataFrame()
    
    if args.prioritize_by_expression:
        print(f"\nPerforming expression-based module prioritization...")
        
        # Look for expression and metadata files (search recursively to include Preprocessing/)
        expression_file = args.expression_file
        metadata_file = args.metadata_file
        
        # If files not provided via command line, search for them
        if expression_file is None or metadata_file is None:
            # priority candidate names (preprocessed outputs)
            expr_candidates = ['lemonpreprocessed_expression.txt', 'normalized_counts.txt', 'expression_matrix.txt', 'counts.txt', 'rna_seq.txt']
            meta_candidates = ['deseq_groups.txt', 'DESeq_groups.txt', 'metadata.txt', 'sample_info.txt', 'groups.txt']

            # walk directory to find any candidate file (prefer Preprocessing/ or top-level)
            for root, dirs, files in os.walk(args.input_dir):
                for f in files:
                    lf = f.lower()
                    if expression_file is None and any(lf == c.lower() for c in expr_candidates):
                        expression_file = os.path.join(root, f)
                    if metadata_file is None and any(lf == c.lower() for c in meta_candidates):
                        metadata_file = os.path.join(root, f)
                    if expression_file and metadata_file:
                        break
                if expression_file and metadata_file:
                    break
        
        if expression_file and metadata_file:
            print(f"Found expression file: {expression_file}")
            print(f"Found metadata file: {metadata_file}")
            
            expression_priority, expression_results_df = auto_prioritize_modules_expression(
                group_column=args.group_column,
                modules_dict=module_genes,
                expression_file=expression_file,
                deseq_groups_file=metadata_file
            )
            
            if not expression_results_df.empty:
                # Save expression results
                expr_results_file = os.path.join(output_dir, 'module_expression_analysis.csv')
                expression_results_df.to_csv(expr_results_file, index=False)
                print(f"Expression analysis results saved to: {expr_results_file}")
        else:
            print("Warning: Could not find expression or metadata files for prioritization")
            print("Available files in input directory:")
            for f in os.listdir(args.input_dir):
                if f.endswith(('.txt', '.csv')):
                    print(f"  - {f}")
    
    # Perform functional clustering using megago
    print(f"\nPerforming functional clustering...")
    module_clusters, go_similarity_matrix = megago_cluster_modules(module_pathways, args.n_clusters, args.use_megago, 
                                           enrichment_data, output_dir)
    
    # Create interactive visualizations
    print(f"\nCreating interactive visualizations...")
    
    # Create network visualization with MegaGO clustering
    network_fig = create_interactive_network_visualization(module_data, module_clusters, output_dir, enrichment_data, go_similarity_matrix)
    
    # COMMENTED OUT: Create cluster heatmap - not needed in output
    # heatmap_fig, cluster_stats = create_cluster_heatmap(module_data, module_clusters, output_dir)
    cluster_stats = pd.DataFrame()  # Empty DataFrame to avoid errors
    
    # Add cluster information to module data
    for module_entry in module_data:
        module_id = str(module_entry['Module'])
        module_entry['Functional_Cluster'] = module_clusters.get(module_id, 0)
        
        # Add expression rank if available
        if expression_priority and module_id in expression_priority:
            module_entry['Expression_rank'] = expression_priority[module_id]
        else:
            module_entry['Expression_rank'] = 'NA'
        
        # Add expression p-value if available
        if not expression_results_df.empty:
            matching_rows = expression_results_df[expression_results_df['Module'].astype(str) == module_id]
            if not matching_rows.empty:
                module_entry['Expression_adjusted_pval'] = matching_rows.iloc[0]['P_adjusted']
            else:
                module_entry['Expression_adjusted_pval'] = 'NA'
        else:
            module_entry['Expression_adjusted_pval'] = 'NA'
    
    # Convert to DataFrame and sort appropriately
    module_df = pd.DataFrame(module_data)
    
    # Sort by expression rank if available, otherwise by cluster and module
    if args.prioritize_by_expression and not expression_results_df.empty:
        # Sort by expression rank (prioritizing significant modules)
        module_df['sort_key'] = module_df['Expression_rank'].apply(
            lambda x: int(x) if x != 'NA' else float('inf')
        )
        module_df = module_df.sort_values(['sort_key', 'Functional_Cluster', 'Module']).reset_index(drop=True)
        module_df = module_df.drop('sort_key', axis=1)
        print("Modules sorted by expression significance")
    else:
        # Sort by cluster, then module
        module_df = module_df.sort_values(['Functional_Cluster', 'Module']).reset_index(drop=True)
    
    # Save enhanced module overview
    output_file = os.path.join(output_dir, 'Module_Overview.csv')
    module_df.to_csv(output_file, sep='\t', index=False)
    print(f"Enhanced module overview saved to: {output_file}")
    
    # Create Module Expression Heatmap if differential expression was performed
    if args.prioritize_by_expression and expression_file and metadata_file:
        # Prepare module p-values dictionary
        if not expression_results_df.empty:
            module_pvalues = dict(zip(
                expression_results_df['Module'].astype(str), 
                expression_results_df['P_adjusted']
            ))
        else:
            module_pvalues = {}
        
        # Generate the heatmap
        create_module_expression_heatmap(
            module_genes=module_genes,
            modules_to_process=modules_to_process,
            expression_file=expression_file,
            metadata_file=metadata_file,
            module_pvalues=module_pvalues,
            output_dir=args.output_dir,  # Save to main output dir, not Module_Overview subdir
            group_column=args.group_column
        )
    
    # Save cluster assignments
    # COMMENTED OUT: Save module functional clusters - not needed in output
    # cluster_file = os.path.join(output_dir, 'module_functional_clusters.csv')
    # cluster_df = pd.DataFrame(list(module_clusters.items()), columns=['Module', 'Functional_Cluster'])
    # cluster_df.to_csv(cluster_file, index=False)
    # print(f"Module cluster assignments saved to: {cluster_file}")
    
    print("Module functional clusters file creation skipped (not needed in output)")
    
    # Generate summary
    print(f"\n" + "="*80)
    print("                    INTERACTIVE MODULE OVERVIEW COMPLETE!")
    print("                           Lemonite Analysis")
    print("="*80)
    print(f"Total modules processed: {len(module_data)}")
    
    # Add filtering information
    if filtered_modules_from_networks is not None:
        print(f"Module filtering: Networks/specific_modules.txt (coherence-filtered by lemontree_to_network.py)")
    elif 'filtered_module_ids' in locals() and filtered_module_ids is not None:
        print(f"Module filtering: Module_coherence_scores.txt (threshold = {args.coherence_threshold})")
    else:
        print(f"Module filtering: None (all modules included)")
    
    print(f"Functional clusters created: {len(set(module_clusters.values()))}")
    print(f"Modules with pathway data: {sum(1 for pathways in module_pathways.values() if pathways)}")
    
    if args.prioritize_by_expression and not expression_results_df.empty:
        significant_count = len(expression_results_df[expression_results_df['P_adjusted'] < 0.05])
        print(f"Expression-prioritized modules: {len(expression_results_df)}")
        print(f"Significantly different modules: {significant_count}")
    
    # Print cluster summary
    cluster_counts = Counter(module_clusters.values())
    print(f"\nCluster distribution:")
    for cluster_id in sorted(cluster_counts.keys()):
        count = cluster_counts[cluster_id]
        print(f"  Cluster {cluster_id}: {count} modules")
    
    print(f"\nFiles created:")
    print(f"  - {output_file}")
    # COMMENTED OUT: These files are no longer generated
    # print(f"  - {cluster_file}")
    # print(f"  - {os.path.join(output_dir, 'cluster_characteristics_heatmap.html')}")
    # print(f"  - {os.path.join(output_dir, 'cluster_statistics.csv')}")
    print(f"  - {os.path.join(output_dir, 'interactive_module_network.html')}")
    
    if args.prioritize_by_expression and not expression_results_df.empty:
        print(f"  - {os.path.join(output_dir, 'module_expression_analysis.csv')}")
    
    print("="*80)
    print("                     SUCCESS! Analysis Complete")
    print("              Open HTML files in web browser for visualizations!")
    print("="*80)

if __name__ == "__main__":
    main()
