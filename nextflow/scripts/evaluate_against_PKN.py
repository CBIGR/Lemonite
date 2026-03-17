#!/usr/bin/env python3

"""
This script evaluates the performance of the Lemonite network against PKN
This will be done in two different approaches:
    1. Classic evaluation metrics (accuracy, precision, recall, F1-score) of Lemonite network against the ground truth data
    2. Evaluation of average shortest paths for metabolite-gene pairs from Lemonite network in the PKN
"""

import os
import sys
import pandas as pd
import numpy as np
import networkx as nx
from multiprocessing import Pool
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import argparse
import shutil
import glob

def get_regulators(regfile):
    """Function to create a dictionary with all regulators per module"""
    regs = {}
    if not os.path.exists(regfile):
        print(f"Warning: {regfile} not found, returning empty dictionary")
        return regs
    
    with open(regfile) as f:
        for line in f:
            if line.strip():
                parts = line.rstrip().split('\t')
                if len(parts) >= 2:
                    module = parts[0]
                    regulators = parts[1].split('|')
                    regs[module] = regulators
    return regs

def load_metabolite_mapping(mapping_file):
    """Load metabolite mapping from Query to HMDB"""
    if not os.path.exists(mapping_file):
        print(f"Warning: Metabolite mapping file {mapping_file} not found")
        return {}
    
    try:
        metabolite_mapping = pd.read_csv(mapping_file, sep=',')
        print(f"Metabolite mapping file columns: {list(metabolite_mapping.columns)}")
        
        # Check if required columns exist
        if 'Query' not in metabolite_mapping.columns:
            print(f"Error: 'Query' column not found in metabolite mapping file. Available columns: {list(metabolite_mapping.columns)}")
            return {}
        if 'HMDB' not in metabolite_mapping.columns:
            print(f"Error: 'HMDB' column not found in metabolite mapping file. Available columns: {list(metabolite_mapping.columns)}")
            return {}
        
        metabolite_mapping = metabolite_mapping.set_index('Query')['HMDB'].dropna().to_dict()
        return metabolite_mapping
    except Exception as e:
        print(f"Error loading metabolite mapping: {e}")
        return {}

def draw_subnetwork(module, target_genes, regulators_dict, PKN, name_to_hmdb):
    """
    Draw subnetwork visualization for a module showing connections in the PKN
    
    Parameters:
    module (str): Module identifier
    target_genes (list): List of target genes for this module
    regulators_dict (dict): Dictionary mapping regulator types to lists of regulators
    PKN (nx.Graph): NetworkX graph representing the PKN
    name_to_hmdb (dict): Dictionary mapping metabolite names to HMDB IDs
    """
    # Create reverse mapping from HMDB ID to metabolite name for display
    hmdb_to_name = {v: k for k, v in name_to_hmdb.items()}

    # Extract regulator lists from the dictionary
    TF_regulators = regulators_dict.get('TF', [])
    metabolite_regulators = regulators_dict.get('Metabolite', [])
    lipids_regulators = regulators_dict.get('Lipid', [])

    # Get metabolite regulators, TF regulators, and target genes for this module
    # These are already passed as parameters, no need to index dictionaries
    print(f'Module {module} has {len(metabolite_regulators)} metabolite regulators, {len(TF_regulators)} TF regulators, {len(lipids_regulators)} lipid regulators and {len(target_genes)} target genes')

    # Initialize lists to store edges and nodes for visualization
    direct_edges = []        # edges directly connecting metabolite -> target gene
    only_PKN_edge = []       # other edges in PKN between nodes
    reachable_targets = set() # target genes reachable directly
    met_regs = []             # list of metabolite regulators (mapped to HMDB IDs)

    # Initialize an empty graph for drawing
    to_draw = nx.Graph()

    # ---- STEP 1: Process metabolite regulators ----
    for reg in metabolite_regulators:
        if reg not in name_to_hmdb:
            print(f'{reg} does not have a valid HMDB ID in the list you provided')
            continue
        
        regulator = name_to_hmdb[reg]  # map name to HMDB ID
        met_regs.append(regulator)

        if regulator in PKN.nodes():  # if the regulator exists in the PKN
            to_draw.add_node(regulator)

            # Check for direct connections from regulator to each target gene
            for target in target_genes:
                try:
                    shortest_path = nx.shortest_path_length(PKN, source=regulator, target=target)
                    if shortest_path == 1:  # direct connection
                        reachable_targets.add(target)
                        to_draw.add_node(target)
                        to_draw.add_edge(regulator, target, weight=1)  # direct edge → weight=1
                        direct_edges.append((regulator, target))
                except nx.NetworkXNoPath:
                    pass
                except nx.NodeNotFound:
                    pass

    # ---- STEP 2: Process TF regulators ----
    reachable_TFs = set()
    for TF in TF_regulators:
        for met in met_regs:
            try:
                shortest_path = nx.shortest_path_length(PKN, source=TF, target=met)

                if shortest_path == 1:  # TF connected directly to metabolite
                    reachable_TFs.add(TF)
                    print(f'{TF} can be reached from regulator {met} in 1 step')
                    to_draw.add_node(TF)
                    to_draw.add_edge(TF, met, weight=1)  # weight=1 for direct

                elif shortest_path == 2:  # TF connected via an intermediate node
                    print(f'{TF} can be reached from regulator {met} in 2 steps')
                    reachable_TFs.add(TF)
                    path = nx.shortest_path(PKN, source=TF, target=met)  # fixed variable name

                    # add all nodes and edges in the path
                    for i in range(len(path)-1):
                        to_draw.add_node(path[i])
                        to_draw.add_node(path[i+1])
                        to_draw.add_edge(path[i], path[i+1], weight=2)  # indirect → weight=2

            except nx.NetworkXNoPath:
                pass
            except nx.NodeNotFound:
                pass

    # ---- STEP 2.5: Process lipid regulators ----
    lipid_regs = []
    for lipid in lipids_regulators:
        lipid_regs.append(lipid)
        if lipid in PKN.nodes():  # if the lipid regulator exists in the PKN
            to_draw.add_node(lipid)

            # Check for direct connections from lipid regulator to each target gene
            for target in target_genes:
                try:
                    shortest_path = nx.shortest_path_length(PKN, source=lipid, target=target)
                    if shortest_path == 1:  # direct connection
                        reachable_targets.add(target)
                        to_draw.add_node(target)
                        to_draw.add_edge(lipid, target, weight=1)  # direct edge → weight=1
                        direct_edges.append((lipid, target))
                except nx.NetworkXNoPath:
                    pass
                except nx.NodeNotFound:
                    pass

    # ---- STEP 3: Annotate nodes with their types ----
    node_attributes = pd.DataFrame(columns=['node', 'type'])
    node_colors = []
    for node in to_draw.nodes():
        if node in name_to_hmdb.values():
            node_colors.append('red')  # metabolite regulator
            node_attributes = pd.concat([node_attributes, pd.DataFrame([{'node': node, 'type': 'metabolite'}])], ignore_index=True)

        elif node in TF_regulators:
            node_colors.append('green')  # TF regulator
            node_attributes = pd.concat([node_attributes, pd.DataFrame([{'node': node, 'type': 'TF'}])], ignore_index=True)

        elif node in lipids_regulators:
            node_colors.append('purple')  # lipid regulator
            node_attributes = pd.concat([node_attributes, pd.DataFrame([{'node': node, 'type': 'lipid'}])], ignore_index=True)

        elif node in target_genes:
            node_colors.append('orange')  # target gene
            node_attributes = pd.concat([node_attributes, pd.DataFrame([{'node': node, 'type': 'target'}])], ignore_index=True)

        else:
            node_colors.append('grey')  # intermediate (bridge) node
            node_attributes = pd.concat([node_attributes, pd.DataFrame([{'node': node, 'type': 'bridge'}])], ignore_index=True)

    # ---- STEP 4: Add additional edges from PKN (between included nodes) ----
    subgraph = PKN.subgraph(to_draw.nodes())
    for edge in subgraph.edges():
        if edge not in to_draw.edges():
            only_PKN_edge.append(edge)
            to_draw.add_edge(edge[0], edge[1], weight=3)  # additional edges → weight=3

    # ---- STEP 5: Clean up graph ----
    to_draw.remove_edges_from(nx.selfloop_edges(to_draw))  # remove self-loops

    # ---- STEP 8: Draw the network ----
    # Dynamically adjust figure size based on max label length
    if len(to_draw.nodes()) == 0:
        print(f'No valid connections found for module {module}')
        return
        
    max_label_length = max([len(str(label)) for label in to_draw.nodes()])
    width = 3 + max_label_length * 0.2
    height = 3 + max_label_length * 0.1  # optioneel iets hoger
    plt.figure(figsize=(width, height))

    ax = plt.gca()
    pos = nx.spring_layout(to_draw, k=3.0, iterations=100, seed=42)  # meer spacing

    # Separate nodes by type
    hmdb_nodes = [n for n in to_draw.nodes() if n in name_to_hmdb.values()]
    tf_nodes = [n for n in to_draw.nodes() if n in TF_regulators]
    lipid_nodes = [n for n in to_draw.nodes() if n in lipids_regulators]
    target_nodes = [n for n in to_draw.nodes() if n in target_genes]
    bridge_nodes = [n for n in to_draw.nodes() if n not in hmdb_nodes + tf_nodes + lipid_nodes + target_nodes]

    # Dynamically scale node sizes based on label length
    sizes_tf = [500 + len(str(node)) * 80 for node in tf_nodes]
    sizes_lipid = [500 + len(str(node)) * 80 for node in lipid_nodes]
    sizes_target = [500 + len(str(node)) * 80 for node in target_nodes]
    sizes_bridge = [500 + len(str(node)) * 80 for node in bridge_nodes]

    # Draw TF nodes as circles (green)
    nx.draw_networkx_nodes(
        to_draw, pos,
        nodelist=tf_nodes,
        node_color='lightgreen',
        node_shape='o',
        node_size=sizes_tf
    )

    # Draw lipid nodes as circles (purple)
    nx.draw_networkx_nodes(
        to_draw, pos,
        nodelist=lipid_nodes,
        node_color='purple',
        node_shape='o',
        node_size=sizes_lipid
    )

    # Draw target nodes as circles (orange)
    nx.draw_networkx_nodes(
        to_draw, pos,
        nodelist=target_nodes,
        node_color='orange',
        node_shape='o',
        node_size=sizes_target
    )

    # Draw bridge nodes as circles (gray)
    nx.draw_networkx_nodes(
        to_draw, pos,
        nodelist=bridge_nodes,
        node_color='lightgrey',
        node_shape='o',
        node_size=sizes_bridge
    )

    # Draw labels for HMDB nodes using metabolite names instead of HMDB IDs
    metabolite_labels = {}
    for node in hmdb_nodes:
        if node in hmdb_to_name:
            metabolite_labels[node] = hmdb_to_name[node]
        else:
            metabolite_labels[node] = node  # fallback to HMDB ID if name not found

    nx.draw_networkx_labels(
        to_draw, pos,
        labels=metabolite_labels,
        font_size=6,
        font_weight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='red', edgecolor='none')
    )

    # Draw edges
    nx.draw_networkx_edges(to_draw, pos, edge_color='black')

    # Draw labels for all other nodes (excluding HMDB nodes which were labeled manually)
    for node in tf_nodes + target_nodes + bridge_nodes:
        x, y = pos[node]
        ax.text(x, y, node, ha='center', va='center', fontsize=5, fontweight='bold')

    plt.axis('off') 
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    plt.title(f'Subnetwork of module {module} in the Lemonite PKN', fontsize=10)

    # Create output directory
    output_dir = os.path.join(os.getcwd(), 'Networks', 'subnetworks')
    os.makedirs(output_dir, exist_ok=True)
    
    # Save figure to file
    plt.savefig(os.path.join(output_dir, f'graph_{module}.png'))
    plt.close()
    
    print(f'Graph with {to_draw.number_of_nodes()} nodes and {to_draw.number_of_edges()} edges')

def calculate_ppi_enrichment(module2genes, PKN_graph, all_module_genes):
    """
    Calculate PPI enrichment for each module using hypergeometric test.
    
    Parameters:
    -----------
    module2genes : dict
        Dictionary mapping module IDs to lists of genes
    PKN_graph : nx.Graph
        NetworkX graph representing the PKN
    all_module_genes : set
        Set of all genes across all modules
    
    Returns:
    --------
    pd.DataFrame : DataFrame with enrichment results containing columns:
        Module, N_genes, N_PPIs_observed, N_PPIs_expected, Fold_enrichment, P_value, FDR, PPI_density
    """
    print("\n" + "=" * 60)
    print("Calculating PPI enrichment for modules...")
    print("=" * 60)
    
    # Get all genes that are in any module
    all_genes_in_modules = list(all_module_genes)
    N_total_genes = len(all_genes_in_modules)
    
    # Count total possible PPIs among genes in modules
    total_possible_ppis = (N_total_genes * (N_total_genes - 1)) // 2
    
    # Create a set of PPI pairs for O(1) lookup (store both directions)
    ppi_set = set()
    for edge in PKN_graph.edges():
        node1, node2 = edge[0], edge[1]
        # Only include if both nodes are genes in modules (not metabolites/lipids)
        # Skip edges containing HMDB IDs (metabolites)
        if ('HMDB' not in str(node1) and 'HMDB' not in str(node2) and 
            node1 in all_module_genes and node2 in all_module_genes):
            # Add both directions (undirected graph)
            ppi_set.add((node1, node2))
            ppi_set.add((node2, node1))
    
    # Count total observed PPIs in the PKN among module genes
    total_observed_ppis = len(ppi_set) // 2  # Divide by 2 because we stored both directions
    
    print(f"Total genes in modules: {N_total_genes}")
    print(f"Total possible gene pairs: {total_possible_ppis}")
    print(f"Total PPIs in PKN among module genes: {total_observed_ppis}")
    
    if total_possible_ppis > 0:
        print(f"Overall PPI density: {total_observed_ppis / total_possible_ppis:.4f}")
    else:
        print("Overall PPI density: N/A (no gene pairs)")
    
    # Calculate enrichment for each module
    enrichment_results = []
    
    for module in module2genes.keys():
        module_genes = module2genes[module]
        n_genes = len(module_genes)
        
        if n_genes < 2:
            continue  # Skip modules with less than 2 genes
        
        # Count PPIs in this module
        n_ppis_in_module = 0
        for i, gene1 in enumerate(module_genes):
            for gene2 in module_genes[i+1:]:
                if (gene1, gene2) in ppi_set:
                    n_ppis_in_module += 1
        
        # Hypergeometric test:
        # Population: all possible gene pairs in all modules
        # Success in population: gene pairs that have a PPI
        # Sample: all possible pairs in this module
        # Observed successes: PPIs in this module
        
        n_possible_pairs_in_module = (n_genes * (n_genes - 1)) // 2
        
        # Calculate expected number of PPIs
        expected_ppis = n_possible_pairs_in_module * (total_observed_ppis / total_possible_ppis) if total_possible_ppis > 0 else 0
        
        # Fold enrichment
        fold_enrichment = n_ppis_in_module / expected_ppis if expected_ppis > 0 else 0
        
        # Hypergeometric p-value (probability of seeing this many or more PPIs by chance)
        # Using survival function (1 - CDF) to get upper tail probability
        if total_possible_ppis > 0 and n_possible_pairs_in_module > 0:
            p_value = hypergeom.sf(
                n_ppis_in_module - 1,  # -1 because sf is P(X >= k), we want P(X > k-1) = P(X >= k)
                M=total_possible_ppis,  # Total possible pairs
                n=total_observed_ppis,  # Total pairs with PPIs
                N=n_possible_pairs_in_module  # Pairs in module
            )
        else:
            p_value = 1.0
        
        enrichment_results.append({
            'Module': module,
            'N_genes': n_genes,
            'N_PPIs_observed': n_ppis_in_module,
            'N_PPIs_expected': expected_ppis,
            'Fold_enrichment': fold_enrichment,
            'P_value': p_value,
            'PPI_density': n_ppis_in_module / n_possible_pairs_in_module if n_possible_pairs_in_module > 0 else 0
        })
    
    # Convert to DataFrame
    enrichment_df = pd.DataFrame(enrichment_results)
    
    if len(enrichment_df) == 0:
        print("No modules with >= 2 genes found for enrichment analysis")
        return enrichment_df
    
    # Add FDR correction using Benjamini-Hochberg
    try:
        # Try scipy's newer false_discovery_control function (scipy >= 1.9.0)
        from scipy.stats import false_discovery_control
        enrichment_df['FDR'] = false_discovery_control(enrichment_df['P_value'].values, method='bh')
    except ImportError:
        # Fallback to manual Benjamini-Hochberg correction
        from scipy.stats import rankdata
        p_values = enrichment_df['P_value'].values
        ranked = rankdata(p_values)
        n = len(p_values)
        enrichment_df['FDR'] = np.minimum(p_values * n / ranked, 1.0)
    
    # Sort by p-value
    enrichment_df = enrichment_df.sort_values('P_value')
    
    print(f"\nPPI Enrichment Results (top 10 most enriched modules):")
    print(enrichment_df.head(10).to_string())
    
    # Summary statistics
    n_significant = (enrichment_df['FDR'] < 0.05).sum()
    print(f"\nNumber of modules significantly enriched for PPIs (FDR < 0.05): {n_significant} / {len(enrichment_df)}")
    if len(enrichment_df) > 0:
        print(f"Mean fold enrichment: {enrichment_df['Fold_enrichment'].mean():.2f}")
        print(f"Median fold enrichment: {enrichment_df['Fold_enrichment'].median():.2f}")
    
    return enrichment_df

def calculate_metabolite_gene_enrichment(module2genes, module2mets, module2lipids, interactions, metabolite_mapping, PKN_graph, all_module_genes):
    """
    Calculate metabolite-gene interaction enrichment for each module using hypergeometric test.
    
    Tests whether a module's metabolite regulators have more known interactions with the
    module's target genes than expected by chance.
    
    Parameters:
    -----------
    module2genes : dict
        Dictionary mapping module IDs to lists of genes
    module2mets : dict
        Dictionary mapping module IDs to lists of metabolite regulators
    module2lipids : dict
        Dictionary mapping module IDs to lists of lipid regulators
    interactions : pd.DataFrame
        DataFrame with columns 'HMDB' and 'All_interactions' (pipe-separated gene names)
    metabolite_mapping : dict
        Dictionary mapping metabolite names to HMDB IDs
    PKN_graph : nx.Graph
        NetworkX graph representing the PKN (used for lipid-gene interactions)
    all_module_genes : set
        Set of all genes across all modules
    
    Returns:
    --------
    pd.DataFrame : DataFrame with enrichment results
    """
    print("\n" + "=" * 60)
    print("Calculating metabolite-gene interaction enrichment for modules...")
    print("=" * 60)
    
    # Build set of known metabolite/lipid-gene interactions from PKN
    known_interactions = set()  # (metabolite_hmdb_or_lipid, gene) pairs
    
    # From interactions file (HMDB metabolites)
    if len(interactions) > 0 and 'HMDB' in interactions.columns and 'All_interactions' in interactions.columns:
        for _, row in interactions.iterrows():
            hmdb = row['HMDB']
            if pd.notna(row.get('All_interactions')):
                for gene in str(row['All_interactions']).split('|'):
                    gene = gene.strip()
                    if gene:
                        known_interactions.add((hmdb, gene))
    
    # From PKN graph: metabolite nodes (HMDB IDs) and their gene neighbors
    for node in PKN_graph.nodes():
        if 'HMDB' in str(node):
            for neighbor in PKN_graph.neighbors(node):
                if 'HMDB' not in str(neighbor):
                    known_interactions.add((str(node), str(neighbor)))
    
    # Also include lipid-gene interactions from PKN
    all_lipids_in_modules = set()
    for lipids in module2lipids.values():
        all_lipids_in_modules.update(lipids)
    
    for lipid in all_lipids_in_modules:
        if lipid in PKN_graph.nodes():
            for neighbor in PKN_graph.neighbors(lipid):
                if lipid != neighbor:
                    known_interactions.add((lipid, str(neighbor)))
    
    print(f"Total known metabolite/lipid-gene interactions in PKN: {len(known_interactions)}")
    
    # Collect all metabolite/lipid IDs across all modules
    all_regulator_ids = set()
    module_regulator_ids = {}  # module -> list of metabolite/lipid IDs
    
    for module in set(list(module2mets.keys()) + list(module2lipids.keys())):
        mod_ids = []
        
        # Metabolite regulators -> map to HMDB
        for met in module2mets.get(module, []):
            hmdb = metabolite_mapping.get(met)
            if hmdb:
                mod_ids.append(hmdb)
                all_regulator_ids.add(hmdb)
        
        # Lipid regulators -> use directly
        for lipid in module2lipids.get(module, []):
            mod_ids.append(lipid)
            all_regulator_ids.add(lipid)
        
        module_regulator_ids[module] = mod_ids
    
    total_regulators = len(all_regulator_ids)
    total_genes = len(all_module_genes)
    total_possible_pairs = total_regulators * total_genes
    
    # Count total known interactions among module regulators and module genes
    total_known = 0
    for reg_id in all_regulator_ids:
        for gene in all_module_genes:
            if (reg_id, gene) in known_interactions:
                total_known += 1
    
    print(f"Total metabolite/lipid regulators across modules: {total_regulators}")
    print(f"Total genes across modules: {total_genes}")
    print(f"Total possible regulator-gene pairs: {total_possible_pairs}")
    print(f"Total known interactions among module regulators and genes: {total_known}")
    
    if total_possible_pairs > 0:
        print(f"Overall interaction density: {total_known / total_possible_pairs:.6f}")
    
    # Calculate enrichment for each module
    enrichment_results = []
    
    for module in module2genes.keys():
        module_genes = module2genes[module]
        mod_regs = module_regulator_ids.get(module, [])
        
        if not mod_regs or not module_genes:
            continue
        
        # Count observed interactions
        n_observed = 0
        for reg_id in mod_regs:
            for gene in module_genes:
                if (reg_id, gene) in known_interactions:
                    n_observed += 1
        
        # Possible pairs in this module
        n_possible = len(mod_regs) * len(module_genes)
        
        if n_possible == 0 or total_possible_pairs == 0:
            continue
        
        # Expected number of interactions
        expected = n_possible * (total_known / total_possible_pairs) if total_possible_pairs > 0 else 0
        
        # Fold enrichment
        fold_enrichment = n_observed / expected if expected > 0 else 0
        
        # Hypergeometric p-value
        if total_possible_pairs > 0 and n_possible > 0 and total_known > 0:
            p_value = hypergeom.sf(
                n_observed - 1,
                M=total_possible_pairs,
                n=total_known,
                N=n_possible
            )
        else:
            p_value = 1.0
        
        enrichment_results.append({
            'Module': module,
            'N_metabolites': len(mod_regs),
            'N_genes': len(module_genes),
            'N_interactions_observed': n_observed,
            'N_interactions_expected': expected,
            'Fold_enrichment': fold_enrichment,
            'P_value': p_value,
            'Interaction_density': n_observed / n_possible if n_possible > 0 else 0
        })
    
    enrichment_df = pd.DataFrame(enrichment_results)
    
    if len(enrichment_df) == 0:
        print("No modules with metabolite regulators and genes found for enrichment analysis")
        return enrichment_df
    
    # Add FDR correction using Benjamini-Hochberg
    try:
        from scipy.stats import false_discovery_control
        enrichment_df['FDR'] = false_discovery_control(enrichment_df['P_value'].values, method='bh')
    except ImportError:
        from scipy.stats import rankdata
        p_values = enrichment_df['P_value'].values
        ranked = rankdata(p_values)
        n = len(p_values)
        enrichment_df['FDR'] = np.minimum(p_values * n / ranked, 1.0)
    
    # Sort by p-value
    enrichment_df = enrichment_df.sort_values('P_value')
    
    print(f"\nMetabolite-Gene Interaction Enrichment Results (top 10):")
    print(enrichment_df.head(10).to_string())
    
    n_significant = (enrichment_df['FDR'] < 0.05).sum()
    print(f"\nModules significantly enriched (FDR < 0.05): {n_significant} / {len(enrichment_df)}")
    if len(enrichment_df) > 0:
        print(f"Mean fold enrichment: {enrichment_df['Fold_enrichment'].mean():.2f}")
        print(f"Median fold enrichment: {enrichment_df['Fold_enrichment'].median():.2f}")
    
    return enrichment_df

def evaluate_network_against_pkn(args):
    """Main evaluation function"""
    
    # Set up paths
    workdir = args.workdir
    os.chdir(workdir)
    
    lemonnetwork_file = args.network_file
    PKNnetwork_file = args.pkn_file
    metabolite_interactions_file = args.metabolite_interactions_file
    annotated_mets = args.metabolite_mapping_file
    
    # Parse regulator files dynamically
    regulator_files = {}
    if args.regulator_files:
        for reg_spec in args.regulator_files.split(','):
            reg_type, reg_path = reg_spec.split(':', 1)
            reg_type = reg_type.strip()
            reg_path = reg_path.strip()
            regulator_files[reg_type] = reg_path
    
    print(f"Configured regulator types: {list(regulator_files.keys())}")
    
    # Module and regulator files
    clusterfile = args.cluster_file
    
    print("Loading LemonTree network...")
    
    # Try to find the main network file
    network_file = None
    possible_patterns = [
        "LemonNetwork_*.txt",
        "Metabolites2targets_*.txt",
        "TFs2targets_*.txt"
    ]
    
    import glob
    for pattern in possible_patterns:
        files = glob.glob(pattern)
        if files:
            network_file = files[0]
            print(f"Found network file: {network_file}")
            break
    
    if not network_file:
        print("Error: No network file found")
        print("Available files:")
        for f in os.listdir('.'):
            if f.endswith('.txt'):
                print(f"  {f}")
        return
    
    # Read the network file and create a comprehensive network
    if "LemonNetwork" in network_file:
        # Full network file
        lemonnetwork = pd.read_csv(network_file, sep='\t')
        print(f"Loaded full network with {len(lemonnetwork)} interactions")
    else:
        # Rebuild network from component files
        print("Rebuilding network from component files...")
        network_parts = []
        
        # Find all target files
        metabolite_files = glob.glob("Metabolites2targets_*.txt")
        tf_files = glob.glob("TFs2targets_*.txt")
        
        for met_file in metabolite_files:
            if os.path.exists(met_file):
                df = pd.read_csv(met_file, sep='\t')
                df['Type'] = 'Metabolite-gene'
                network_parts.append(df)
                
        for tf_file in tf_files:
            if os.path.exists(tf_file):
                df = pd.read_csv(tf_file, sep='\t')
                df['Type'] = 'TF-gene'
                network_parts.append(df)
        
        if network_parts:
            lemonnetwork = pd.concat(network_parts, ignore_index=True)
            print(f"Rebuilt network with {len(lemonnetwork)} interactions from {len(network_parts)} files")
        else:
            print("Error: Could not rebuild network from available files")
            return
    
    # Dynamically collect regulator-to-gene interactions for all regulator types
    # Build Type patterns from regulator configuration
    regulator_type_patterns = {}
    for reg_type in regulator_files.keys():
        # Handle both "TypeName-gene" and "TypeNames-gene" (singular and plural)
        regulator_type_patterns[reg_type] = [f"{reg_type}-gene", f"{reg_type}s-gene"]
    
    # Collect all regulator-gene interactions
    regulators2genes = {}  # Will store {regulator_type: {regulator: [genes]}}
    all_genes = set()
    
    # Check which format we have
    has_targets_column = 'Targets' in lemonnetwork.columns
    has_target_column = 'Target' in lemonnetwork.columns
    
    if not has_targets_column and not has_target_column:
        print("Error: Neither 'Targets' nor 'Target' column found in network")
        print(f"Available columns: {list(lemonnetwork.columns)}")
        return
    
    print("\nExtracting regulator-gene interactions by type:")
    for reg_type, patterns in regulator_type_patterns.items():
        # Filter for this regulator type
        type_network = lemonnetwork[lemonnetwork['Type'].isin(patterns)]
        print(f"  {reg_type}: Found {len(type_network)} interactions")
        
        if len(type_network) == 0:
            regulators2genes[reg_type] = {}
            continue
        
        # Build regulator->genes mapping for this type
        reg_dict = {}
        for index, row in type_network.iterrows():
            regulator = row['Regulator']
            
            try:
                if has_targets_column and pd.notna(row['Targets']):
                    targets = row['Targets'].split('|')
                elif has_target_column and pd.notna(row['Target']):
                    targets = [row['Target']]
                else:
                    targets = []
            except KeyError as e:
                print(f"    KeyError accessing target column: {e}")
                targets = []
            
            if regulator not in reg_dict:
                reg_dict[regulator] = targets
            else:
                reg_dict[regulator].extend(targets)
            
            all_genes.update(targets)
        
        regulators2genes[reg_type] = reg_dict
    
    # For backward compatibility, extract metabolites2genes if it exists
    metabolites2genes = regulators2genes.get('Metabolites', {})
    
    genes_in_dataset = list(all_genes)
    print(f'There are {len(genes_in_dataset)} genes in the dataset')
    
    # Load regulators per module dynamically
    print("Loading module regulators...")
    module_regulators = {}
    for reg_type, reg_file in regulator_files.items():
        if os.path.exists(reg_file):
            module_regulators[reg_type] = get_regulators(reg_file)
            print(f'Loaded {reg_type} regulators for {len(module_regulators[reg_type])} modules')
        else:
            print(f'Warning: {reg_type} regulator file {reg_file} not found')
            module_regulators[reg_type] = {}
    
    module2genes = get_regulators(clusterfile)
    print(f'Loaded {len(module2genes)} modules')
    
    # For backward compatibility, create the old variable names
    # Try both singular and plural forms to match configuration
    module2TFs = module_regulators.get('TFs', module_regulators.get('TF', {}))
    module2mets = module_regulators.get('Metabolites', module_regulators.get('Metabolite', {}))
    module2lipids = module_regulators.get('Lipids', module_regulators.get('Lipid', {}))
    module2proteins = module_regulators.get('Proteins', module_regulators.get('Protein', {}))
    
    # Load metabolite-gene interactions from PKN
    print("Loading PKN interactions...")
    if not os.path.exists(metabolite_interactions_file):
        print(f"Warning: PKN interactions file {metabolite_interactions_file} not found")
        interactions = pd.DataFrame()
    else:
        interactions = pd.read_csv(metabolite_interactions_file, sep='\t')
        print(f"Loaded {len(interactions)} PKN interactions")
        print(f"Metabolite interactions file columns: {list(interactions.columns)}")
        
        # Handle different file formats
        if 'HMDB' not in interactions.columns and 'Metabolite' in interactions.columns:
            # Format: Metabolite | Gene | Source
            # Convert to: HMDB | All_interactions
            print("Converting metabolite interactions format from Metabolite-Gene to HMDB-based format...")
            
            # Extract HMDB ID from metabolite name (e.g., "1-Methylhistidine_HMDB0000001" -> "HMDB0000001")
            interactions['HMDB'] = interactions['Metabolite'].str.extract(r'(HMDB\d+)')
            
            # Group by HMDB (or Metabolite if extraction failed) and collect all genes
            grouped = interactions.groupby('HMDB')['Gene'].apply(lambda x: '|'.join(x)).reset_index()
            grouped.columns = ['HMDB', 'All_interactions']
            interactions = grouped
            print(f"Converted to {len(interactions)} unique HMDB entries with gene interactions")
        elif 'HMDB' not in interactions.columns and 'All_interactions' not in interactions.columns:
            print(f"Warning: Unexpected column format. Expected 'HMDB' or 'Metabolite' columns. Found: {list(interactions.columns)}")
    
    # Load metabolite mapping
    metabolite_mapping = load_metabolite_mapping(annotated_mets)
    print(f"Loaded {len(metabolite_mapping)} metabolite mappings")
    
    # Load PKN network as NetworkX graph
    print("Loading PKN network...")
    if not os.path.exists(PKNnetwork_file):
        print(f"Error: PKN network file {PKNnetwork_file} not found")
        PKN_graph = nx.Graph()
    else:
        try:
            # Try to load as tab-separated file first (common format)
            PKN_df = pd.read_csv(PKNnetwork_file, sep='\t')
            print(f"Loaded PKN network with {len(PKN_df)} interactions")

            # Create NetworkX graph from the interactions
            PKN_graph = nx.Graph()
            for _, row in PKN_df.iterrows():
                if len(row) >= 2:
                    source = str(row.iloc[0]).strip()
                    target = str(row.iloc[1]).strip()

                    # Parse metabolite names to extract HMDB IDs
                    if 'HMDB' in source:
                        # Extract HMDB ID from metabolite name (format: Name_HMDBXXXXXXX)
                        source_parts = source.split('_')
                        for part in reversed(source_parts):
                            if part.startswith('HMDB'):
                                source = part
                                break

                    if 'HMDB' in target:
                        # Extract HMDB ID from metabolite name (format: Name_HMDBXXXXXXX)
                        target_parts = target.split('_')
                        for part in reversed(target_parts):
                            if part.startswith('HMDB'):
                                target = part
                                break

                    PKN_graph.add_edge(source, target)

            print(f"Created PKN graph with {PKN_graph.number_of_nodes()} nodes and {PKN_graph.number_of_edges()} edges")
            
        except Exception as e:
            print(f"Error loading PKN network: {e}")
            print("Creating empty PKN graph")
            PKN_graph = nx.Graph()
    
    if len(interactions) > 0 and len(metabolite_mapping) > 0:
        # Contextualize metabolite-gene interactions to dataset
        metabolites_hmdb = list(metabolite_mapping.values())
        print(f'{len(metabolites_hmdb)} metabolites in this dataset are mapped to a HMDB id')
        
        # Check if 'HMDB' column exists in interactions DataFrame
        if 'HMDB' in interactions.columns:
            # How many metabolites_hmdb are in the interactions file?
            metabolites_in_interactions = interactions['HMDB'].isin(metabolites_hmdb)
            print(f'{metabolites_in_interactions.sum()} metabolites in the interactions file are in the dataset')
            
            # How many metabolites in the dataset are in the interactions file?
            metabolites_in_dataset = pd.Series(list(metabolites_hmdb)).isin(interactions['HMDB'])
            print(f'{metabolites_in_dataset.sum()} metabolites in the dataset are in the interactions file')
        else:
            print(f"Warning: 'HMDB' column not found in interactions file. Available columns: {list(interactions.columns)}")
            print("Skipping metabolite interaction contextualization check")
    
    # Create output directory in ./ModuleViewer_files
    moduleviewer_dir = os.path.join(workdir, 'ModuleViewer_files')
    os.makedirs(moduleviewer_dir, exist_ok=True)
    
    # Generate ModuleViewer file for metabolite-KG interactions
    mvf_file = os.path.join(moduleviewer_dir, 'metabolite_LemoniteKG_interactions.mvf')
    
    print(f"Generating ModuleViewer file: {mvf_file}")
    interactions_found = 0
    modules_with_metabolites = 0
    total_metabolites_checked = 0
    
    with open(mvf_file, 'w') as handle:
        # Write header lines
        handle.write('::TYPE=Lemonite_KG\n')
        handle.write('::TITLE:Lemonite_KG\n')
        handle.write('::OBJECT=GENES\n')
        handle.write('::COLOR=YELLOW\n')
        
        print(f"Checking {len(module2mets)} modules with metabolite regulators and {len(module2lipids)} modules with lipid regulators...")
        
        # Loop over every module, and for every metabolite in the module, 
        # check if it has interactions with genes in the module
        for module in module2mets.keys():
            if module not in module2genes:
                print(f"  Module {module}: No genes found, skipping")
                continue
                
            module_genes = module2genes[module]
            module_mets = module2mets[module]
            modules_with_metabolites += 1
            
            print(f"  Module {module}: {len(module_mets)} metabolites, {len(module_genes)} genes")
            
            for met in module_mets:
                total_metabolites_checked += 1
                met_hmdb = metabolite_mapping.get(met)
                if met_hmdb is None:
                    print(f"    Metabolite {met}: No HMDB mapping found")
                    continue
                    
                if len(interactions) == 0:
                    print(f"    Metabolite {met}: No PKN interactions available")
                    continue
                
                # Check if 'HMDB' column exists in interactions
                if 'HMDB' not in interactions.columns:
                    print(f"    Metabolite {met}: 'HMDB' column not found in PKN interactions file (columns: {list(interactions.columns)})")
                    continue
                    
                met_interactions = interactions[interactions['HMDB'] == met_hmdb]
                if met_interactions.empty:
                    print(f"    Metabolite {met} ({met_hmdb}): No PKN interactions found")
                    continue
                    
                # Get gene interactions for this metabolite
                if 'All_interactions' in met_interactions.columns:
                    interaction_col = met_interactions['All_interactions'].values[0]
                    if pd.isna(interaction_col):
                        print(f"    Metabolite {met} ({met_hmdb}): Empty interactions column")
                        continue
                    met_genes = interaction_col.split('|')
                else:
                    print(f"    Metabolite {met} ({met_hmdb}): No 'All_interactions' column in PKN")
                    continue
                
                shared_genes = set(module_genes).intersection(set(met_genes))
                print(f'    {met} ({met_hmdb}) has {len(shared_genes)} interactions with genes in module {module}')
                
                if len(shared_genes) > 0:
                    # write genes in a deterministic order
                    handle.write(str(int(module)) + '\t' + '|'.join(sorted(shared_genes)) + '\t' + str(met) + '\n')
                    interactions_found += 1
        
        # Process lipid regulators
        for module in module2lipids.keys():
            if module not in module2genes:
                print(f"  Module {module}: No genes found, skipping lipids")
                continue
                
            module_genes = module2genes[module]
            module_lipids = module2lipids[module]
            
            print(f"  Module {module}: {len(module_lipids)} lipids, {len(module_genes)} genes")
            
            for lipid in module_lipids:
                # For lipids, we assume they are already in the correct format for PKN
                # Check if lipid has interactions with genes in the module
                if lipid not in PKN_graph.nodes():
                    print(f"    Lipid {lipid}: Not found in PKN")
                    continue
                    
                # Get neighbors of the lipid in PKN
                lipid_neighbors = set(PKN_graph.neighbors(lipid))
                shared_genes = set(module_genes).intersection(lipid_neighbors)
                
                print(f'    {lipid} has {len(shared_genes)} interactions with genes in module {module}')
                
                if len(shared_genes) > 0:
                    # write genes in a deterministic order
                    handle.write(str(int(module)) + '\t' + '|'.join(sorted(shared_genes)) + '\t' + str(lipid) + '\n')
                    interactions_found += 1
        
        # If no interactions found, write a comment line to indicate the file was processed
        if interactions_found == 0:
            handle.write('# No metabolite-gene interactions found that are supported by PKN\n')
    
    print(f"ModuleViewer file generation summary:")
    print(f"  - File created: {mvf_file}")
    print(f"  - Modules with metabolites checked: {modules_with_metabolites}")
    print(f"  - Modules with lipids checked: {len(module2lipids)}")
    print(f"  - Total metabolites checked: {total_metabolites_checked}")
    print(f"  - PKN-supported interactions found: {interactions_found}")
    
    # ========================================
    # Generate PPI interactions MVF file
    # ========================================
    print("\n" + "=" * 60)
    print("Generating PPI interactions MVF file for ModuleViewer...")
    print("=" * 60)
    
    ppi_mvf_file = os.path.join(moduleviewer_dir, 'PPI_interactions.mvf')
    
    # Pre-process: Get all genes across all modules
    all_module_genes = set()
    for genes_list in module2genes.values():
        all_module_genes.update(genes_list)
    
    print(f"Total unique genes across all modules: {len(all_module_genes)}")
    
    # Filter PKN to only include genes that are in modules (both nodes must be in module genes)
    # PKN_graph already contains the network, so we can use it directly
    ppi_set = set()
    ppi_count_in_pkn = 0
    
    if PKN_graph.number_of_nodes() > 0:
        # Create a set of PPI pairs for O(1) lookup (store both directions)
        for edge in PKN_graph.edges():
            node1, node2 = edge[0], edge[1]
            # Only include if both nodes are genes in modules (not metabolites/lipids)
            # Skip edges containing HMDB IDs (metabolites)
            if ('HMDB' not in str(node1) and 'HMDB' not in str(node2) and 
                node1 in all_module_genes and node2 in all_module_genes):
                # Add both directions (undirected graph)
                ppi_set.add((node1, node2))
                ppi_set.add((node2, node1))
                ppi_count_in_pkn += 1
        
        print(f"Total PPI pairs in PKN (filtered for module genes): {ppi_count_in_pkn}")
        print(f"Total PPI pairs (bidirectional): {len(ppi_set)}")
    else:
        print("Warning: No PKN graph loaded, PPI file will be empty")
    
    # Now create the MVF file
    total_ppis_written = 0
    modules_with_ppis = 0
    
    with open(ppi_mvf_file, 'w') as handle:
        # Write header lines
        handle.write('::TYPE=PPI\n')
        handle.write('::TITLE=Protein-Protein Interactions\n')
        handle.write('::OBJECT=GENE_PAIRS\n')
        handle.write('::COLOR=BLUE\n')
        
        # Loop over every module and find PPIs between genes in the module
        for module in module2genes.keys():
            module_genes = module2genes[module]
            module_ppis = []
            
            # Check all pairs of genes in the module for PPIs in the PKN
            for i, gene1 in enumerate(module_genes):
                for gene2 in module_genes[i+1:]:  # Only check each pair once
                    # Fast O(1) lookup in set
                    if (gene1, gene2) in ppi_set:
                        module_ppis.append(f"{gene1}|{gene2}")
            
            # Write PPIs for this module
            for ppi in module_ppis:
                handle.write(f"{int(module)}\t{ppi}\n")
                total_ppis_written += 1
            
            if len(module_ppis) > 0:
                modules_with_ppis += 1
                print(f"Module {module}: Found {len(module_ppis)} PPIs between module genes")
        
        # If no PPIs found, write a comment line to indicate the file was processed
        if total_ppis_written == 0:
            handle.write('# No protein-protein interactions found in PKN for module genes\n')
    
    print(f"\nPPI MVF file generation summary:")
    print(f"  - File created: {ppi_mvf_file}")
    print(f"  - Total PPIs written: {total_ppis_written}")
    print(f"  - Modules with PPIs: {modules_with_ppis}/{len(module2genes)}")
    
    # Verify the PPI file was created
    if not os.path.exists(ppi_mvf_file):
        print(f"Warning: {ppi_mvf_file} was not created")
    elif os.path.getsize(ppi_mvf_file) == 0:
        print(f"Warning: {ppi_mvf_file} exists but is empty")
    else:
        print(f"PPI MVF file successfully created with size: {os.path.getsize(ppi_mvf_file)} bytes")
    
    # ========================================
    # Calculate PPI enrichment for modules
    # ========================================
    ppi_enrichment_df = calculate_ppi_enrichment(module2genes, PKN_graph, all_module_genes)
    
    # Save PPI enrichment results
    if len(ppi_enrichment_df) > 0:
        ppi_enrichment_file = os.path.join(moduleviewer_dir, 'PPI_enrichment_results.csv')
        ppi_enrichment_df.to_csv(ppi_enrichment_file, index=False)
        print(f"\nPPI enrichment results saved to: {ppi_enrichment_file}")
        print(f"  - Total modules analyzed: {len(ppi_enrichment_df)}")
        n_significant = (ppi_enrichment_df['FDR'] < 0.05).sum()
        print(f"  - Significantly enriched modules (FDR < 0.05): {n_significant}")
    else:
        print("\nNo PPI enrichment results to save (no modules with >= 2 genes)")
    
    # ========================================
    # Calculate metabolite-gene interaction enrichment for modules
    # ========================================
    metgene_enrichment_df = calculate_metabolite_gene_enrichment(
        module2genes, module2mets, module2lipids, interactions,
        metabolite_mapping, PKN_graph, all_module_genes
    )
    
    # Save metabolite-gene enrichment results
    if len(metgene_enrichment_df) > 0:
        metgene_enrichment_file = os.path.join(moduleviewer_dir, 'Metabolite_Gene_enrichment_results.csv')
        metgene_enrichment_df.to_csv(metgene_enrichment_file, index=False)
        print(f"\nMetabolite-gene enrichment results saved to: {metgene_enrichment_file}")
        print(f"  - Total modules analyzed: {len(metgene_enrichment_df)}")
        n_significant = (metgene_enrichment_df['FDR'] < 0.05).sum()
        print(f"  - Significantly enriched modules (FDR < 0.05): {n_significant}")
    else:
        print("\nNo metabolite-gene enrichment results to save")
    
    # Ensure sample_mapping.mvf is available in ModuleViewer_files for downstream viewer
    def find_and_copy_sample_mapping(dest_dir):
        # common locations to search: current directory (staged files), work/*/results/ModuleViewer_files and results/ModuleViewer_files
        candidates = glob.glob(os.path.join(workdir, 'sample_mapping.mvf'))
        candidates += glob.glob(os.path.join(workdir, 'work', '*', '*', 'results', 'ModuleViewer_files', 'sample_mapping.mvf'))
        candidates += glob.glob(os.path.join(workdir, 'results', 'ModuleViewer_files', 'sample_mapping.mvf'))
        candidates += glob.glob(os.path.join(workdir, 'work', '*', '*', 'ModuleViewer_files', 'sample_mapping.mvf'))
        candidates += glob.glob(os.path.join(workdir, 'ModuleViewer_files', 'sample_mapping.mvf'))
        for src in candidates:
            try:
                shutil.copy(src, dest_dir)
                print(f"Copied sample_mapping.mvf from {src} to {dest_dir}")
                return True
            except Exception as e:
                print(f"Failed to copy sample_mapping from {src}: {e}")
        print("No sample_mapping.mvf found in common locations")
        return False

    copied = find_and_copy_sample_mapping(moduleviewer_dir)
    # verify both MVF files exist and are non-empty
    sample_map_file = os.path.join(moduleviewer_dir, 'sample_mapping.mvf')
    if not os.path.exists(sample_map_file):
        print(f"Warning: {sample_map_file} not found after copy attempt")
    else:
        if os.path.getsize(sample_map_file) == 0:
            print(f"Warning: {sample_map_file} exists but is empty")

    if not os.path.exists(mvf_file) or os.path.getsize(mvf_file) == 0:
        print(f"Warning: {mvf_file} missing or empty")
    
    # Generate subnetwork visualizations for each module
    print("\nGenerating subnetwork visualizations...")
    subnetworks_dir = os.path.join(workdir, 'Networks', 'subnetworks')
    os.makedirs(subnetworks_dir, exist_ok=True)
    
    modules_processed = 0
    modules_with_subnetworks = 0
    
    for module in module2genes.keys():
        try:
            if module not in module2TFs and module not in module2mets and module not in module2lipids:
                print(f'Module {module}: No regulators found, skipping subnetwork visualization')
                continue
            
            target_genes = module2genes[module]
            
            # Create regulators dictionary for this module
            regulators_dict = {}
            for reg_type, module_regs in module_regulators.items():
                regulators_dict[reg_type] = module_regs.get(module, [])
            
            # Check if any regulators exist for this module
            has_regulators = any(regulators_dict.values())
            if not has_regulators:
                print(f'Module {module}: No regulators found, skipping subnetwork visualization')
                continue
            
            # Draw subnetwork for this module (DISABLED - using categorized graphs from SUBNETWORK_GRAPHS process instead)
            # draw_subnetwork(module, target_genes, regulators_dict, PKN_graph, metabolite_mapping)
            
            modules_processed += 1
            # modules_with_subnetworks += 1  # Commented out since we're not creating subnetworks here anymore
            
        except KeyError:
            print(f'Module {module}: Key error during subnetwork visualization')
            modules_processed += 1
            continue
        except Exception as e:
            print(f'Error during subnetwork visualization for module {module}: {e}')
            modules_processed += 1
            continue
    
    print(f"Subnetwork visualization summary:")
    print(f"  - Modules processed: {modules_processed}")
    print(f"  - Subnetworks generated: {modules_with_subnetworks}")
    print(f"  - Output directory: {subnetworks_dir}")
    
    # Create evaluation summary
    summary_file = os.path.join(moduleviewer_dir, 'evaluation_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Lemonite Network PKN Evaluation Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Network file: {lemonnetwork_file}\n")
        f.write(f"Total network interactions: {len(lemonnetwork)}\n")
        # Write stats for all regulator types
        for reg_type in regulator_files.keys():
            if reg_type in regulators2genes:
                num_interactions = sum(len(targets) for targets in regulators2genes[reg_type].values())
                f.write(f"{reg_type}-gene interactions: {num_interactions}\n")
        f.write(f"Genes in dataset: {len(genes_in_dataset)}\n")
        f.write(f"Modules analyzed: {len(module2genes)}\n")
        f.write(f"Modules with TF regulators: {len(module2TFs)}\n")
        f.write(f"Modules with metabolite regulators: {len(module2mets)}\n")
        f.write(f"Modules with lipid regulators: {len(module2lipids)}\n")
        if 'Proteins' in regulator_files:
            f.write(f"Modules with protein regulators: {len(module2proteins)}\n")
        f.write(f"Modules with metabolites checked: {modules_with_metabolites}\n")
        f.write(f"Total metabolites checked: {total_metabolites_checked}\n")
        f.write(f"PKN-supported interactions: {interactions_found}\n")
        f.write(f"ModuleViewer file created: {mvf_file}\n")
        f.write(f"PPI MVF file created: {ppi_mvf_file}\n")
        f.write(f"Total PPIs in modules: {total_ppis_written}\n")
        f.write(f"Modules with PPIs: {modules_with_ppis}/{len(module2genes)}\n")
        if len(ppi_enrichment_df) > 0:
            n_ppi_sig = int((ppi_enrichment_df['FDR'] < 0.05).sum())
            f.write(f"PPI enrichment: {n_ppi_sig} modules significantly enriched (FDR < 0.05)\n")
        if len(metgene_enrichment_df) > 0:
            n_metgene_sig = int((metgene_enrichment_df['FDR'] < 0.05).sum())
            f.write(f"Metabolite-gene enrichment: {n_metgene_sig} modules significantly enriched (FDR < 0.05)\n")
        f.write(f"Subnetworks generated: {modules_with_subnetworks}/{modules_processed}\n")
        f.write(f"Subnetworks directory: {subnetworks_dir}\n")
        if len(metabolite_mapping) > 0:
            f.write(f"Metabolites with HMDB mapping: {len(metabolite_mapping)}\n")
        if PKN_graph.number_of_nodes() > 0:
            f.write(f"PKN graph nodes: {PKN_graph.number_of_nodes()}\n")
            f.write(f"PKN graph edges: {PKN_graph.number_of_edges()}\n")
        else:
            f.write("PKN graph: None loaded\n")
        if len(interactions) > 0:
            f.write(f"PKN interactions loaded: {len(interactions)}\n")
        else:
            f.write("PKN interactions: None loaded\n")
    
    print(f"Evaluation summary written to: {summary_file}")
    print("PKN evaluation completed successfully!")

def main():
    parser = argparse.ArgumentParser(description='Evaluate Lemonite network against PKN')
    
    # Required arguments
    parser.add_argument('--workdir', required=True, help='Working directory')
    parser.add_argument('--network_file', required=True, help='LemonTree network file')
    parser.add_argument('--pkn_file', required=True, help='PKN network file')
    parser.add_argument('--metabolite_interactions_file', required=True, help='Metabolite-gene interactions file')
    parser.add_argument('--metabolite_mapping_file', required=True, help='Metabolite mapping file (Query to HMDB)')
    parser.add_argument('--regulator_files', required=True, 
                        help='Regulator files in format "Type:Path,Type:Path" (e.g., "TF:/path/to/tf.txt,Metabolite:/path/to/met.txt,Lipid:/path/to/lipids.txt")')
    parser.add_argument('--cluster_file', required=True, help='Cluster file')
    
    # Optional arguments
    parser.add_argument('--top_n_percent_regulators', type=float, default=2.0, help='Top N percent of regulators to select')
    parser.add_argument('--output_dir', help='Output directory (default: workdir/PKN_Evaluation)')
    parser.add_argument('--cores', type=int, default=1, help='Number of CPU cores to use')
    
    args = parser.parse_args()
    
    # Set number of cores for potential future parallelization
    if args.cores > 1:
        print(f"PKN evaluation configured to use {args.cores} cores (multiprocessing ready)")
    
    try:
        evaluate_network_against_pkn(args)
    except Exception as e:
        print(f"Error during PKN evaluation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
