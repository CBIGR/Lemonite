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
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')  # headless backend required in Singularity/HPC environments
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.lines import Line2D
import seaborn as sns
import argparse
import shutil
import glob

try:
    import pygraphviz as pgv
    PYGRAPHVIZ_AVAILABLE = True
except ImportError:
    PYGRAPHVIZ_AVAILABLE = False
    print('Warning: pygraphviz not available; subnetwork plots will be skipped.')

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

# ---------------------------------------------------------------------------
# Edge-styling helpers (match Wang_GBM/Lemonite/Investigate_shortest_paths.ipynb)
# ---------------------------------------------------------------------------

def get_edge_type_category(source, edge_type):
    """Categorize a PKN edge into Causal / Metabolic_pathway / Ambiguous / PPI."""
    if edge_type == 'metabolite-gene':
        if source in ('LINCS', 'chEMBL'):
            return 'Causal'
        elif source in ('Human1_GEM_dist1', 'Human1_GEM_dist2'):
            return 'Metabolic_pathway'
        else:
            return 'Ambiguous'
    return 'PPI'


def get_edge_style_mapping():
    """Visual styles for each interaction category."""
    return {
        'Causal': {
            'color': '#555555', 'style': 'solid', 'width': 3.0, 'alpha': 0.9,
            'label': 'Causal (LINCS/chEMBL)',
        },
        'Metabolic_pathway': {
            'color': '#555555', 'style': 'solid', 'width': 2.5, 'alpha': 0.85,
            'label': 'Metabolic pathway (GEM)',
        },
        'Ambiguous': {
            'color': '#555555', 'style': 'solid', 'width': 2.0, 'alpha': 0.7,
            'label': 'Ambiguous (other)',
        },
        'PPI': {
            'color': '#A9A9A9', 'style': 'solid', 'width': 1.5, 'alpha': 0.5,
            'label': 'PPI',
        },
        'HumanNet': {
            'color': 'orange', 'style': 'solid', 'width': 1.5, 'alpha': 0.6,
            'label': 'HumanNet (functional)',
        },
    }


def get_edge_source_and_type(PKN_df, node1, node2):
    """Return (source, edge_type, category) for an edge, using PKN_df if available."""
    if PKN_df is None or PKN_df.empty:
        return 'unknown', 'unknown', 'Ambiguous'

    # Detect column names: prefer Node1/Node2, fall back to first two columns
    cols = list(PKN_df.columns)
    if 'Node1' in cols and 'Node2' in cols:
        c1, c2 = 'Node1', 'Node2'
    else:
        c1, c2 = cols[0], cols[1]

    mask = (
        ((PKN_df[c1] == node1) & (PKN_df[c2] == node2)) |
        ((PKN_df[c1] == node2) & (PKN_df[c2] == node1))
    )
    match = PKN_df[mask]
    if match.empty:
        return 'unknown', 'unknown', 'Ambiguous'

    source = match.iloc[0]['Source'] if 'Source' in cols else 'unknown'
    edge_type = match.iloc[0]['Type'] if 'Type' in cols else 'unknown'
    category = get_edge_type_category(source, edge_type)
    return source, edge_type, category


def create_edge_legend_by_category(edge_styles, edge_info):
    """Return legend Line2D handles for edge categories present in the network."""
    present = {info.get('category', 'PPI') for info in edge_info.values()}
    handles = []
    for category in sorted(edge_styles.keys()):
        style = edge_styles[category]
        if category == 'Metabolic_pathway':
            handles.append(Line2D([0, 1], [0, 0],
                marker='o', markersize=6, markerfacecolor=style['color'],
                markeredgecolor='black', markeredgewidth=0.5,
                color=style['color'], linewidth=style['width'],
                linestyle=style['style'], alpha=style['alpha'],
                label=style['label']))
        elif category == 'Causal':
            handles.append(Line2D([0, 1], [0, 0],
                marker='>', markersize=8, markerfacecolor=style['color'],
                markeredgecolor='black', markeredgewidth=0.5,
                color=style['color'], linewidth=style['width'],
                linestyle=style['style'], alpha=style['alpha'],
                label=style['label']))
        else:
            handles.append(Line2D([0, 1], [0, 0],
                color=style['color'], linewidth=style['width'],
                linestyle=style['style'], alpha=style['alpha'],
                label=style['label']))
    return handles


def create_node_legend(has_metabolites, has_tfs, has_targets, has_bridge, has_lipids=False):
    """Return Patch legend handles for node types present in the network."""
    handles = []
    if has_metabolites:
        handles.append(Patch(facecolor='red', edgecolor='black', label='Metabolites'))
    if has_lipids:
        handles.append(Patch(facecolor='purple', edgecolor='black', label='Lipids'))
    if has_tfs:
        handles.append(Patch(facecolor='lightgreen', edgecolor='black', label='TFs'))
    if has_targets:
        handles.append(Patch(facecolor='orange', edgecolor='black', label='Targets'))
    if has_bridge:
        handles.append(Patch(facecolor='lightgrey', edgecolor='black', label='Bridge'))
    return handles


def print_edge_statistics(edge_info):
    """Print counts of edge categories, types, and sources."""
    from collections import Counter
    categories = [i.get('category', 'PPI') for i in edge_info.values()]
    types      = [i.get('type', 'unknown') for i in edge_info.values()]
    sources    = [i.get('source', 'unknown') for i in edge_info.values()]
    print("\n=== Edge Statistics ===")
    print("By Category:", dict(Counter(categories)))
    print("By Type:",     dict(Counter(types)))
    print("By Source:",   dict(Counter(sources)))


# ---------------------------------------------------------------------------
# Main subnetwork drawing function
# ---------------------------------------------------------------------------

def draw_subnetwork(module, target_genes, regulators_dict, PKN_graph, name_to_hmdb,
                    PKN_df=None, humannet_set=None, layout_engine='neato'):
    """
    Draw subnetwork visualization for a module using Graphviz/pygraphviz.
    Visual style matches Wang_GBM/Lemonite/Investigate_shortest_paths.ipynb.

    layout_engine options: 'neato' (force-directed), 'dot' (hierarchical), 'fdp', 'sfdp'
    """
    if not PYGRAPHVIZ_AVAILABLE:
        print(f'Skipping subnetwork plot for module {module}: pygraphviz not available.')
        return

    hmdb_to_name = {v: k for k, v in name_to_hmdb.items()}

    TF_regulators         = regulators_dict.get('TF', [])
    metabolite_regulators = regulators_dict.get('Metabolite', [])
    lipid_regulators      = regulators_dict.get('Lipid', [])

    print(f'Module {module}: {len(metabolite_regulators)} metabolites, '
          f'{len(TF_regulators)} TFs, {len(lipid_regulators)} lipids, '
          f'{len(target_genes)} targets.')

    to_draw   = nx.Graph()
    edge_info = {}
    met_regs  = []

    # --- STEP 1: Metabolite regulators (direct edges to target genes) ---
    for reg in metabolite_regulators:
        if reg not in name_to_hmdb:
            continue
        regulator = name_to_hmdb[reg]
        met_regs.append(regulator)
        if regulator not in PKN_graph.nodes():
            continue
        to_draw.add_node(regulator)
        for target in target_genes:
            try:
                if nx.shortest_path_length(PKN_graph, source=regulator, target=target) == 1:
                    to_draw.add_nodes_from([regulator, target])
                    to_draw.add_edge(regulator, target)
                    src, etype, cat = get_edge_source_and_type(PKN_df, regulator, target)
                    edge_info[(regulator, target)] = {'source': src, 'type': etype, 'category': cat}
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                pass

    # --- STEP 1.5: Lipid regulators (direct edges to target genes) ---
    for lipid in lipid_regulators:
        if lipid not in PKN_graph.nodes():
            continue
        to_draw.add_node(lipid)
        for target in target_genes:
            try:
                if nx.shortest_path_length(PKN_graph, source=lipid, target=target) == 1:
                    to_draw.add_nodes_from([lipid, target])
                    to_draw.add_edge(lipid, target)
                    src, etype, cat = get_edge_source_and_type(PKN_df, lipid, target)
                    edge_info[(lipid, target)] = {'source': src, 'type': etype, 'category': cat}
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                pass

    # --- STEP 2: TF regulators (via PPI path to metabolite, ≤2 steps) ---
    for TF in TF_regulators:
        for regulator in met_regs:
            try:
                if nx.shortest_path_length(PKN_graph, source=TF, target=regulator) <= 2:
                    path = nx.shortest_path(PKN_graph, source=TF, target=regulator)
                    for i in range(len(path) - 1):
                        u, v = path[i], path[i + 1]
                        to_draw.add_edge(u, v)
                        src, etype, cat = get_edge_source_and_type(PKN_df, u, v)
                        edge_info[(u, v)] = {'source': src, 'type': etype, 'category': cat}
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                pass

    # --- STEP 3: Additional PKN edges among already-included nodes ---
    subgraph = PKN_graph.subgraph(to_draw.nodes())
    for u, v in subgraph.edges():
        if not to_draw.has_edge(u, v):
            to_draw.add_edge(u, v)
            src, etype, cat = get_edge_source_and_type(PKN_df, u, v)
            edge_info[(u, v)] = {'source': src, 'type': etype, 'category': cat}

    # --- STEP 4: HumanNet edges among included nodes ---
    if humannet_set is not None:
        nodes = list(to_draw.nodes())
        for i, u in enumerate(nodes):
            for v in nodes[i + 1:]:
                if not to_draw.has_edge(u, v) and (
                        (u, v) in humannet_set or (v, u) in humannet_set):
                    to_draw.add_edge(u, v)
                    edge_info[(u, v)] = {'source': 'HumanNet', 'type': 'functional',
                                         'category': 'HumanNet'}

    to_draw.remove_edges_from(nx.selfloop_edges(to_draw))
    to_draw.remove_nodes_from([n for n in list(to_draw) if to_draw.degree(n) == 0])

    if not to_draw.nodes:
        print(f'No valid connections found for module {module}')
        return

    # === Node categories ===
    hmdb_nodes   = set(n for n in to_draw if n in name_to_hmdb.values())
    lipid_nodes  = set(n for n in to_draw if n in lipid_regulators)
    tf_nodes     = set(n for n in to_draw if n in TF_regulators)
    target_nodes = set(n for n in to_draw if n in target_genes)
    bridge_nodes = set(n for n in to_draw
                       if n not in hmdb_nodes | lipid_nodes | tf_nodes | target_nodes)

    # === Build Graphviz graph ===
    G = pgv.AGraph(directed=False, strict=False)
    G.graph_attr['overlap'] = 'false'
    G.graph_attr['splines'] = 'curved'
    G.graph_attr['sep']     = '+0.5'
    G.graph_attr['pad']     = '0.5'

    for node in to_draw.nodes():
        label = hmdb_to_name.get(node, node)
        if node in hmdb_nodes:
            fillcolor, edgecolor = 'red', 'darkred'
        elif node in lipid_nodes:
            fillcolor, edgecolor = 'mediumpurple', 'purple'
        elif node in tf_nodes:
            fillcolor, edgecolor = 'lightgreen', 'darkgreen'
        elif node in target_nodes:
            fillcolor, edgecolor = 'orange', 'darkorange'
        else:
            fillcolor, edgecolor = 'lightgrey', 'darkgrey'
        G.add_node(node, label=label, shape='box', style='rounded,filled',
                   fillcolor=fillcolor, color=edgecolor,
                   fontname='Arial', fontsize='10', penwidth='2.0',
                   height='0.55', width=str(0.6 + len(label) * 0.07))

    for u, v in to_draw.edges():
        info     = edge_info.get((u, v)) or edge_info.get((v, u), {})
        category = info.get('category', 'PPI')
        if category == 'Causal':
            G.add_edge(u, v, color='#555555', penwidth='3.0',
                       arrowhead='vee', arrowtail='none', dir='forward', arrowsize='1.2')
        elif category == 'Metabolic_pathway':
            G.add_edge(u, v, color='#555555', penwidth='2.5',
                       arrowhead='dot', arrowtail='none', dir='forward', arrowsize='1.8')
        elif category == 'Ambiguous':
            G.add_edge(u, v, color='#555555', penwidth='2.0',
                       arrowhead='none', arrowtail='none', dir='none')
        elif category == 'HumanNet':
            G.add_edge(u, v, color='orange', penwidth='1.5',
                       arrowhead='none', arrowtail='none', dir='none')
        else:  # PPI
            G.add_edge(u, v, color='#A9A9A9', penwidth='1.5',
                       arrowhead='none', arrowtail='none', dir='none')

    # === Layout and save ===
    G.layout(prog=layout_engine)
    output_dir = os.path.join(os.getcwd(), 'Networks', 'subnetworks')
    os.makedirs(output_dir, exist_ok=True)
    out_png = os.path.join(output_dir, f'graph_{module}.png')
    G.draw(out_png, prog=layout_engine, format='png', args='-Gdpi=200')

    print_edge_statistics(edge_info)
    print(f'Saved subnetwork: {out_png} '
          f'({to_draw.number_of_nodes()} nodes, {to_draw.number_of_edges()} edges)')

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
    _, fdr_values, _, _ = multipletests(enrichment_df['P_value'].values, method='fdr_bh')
    enrichment_df['FDR'] = fdr_values
    
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
    _, fdr_values, _, _ = multipletests(enrichment_df['P_value'].values, method='fdr_bh')
    enrichment_df['FDR'] = fdr_values
    
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
    # Generate HumanNet interactions MVF file
    # ========================================
    print("\n" + "=" * 60)
    print("Generating HumanNet interactions MVF file for ModuleViewer...")
    print("=" * 60)

    humannet_mvf_file = os.path.join(moduleviewer_dir, 'HumanNet_interactions.mvf')

    # Locate HumanNet CSV file
    humannet_csv = None
    for candidate in [
        os.path.join(os.path.dirname(args.pkn_file), 'HumanNet_interactions.csv'),
        '/opt/PKN/HumanNet_interactions.csv',
    ]:
        if os.path.exists(candidate):
            humannet_csv = candidate
            print(f"Found HumanNet file: {humannet_csv}")
            break

    if humannet_csv is None:
        print("Warning: HumanNet_interactions.csv not found - HumanNet MVF will be empty")

    humannet_set = set()
    if humannet_csv is not None:
        hn_df = pd.read_csv(humannet_csv, header=0)
        for _, row in hn_df.iterrows():
            g1, g2 = str(row.iloc[0]), str(row.iloc[1])
            if g1 in all_module_genes and g2 in all_module_genes:
                humannet_set.add((g1, g2))
                humannet_set.add((g2, g1))
        print(f"HumanNet pairs filtered for module genes (bidirectional): {len(humannet_set)}")

    total_hn_written = 0
    modules_with_hn = 0

    with open(humannet_mvf_file, 'w') as handle:
        handle.write('::TYPE=HumanNet\n')
        handle.write('::TITLE=HumanNet Functional Interactions\n')
        handle.write('::OBJECT=GENE_PAIRS\n')
        handle.write('::COLOR=saddlebrown\n')

        for module in module2genes.keys():
            module_genes = module2genes[module]
            module_hn = []
            for i, gene1 in enumerate(module_genes):
                for gene2 in module_genes[i+1:]:
                    if (gene1, gene2) in humannet_set:
                        module_hn.append(f"{gene1}|{gene2}")
            for pair in module_hn:
                handle.write(f"{int(module)}\t{pair}\n")
                total_hn_written += 1
            if len(module_hn) > 0:
                modules_with_hn += 1
                print(f"Module {module}: Found {len(module_hn)} HumanNet interactions")

        if total_hn_written == 0:
            handle.write('# No HumanNet interactions found for module genes\n')

    print(f"\nHumanNet MVF file generation summary:")
    print(f"  - File created: {humannet_mvf_file}")
    print(f"  - Total HumanNet pairs written: {total_hn_written}")
    print(f"  - Modules with HumanNet interactions: {modules_with_hn}/{len(module2genes)}")

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
            
            # Draw subnetwork for this module
            draw_subnetwork(module, target_genes, regulators_dict, PKN_graph, metabolite_mapping,
                            PKN_df=PKN_df, humannet_set=humannet_set)

            modules_processed += 1
            modules_with_subnetworks += 1
            
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
