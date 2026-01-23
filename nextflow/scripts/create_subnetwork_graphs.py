#!/usr/bin/env python3
"""
Create subnetwork visualization graphs showing how metabolite regulators
connect to module genes through the Lemonite PKN.

Generates graph_{module}.png files in the Networks/subnetworks/ directory.
Enhanced with edge categorization and sophisticated visualization.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def get_edge_source_and_type(PKN_df, node1, node2):
    """
    Get the source and type information for an edge from the PKN DataFrame.
    
    Parameters:
    -----------
    PKN_df : pd.DataFrame
        The PKN network DataFrame
    node1, node2 : str
        The two nodes forming the edge
        
    Returns:
    --------
    tuple
        (source, edge_type, category)
    """
    # Find the edge in either direction
    edge = PKN_df[((PKN_df['Node1'] == node1) & (PKN_df['Node2'] == node2)) |
                  ((PKN_df['Node1'] == node2) & (PKN_df['Node2'] == node1))]
    
    if len(edge) == 0:
        return 'unknown', 'unknown', 'PPI'
    
    row = edge.iloc[0]
    source = row.get('Source', 'unknown')
    edge_type = row.get('Type', 'unknown')
    
    # Categorize based on source and type
    if source in ['LINCS', 'chEMBL'] or 'LINCS' in source or 'chEMBL' in source:
        category = 'Causal'
    elif source == 'GEM' or 'GEM' in source:
        category = 'Metabolic_pathway'
    elif source in ['STITCH', 'BioGRID', 'IntAct', 'STRING'] or any(db in source for db in ['STITCH', 'BioGRID', 'IntAct', 'STRING']):
        category = 'PPI'
    else:
        category = 'Ambiguous'
    
    return source, edge_type, category


def create_edge_legend_by_category(ax, edge_styles, edge_info):
    """Create a legend showing edge categories present in the network"""
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch
    
    # Get unique categories in this network
    present_categories = set()
    for info in edge_info.values():
        present_categories.add(info.get('category', 'PPI'))
    
    # Create legend elements for all categories in edge_styles
    legend_elements = []
    for category in sorted(edge_styles.keys()):
        style = edge_styles[category]
        
        if category == 'Metabolic_pathway':
            # For Metabolic_pathway, show a line with a dot marker at the end
            legend_elements.append(
                Line2D([0, 1], [0, 0], 
                      marker='o', markersize=6, markerfacecolor=style['color'],
                      markeredgecolor='black', markeredgewidth=0.5,
                      color=style['color'], linewidth=style['width'],
                      linestyle=style['style'], alpha=style['alpha'],
                      label=style['label'])
            )
        elif category == 'Causal':
            # For Causal, show a line with an arrow marker
            legend_elements.append(
                Line2D([0, 1], [0, 0], 
                      marker='>', markersize=8, markerfacecolor=style['color'],
                      markeredgecolor='black', markeredgewidth=0.5,
                      color=style['color'], linewidth=style['width'],
                      linestyle=style['style'], alpha=style['alpha'],
                      label=style['label'])
            )
        else:
            # Normal line for others (PPI, Ambiguous)
            legend_elements.append(
                Line2D([0, 1], [0, 0], 
                      color=style['color'], 
                      linewidth=style['width'],
                      linestyle=style['style'],
                      alpha=style['alpha'],
                      label=style['label'])
            )
    
    # Add legend if there are elements
    if legend_elements:
        return legend_elements
    return []


def create_node_legend(ax, has_metabolites=False, has_tfs=False, has_targets=False, has_bridge=False):
    """Create a legend for node types"""
    legend_elements = []
    
    if has_metabolites:
        legend_elements.append(
            Patch(facecolor='red', edgecolor='none', alpha=0.7, 
                  label='Metabolite regulators')
        )
    
    if has_tfs:
        legend_elements.append(
            plt.Circle((0, 0), 0.05, facecolor='lightgreen', edgecolor='black', 
                      linewidth=1, label='TF regulators')
        )
    
    if has_targets:
        legend_elements.append(
            plt.Circle((0, 0), 0.05, facecolor='orange', edgecolor='black', 
                      linewidth=1, label='Module target genes')
        )
    
    if has_bridge:
        legend_elements.append(
            plt.Circle((0, 0), 0.05, facecolor='lightgrey', edgecolor='black', 
                      linewidth=1, label='Bridge nodes')
        )
    
    return legend_elements


def export_to_cytoscape_with_categories(graph, module, hmdb_to_name, TF_regulators, target_genes, edge_info, PKN_df):
    """Export network with edge category information for Cytoscape"""
    import os
    
    cytoscape_dir = "./Networks/cytoscape"
    os.makedirs(cytoscape_dir, exist_ok=True)
    
    # Export edges with full information
    edges_data = []
    for u, v in graph.edges():
        node1 = hmdb_to_name.get(u, u)
        node2 = hmdb_to_name.get(v, v)
        
        info = edge_info.get((u, v)) or edge_info.get((v, u), {})
        source = info.get('source', 'unknown')
        edge_type = info.get('type', 'unknown')
        category = info.get('category', 'Unknown')
        
        edges_data.append({
            'Node1': node1,
            'Node2': node2,
            'Source': source,
            'Type': edge_type,
            'Category': category
        })
    
    edges_df = pd.DataFrame(edges_data)
    edges_file = f"{cytoscape_dir}/module_{module}_edges_categorized.tsv"
    edges_df.to_csv(edges_file, sep='\t', index=False)
    
    # Export node attributes
    attributes_data = []
    for node in graph.nodes():
        if node in hmdb_to_name:
            node_type = 'metabolite'
            label = hmdb_to_name[node]
        elif node in TF_regulators:
            node_type = 'TF'
            label = node
        elif node in target_genes:
            node_type = 'module_gene'
            label = node
        else:
            node_type = 'bridge_node'
            label = node
        
        attributes_data.append({
            'Node': label,
            'Label': label,
            'Type': node_type
        })
    
    attributes_df = pd.DataFrame(attributes_data)
    attributes_file = f"{cytoscape_dir}/module_{module}_attributes.tsv"
    attributes_df.to_csv(attributes_file, sep='\t', index=False)
    
    print(f"\nCytoscape files with categorized edges saved for module {module}:")
    print(f"  Edges: {edges_file}")
    print(f"  Attributes: {attributes_file}")


def get_edge_style_mapping():
    """
    Define visual styles for each interaction category
    Metabolite-gene interactions have 3 types (Causal, Metabolic_pathway, Ambiguous)
    PPI interactions are all shown as gray lines
    """
    return {
        # Metabolite-gene interaction styles - all dark grey, distinguished by markers
        'Causal': {
            'color': '#555555',  # Dark grey - causal relationships
            'style': 'solid',
            'width': 3.0,
            'alpha': 0.9,
            'label': 'Causal (LINCS/chEMBL)'
        },
        'Metabolic_pathway': {
            'color': '#555555',  # Dark grey - metabolic pathways
            'style': 'solid',
            'width': 2.5,
            'alpha': 0.85,
            'label': 'Metabolic pathway (GEM)'
        },
        'Ambiguous': {
            'color': '#555555',  # Dark grey - ambiguous
            'style': 'solid',
            'width': 2.0,
            'alpha': 0.7,
            'label': 'Ambiguous (other)'
        },
        
        # PPI interaction style (light grey)
        'PPI': {
            'color': '#A9A9A9',  # Darker light grey - protein-protein interactions
            'style': 'solid',
            'width': 1.5,
            'alpha': 0.5,
            'label': 'PPI'
        }
    }


def print_edge_statistics(edge_info):
    """Print statistics about edge categories in the network"""
    categories = {}
    for info in edge_info.values():
        cat = info.get('category', 'Unknown')
        categories[cat] = categories.get(cat, 0) + 1
    
    print(f"Edge statistics: {categories}")


def get_regulators(regfile):
    """
    Read regulators file and create a dictionary mapping modules to regulators.
    
    Parameters:
    -----------
    regfile : str
        Path to the regulators file (format: module\tregulator1|regulator2|...)
    
    Returns:
    --------
    dict
        Dictionary mapping module IDs to lists of regulators
    """
    regs = {}
    with open(regfile) as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) >= 2:
                module = parts[0]
                regulators = parts[1].split('|')
                regs[module] = regulators
    return regs


def draw_subnetwork(module, target_genes, TF_regulators, metabolite_regulators,
                    PKN_df, PKN_nx, name_to_hmdb, output_dir):
    """
    Draw subnetwork with clearer arrows, dots, and spacing.
    Uses edge categories and node types for styling.
    """

    # === Setup ===
    hmdb_to_name = {v: k for k, v in name_to_hmdb.items()}
    print(f'Module {module}: {len(metabolite_regulators)} metabolites, '
          f'{len(TF_regulators)} TFs, {len(target_genes)} targets.')

    to_draw = nx.Graph()
    edge_info = {}
    met_regs = []

    # --- STEP 1: Metabolite regulators ---
    for reg in metabolite_regulators:
        if reg not in name_to_hmdb:
            continue
        regulator = name_to_hmdb[reg]
        met_regs.append(regulator)
        if regulator not in PKN_nx:
            continue
        to_draw.add_node(regulator)
        for target in target_genes:
            try:
                if nx.shortest_path_length(PKN_nx, source=regulator, target=target) == 1:
                    to_draw.add_nodes_from([regulator, target])
                    to_draw.add_edge(regulator, target)
                    source, edge_type, category = get_edge_source_and_type(PKN_df, regulator, target)
                    edge_info[(regulator, target)] = {'source': source, 'type': edge_type, 'category': category}
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                continue

    # --- STEP 2: TF regulators (PPI edges) ---
    for TF in TF_regulators:
        for regulator in met_regs:
            try:
                length = nx.shortest_path_length(PKN_nx, source=TF, target=regulator)
                if length <= 2:
                    path = nx.shortest_path(PKN_nx, source=TF, target=regulator)
                    for i in range(len(path) - 1):
                        u, v = path[i], path[i + 1]
                        to_draw.add_edge(u, v)
                        source, edge_type, category = get_edge_source_and_type(PKN_df, u, v)
                        edge_info[(u, v)] = {'source': source, 'type': edge_type, 'category': category}
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                continue

    # --- STEP 3: Add edges among included nodes ---
    subgraph = PKN_nx.subgraph(to_draw.nodes())
    for u, v in subgraph.edges():
        if not to_draw.has_edge(u, v):
            to_draw.add_edge(u, v)
            source, edge_type, category = get_edge_source_and_type(PKN_df, u, v)
            edge_info[(u, v)] = {'source': source, 'type': edge_type, 'category': category}
    
    # Clean up isolated nodes
    to_draw.remove_edges_from(nx.selfloop_edges(to_draw))
    disconnected = [n for n in to_draw if to_draw.degree(n) == 0]
    to_draw.remove_nodes_from(disconnected)

    if not to_draw.nodes:
        print(f"No valid connections for module {module}.")
        return

    # === Layout ===
    pos = nx.spring_layout(to_draw, k=1.8, iterations=200, seed=42)

    # === Node categories ===
    hmdb_nodes = [n for n in to_draw if n in name_to_hmdb.values()]
    tf_nodes = [n for n in to_draw if n in TF_regulators]
    target_nodes = [n for n in to_draw if n in target_genes]
    bridge_nodes = [n for n in to_draw if n not in hmdb_nodes + tf_nodes + target_nodes]

    # === Figure ===
    fig, ax = plt.subplots(figsize=(8, 6))
    edge_styles = get_edge_style_mapping()

    # === Draw nodes ===
    nx.draw_networkx_nodes(
        to_draw, pos,
        nodelist=tf_nodes,
        node_color="lightgreen",
        node_size=800,
        node_shape='o',
        ax=ax
    )
    nx.draw_networkx_nodes(
        to_draw, pos,
        nodelist=target_nodes,
        node_color="orange",
        node_size=800,
        node_shape='o',
        ax=ax
    )
    nx.draw_networkx_nodes(
        to_draw, pos,
        nodelist=bridge_nodes,
        node_color="lightgrey",
        node_size=700,
        node_shape='o',
        ax=ax
    )
    # Metabolites are represented only by their label boxes, not as separate nodes
    
    # === Draw edges ===
    for category, style in edge_styles.items():
        edges_of_category = [
            (u, v) for (u, v) in to_draw.edges()
            if edge_info.get((u, v), {}).get("category") == category
            or edge_info.get((v, u), {}).get("category") == category
        ]
        if not edges_of_category:
            continue

        if category == "Metabolic_pathway":
            # Draw line ending in a dot at the edge of the target node
            # Orient: metabolite -> target
            for u, v in edges_of_category:
                # Determine direction: metabolite should be source
                if u in hmdb_nodes:
                    source, target = u, v
                elif v in hmdb_nodes:
                    source, target = v, u
                else:
                    # If neither is metabolite, just use u,v
                    source, target = u, v
                
                p1 = np.array(pos[source])
                p2 = np.array(pos[target])
                vec = p2 - p1
                dist = np.linalg.norm(vec)
                if dist == 0:
                    continue
                vec /= dist
                
                # Adjust node radius based on node type
                # Metabolites are squares with size 900, others are circles with size 700-800
                if target in hmdb_nodes:
                    node_radius = 0.065  # Slightly larger for metabolites
                else:
                    node_radius = 0.055
                
                # Stop line at node boundary
                p2_adj = p2 - vec * node_radius
                
                # Draw the line
                ax.plot([p1[0], p2_adj[0]], [p1[1], p2_adj[1]],
                        color=style["color"], linewidth=style["width"], 
                        alpha=style["alpha"], zorder=1)
                
                # Place dot exactly at the node boundary
                dot_point = p2 - vec * node_radius
                ax.plot(dot_point[0], dot_point[1], "o",
                        color=style["color"], markersize=8, 
                        markerfacecolor=style["color"],
                        markeredgecolor='black', markeredgewidth=0.5,
                        alpha=style["alpha"], zorder=3)

        elif category == "Causal":
            # Draw arrows from metabolite to target
            for u, v in edges_of_category:
                # Determine direction: arrow should point from metabolite to target
                start, end = (u, v) if u in hmdb_nodes else (v, u)
                p1 = np.array(pos[start])
                p2 = np.array(pos[end])
                vec = p2 - p1
                dist = np.linalg.norm(vec)
                if dist == 0:
                    continue
                vec /= dist
                
                # Adjust for node sizes
                start_radius = 0.065 if start in hmdb_nodes else 0.055
                end_radius = 0.065 if end in hmdb_nodes else 0.055
                
                p1_adj = p1 + vec * start_radius
                p2_adj = p2 - vec * end_radius
                
                arrow = FancyArrowPatch(
                    p1_adj, p2_adj,
                    arrowstyle='-|>',
                    color=style["color"],
                    linewidth=style["width"],
                    alpha=style["alpha"],
                    mutation_scale=15,
                    connectionstyle='arc3,rad=0',
                    shrinkA=0, shrinkB=0,
                    zorder=2
                )
                ax.add_patch(arrow)

        else:
            # Draw simple edges for PPI, Ambiguous, and Unknown
            nx.draw_networkx_edges(
                to_draw, pos,
                edgelist=edges_of_category,
                edge_color=style["color"],
                style=style["style"],
                width=style["width"],
                alpha=style["alpha"],
                arrows=False,
                ax=ax
            )
    
    # === Labels ===
    metabolite_labels = {n: hmdb_to_name.get(n, n) for n in hmdb_nodes}
    nx.draw_networkx_labels(to_draw, pos, labels=metabolite_labels,
                            font_size=6, font_weight='bold',
                            bbox=dict(boxstyle='round,pad=0.4', fc='red', ec='none', alpha=0.7), ax=ax)
    other_labels = {n: n for n in to_draw if n not in hmdb_nodes}
    nx.draw_networkx_labels(to_draw, pos, labels=other_labels, font_size=5, font_weight='bold', ax=ax)

    # === Legends ===
    edge_legend_elements = create_edge_legend_by_category(ax, edge_styles, edge_info)
    node_legend_elements = create_node_legend(ax,
                                              has_metabolites=bool(hmdb_nodes),
                                              has_tfs=bool(tf_nodes),
                                              has_targets=bool(target_nodes),
                                              has_bridge=bool(bridge_nodes))
    legend_handles = edge_legend_elements + [Line2D([], [], color='none', label='')] + node_legend_elements
    ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=6,
              framealpha=0.9, borderpad=0.5, labelspacing=0.3, handlelength=1.5)

    # === Finalize ===
    plt.title(f'Module {module} - Subnetwork in Lemonite PKN\n'
              f'({len(to_draw.nodes())} nodes, {len(to_draw.edges())} edges)', fontsize=10)
    plt.axis('off')
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f'{output_dir}/graph_{module}_categorized_clean.png', dpi=200, bbox_inches='tight')
    plt.close()

    print_edge_statistics(edge_info)
    export_to_cytoscape_with_categories(to_draw, module, hmdb_to_name, TF_regulators, target_genes, edge_info, PKN_df)


def main():
    parser = argparse.ArgumentParser(
        description='Create subnetwork visualization graphs for LemonTree modules'
    )
    parser.add_argument('--regulator_files', required=True,
                       help='Comma-separated list of regulator files (format: Type:Path,Type:Path)')
    parser.add_argument('--clusters', required=True,
                       help='Path to clusters file (module genes)')
    parser.add_argument('--pkn', required=True,
                       help='Path to PKN network file (TSV format)')
    parser.add_argument('--metabolite_mapping', required=False, default=None,
                       help='Path to metabolite name to HMDB ID mapping file (optional, only needed for metabolite regulators)')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for subnetwork graphs')
    
    args = parser.parse_args()
    
    print("="*60)
    print("Creating subnetwork visualizations")
    print("="*60)
    
    # Parse regulator files
    print(f"Parsing regulator files: {args.regulator_files}")
    regulator_configs = {}
    for reg_config in args.regulator_files.split(','):
        reg_type, reg_path = reg_config.split(':', 1)
        regulator_configs[reg_type] = reg_path.strip()
    
    print(f"Found regulator types: {list(regulator_configs.keys())}")
    
    # Read regulator files dynamically
    module2regulators = {}
    for reg_type, reg_path in regulator_configs.items():
        print(f"Reading {reg_type} regulators from {reg_path}")
        module2regulators[reg_type] = get_regulators(reg_path)
    
    print(f"Reading module genes from {args.clusters}")
    module2genes = get_regulators(args.clusters)
    
    # Read PKN network
    print(f"Reading PKN from {args.pkn}")
    PKN = pd.read_csv(args.pkn, sep='\t', header=0)
    
    # Process PKN: split Node1 on '_' and keep the last part (HMDB ID extraction)
    if 'Node1' in PKN.columns and any('_' in str(x) for x in PKN['Node1'].dropna()):
        PKN['Node1'] = PKN['Node1'].str.split('_').str[-1]
    
    # Create NetworkX graph from PKN
    PKN_nx = nx.from_pandas_edgelist(PKN, 'Node1', 'Node2')
    print(f'PKN contains {PKN_nx.number_of_nodes()} nodes and {PKN_nx.number_of_edges()} edges')
    
    # Read metabolite name to HMDB mapping if provided
    name_to_hmdb = {}
    if args.metabolite_mapping:
        print(f"Reading metabolite mapping from {args.metabolite_mapping}")
        metabolite_mapping = pd.read_csv(args.metabolite_mapping, sep=',')
        
        # Create a dictionary that maps metabolite names (Query) to HMDB IDs
        if 'Query' in metabolite_mapping.columns and 'HMDB' in metabolite_mapping.columns:
            name_to_hmdb = metabolite_mapping.set_index('Query')['HMDB'].dropna().to_dict()
        else:
            # Try alternative column names
            cols = metabolite_mapping.columns.tolist()
            if len(cols) >= 2:
                name_to_hmdb = metabolite_mapping.set_index(cols[0])[cols[1]].dropna().to_dict()
            else:
                print("ERROR: Could not parse metabolite mapping file")
                sys.exit(1)
        
        print(f'Found {len(name_to_hmdb)} metabolite mappings')
    else:
        print("No metabolite mapping provided - metabolite names will not be converted to HMDB IDs")
    
    # Create subnetwork graphs for each module
    # Use the first regulator type to determine available modules
    first_reg_type = list(regulator_configs.keys())[0]
    modules = list(module2regulators[first_reg_type].keys())
    print(f"\nCreating subnetwork graphs for {len(modules)} modules")
    successful = 0
    skipped = 0
    
    for module in modules:
        try:
            # Get module data
            target_genes = module2genes.get(module, [])
            
            # Collect regulators for this module from all types
            regulators_dict = {}
            for reg_type in regulator_configs.keys():
                regulators_dict[reg_type] = module2regulators[reg_type].get(module, [])
            
            if not target_genes:
                print(f'Module {module}: No target genes found, skipping')
                skipped += 1
                continue
            
            # Separate TF and metabolite regulators
            TF_regulators = []
            metabolite_regulators = []
            
            for reg_type, regulators in regulators_dict.items():
                if reg_type.lower() in ['tf', 'transcription_factor', 'tfs']:
                    TF_regulators.extend(regulators)
                elif reg_type.lower() in ['metabolite', 'lipid', 'metabolites', 'lipids']:
                    metabolite_regulators.extend(regulators)
                else:
                    # For other types, assume they are TF-like
                    TF_regulators.extend(regulators)
            
            if not TF_regulators and not metabolite_regulators:
                print(f'Module {module}: No TF or metabolite regulators found, skipping')
                skipped += 1
                continue
            
            # Create subnetwork visualization
            draw_subnetwork(
                module,
                target_genes,
                TF_regulators,
                metabolite_regulators,
                PKN,  # Pass the PKN DataFrame
                PKN_nx,  # Pass the NetworkX graph
                name_to_hmdb,
                args.output_dir
            )
            
            successful += 1
            
        except Exception as e:
            print(f'ERROR processing module {module}: {e}')
            import traceback
            traceback.print_exc()
            skipped += 1
            continue


if __name__ == '__main__':
    main()
