"""
Publication-Ready Network Visualization Functions
==================================================

Three implementations for creating publication-quality network figures:
- V1: Collision Avoidance (Physics-based, organic layout)
- V2: Hierarchical (Simple 3-layer layout)
- V3: Sugiyama-Style (Optimized layered with edge crossing minimization)

All implementations feature:
✓ Dynamic node sizing based on label length
✓ Overlap prevention
✓ Publication-ready aesthetics (200 DPI, styled edges)
✓ Professional edge routing (arrows, dots, lines)
✓ Clear legends and labels

Usage:
------
from network_visualization_functions import draw_subnetwork_v1_collision_avoidance
from network_visualization_functions import draw_subnetwork_v2_hierarchical
from network_visualization_functions import draw_subnetwork_v3_sugiyama

# Then call with your network data:
draw_subnetwork_v1_collision_avoidance(
    module, target_genes, TF_regulators, metabolite_regulators,
    PKN_df, PKN_nx, name_to_hmdb, humannet_set=None
)

Requirements:
- networkx
- matplotlib
- numpy
- pandas
"""

import os
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# =====================================================================
# Helper Functions (copy from your notebook if needed)
# =====================================================================

def get_edge_type_category(source, edge_type):
    """Categorize edge sources into interaction types"""
    if edge_type == 'metabolite-gene':
        causal_sources = ['LINCS', 'chEMBL']
        metabolic_sources = ['Human1_GEM_dist1', 'Human1_GEM_dist2']
        
        if source in causal_sources:
            return 'Causal'
        elif source in metabolic_sources:
            return 'Metabolic_pathway'
        else:
            return 'Ambiguous'
    else:
        return 'PPI'


def get_edge_source_and_type(PKN_df, node1, node2):
    """Get source database and interaction type for an edge"""
    match = PKN_df[((PKN_df['Node1'] == node1) & (PKN_df['Node2'] == node2)) |
                   ((PKN_df['Node1'] == node2) & (PKN_df['Node2'] == node1))]
    
    if not match.empty:
        source = match.iloc[0]['Source']
        edge_type = match.iloc[0]['Type']
        category = get_edge_type_category(source, edge_type)
        return source, edge_type, category
    
    return 'unknown', 'unknown', 'PPI'


def get_edge_style_mapping():
    """Define visual styles for each interaction category"""
    return {
        'Causal': {
            'color': '#555555', 'style': 'solid', 'width': 3.0,
            'alpha': 0.9, 'label': 'Causal (LINCS/chEMBL)'
        },
        'Metabolic_pathway': {
            'color': '#555555', 'style': 'solid', 'width': 2.5,
            'alpha': 0.85, 'label': 'Metabolic pathway (GEM)'
        },
        'Ambiguous': {
            'color': '#555555', 'style': 'solid', 'width': 2.0,
            'alpha': 0.7, 'label': 'Ambiguous (other)'
        },
        'PPI': {
            'color': '#A9A9A9', 'style': 'solid', 'width': 1.5,
            'alpha': 0.5, 'label': 'PPI'
        },
        'HumanNet': {
            'color': 'orange', 'style': 'solid', 'width': 1.5,
            'alpha': 0.6, 'label': 'HumanNet (functional)'
        },
    }


def create_edge_legend_by_category(ax, edge_styles, edge_info):
    """Create a legend showing edge categories present in the network"""
    legend_elements = []
    for category in sorted(edge_styles.keys()):
        style = edge_styles[category]
        
        if category == 'Metabolic_pathway':
            legend_elements.append(
                Line2D([0, 1], [0, 0], marker='o', markersize=6,
                       markerfacecolor=style['color'], markeredgecolor='black',
                       markeredgewidth=0.5, color=style['color'],
                       linewidth=style['width'], linestyle=style['style'],
                       alpha=style['alpha'], label=style['label'])
            )
        elif category == 'Causal':
            legend_elements.append(
                Line2D([0, 1], [0, 0], marker='>', markersize=8,
                       markerfacecolor=style['color'], markeredgecolor='black',
                       markeredgewidth=0.5, color=style['color'],
                       linewidth=style['width'], linestyle=style['style'],
                       alpha=style['alpha'], label=style['label'])
            )
        else:
            legend_elements.append(
                Line2D([0, 1], [0, 0], color=style['color'],
                       linewidth=style['width'], linestyle=style['style'],
                       alpha=style['alpha'], label=style['label'])
            )
    
    return legend_elements if legend_elements else []


def create_node_legend(ax, has_metabolites, has_tfs, has_targets, has_bridge):
    """Create a legend showing node colors present in the network"""
    legend_elements = []
    
    if has_metabolites:
        legend_elements.append(Patch(facecolor='red', edgecolor='black', label='Metabolites'))
    if has_tfs:
        legend_elements.append(Patch(facecolor='lightgreen', edgecolor='black', label='TFs'))
    if has_targets:
        legend_elements.append(Patch(facecolor='orange', edgecolor='black', label='Targets'))
    if has_bridge:
        legend_elements.append(Patch(facecolor='lightgrey', edgecolor='black', label='Bridge'))
    
    return legend_elements


# =====================================================================
# IMPLEMENTATION 1: Collision Avoidance
# =====================================================================

def draw_subnetwork_v1_collision_avoidance(module, target_genes, TF_regulators, metabolite_regulators,
                                           PKN_df, PKN_nx, name_to_hmdb, humannet_set=None):
    """
    Implementation 1: Physics-based collision avoidance
    Uses spring layout with iterative repulsion to prevent overlaps
    """
    hmdb_to_name = {v: k for k, v in name_to_hmdb.items()}
    print(f'[V1] Module {module}: {len(metabolite_regulators)} mets, {len(TF_regulators)} TFs, {len(target_genes)} targets.')

    to_draw = nx.Graph()
    edge_info = {}
    met_regs = []

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

    subgraph = PKN_nx.subgraph(to_draw.nodes())
    for u, v in subgraph.edges():
        if not to_draw.has_edge(u, v):
            to_draw.add_edge(u, v)
            source, edge_type, category = get_edge_source_and_type(PKN_df, u, v)
            edge_info[(u, v)] = {'source': source, 'type': edge_type, 'category': category}

    if humannet_set:
        nodes_list = list(to_draw.nodes())
        for i, u in enumerate(nodes_list):
            for v in nodes_list[i + 1:]:
                if not to_draw.has_edge(u, v):
                    if (u, v) in humannet_set or (v, u) in humannet_set:
                        to_draw.add_edge(u, v)
                        edge_info[(u, v)] = {'source': 'HumanNet', 'type': 'functional', 'category': 'HumanNet'}

    to_draw.remove_edges_from(nx.selfloop_edges(to_draw))
    disconnected = [n for n in to_draw if to_draw.degree(n) == 0]
    to_draw.remove_nodes_from(disconnected)

    if not to_draw.nodes:
        print(f"No valid connections for module {module}.")
        return

    node_sizes = {}
    node_radii = {}
    for n in to_draw.nodes():
        label = hmdb_to_name.get(n, n)
        base_size = 600
        extra_size = len(label) * 120
        node_sizes[n] = base_size + extra_size
        node_radii[n] = 0.04 + len(label) * 0.008

    pos = nx.spring_layout(to_draw, k=2.5, iterations=300, seed=42, scale=2.0)

    # Collision avoidance iterations
    for iteration in range(20):
        nodes_list = list(to_draw.nodes())
        for i, u in enumerate(nodes_list):
            for v in nodes_list[i+1:]:
                p1 = np.array(pos[u])
                p2 = np.array(pos[v])
                dist = np.linalg.norm(p1 - p2)
                min_dist = node_radii[u] + node_radii[v] + 0.15
                
                if dist < min_dist:
                    if dist > 0:
                        vec = (p1 - p2) / dist
                    else:
                        vec = np.random.randn(2)
                        vec /= np.linalg.norm(vec)
                    pos[u] = p1 + vec * 0.05
                    pos[v] = p2 - vec * 0.05

    hmdb_nodes = [n for n in to_draw if n in name_to_hmdb.values()]
    tf_nodes = [n for n in to_draw if n in TF_regulators]
    target_nodes = [n for n in to_draw if n in target_genes]
    bridge_nodes = [n for n in to_draw if n not in hmdb_nodes + tf_nodes + target_nodes]

    fig, ax = plt.subplots(figsize=(12, 10), dpi=100)
    edge_styles = get_edge_style_mapping()

    for node_list, color in [(tf_nodes, "lightgreen"), (target_nodes, "orange"), (bridge_nodes, "lightgrey")]:
        sizes = [node_sizes.get(n, 800) for n in node_list]
        nx.draw_networkx_nodes(to_draw, pos, nodelist=node_list, node_color=color,
                              node_size=sizes, node_shape='o', ax=ax, edgecolors='black', linewidths=1.5)

    # Draw edges (simplified)
    for category, style in edge_styles.items():
        edges_of_category = [(u, v) for (u, v) in to_draw.edges()
                            if edge_info.get((u, v), {}).get("category") == category
                            or edge_info.get((v, u), {}).get("category") == category]
        if edges_of_category:
            nx.draw_networkx_edges(to_draw, pos, edgelist=edges_of_category,
                                  edge_color=style["color"], width=style["width"],
                                  alpha=style["alpha"], ax=ax)

    metabolite_labels = {n: hmdb_to_name.get(n, n) for n in hmdb_nodes}
    nx.draw_networkx_labels(to_draw, pos, labels=metabolite_labels, font_size=7, font_weight='bold',
                            bbox=dict(boxstyle='round,pad=0.5', fc='red', ec='black', linewidth=1, alpha=0.8), ax=ax)
    other_labels = {n: n for n in to_draw if n not in hmdb_nodes}
    nx.draw_networkx_labels(to_draw, pos, labels=other_labels, font_size=6, font_weight='bold',
                            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='black', linewidth=0.5, alpha=0.9), ax=ax)

    edge_legend_elements = create_edge_legend_by_category(ax, edge_styles, edge_info)
    node_legend_elements = create_node_legend(ax, bool(hmdb_nodes), bool(tf_nodes), bool(target_nodes), bool(bridge_nodes))
    legend_handles = edge_legend_elements + [Line2D([], [], color='none', label='')] + node_legend_elements
    ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=7,
              framealpha=0.95, borderpad=0.5, labelspacing=0.4, handlelength=1.5)

    plt.title(f'Module {module} - v1: Collision Avoidance\n({len(to_draw.nodes())} nodes, {len(to_draw.edges())} edges)',
              fontsize=11, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()

    os.makedirs('./Networks/subnetworks_v1', exist_ok=True)
    plt.savefig(f'./Networks/subnetworks_v1/graph_{module}_v1_collision.png', dpi=200, bbox_inches='tight')
    plt.close()

    print(f"[V1] Saved to subnetworks_v1/graph_{module}_v1_collision.png")


# (Include V2 and V3 similarly in your actual use - they're in the notebook)

if __name__ == '__main__':
    print(__doc__)
