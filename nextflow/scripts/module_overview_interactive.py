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
import shutil
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

CANONICAL_CLUSTER_COLUMN = 'top_30'
CANONICAL_ONTOLOGY = 'BP'
CANONICAL_TOP_N = 30
RRVGO_THRESHOLD = 0.7
RRVGO_METHOD = 'Rel'
ENRICHMENT_SOURCE_PRIORITY = ('GSEA', 'EnrichR')

# ── Cytoscape.JS network constants ───────────────────────────────────────────
import colorsys

CLUSTER_PALETTE = ['#F44336', '#2196F3', '#4CAF50', '#FF9800', '#00E5FF']
UNASSIGNED_COLOR = '#B0B0B0'

# Per-regulator-type gradient anchors (light low-out-degree -> dark high). Known
# regulator types map to a curated colour; any other type present in the data is
# assigned the next free anchor from REGULATOR_GRADIENT_POOL. This keeps the
# config-driven regulator-types design working without code changes.
REGULATOR_TYPE_COLORS = {
    'TF':          ('#D7C3F0', '#5E35B1'),   # lavender   -> deep purple
    'TFs':         ('#D7C3F0', '#5E35B1'),
    'metabolite':  ('#FFE0B2', '#E65100'),   # pale amber -> burnt orange
    'metabolites': ('#FFE0B2', '#E65100'),
    'Metabolites': ('#FFE0B2', '#E65100'),
    'Metabolomics':('#FFE0B2', '#E65100'),
    'lipid':       ('#B2EBF2', '#006064'),   # pale cyan  -> deep teal
    'lipids':      ('#B2EBF2', '#006064'),
    'Lipids':      ('#B2EBF2', '#006064'),
}
# Fallback anchors handed out (in order) to any unrecognised regulator type.
REGULATOR_GRADIENT_POOL = [
    ('#C8E6C9', '#1B5E20'),   # green
    ('#F8BBD0', '#880E4F'),   # pink
    ('#FFF9C4', '#F57F17'),   # yellow
    ('#CFD8DC', '#37474F'),   # blue-grey
]
REGULATOR_FALLBACK_COLORS = ('#D0D0D0', '#606060')


def _hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def _rgb_to_hex(rgb):
    return '#{:02X}{:02X}{:02X}'.format(*(max(0, min(255, int(round(c)))) for c in rgb))


def _interp_color(light_hex, dark_hex, frac):
    """Interpolate light->dark in HSV so the gradient stays vivid (frac in [0,1])."""
    l = colorsys.rgb_to_hsv(*[c / 255 for c in _hex_to_rgb(light_hex)])
    d = colorsys.rgb_to_hsv(*[c / 255 for c in _hex_to_rgb(dark_hex)])
    hsv = tuple(l[i] + (d[i] - l[i]) * frac for i in range(3))
    rgb = tuple(c * 255 for c in colorsys.hsv_to_rgb(*hsv))
    return _rgb_to_hex(rgb)


def resolve_regulator_type_anchors(reg_types):
    """Return {reg_type: (light, dark)} for every regulator type present.

    Known types use their curated colour; unknown types are assigned, in sorted
    order, from REGULATOR_GRADIENT_POOL so each type gets a distinct hue.
    """
    anchors = {}
    pool_idx = 0
    for rt in sorted(reg_types):
        if rt in REGULATOR_TYPE_COLORS:
            anchors[rt] = REGULATOR_TYPE_COLORS[rt]
        else:
            anchors[rt] = REGULATOR_GRADIENT_POOL[pool_idx % len(REGULATOR_GRADIENT_POOL)]
            pool_idx += 1
    return anchors


def compute_regulator_colors(nodes, edges):
    """Return ({regulator_id: hex_color}, {regulator_id: out_degree}).

    Colour is shaded light->dark by out-degree, normalised within each regulator
    type so every type uses its full gradient range independently.
    """
    out_degree = {}
    for e in edges:
        src = str(e['source'])
        out_degree[src] = out_degree.get(src, 0) + 1

    by_type = {}
    for node in nodes:
        ntype = str(node.get('type', ''))
        if ntype == 'module':
            continue
        nid = str(node['id'])
        by_type.setdefault(ntype, []).append((nid, out_degree.get(nid, 0)))

    anchors = resolve_regulator_type_anchors(by_type.keys())
    colors = {}
    for ntype, items in by_type.items():
        light, dark = anchors.get(ntype, REGULATOR_FALLBACK_COLORS)
        degs = [d for _, d in items]
        lo, hi = min(degs), max(degs)
        for nid, deg in items:
            frac = (deg - lo) / (hi - lo) if hi > lo else 0.6
            colors[nid] = _interp_color(light, dark, frac)
    return colors, out_degree


def build_cluster_color_map(module_clusters):
    all_clusters = sorted(set(v for v in module_clusters.values() if v != 'Unassigned'))
    color_map = {cl: CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)] for i, cl in enumerate(all_clusters)}
    color_map['Unassigned'] = UNASSIGNED_COLOR
    return color_map


def _build_cytoscape_hover_text(module_id, module_overview, module_enrichment_all, edges, name_lookup=None):
    hover_text = f"<b>Module {module_id}</b><br>"
    if module_overview is not None and not module_overview.empty:
        rows = module_overview[module_overview['Module'].astype(str) == str(module_id)]
        if len(rows) > 0:
            row = rows.iloc[0]
            expr_p = row.get('Expression_adjusted_pval', 'NA')
            if expr_p != 'NA' and pd.notna(expr_p):
                try:
                    hover_text += f"<b>Expression Analysis:</b><br>"
                    hover_text += f"  \u2022 adj. p-value: {float(expr_p):.2e}<br>"
                    hover_text += f"  \u2022 Rank: {row.get('Expression_rank', 'NA')}<br><br>"
                except (ValueError, TypeError):
                    pass
            genes_val = row.get('Module_genes', 'NA')
            if genes_val != 'NA' and pd.notna(genes_val):
                hover_text += f"<b>Genes:</b> {len(str(genes_val).split('|'))} genes<br><br>"
    module_target = f"Module_{module_id}"
    reg_type_labels = {}
    for e in edges:
        if e['target'] == module_target:
            rtype = e.get('type', '').replace('_regulation', '')
            reg_type_labels.setdefault(rtype, set()).add(e['source'])
    for reg_type in sorted(reg_type_labels.keys()):
        regs = sorted(reg_type_labels[reg_type])
        if regs:
            display_regs = [name_lookup.get(r, r) for r in regs] if name_lookup else regs
            hover_text += f"<b>{reg_type.capitalize()} ({len(regs)}):</b> "
            hover_text += ', '.join(display_regs[:10])
            if len(regs) > 10:
                hover_text += f", ... (+{len(regs)-10} more)"
            hover_text += '<br>'
    hover_text += '<br>'
    if module_enrichment_all is not None and not module_enrichment_all.empty:
        mod_enrich = module_enrichment_all[module_enrichment_all['Module'].astype(str) == str(module_id)]
        for db in ['BP', 'MF', 'CC', 'KEGG', 'Reactome']:
            db_enrich = mod_enrich[mod_enrich['Database'].str.upper() == db.upper()].sort_values('p.adjust').head(3)
            if len(db_enrich) > 0:
                hover_text += f"<b>Top {db}:</b><br>"
                for _, erow in db_enrich.iterrows():
                    term = str(erow['Term'])
                    if len(term) > 55:
                        term = term[:52] + '...'
                    try:
                        hover_text += f"  \u2022 {term} (p={float(erow['p.adjust']):.1e})<br>"
                    except (ValueError, TypeError):
                        hover_text += f"  \u2022 {term}<br>"
    return hover_text


def create_cytoscape_html_network(nodes, edges, module_clusters, config_name,
                                   variant_name, output_file,
                                   module_overview=None, module_enrichment_all=None,
                                   name_lookup=None, filter_description=None):
    """Create a self-contained Cytoscape.JS interactive network HTML."""
    # The pipeline includes every regulator-module connection from module_data
    # (no top-N/connectivity filtering), so the header always says so -
    # regardless of variant_name, which is currently always 'filtered' here.
    if filter_description is None:
        filter_description = "Showing all regulator–module connections (unfiltered)"

    cluster_color_map = build_cluster_color_map(module_clusters)
    regulator_colors, regulator_out_degree = compute_regulator_colors(nodes, edges)

    # Map each MegaGO cluster -> its rrvgo descriptive label (for the legend).
    cluster_label_map = {}
    if module_overview is not None and not module_overview.empty and 'MegaGO_Label' in module_overview.columns:
        _ccol = 'MegaGO_Cluster' if 'MegaGO_Cluster' in module_overview.columns else 'Functional_Cluster'
        if _ccol in module_overview.columns:
            for _cl, _grp in module_overview.dropna(subset=[_ccol]).groupby(_ccol):
                _lbls = _grp['MegaGO_Label'].dropna()
                if len(_lbls):
                    cluster_label_map[str(_cl)] = str(_lbls.iloc[0])

    cy_nodes = []
    for node in nodes:
        node_id = str(node['id'])
        node_type = str(node.get('type', ''))
        if node_type == 'module':
            module_num = node_id.replace('Module_', '')
            label = f"M{module_num}"
            cluster = module_clusters.get(module_num, 'Unassigned')
            color = cluster_color_map.get(cluster, UNASSIGNED_COLOR)
            node_class = 'module'
            hover_info = _build_cytoscape_hover_text(
                module_num, module_overview, module_enrichment_all, edges, name_lookup
            )
            hover_info += f"<b>MegaGO Cluster:</b> {cluster}"
            # Also show the rrvgo descriptive label for the cluster, if available.
            if module_overview is not None and not module_overview.empty:
                _lrows = module_overview[module_overview['Module'].astype(str) == str(module_num)]
                if len(_lrows) > 0:
                    _lbl = _lrows.iloc[0].get('MegaGO_Label', '')
                    if pd.notna(_lbl) and str(_lbl).strip() not in ('', 'NA', 'nan'):
                        hover_info += f"<br><b>MegaGO Label:</b> {_lbl}"
        else:
            # Full label (no truncation); shaded light->dark by out-degree.
            label = str(node.get('label', node_id))
            cluster = None
            color = regulator_colors.get(node_id, REGULATOR_FALLBACK_COLORS[1])
            node_class = 'regulator'
            target_modules = sorted(set(
                e['target'].replace('Module_', 'M')
                for e in edges if str(e['source']) == node_id
            ))
            hover_info = f"<b>{node.get('label', node_id)}</b><br>"
            hover_info += f"<b>Type:</b> {node_type}<br>"
            hover_info += f"<b>Out-degree:</b> {regulator_out_degree.get(node_id, 0)} module(s)<br>"
            if target_modules:
                hover_info += f"<b>Targets ({len(target_modules)}):</b> {', '.join(target_modules)}"
        cy_nodes.append({
            'data': {'id': node_id, 'label': label, 'type': node_type,
                     'cluster': cluster if cluster else 'N/A', 'hover_info': hover_info},
            'classes': node_class,
            'style': {'background-color': color},
        })

    cy_edges = []
    for edge in edges:
        category = str(edge.get('category', 'Other'))
        corr = edge.get('correlation', 0)
        color = '#27AE60' if corr > 0 else ('#E74C3C' if corr < 0 else '#5D6D7E')
        source = str(edge['source'])
        target = str(edge['target'])
        cy_edges.append({'data': {
            'id': f"{source}-{target}", 'source': source, 'target': target,
            'category': category, 'correlation': round(float(corr), 4), 'color': color,
        }})

    # Group module ids by cluster so the fcose layout can seed each cluster
    # around its own anchor (keeps same-cluster modules clumped together).
    cluster_groups = {}
    for n in cy_nodes:
        if n['classes'] == 'module':
            cl = n['data']['cluster']
            if cl and cl != 'Unassigned' and cl != 'N/A':
                cluster_groups.setdefault(cl, []).append(n['data']['id'])
    ordered_clusters = sorted(cluster_groups.keys())
    cluster_members = {cl: cluster_groups[cl] for cl in ordered_clusters}

    elements_json = json.dumps(cy_nodes + cy_edges, indent=2)
    cluster_members_json = json.dumps(cluster_members)

    # Legends ----------------------------------------------------------
    # Cluster legend (module colours)
    cluster_legend_html = ""
    for cl_name, cl_color in cluster_color_map.items():
        cluster_legend_html += (
            f'        <div class="legend-item">\n'
            f'            <div class="legend-color" style="background: {cl_color};"></div>\n'
            f'            <span class="legend-label">{cl_name}</span>\n'
            f'        </div>\n'
        )
    # Regulator legend (one gradient swatch per type that is actually present)
    present_reg_types = sorted({str(n.get('type', '')) for n in nodes if str(n.get('type', '')) != 'module'})
    reg_anchors = resolve_regulator_type_anchors(present_reg_types)
    regulator_legend_html = ""
    for rt in present_reg_types:
        light, dark = reg_anchors.get(rt, REGULATOR_FALLBACK_COLORS)
        rt_label = (rt[:1].upper() + rt[1:]) if rt else rt
        regulator_legend_html += (
            f'        <div class="legend-item">\n'
            f'            <div class="legend-color" '
            f'style="background: linear-gradient(135deg, {light} 0%, {dark} 100%);"></div>\n'
            f'            <span class="legend-label">{rt_label} (light→dark = out-degree)</span>\n'
            f'        </div>\n'
        )

    # Structured legend data (mirrors the on-screen legend) so the PNG
    # export can redraw the legend directly onto the exported canvas.
    legend_data = {
        'sections': [
            {
                'title': 'Regulator types',
                'items': [
                    {
                        'label': f'{(rt[:1].upper() + rt[1:]) if rt else rt} (light→dark = out-degree)',
                        'kind': 'gradient',
                        'light': reg_anchors.get(rt, REGULATOR_FALLBACK_COLORS)[0],
                        'dark': reg_anchors.get(rt, REGULATOR_FALLBACK_COLORS)[1],
                    }
                    for rt in present_reg_types
                ],
            },
            {
                'title': 'Module clusters',
                'items': [
                    {'label': (f"{cl_name}: {cluster_label_map[cl_name]}"
                               if cluster_label_map.get(cl_name) else cl_name),
                     'kind': 'dot', 'color': cl_color}
                    for cl_name, cl_color in cluster_color_map.items()
                ],
            },
            {
                'title': 'Edge types',
                'items': [
                    {'label': 'Causal (arrow)', 'kind': 'line', 'color': '#555555', 'width': 3, 'arrow': 'triangle'},
                    {'label': 'Metabolic (dot)', 'kind': 'line', 'color': '#555555', 'width': 2, 'arrow': 'dot'},
                    {'label': 'Other (line)', 'kind': 'line', 'color': '#555555', 'width': 1, 'arrow': 'none'},
                ],
            },
            {
                'title': 'Correlation',
                'items': [
                    {'label': 'Positive', 'kind': 'line', 'color': '#27AE60', 'width': 3, 'arrow': 'none'},
                    {'label': 'Negative', 'kind': 'line', 'color': '#E74C3C', 'width': 3, 'arrow': 'none'},
                    {'label': 'N/A', 'kind': 'line', 'color': '#5D6D7E', 'width': 2, 'arrow': 'none'},
                ],
            },
        ]
    }
    legend_data_json = json.dumps(legend_data)

    # Generate HTML ----------------------------------------------------
    html_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Module-Regulator Network: {variant_name}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.23.0/cytoscape.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/layout-base@2.0.1/layout-base.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cose-base@2.2.0/cose-base.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cytoscape-fcose@2.2.0/cytoscape-fcose.js"></script>
<script src="https://cdn.jsdelivr.net/npm/cytoscape-svg@0.4.0/cytoscape-svg.min.js"></script>
<style>
    * {{
        margin: 0;
        padding: 0;
        box-sizing: border-box;
    }}

    body {{
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        overflow: hidden;
    }}

    #cy {{
        width: 100%;
        height: calc(100vh - 100px);
        display: block;
        background: linear-gradient(to bottom, #f5f5f5 0%, #ffffff 100%);
    }}

    .header {{
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 15px 20px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    }}

    .header h1 {{
        margin: 0;
        font-size: 20px;
    }}

    .header p {{
        margin: 5px 0 0 0;
        font-size: 12px;
        opacity: 0.9;
    }}

    .controls {{
        position: absolute;
        top: 110px;
        left: 20px;
        background: white;
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 15px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        z-index: 100;
        font-size: 12px;
    }}

    .control-group {{
        margin-bottom: 12px;
    }}

    .control-group label {{
        display: block;
        font-weight: bold;
        margin-bottom: 5px;
        color: #333;
    }}

    button {{
        padding: 8px 12px;
        margin: 3px 0;
        cursor: pointer;
        background: #667eea;
        color: white;
        border: none;
        border-radius: 3px;
        font-size: 11px;
        font-weight: 600;
        transition: background 0.2s;
        display: block;
        width: 100%;
    }}

    button:hover {{
        background: #764ba2;
    }}

    .legend {{
        position: absolute;
        bottom: 20px;
        right: 20px;
        background: white;
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 15px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        font-size: 12px;
        max-width: 280px;
        max-height: calc(100vh - 200px);
        overflow-y: auto;
    }}

    .legend-item {{
        display: flex;
        align-items: center;
        margin-bottom: 8px;
    }}

    .legend-color {{
        width: 20px;
        height: 20px;
        border-radius: 50%;
        margin-right: 8px;
        border: 2px solid #fff;
        flex-shrink: 0;
    }}

    .legend-line {{
        width: 30px;
        height: 0;
        margin-right: 8px;
        flex-shrink: 0;
    }}

    .legend-label {{
        color: #333;
    }}

    .info-box {{
        position: absolute;
        bottom: 20px;
        left: 20px;
        background: white;
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 10px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        font-size: 11px;
        color: #666;
        max-width: 250px;
    }}

    .search-box {{
        position: relative;
    }}

    .search-box input {{
        width: 100%;
        padding: 6px 8px;
        border: 1px solid #ccc;
        border-radius: 3px;
        font-size: 11px;
        outline: none;
    }}

    .search-box input:focus {{
        border-color: #667eea;
    }}

    .autocomplete-list {{
        position: absolute;
        top: 100%;
        left: 0;
        right: 0;
        background: white;
        border: 1px solid #ccc;
        border-top: none;
        border-radius: 0 0 3px 3px;
        max-height: 150px;
        overflow-y: auto;
        z-index: 200;
        display: none;
    }}

    .autocomplete-list div {{
        padding: 4px 8px;
        cursor: pointer;
        font-size: 11px;
    }}

    .autocomplete-list div:hover,
    .autocomplete-list div.active {{
        background: #667eea;
        color: white;
    }}

    /* Tooltip for node hover */
    #tooltip {{
        position: absolute;
        display: none;
        background: rgba(255, 255, 255, 0.97);
        border: 1px solid #aaa;
        border-radius: 6px;
        padding: 10px 14px;
        box-shadow: 0 4px 16px rgba(0,0,0,0.18);
        font-size: 12px;
        line-height: 1.5;
        color: #333;
        max-width: 500px;
        z-index: 1000;
        pointer-events: none;
    }}
</style>
</head>
<body>
<div class="header">
    <h1>Module-Regulator Network</h1>
    <p>{filter_description} | Nodes: {len(cy_nodes)} | Edges: {len(cy_edges)}</p>
</div>

<div id="cy"></div>

<!-- Tooltip element for node hover info -->
<div id="tooltip"></div>

<div class="controls">
    <div class="control-group">
        <label>Layout:</label>
        <button onclick="runLayout('fcose')">Cluster grouped (fCoSE)</button>
        <button onclick="runLayout('concentric')">Concentric</button>
        <button onclick="runLayout('breadthfirst')">Breadth-First</button>
    </div>

    <div class="control-group">
        <label>Actions:</label>
        <button onclick="fitView()">Fit View</button>
        <button onclick="resetView()">Reset</button>
        <button onclick="exportPositions()">Export Positions</button>
        <button onclick="exportPNG()">&#128247; Export PNG (hi-res)</button>
        <button onclick="exportSVG()">&#9998; Export SVG (vector)</button>
    </div>
    <div class="control-group">
        <label>Search:</label>
        <div class="search-box">
            <input type="text" id="searchInput" placeholder="Type node name..." autocomplete="off" />
            <div class="autocomplete-list" id="autocompleteList"></div>
        </div>
        <button onclick="clearHighlight()" style="margin-top:6px">Clear Highlight</button>
    </div>
</div>

<div class="legend">
    <strong>Regulator types:</strong>
{regulator_legend_html}
    <hr style="margin: 8px 0; border: none; border-top: 1px solid #ddd;">
    <strong>Module clusters:</strong>
{cluster_legend_html}
    <hr style="margin: 8px 0; border: none; border-top: 1px solid #ddd;">
    <strong>Edge types:</strong>
    <div class="legend-item">
        <div class="legend-line" style="border-top: 3px solid #555;"><span style="float:right;font-size:12px;margin-top:-9px;color:#555">&#9654;</span></div>
        <span class="legend-label">Causal (arrow)</span>
    </div>
    <div class="legend-item">
        <div class="legend-line" style="border-top: 2px solid #555;"><span style="float:right;font-size:8px;margin-top:-6px;color:#555">&#9679;</span></div>
        <span class="legend-label">Metabolic (dot)</span>
    </div>
    <div class="legend-item">
        <div class="legend-line" style="border-top: 1px solid #555;"></div>
        <span class="legend-label">Other (line)</span>
    </div>
    <hr style="margin: 8px 0; border: none; border-top: 1px solid #ddd;">
    <strong>Correlation:</strong>
    <div class="legend-item">
        <div class="legend-line" style="border-top: 3px solid #27AE60;"></div>
        <span class="legend-label">Positive</span>
    </div>
    <div class="legend-item">
        <div class="legend-line" style="border-top: 3px solid #E74C3C;"></div>
        <span class="legend-label">Negative</span>
    </div>
    <div class="legend-item">
        <div class="legend-line" style="border-top: 2px solid #5D6D7E;"></div>
        <span class="legend-label">N/A</span>
    </div>
</div>

<div class="info-box">
    <strong>Controls:</strong><br>
    &bull; Drag nodes to move<br>
    &bull; Scroll to zoom<br>
    &bull; Pan with right-click drag<br>
    &bull; Hover over nodes for details<br>
    &bull; PNG/SVG export keeps your manual layout
</div>

<script>
    // Register fcose extension (cluster-grouping force layout)
    if (typeof cytoscapeFcose !== 'undefined') {{
        cytoscape.use(cytoscapeFcose);
    }}
    // Register cytoscape-svg extension
    if (typeof cytoscapeSvg !== 'undefined') {{
        cytoscape.use(cytoscapeSvg);
    }}

    // Build elements data
    var elements = {elements_json};

    // Module ids grouped by cluster -> used to seed each cluster around a
    // big circle so fcose relaxes into visibly separated clumps.
    var clusterMembers = {cluster_members_json};

    // Legend data for compositing onto the exported PNG.
    var legendData = {legend_data_json};

    // Compute a fixed anchor point per cluster on a large circle, then place
    // every module of that cluster in a small jittered disc around its
    // anchor. fcose then relaxes locally from these seeds (it does not lock
    // them), so clusters stay clumped while nodes remain draggable.
    function clusterAnchors() {{
        var clusters = Object.keys(clusterMembers);
        var R = 350 + clusters.length * 60;   // ring radius scales with #clusters
        var anchors = {{}};
        clusters.forEach(function(cl, i) {{
            var ang = (2 * Math.PI * i) / clusters.length;
            anchors[cl] = {{ x: R * Math.cos(ang), y: R * Math.sin(ang) }};
        }});
        return anchors;
    }}

    // Seed module positions near their cluster anchor before running fcose.
    function seedClusterPositions() {{
        var anchors = clusterAnchors();
        Object.keys(clusterMembers).forEach(function(cl) {{
            var a = anchors[cl];
            var members = clusterMembers[cl];
            var spread = 40 + members.length * 12;
            members.forEach(function(id, k) {{
                var n = cy.getElementById(id);
                if (n.length === 0) return;
                var t = (2 * Math.PI * k) / Math.max(1, members.length);
                n.position({{
                    x: a.x + spread * Math.cos(t) * Math.random(),
                    y: a.y + spread * Math.sin(t) * Math.random()
                }});
            }});
        }});
    }}

    // fcose layout. randomize:false makes it relax FROM the seeded
    // positions instead of starting over, preserving the cluster clumps.
    function fcoseOptions() {{
        return {{
            name: 'fcose',
            quality: 'proof',
            randomize: false,        // start from seeded cluster positions
            animate: true,
            animationDuration: 800,
            fit: true,
            padding: 60,
            nodeSeparation: 110,
            idealEdgeLength: 80,
            nodeRepulsion: 9000,
            gravity: 0.15,           // weak global gravity so clumps don't merge
            gravityRange: 2.0,
            numIter: 2500,
            tile: false
        }};
    }}

    // Initialize Cytoscape
    var cy = cytoscape({{
        container: document.getElementById('cy'),
        elements: elements,

        style: [
            {{
                selector: 'node.module',
                style: {{
                    'label': 'data(label)',
                    'width': '65px',
                    'height': '65px',
                    'border-width': '3px',
                    'border-color': 'rgba(255, 255, 255, 0.8)',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'font-size': '13px',
                    'font-weight': 'bold',
                    'color': '#000',
                    'text-outline-width': '0px',
                    'cursor': 'grab'
                }}
            }},
            {{
                selector: 'node.regulator',
                style: {{
                    'label': 'data(label)',
                    'text-wrap': 'none',
                    'text-events': 'no',
                    'shape': 'round-rectangle',
                    'width': '50px',
                    'height': '50px',
                    'border-width': '2px',
                    'border-color': 'rgba(255, 255, 255, 0.7)',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'font-size': '11px',
                    'font-weight': 'bold',
                    'color': '#000',
                    'text-outline-width': '0px',
                    'cursor': 'grab'
                }}
            }},
            {{
                selector: 'node:selected',
                style: {{
                    'border-width': '4px',
                    'border-color': '#FFD700',
                    'shadow-blur': '10px',
                    'shadow-color': 'rgba(0, 0, 0, 0.5)',
                    'shadow-offset-x': '0px',
                    'shadow-offset-y': '4px'
                }}
            }},
            {{
                selector: 'edge[category="Causal"]',
                style: {{
                    'width': '3px',
                    'line-color': 'data(color)',
                    'target-arrow-color': 'data(color)',
                    'target-arrow-shape': 'triangle',
                    'target-arrow-fill': 'filled',
                    'curve-style': 'bezier',
                    'opacity': 0.8
                }}
            }},
            {{
                selector: 'edge[category="Metabolic_pathway"]',
                style: {{
                    'width': '2px',
                    'line-color': 'data(color)',
                    'target-arrow-color': 'data(color)',
                    'target-arrow-shape': 'circle',
                    'target-arrow-fill': 'filled',
                    'curve-style': 'bezier',
                    'line-dash-pattern': [5, 5],
                    'opacity': 0.7
                }}
            }},
            {{
                selector: 'edge[category="Other"]',
                style: {{
                    'width': '1px',
                    'line-color': 'data(color)',
                    'target-arrow-color': 'data(color)',
                    'curve-style': 'bezier',
                    'opacity': 0.5
                }}
            }},
            {{
                selector: 'edge:selected',
                style: {{
                    'width': '4px',
                    'opacity': 1,
                    'shadow-blur': '8px',
                    'shadow-color': 'rgba(0, 0, 0, 0.4)'
                }}
            }},
            {{
                selector: '.dimmed',
                style: {{
                    'opacity': 0.15
                }}
            }},
            {{
                selector: '.highlighted',
                style: {{
                    'opacity': 1,
                    'border-width': '4px',
                    'border-color': '#FFD700'
                }}
            }}
        ],

        layout: {{ name: 'preset' }},   // positions set by seedClusterPositions(), then fcose runs

        wheelSensitivity: 0.1,
        autoungrabify: false,
        userPanningEnabled: true
    }});

    // Seed cluster positions then run fcose so clusters relax into clumps.
    seedClusterPositions();
    cy.layout(fcoseOptions()).run();

    // Store initial layout positions
    var initialPositions = {{}};
    cy.on('layoutstop', function() {{
        initialPositions = {{}};
        cy.nodes().forEach(function(n) {{
            initialPositions[n.id()] = {{ x: n.position('x'), y: n.position('y') }};
        }});
    }});

    // Fit to view on load
    setTimeout(() => cy.fit(undefined, 50), 500);

    // Mouse tracking for tooltip
    var tooltip = document.getElementById('tooltip');
    cy.on('mouseover', 'node', function(e) {{
        var node = e.target;
        var hover_info = node.data('hover_info');
        tooltip.innerHTML = hover_info;
        tooltip.style.display = 'block';
    }});

    cy.on('mouseout', 'node', function(e) {{
        tooltip.style.display = 'none';
    }});

    cy.on('mousemove', function(e) {{
        if (tooltip.style.display === 'block') {{
            tooltip.style.left = (e.originalEvent.pageX + 10) + 'px';
            tooltip.style.top = (e.originalEvent.pageY + 10) + 'px';
        }}
    }});

    // Layout switching
    function runLayout(layoutName) {{
        if (layoutName === 'fcose') {{
            seedClusterPositions();
            cy.layout(fcoseOptions()).run();
            return;
        }}
        cy.layout({{
            name: layoutName,
            directed: false,
            animate: true,
            animationDuration: 800,
            fit: true,
            padding: 50,
            spacingFactor: 1.2
        }}).run();
    }}

    // Fit to view
    function fitView() {{
        cy.fit(undefined, 50);
    }}

    // Reset view and positions
    function resetView() {{
        if (Object.keys(initialPositions).length > 0) {{
            cy.nodes().forEach(function(n) {{
                if (initialPositions[n.id()]) {{
                    n.position(initialPositions[n.id()]);
                }}
            }});
        }}
        cy.fit(undefined, 50);
    }}

    // Export node positions as JSON
    function exportPositions() {{
        var positions = {{}};
        cy.nodes().forEach(function(n) {{
            positions[n.id()] = {{ x: n.position('x'), y: n.position('y') }};
        }});
        console.log(JSON.stringify(positions, null, 2));
        alert('Positions exported to console (F12)');
    }}

    // Export the WHOLE network at its current node positions as a high-res
    // PNG. `full: true` renders the entire graph model (not just the visible
    // viewport), so any manual node rearrangement is preserved in the image.
    // Draw the legend onto a canvas context at the given offset. Returns
    // the total height consumed, so callers can size the canvas. Pass
    // measureOnly=true to compute height without drawing.
    function drawLegend(ctx, x, y, scale, measureOnly) {{
        var pad = 18 * scale;
        var titleFont = 'bold ' + (15 * scale) + 'px Segoe UI, sans-serif';
        var itemFont = (13 * scale) + 'px Segoe UI, sans-serif';
        var rowH = 26 * scale;
        var sectGap = 14 * scale;
        var swatch = 18 * scale;
        var cy0 = y + pad;

        legendData.sections.forEach(function(section) {{
            if (!measureOnly) {{
                ctx.fillStyle = '#222';
                ctx.font = titleFont;
                ctx.textBaseline = 'middle';
                ctx.fillText(section.title, x + pad, cy0 + rowH / 2);
            }}
            cy0 += rowH;

            section.items.forEach(function(item) {{
                var sx = x + pad;
                var sy = cy0 + rowH / 2;
                if (!measureOnly) {{
                    if (item.kind === 'dot') {{
                        ctx.beginPath();
                        ctx.arc(sx + swatch / 2, sy, swatch / 2, 0, 2 * Math.PI);
                        ctx.fillStyle = item.color;
                        ctx.fill();
                    }} else if (item.kind === 'gradient') {{
                        var g = ctx.createLinearGradient(sx, sy, sx + swatch, sy);
                        g.addColorStop(0, item.light);
                        g.addColorStop(1, item.dark);
                        ctx.beginPath();
                        ctx.arc(sx + swatch / 2, sy, swatch / 2, 0, 2 * Math.PI);
                        ctx.fillStyle = g;
                        ctx.fill();
                    }} else if (item.kind === 'line') {{
                        ctx.strokeStyle = item.color;
                        ctx.lineWidth = (item.width || 2) * scale;
                        ctx.beginPath();
                        ctx.moveTo(sx, sy);
                        ctx.lineTo(sx + swatch + 8 * scale, sy);
                        ctx.stroke();
                        if (item.arrow === 'triangle') {{
                            var ax = sx + swatch + 8 * scale;
                            ctx.beginPath();
                            ctx.moveTo(ax, sy);
                            ctx.lineTo(ax - 6 * scale, sy - 4 * scale);
                            ctx.lineTo(ax - 6 * scale, sy + 4 * scale);
                            ctx.closePath();
                            ctx.fillStyle = item.color;
                            ctx.fill();
                        }} else if (item.arrow === 'dot') {{
                            ctx.beginPath();
                            ctx.arc(sx + swatch + 8 * scale, sy, 3 * scale, 0, 2 * Math.PI);
                            ctx.fillStyle = item.color;
                            ctx.fill();
                        }}
                    }}
                    ctx.fillStyle = '#333';
                    ctx.font = itemFont;
                    ctx.textBaseline = 'middle';
                    ctx.fillText(item.label, sx + swatch + 16 * scale, sy);
                }}
                cy0 += rowH;
            }});
            cy0 += sectGap;
        }});
        return (cy0 - y);
    }}

    // Export the WHOLE network at its current node positions as a high-res
    // PNG, with the legend composited onto a white panel on the right.
    function exportPNG() {{
        var scale = 4;   // print-quality raster; legend is matched to this scale

        // WYSIWYG export: render exactly what is currently on screen (the
        // visible viewport, current pan/zoom) so the legend can be placed at
        // the *same position and size* it has in the HTML at this moment.
        var cyEl = document.getElementById('cy');
        var legendEl = document.querySelector('.legend');
        var cyRect = cyEl.getBoundingClientRect();
        var legRect = legendEl ? legendEl.getBoundingClientRect() : null;

        var netDataUrl = cy.png({{
            output: 'base64uri',
            scale: scale,
            full: false,        // current viewport only -> matches the screen
            bg: 'white'
        }});
        var netImg = new Image();
        netImg.onload = function() {{
            var canvas = document.createElement('canvas');
            canvas.width = netImg.width;
            canvas.height = netImg.height;
            var ctx = canvas.getContext('2d');
            ctx.fillStyle = '#ffffff';
            ctx.fillRect(0, 0, canvas.width, canvas.height);
            ctx.drawImage(netImg, 0, 0);

            if (legRect) {{
                // Legend position/size relative to the #cy container, in CSS
                // px, scaled to output px. (0,0) of the rendered viewport PNG
                // is the top-left of #cy, so this reproduces the on-screen
                // placement and footprint exactly.
                var lx = (legRect.left - cyRect.left) * scale;
                var ly = (legRect.top  - cyRect.top ) * scale;
                var lw = legRect.width  * scale;
                var lh = legRect.height * scale;

                // Pick the legend draw-scale so its rendered height fills the
                // on-screen legend box (height parity with the HTML legend).
                var natH = drawLegend(ctx, 0, 0, 1, true);
                var legendScale = lh / natH;

                // White rounded panel + border + soft shadow, mirroring the
                // .legend CSS so the exported legend looks identical.
                var r = 4 * scale;
                ctx.save();
                ctx.shadowColor = 'rgba(0,0,0,0.10)';
                ctx.shadowBlur = 8 * scale;
                ctx.shadowOffsetY = 2 * scale;
                roundRect(ctx, lx, ly, lw, lh, r);
                ctx.fillStyle = '#ffffff';
                ctx.fill();
                ctx.restore();
                roundRect(ctx, lx, ly, lw, lh, r);
                ctx.strokeStyle = '#dddddd';
                ctx.lineWidth = 1 * scale;
                ctx.stroke();

                // Clip to the panel and draw the legend content.
                ctx.save();
                roundRect(ctx, lx, ly, lw, lh, r);
                ctx.clip();
                drawLegend(ctx, lx, ly, legendScale, false);
                ctx.restore();
            }}

            canvas.toBlob(function(blob) {{
                var filename = '{config_name}_{variant_name}.png';
                var url = URL.createObjectURL(blob);
                var a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(url);
            }}, 'image/png');
        }};
        netImg.src = netDataUrl;
    }}

    // Rounded-rectangle path helper for the legend panel.
    function roundRect(ctx, x, y, w, h, r) {{
        r = Math.min(r, w / 2, h / 2);
        ctx.beginPath();
        ctx.moveTo(x + r, y);
        ctx.arcTo(x + w, y,     x + w, y + h, r);
        ctx.arcTo(x + w, y + h, x,     y + h, r);
        ctx.arcTo(x,     y + h, x,     y,     r);
        ctx.arcTo(x,     y,     x + w, y,     r);
        ctx.closePath();
    }}

    // Export the whole network as SVG (vector — best for publication).
    // Also uses current node positions.
    function exportSVG() {{
        if (typeof cy.svg !== 'function') {{
            alert('SVG export plugin not loaded — check your internet connection and reload.');
            return;
        }}
        var svgContent = cy.svg({{scale: 1, full: true, bg: 'white'}});
        var blob = new Blob([svgContent], {{type: 'image/svg+xml'}});
        var filename = '{config_name}_{variant_name}.svg';
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }}

    // Search / filter functionality
    var nodeList = [];
    cy.nodes().forEach(function(n) {{
        nodeList.push(n.data('label'));
    }});

    var searchInput = document.getElementById('searchInput');
    var autocompleteList = document.getElementById('autocompleteList');

    searchInput.addEventListener('input', function() {{
        var query = searchInput.value.toLowerCase();
        var matches = nodeList.filter(n => n.toLowerCase().includes(query));

        autocompleteList.innerHTML = '';
        if (query && matches.length > 0) {{
            matches.slice(0, 10).forEach(function(match) {{
                var div = document.createElement('div');
                div.textContent = match;
                div.onclick = function() {{
                    searchInput.value = match;
                    autocompleteList.innerHTML = '';
                    highlightNode(match);
                }};
                autocompleteList.appendChild(div);
            }});
            autocompleteList.style.display = 'block';
        }} else {{
            autocompleteList.style.display = 'none';
        }}
    }});

    searchInput.addEventListener('keypress', function(e) {{
        if (e.key === 'Enter') {{
            var query = searchInput.value;
            highlightNode(query);
            autocompleteList.innerHTML = '';
            autocompleteList.style.display = 'none';
        }}
    }});

    function highlightNode(nodeLabel) {{
        cy.nodes().removeClass('dimmed').removeClass('highlighted');
        var node = cy.nodes().filter(n => n.data('label') === nodeLabel)[0];

        if (node) {{
            node.addClass('highlighted');
            cy.nodes().not(node).addClass('dimmed');
            cy.center(node);
        }}
    }}

    function clearHighlight() {{
        cy.nodes().removeClass('dimmed').removeClass('highlighted');
        searchInput.value = '';
        autocompleteList.innerHTML = '';
    }}
</script>
</body>
</html>
"""

    # Save HTML
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"✓ Cytoscape network saved: {output_file}")

    return output_file



def module_sort_key(module_id):
    text = str(module_id)
    try:
        return (0, float(text), text)
    except ValueError:
        return (1, text, text)


def empty_enrichment_frame():
    return pd.DataFrame(columns=['Module', 'Database', 'Term', 'p.adjust', '__direction__', '__source__'])


def standardize_enrichment_df(df, direction=None, source=None):
    if df is None or df.empty:
        return empty_enrichment_frame()

    standardized = df.copy()

    mod_col = next((c for c in standardized.columns if c.lower() in ('module', 'cluster', 'moduleid', 'mod')), None)
    if mod_col and mod_col != 'Module':
        standardized = standardized.rename(columns={mod_col: 'Module'})

    term_col = next((c for c in standardized.columns if c.lower() in ('term', 'description', 'pathway', 'name')), None)
    if term_col and term_col != 'Term':
        standardized = standardized.rename(columns={term_col: 'Term'})

    padj_col = next((c for c in standardized.columns if c.lower() in ('padj', 'p.adjust', 'fdr', 'qvalue', 'p_adj', 'adjusted.p.value', 'p.value')), None)
    if padj_col and padj_col != 'p.adjust':
        standardized = standardized.rename(columns={padj_col: 'p.adjust'})

    db_col = next((c for c in standardized.columns if c.lower() in ('database', 'db', 'source')), None)
    if db_col and db_col != 'Database':
        standardized = standardized.rename(columns={db_col: 'Database'})

    if 'Module' not in standardized.columns or 'Term' not in standardized.columns:
        return empty_enrichment_frame()

    if 'Database' not in standardized.columns:
        standardized['Database'] = ''
    if 'p.adjust' not in standardized.columns:
        standardized['p.adjust'] = np.nan

    standardized['Module'] = standardized['Module'].astype(str).str.strip()
    standardized['Database'] = standardized['Database'].astype(str).str.strip()
    standardized['Term'] = standardized['Term'].astype(str).str.strip()
    standardized['p.adjust'] = pd.to_numeric(standardized['p.adjust'], errors='coerce')
    standardized['__direction__'] = direction if direction is not None else standardized.get('__direction__', 'Up')
    standardized['__source__'] = source if source is not None else standardized.get('__source__', '')

    standardized = standardized[['Module', 'Database', 'Term', 'p.adjust', '__direction__', '__source__']]
    standardized = standardized.dropna(subset=['Module', 'Term', 'p.adjust'])
    if standardized.empty:
        return empty_enrichment_frame()

    return standardized.sort_values(['Module', 'p.adjust', 'Term', '__direction__']).reset_index(drop=True)


def infer_direction_from_filename(filename):
    lowered = filename.lower()
    if re.search(r'(^|_)down(_|\.)', lowered):
        return 'Down'
    return 'Up'


def normalize_enrichment_method(method):
    lowered = str(method).strip().lower()
    if lowered == 'enrichr':
        return 'EnrichR'
    if lowered == 'gsea':
        return 'GSEA'
    return 'auto'


def bucket_enrichment_files(csv_files):
    buckets = {
        source: {
            'all_up': [],
            'all_down': [],
            'top_10_up': [],
            'top_10_down': []
        }
        for source in ENRICHMENT_SOURCE_PRIORITY
    }

    for path in sorted(dict.fromkeys(csv_files)):
        filename = os.path.basename(path).lower()
        if 'enrichr' in filename:
            source = 'EnrichR'
        elif 'gsea' in filename:
            source = 'GSEA'
        else:
            continue

        if 'all_enriched_pathways' in filename:
            granularity = 'all'
        elif 'top_10_enriched_pathways' in filename:
            granularity = 'top_10'
        else:
            continue

        direction = infer_direction_from_filename(filename).lower()
        buckets[source][f'{granularity}_{direction}'].append(path)

    return buckets


def choose_enrichment_source(file_buckets, requested_method):
    requested = normalize_enrichment_method(requested_method)
    available_sources = [
        source for source in ENRICHMENT_SOURCE_PRIORITY
        if any(file_buckets[source][key] for key in file_buckets[source])
    ]

    if not available_sources:
        return None

    if requested in ('EnrichR', 'GSEA'):
        if requested in available_sources:
            return requested
        print(f"Requested enrichment source {requested} not found. Falling back to {available_sources[0]}.")
        return available_sources[0]

    if len(available_sources) > 1:
        print(f"Both EnrichR and GSEA outputs are available; preferring {available_sources[0]} for overview labeling.")
    else:
        print(f"Auto-detected enrichment source: {available_sources[0]}")

    return available_sources[0]


def database_key(database_value):
    value = str(database_value).strip().lower()
    if value in ('bp', 'go_bp') or 'biological_process' in value or 'biological process' in value:
        return 'bp'
    if value in ('mf', 'go_mf') or 'molecular_function' in value or 'molecular function' in value:
        return 'mf'
    if value in ('cc', 'go_cc') or 'cellular_component' in value or 'cellular component' in value:
        return 'cc'
    if 'kegg' in value:
        return 'kegg'
    if 'react' in value:
        return 'reactome'
    return 'other'


def extract_go_id(term):
    match = re.search(r'GO:\d+', str(term))
    return match.group(0) if match else None


def prepare_top30_bp_terms(enrichment_data, output_dir, selected_source, top_n=CANONICAL_TOP_N):
    top30_dir = os.path.join(output_dir, CANONICAL_CLUSTER_COLUMN)
    os.makedirs(top30_dir, exist_ok=True)

    bp_terms_path = os.path.join(top30_dir, f'bp_terms_{CANONICAL_CLUSTER_COLUMN}.csv')
    bp_df = enrichment_data.get('bp', empty_enrichment_frame()).copy()

    if bp_df.empty:
        empty_terms = empty_enrichment_frame().assign(GO_ID=pd.Series(dtype='object'))
        empty_terms.to_csv(bp_terms_path, index=False)
        print(f"No BP enrichment terms available. Wrote empty top_30 file to: {bp_terms_path}")
        return empty_terms, bp_terms_path, top30_dir

    bp_df = standardize_enrichment_df(bp_df, source=selected_source)
    bp_df['GO_ID'] = bp_df['Term'].map(extract_go_id)
    bp_df['__term_key__'] = bp_df['GO_ID'].fillna(bp_df['Term'])
    bp_df = bp_df.drop_duplicates(subset=['Module', '__term_key__'], keep='first')
    bp_df = bp_df.sort_values(['Module', 'p.adjust', 'Term', '__direction__'])
    top_terms = bp_df.groupby('Module', sort=True).head(top_n).copy()
    top_terms = top_terms.drop(columns=['__term_key__'])
    top_terms.to_csv(bp_terms_path, index=False)

    print(f"Prepared canonical BP top_30 terms for {top_terms['Module'].nunique()} modules")
    print(f"Top_30 BP term file: {bp_terms_path}")
    return top_terms, bp_terms_path, top30_dir


def write_cluster_assignments(module_clusters, output_dir, cluster_column=CANONICAL_CLUSTER_COLUMN):
    cluster_path = os.path.join(output_dir, f'cluster_assignments_{cluster_column}.csv')
    cluster_df = pd.DataFrame([
        {'Module': str(module_id), cluster_column: cluster_name}
        for module_id, cluster_name in sorted(module_clusters.items(), key=lambda item: module_sort_key(item[0]))
    ])
    cluster_df.to_csv(cluster_path, index=False)
    print(f"Cluster assignments saved to: {cluster_path}")
    return cluster_df, cluster_path


def load_optional_csv(csv_path, columns):
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return pd.DataFrame(columns=columns)


def run_rrvgo_labeler(bp_terms_path, cluster_path, output_dir, organism, cluster_column=CANONICAL_CLUSTER_COLUMN):
    helper_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rrvgo_cluster_labels.R')
    empty_outputs = {
        'cluster_labels': pd.DataFrame(columns=['MegaGO_cluster', 'MegaGO_label']),
        'module_labels': pd.DataFrame(columns=['Module', 'MegaGO_cluster', 'MegaGO_label', 'representative_go_id', 'label_source']),
        'reduced_terms': pd.DataFrame(columns=['MegaGO_cluster', 'representative_label'])
    }

    if not os.path.exists(helper_path):
        print(f"Warning: rrvgo helper not found at {helper_path}; skipping label generation.")
        return empty_outputs

    command = [
        'Rscript',
        helper_path,
        f'--enrichment={bp_terms_path}',
        f'--clusters={cluster_path}',
        f'--out-dir={output_dir}',
        f'--cluster-column={cluster_column}',
        f'--organism={organism}',
        f'--ontology={CANONICAL_ONTOLOGY}',
        f'--top-n={CANONICAL_TOP_N}',
        f'--threshold={RRVGO_THRESHOLD}',
        f'--method={RRVGO_METHOD}'
    ]

    print("Running rrvgo cluster labeling helper...")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout.strip())
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr.strip())
        print("Warning: rrvgo helper failed; continuing without cluster labels.")
        return empty_outputs

    cluster_labels_path = os.path.join(output_dir, f'rrvgo_cluster_labels_{cluster_column}.csv')
    module_labels_path = os.path.join(output_dir, f'rrvgo_module_labels_{cluster_column}.csv')
    reduced_terms_path = os.path.join(output_dir, f'rrvgo_reduced_terms_{cluster_column}.csv')

    return {
        'cluster_labels': load_optional_csv(cluster_labels_path, ['MegaGO_cluster', 'MegaGO_label']),
        'module_labels': load_optional_csv(module_labels_path, ['Module', 'MegaGO_cluster', 'MegaGO_label', 'representative_go_id', 'label_source']),
        'reduced_terms': load_optional_csv(reduced_terms_path, ['MegaGO_cluster', 'representative_label'])
    }

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


def load_regulator_scores(score_file):
    """
    Load selected regulator score file (TSV with header: Regulator, Target, Score, Overall_rank).
    
    Parameters:
    score_file (str): Path to the selected regulators score file
    
    Returns:
    pd.DataFrame: DataFrame with columns [Regulator, Module, Score, Overall_rank]
    """
    try:
        df = pd.read_csv(score_file, sep='\t')
        # Standardize column names
        rename_map = {}
        if 'Target' in df.columns:
            rename_map['Target'] = 'Module'
        df = df.rename(columns=rename_map)
        df['Module'] = df['Module'].astype(str)
        return df
    except FileNotFoundError:
        print(f"Warning: Score file {score_file} not found.")
        return pd.DataFrame(columns=['Regulator', 'Module', 'Score', 'Overall_rank'])
    except Exception as e:
        print(f"Error reading score file {score_file}: {e}")
        return pd.DataFrame(columns=['Regulator', 'Module', 'Score', 'Overall_rank'])
 

def generate_regulator_tables_html(regulator_scores_dict, output_dir, name_lookup=None):
    """
    Generate an HTML file with regulator ranking tables.
    
    For each regulator type, creates:
    1. A table of regulator-module pairs ranked by score
    2. A summary table per regulator (total score, number of target modules, list of modules)
    
    Parameters:
    regulator_scores_dict (dict): {type_name: pd.DataFrame} with score data
    output_dir (str): Directory to save the HTML file
    
    Returns:
    str: Path to the generated HTML file
    """
    if not regulator_scores_dict:
        print("No regulator score data available — skipping regulator tables.")
        return None

    html_parts = []
    html_parts.append("""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Regulator Rankings</title>
<style>
  body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background: #f8f9fa; color: #333; }
  h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
  h2 { color: #2c3e50; margin-top: 40px; }
  h3 { color: #34495e; margin-top: 25px; }
  .table-container { overflow-x: auto; margin: 15px 0; }
  table { border-collapse: collapse; width: 100%; background: white; box-shadow: 0 1px 3px rgba(0,0,0,0.12); }
  th { background: #3498db; color: white; padding: 10px 14px; text-align: left; font-weight: 600;
       position: sticky; top: 0; cursor: pointer; user-select: none; white-space: nowrap; }
  th:hover { background: #2980b9; }
  th .sort-icon { font-size: 0.7em; margin-left: 4px; opacity: 0.7; }
  td { padding: 8px 14px; border-bottom: 1px solid #ecf0f1; }
  tr:hover { background: #ebf5fb; }
  tr:nth-child(even) { background: #f9f9f9; }
  tr:nth-child(even):hover { background: #ebf5fb; }
  .positive { color: #27ae60; font-weight: 600; }
  .negative { color: #e74c3c; }
  .module-list { font-size: 0.9em; max-width: 500px; word-wrap: break-word; }
  .section { margin-bottom: 50px; padding: 20px; background: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
  .toc { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 30px; }
  .toc a { color: #3498db; text-decoration: none; display: block; padding: 4px 0; }
  .toc a:hover { text-decoration: underline; }
  .search-box { padding: 8px 12px; border: 1px solid #ddd; border-radius: 4px; width: 300px; margin: 10px 0; font-size: 14px; }
  .search-box:focus { outline: none; border-color: #3498db; box-shadow: 0 0 3px rgba(52,152,219,0.3); }
</style>
</head>
<body>
<h1>Regulator Rankings</h1>
""")

    # Table of contents
    html_parts.append('<div class="toc"><strong>Contents</strong>')
    for reg_type in regulator_scores_dict:
        safe_id = reg_type.replace(' ', '_')
        html_parts.append(f'<a href="#{safe_id}_pairs">{reg_type} — Regulator-Module Pairs</a>')
        html_parts.append(f'<a href="#{safe_id}_summary">{reg_type} — Regulator Summary</a>')
    html_parts.append('</div>')

    for reg_type, df in regulator_scores_dict.items():
        if df.empty:
            continue

        safe_id = reg_type.replace(' ', '_')

        html_parts.append(f'<div class="section">')
        html_parts.append(f'<h2>{reg_type}</h2>')

        # --- Table 1: Regulator-Module pairs ranked by score ---
        html_parts.append(f'<h3 id="{safe_id}_pairs">Regulator–Module Pairs (ranked by score)</h3>')
        html_parts.append(f'<input type="text" class="search-box" id="search_{safe_id}_pairs" '
                          f'onkeyup="filterTable(\'{safe_id}_pairs_table\', this.value)" '
                          f'placeholder="Search regulators or modules...">')
        html_parts.append('<div class="table-container">')

        pairs_df = df.sort_values('Score', ascending=False).reset_index(drop=True)
        html_parts.append(f'<table id="{safe_id}_pairs_table">')
        html_parts.append('<thead><tr>'
                          '<th onclick="sortTable(\'{safe_id}_pairs_table\', 0)">Rank <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_pairs_table\', 1)">Regulator <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_pairs_table\', 2)">Module <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_pairs_table\', 3)">Score <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_pairs_table\', 4)">Overall Rank <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          '</tr></thead><tbody>')

        for rank, (_, row) in enumerate(pairs_df.iterrows(), 1):
            score = row['Score']
            score_class = 'positive' if score > 0 else 'negative'
            overall_rank = row.get('Overall_rank', 'NA')
            display_reg = name_lookup.get(row['Regulator'], row['Regulator']) if name_lookup else row['Regulator']
            html_parts.append(
                f'<tr><td>{rank}</td>'
                f'<td>{display_reg}</td>'
                f'<td>{row["Module"]}</td>'
                f'<td class="{score_class}">{score:.4f}</td>'
                f'<td>{overall_rank}</td></tr>'
            )

        html_parts.append('</tbody></table></div>')

        # --- Table 2: Regulator summary ---
        html_parts.append(f'<h3 id="{safe_id}_summary">Regulator Summary</h3>')
        html_parts.append(f'<input type="text" class="search-box" id="search_{safe_id}_summary" '
                          f'onkeyup="filterTable(\'{safe_id}_summary_table\', this.value)" '
                          f'placeholder="Search regulators...">')
        html_parts.append('<div class="table-container">')

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

        html_parts.append(f'<table id="{safe_id}_summary_table">')
        html_parts.append('<thead><tr>'
                          f'<th onclick="sortTable(\'{safe_id}_summary_table\', 0)">Rank <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_summary_table\', 1)">Regulator <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_summary_table\', 2)">Sum of Scores <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_summary_table\', 3)">N Target Modules <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          f'<th onclick="sortTable(\'{safe_id}_summary_table\', 4)">Target Modules <span class="sort-icon">&#x25B2;&#x25BC;</span></th>'
                          '</tr></thead><tbody>')

        for rank, (_, row) in enumerate(summary_df.iterrows(), 1):
            score_class = 'positive' if row['Total_Score'] > 0 else 'negative'
            display_reg = name_lookup.get(row['Regulator'], row['Regulator']) if name_lookup else row['Regulator']
            html_parts.append(
                f'<tr><td>{rank}</td>'
                f'<td>{display_reg}</td>'
                f'<td class="{score_class}">{row["Total_Score"]:.4f}</td>'
                f'<td>{row["N_Modules"]}</td>'
                f'<td class="module-list">{row["Target_Modules"]}</td></tr>'
            )

        html_parts.append('</tbody></table></div>')
        html_parts.append('</div>')  # close section

    # JavaScript for sorting and filtering
    html_parts.append("""
<script>
function sortTable(tableId, colIdx) {
  var table = document.getElementById(tableId);
  var tbody = table.querySelector('tbody');
  var rows = Array.from(tbody.querySelectorAll('tr'));
  var asc = table.getAttribute('data-sort-col') == colIdx && table.getAttribute('data-sort-dir') == 'asc';
  var dir = asc ? 'desc' : 'asc';
  table.setAttribute('data-sort-col', colIdx);
  table.setAttribute('data-sort-dir', dir);
  rows.sort(function(a, b) {
    var aVal = a.cells[colIdx].textContent.trim();
    var bVal = b.cells[colIdx].textContent.trim();
    var aNum = parseFloat(aVal);
    var bNum = parseFloat(bVal);
    if (!isNaN(aNum) && !isNaN(bNum)) {
      return dir === 'asc' ? aNum - bNum : bNum - aNum;
    }
    return dir === 'asc' ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
  });
  rows.forEach(function(row) { tbody.appendChild(row); });
}
function filterTable(tableId, query) {
  var table = document.getElementById(tableId);
  var rows = table.querySelectorAll('tbody tr');
  var q = query.toLowerCase();
  rows.forEach(function(row) {
    row.style.display = row.textContent.toLowerCase().indexOf(q) > -1 ? '' : 'none';
  });
}
</script>
</body>
</html>""")

    output_file = os.path.join(output_dir, 'regulator_rankings.html')
    with open(output_file, 'w') as f:
        f.write('\n'.join(html_parts))
    print(f"Regulator rankings HTML saved to: {output_file}")
    return output_file


def create_comprehensive_html_report(regulator_tables_html_path, network_html_path, output_dir):
    """
    Create comprehensive HTML report combining network visualization and regulator rankings
    
    Parameters:
    regulator_tables_html_path (str): Path to regulator rankings HTML file
    network_html_path (str): Path to network visualization HTML file
    output_dir (str): Output directory
    
    Returns:
    str: Path to the comprehensive HTML report
    """
    print("\nCreating comprehensive HTML report...")
    
    # Read regulator tables HTML
    regulator_content = ""
    if regulator_tables_html_path and os.path.exists(regulator_tables_html_path):
        try:
            with open(regulator_tables_html_path, 'r') as f:
                regulator_html = f.read()
            # Extract body content
            import re
            body_match = re.search(r'<body>(.*?)</body>', regulator_html, re.DOTALL)
            if body_match:
                regulator_content = body_match.group(1)
        except Exception as e:
            print(f"Warning: Could not read regulator rankings file: {e}")
            regulator_content = f'<p>Regulator rankings available in: <a href="{os.path.basename(regulator_tables_html_path)}">{os.path.basename(regulator_tables_html_path)}</a></p>'
    
    # Read network visualization HTML
    network_iframe = ""
    if network_html_path and os.path.exists(network_html_path):
        network_iframe = f'<iframe src="{os.path.basename(network_html_path)}" width="100%" height="800px" frameborder="0"></iframe>'
    
    # Create comprehensive HTML
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Lemonite Module Overview - Comprehensive Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 0;
            background: #f8f9fa;
            color: #333;
        }}
        .header {{
            background: linear-gradient(135deg, #2E7D32 0%, #1B5E20 100%);
            color: white;
            padding: 2rem;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            font-size: 2.5rem;
            margin: 0 0 0.5rem 0;
        }}
        .header .subtitle {{
            font-size: 1.1rem;
            opacity: 0.9;
        }}
        .nav {{
            background: white;
            padding: 1rem;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            display: flex;
            justify-content: center;
            flex-wrap: wrap;
            gap: 0.5rem;
        }}
        .nav a {{
            color: #2E7D32;
            text-decoration: none;
            padding: 0.5rem 1rem;
            border-radius: 20px;
            transition: all 0.3s ease;
            font-weight: 500;
        }}
        .nav a:hover {{
            background: #81C784;
            color: white;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 2rem;
        }}
        .section {{
            background: white;
            border-radius: 12px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }}
        .section h2 {{
            color: #2E7D32;
            border-bottom: 2px solid #81C784;
            padding-bottom: 1rem;
            margin-bottom: 1.5rem;
        }}
        .footer {{
            text-align: center;
            padding: 2rem;
            background: #2E7D32;
            color: white;
            margin-top: 3rem;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🍋🌳 Lemonite Module Overview</h1>
        <div class="subtitle">Comprehensive Multi-Omics Integration & Network Analysis</div>
    </div>
    
    <nav class="nav">
        <a href="#network-section">📊 Network Visualization</a>
        <a href="#regulator-section">🔬 Regulator Rankings</a>
    </nav>
    
    <div class="container">
        <section class="section" id="network-section">
            <h2>📊 Interactive Module Network Visualization</h2>
            <p>Explore the module-regulator network. Hover over nodes for details, zoom and pan to navigate.</p>
            {network_iframe}
        </section>
        
        <section class="section" id="regulator-section">
            <h2>🔬 Regulator-Module Rankings</h2>
            <p>Regulator-module pairs ranked by association scores. Higher scores indicate stronger predicted regulatory relationships.</p>
            {regulator_content}
        </section>
    </div>
    
    <div class="footer">
        <p>🍋🌳 <strong>Lemonite Pipeline</strong> | Module Overview Report</p>
        <p>Interactive analysis of regulatory modules and multi-omics integration</p>
    </div>
</body>
</html>'''
    
    # Save comprehensive report
    output_file = os.path.join(output_dir, 'Module_Overview_Comprehensive.html')
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"✓ Comprehensive HTML report created: {output_file}")
    print(f"  This report includes both network visualization and regulator rankings")
    
    return output_file


def load_enrichment_data(enrichment_dir, method='auto'):
    """
    Load pathway enrichment results from either EnrichR or GSEA methods
    
    Parameters:
    enrichment_dir (str): Directory containing enrichment results
    method (str): 'EnrichR' or 'GSEA'
    
    Returns:
    tuple: (enrichment_data dict, metadata dict)
    """
    enrichment_data = {
        'bp': empty_enrichment_frame(),
        'mf': empty_enrichment_frame(),
        'cc': empty_enrichment_frame(),
        'kegg': empty_enrichment_frame(),
        'reactome': empty_enrichment_frame()
    }
    metadata = {
        'selected_source': None,
        'selected_granularity': None,
        'selected_files': []
    }

    try:
        if not os.path.isdir(enrichment_dir):
            return enrichment_data, metadata

        csv_files = []
        for root, dirs, files in os.walk(enrichment_dir):
            for f in files:
                if f.lower().endswith('.csv'):
                    csv_files.append(os.path.join(root, f))

        if csv_files:
            print(f"Found enrichment CSV files ({len(csv_files)}):")
            for p in csv_files:
                print(f"  - {p}")
        else:
            print(f"No enrichment CSV files found under: {enrichment_dir}")
        if not csv_files:
            return enrichment_data, metadata

        file_buckets = bucket_enrichment_files(csv_files)
        selected_source = choose_enrichment_source(file_buckets, method)
        if selected_source is None:
            print("No canonical EnrichR or GSEA enrichment files were found.")
            return enrichment_data, metadata

        source_bucket = file_buckets[selected_source]
        selected_granularity = 'all' if (source_bucket['all_up'] or source_bucket['all_down']) else 'top_10'
        selected_files = list(dict.fromkeys(
            source_bucket[f'{selected_granularity}_up'] + source_bucket[f'{selected_granularity}_down']
        ))

        metadata = {
            'selected_source': selected_source,
            'selected_granularity': selected_granularity,
            'selected_files': selected_files
        }

        print(f"Using {selected_source} enrichment files ({selected_granularity}) for module overview integration")

        if not selected_files:
            return enrichment_data, metadata

        loaded_frames = []
        for path in selected_files:
            try:
                direction = infer_direction_from_filename(os.path.basename(path))
                frame = pd.read_csv(path)
                standardized = standardize_enrichment_df(frame, direction=direction, source=selected_source)
                if not standardized.empty:
                    loaded_frames.append(standardized)
                    print(f"Loaded enrichment file: {os.path.basename(path)} ({len(standardized)} rows)")
            except Exception as exc:
                print(f"Warning: could not read enrichment file {path}: {exc}")

        if not loaded_frames:
            return enrichment_data, metadata

        combined = pd.concat(loaded_frames, ignore_index=True)
        for db_value, subframe in combined.groupby('Database', dropna=False):
            key = database_key(db_value)
            if key == 'other':
                continue
            enrichment_data[key] = pd.concat([
                enrichment_data[key],
                standardize_enrichment_df(subframe)
            ], ignore_index=True)

        for key in enrichment_data:
            enrichment_data[key] = standardize_enrichment_df(enrichment_data[key])

    except Exception as e:
        print(f"Error loading enrichment data: {e}")

    return enrichment_data, metadata


def load_pkn_data(pkn_file, metabolite_mapping_file=None):
    """
    Load PKN and metabolite mapping to build a lookup dictionary for edge categorization.
    
    Parameters:
    pkn_file (str): Path to PKN TSV file (Lemonite_PKN.tsv)
    metabolite_mapping_file (str): Path to name_map.csv (Query -> HMDB mapping)
    
    Returns:
    tuple: (pkn_lookup dict, name_to_hmdb dict)
    """
    pkn_lookup = {}
    name_to_hmdb = {}

    try:
        pkn_df = pd.read_csv(pkn_file, sep='\t', header=0)
        pkn_df['Node1'] = pkn_df['Node1'].astype(str).str.split('_').str[-1]
        pkn_df['Node2'] = pkn_df['Node2'].astype(str).str.split('_').str[-1]

        for _, r in pkn_df.iterrows():
            n1 = r['Node1']
            n2 = r['Node2']
            source = r.get('Source', '')
            edge_type = r.get('Type', '')

            # Classify edge
            if edge_type == 'metabolite-gene':
                causal_sources = ['LINCS', 'chEMBL']
                metabolic_sources = ['Human1_GEM_dist1', 'Human1_GEM_dist2']
                if source in causal_sources:
                    cat = 'Causal'
                elif source in metabolic_sources:
                    cat = 'Metabolic_pathway'
                else:
                    cat = 'Other'
            else:
                cat = 'PPI'

            pkn_lookup.setdefault((n1, n2), []).append(cat)
            pkn_lookup.setdefault((n2, n1), []).append(cat)

        print(f"Loaded PKN with {len(pkn_df)} edges, {len(pkn_lookup)} lookup entries")
    except Exception as e:
        print(f"Warning: Could not load PKN file: {e}")

    if metabolite_mapping_file and os.path.exists(metabolite_mapping_file):
        try:
            mapping_df = pd.read_csv(metabolite_mapping_file, sep=',')
            name_to_hmdb = mapping_df.set_index('Query')['HMDB'].dropna().to_dict()
            print(f"Loaded metabolite mapping with {len(name_to_hmdb)} entries")
        except Exception as e:
            print(f"Warning: Could not load metabolite mapping: {e}")

    return pkn_lookup, name_to_hmdb


def categorize_regulator_module_edge(regulator, module_id, pkn_lookup, name_to_hmdb, module_genes_map):
    """
    Determine the best interaction category between a metabolite regulator and a module's
    genes via PKN lookup. Priority: Causal > Metabolic_pathway > Other.
    """
    hmdb = name_to_hmdb.get(regulator, regulator)
    best = 'Other'
    for gene in module_genes_map.get(str(module_id), []):
        cats = pkn_lookup.get((hmdb, gene), [])
        if not cats:
            cats = pkn_lookup.get((gene, hmdb), [])
        for c in cats:
            if c == 'Causal':
                return 'Causal'
            elif c == 'Metabolic_pathway' and best != 'Causal':
                best = 'Metabolic_pathway'
    return best


def annotate_edges_with_category(edges, pkn_lookup, name_to_hmdb, module_genes_map):
    """
    Annotate each edge in the edge list with an interaction category
    (Causal / Metabolic_pathway / Other) based on PKN data.
    Only metabolite edges get PKN-based annotation; TF edges are always 'Other'.
    """
    _cache = {}
    for edge in edges:
        if edge['type'].startswith('metabolite') or edge['type'].startswith('Metabolite') or edge['type'].startswith('Lipid') or edge['type'].startswith('lipid'):
            key = (edge['source'], edge['target'])
            if key in _cache:
                edge['category'] = _cache[key]
            else:
                reg = edge['source']
                mod = edge['target'].replace('Module_', '')
                cat = categorize_regulator_module_edge(reg, mod, pkn_lookup, name_to_hmdb, module_genes_map)
                edge['category'] = cat
                _cache[key] = cat
        else:
            edge['category'] = 'Other'
    return edges


def build_enriched_hover_text(module_id, module_overview_df, enrichment_all_df, edges, name_lookup=None):
    """
    Build rich hover text for a module node including expression info,
    regulators (up to 10), and top 3 pathways per database with p-values.
    
    Parameters:
    module_id (str): Module ID
    module_overview_df (pd.DataFrame): Module overview DataFrame
    enrichment_all_df (pd.DataFrame or None): Combined enrichment DataFrame with Database column
    edges (list): Edge list for extracting connected regulators
    
    Returns:
    str: HTML-formatted hover text
    """
    hover_text = f"<b>Module {module_id}</b><br>"

    row = None
    if module_overview_df is not None and not module_overview_df.empty:
        matches = module_overview_df[module_overview_df['Module'].astype(str) == str(module_id)]
        if len(matches) > 0:
            row = matches.iloc[0]

    if row is not None:
        # Expression analysis
        expr_p = row.get('Expression_adjusted_pval', 'NA')
        if expr_p != 'NA' and pd.notna(expr_p):
            try:
                hover_text += f"<b>Expression Analysis:</b><br>"
                hover_text += f"  \u2022 adj. p-value: {float(expr_p):.2e}<br>"
                expr_rank = row.get('Expression_rank', 'NA')
                hover_text += f"  \u2022 Rank: {expr_rank}<br>"
                expr_sig = 'Yes' if isinstance(expr_p, (int, float)) and float(expr_p) < 0.05 else 'No'
                hover_text += f"  \u2022 Significant: {expr_sig}<br><br>"
            except (ValueError, TypeError):
                pass

        # Gene count
        genes_val = row.get('Module_genes', 'NA')
        if genes_val != 'NA' and pd.notna(genes_val):
            gene_count = len(str(genes_val).split('|'))
            hover_text += f"<b>Genes:</b> {gene_count} genes<br>"

        megago_cluster = row.get('MegaGO_Cluster', row.get('Functional_Cluster', 'NA'))
        if megago_cluster != 'NA' and pd.notna(megago_cluster):
            hover_text += f"<b>MegaGO Cluster:</b> {megago_cluster}<br>"

        megago_label = row.get('MegaGO_Label', 'NA')
        if megago_label != 'NA' and pd.notna(megago_label) and str(megago_label).strip() != '':
            hover_text += f"<b>MegaGO Label:</b> {megago_label}<br>"
            representative_go_id = row.get('MegaGO_Representative_GO_ID', 'NA')
            if representative_go_id != 'NA' and pd.notna(representative_go_id) and str(representative_go_id).strip() != '':
                hover_text += f"<b>Representative GO:</b> {representative_go_id}<br>"

        hover_text += '<br>'

    # Regulators per type from edges
    module_target = f"Module_{module_id}"
    reg_type_labels = {}
    for e in edges:
        if e['target'] == module_target:
            rtype = e.get('type', '').replace('_regulation', '').replace('_to_module', '')
            reg_type_labels.setdefault(rtype, set()).add(e['source'])

    for reg_type in sorted(reg_type_labels.keys()):
        regs = sorted(reg_type_labels[reg_type])
        if regs:
            display_regs = [name_lookup.get(r, r) for r in regs] if name_lookup else regs
            hover_text += f"<b>{reg_type.capitalize()} ({len(regs)}):</b> "
            hover_text += ', '.join(display_regs[:10])
            if len(regs) > 10:
                hover_text += f", ... (+{len(regs) - 10} more)"
            hover_text += '<br>'
    hover_text += '<br>'

    # Top enriched pathways from ALL databases (BP, MF, CC, KEGG, Reactome)
    if enrichment_all_df is not None and not enrichment_all_df.empty:
        mod_col = 'Module' if 'Module' in enrichment_all_df.columns else None
        db_col = next((c for c in enrichment_all_df.columns if c.lower() in ('database', 'db', 'source')), None)
        padj_col = 'p.adjust' if 'p.adjust' in enrichment_all_df.columns else None

        if mod_col and db_col and padj_col:
            mod_enrich = enrichment_all_df[enrichment_all_df[mod_col].astype(str) == str(module_id)]
            for db in ['BP', 'MF', 'CC', 'KEGG', 'Reactome']:
                db_enrich = mod_enrich[mod_enrich[db_col].str.upper() == db.upper()]
                if db_enrich.empty:
                    # Try case-insensitive partial match
                    db_enrich = mod_enrich[mod_enrich[db_col].str.lower().str.contains(db.lower(), na=False)]
                db_enrich = db_enrich.sort_values(padj_col).head(3)
                if len(db_enrich) > 0:
                    hover_text += f"<b>Top {db}:</b><br>"
                    for _, erow in db_enrich.iterrows():
                        term = str(erow.get('Term', ''))
                        if len(term) > 55:
                            term = term[:52] + '...'
                        p_val = erow[padj_col]
                        try:
                            hover_text += f"  \u2022 {term} (p={float(p_val):.1e})<br>"
                        except (ValueError, TypeError):
                            hover_text += f"  \u2022 {term}<br>"

    return hover_text


def adjust_positions_to_avoid_overlap(pos, min_dist=60):
    """
    Iteratively push module nodes apart until every pair is at least min_dist apart.
    Only adjusts nodes whose ID starts with 'Module_'.
    """
    import random
    moved = True
    max_iterations = 100
    iteration = 0
    while moved and iteration < max_iterations:
        moved = False
        iteration += 1
        items = list(pos.items())
        for i, (n1, (x1, y1)) in enumerate(items):
            if not n1.startswith('Module_'):
                continue
            for n2, (x2, y2) in items[i + 1:]:
                if not n2.startswith('Module_'):
                    continue
                dx = x2 - x1
                dy = y2 - y1
                dist = np.hypot(dx, dy)
                if dist == 0:
                    dx = random.uniform(-1, 1)
                    dy = random.uniform(-1, 1)
                    dist = np.hypot(dx, dy)
                if dist < min_dist:
                    shift = (min_dist - dist) / 2.0
                    ang = np.arctan2(dy, dx)
                    pos[n1] = (x1 - shift * np.cos(ang), y1 - shift * np.sin(ang))
                    pos[n2] = (x2 + shift * np.cos(ang), y2 + shift * np.sin(ang))
                    moved = True
    return pos


def export_cytoscape_files(nodes, edges, module_clusters, module_overview_df, output_dir):
    """
    Export Cytoscape-compatible edge list and node attribute TSV files.
    
    Parameters:
    nodes (list): Node dicts with 'id', 'type', etc.
    edges (list): Edge dicts with 'source', 'target', 'category'
    module_clusters (dict): module_id -> 'Cluster_N' string
    module_overview_df (pd.DataFrame): Module overview for attribute lookup
    output_dir (str): Output directory
    """
    print("\nExporting Cytoscape-compatible files...")

    # Edge file
    edge_file = os.path.join(output_dir, 'module_network_edges.txt')
    with open(edge_file, 'w') as f:
        f.write("source\ttarget\tCategory\tArrowShape\n")
        for edge in edges:
            cat = edge.get('category', 'Other')
            if cat == 'Causal':
                arrow = 'DELTA'
            elif cat == 'Metabolic_pathway':
                arrow = 'DOT'
            else:
                arrow = ''
            f.write(f"{edge['source']}\t{edge['target']}\t{cat}\t{arrow}\n")
    print(f"  Saved edge list: {edge_file}")

    # Node attributes file
    node_file = os.path.join(output_dir, 'module_network_node_attributes.txt')
    with open(node_file, 'w') as f:
        header = "Node\tMegaGO_Cluster\tMegaGO_Label\tMegaGO_Representative_GO_ID\tMegaGO_Label_Source\tNode_Type\tExpression_significant\tPPI_significant\tModule_genes_count\n"
        f.write(header)
        for node in nodes:
            node_id = node['id']
            node_type = node['type']
            if node_type == 'module':
                module_num = node_id.replace('Module_', '')
                cluster = module_clusters.get(module_num, 'Unassigned')
                megago_label = ''
                representative_go_id = ''
                label_source = ''

                expr_flag = ''
                ppi_flag = ''
                gene_count = ''
                if module_overview_df is not None and not module_overview_df.empty:
                    mod_row = module_overview_df[module_overview_df['Module'].astype(str) == module_num]
                    if len(mod_row) > 0:
                        try:
                            p_adj = mod_row.iloc[0].get('Expression_adjusted_pval', 'NA')
                            if p_adj != 'NA' and pd.notna(p_adj):
                                expr_flag = 'Yes' if float(p_adj) < 0.05 else 'No'
                        except Exception:
                            pass
                        try:
                            ppi_val = mod_row.iloc[0].get('PPI_FDR', 'NA')
                            if ppi_val != 'NA' and pd.notna(ppi_val):
                                ppi_flag = 'Yes' if float(ppi_val) < 0.05 else 'No'
                        except Exception:
                            pass
                        try:
                            megago_label = mod_row.iloc[0].get('MegaGO_Label', '')
                            representative_go_id = mod_row.iloc[0].get('MegaGO_Representative_GO_ID', '')
                            label_source = mod_row.iloc[0].get('MegaGO_Label_Source', '')
                        except Exception:
                            pass
                        try:
                            genes = mod_row.iloc[0].get('Module_genes', '')
                            if genes != 'NA' and pd.notna(genes):
                                gene_count = str(len(str(genes).split('|')))
                        except Exception:
                            pass
                f.write(f"{node_id}\t{cluster}\t{megago_label}\t{representative_go_id}\t{label_source}\t{node_type}\t{expr_flag}\t{ppi_flag}\t{gene_count}\n")
            else:
                f.write(f"{node_id}\t\t\t\t\t{node_type}\t\t\t\n")
    print(f"  Saved node attributes: {node_file}")


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

def create_megago_files(bp_terms_df, output_dir):
    """
    Create MegaGO files with canonical BP top_30 terms for each module.
    
    Parameters:
    bp_terms_df (pd.DataFrame): Canonical BP top_30 terms per module
    output_dir (str): Output directory for megaGO files
    
    Returns:
    str: Directory path where megaGO files are created
    """
    megago_dir = os.path.join(output_dir, 'megaGO_files')
    os.makedirs(megago_dir, exist_ok=True)
    
    files_created = 0
    if bp_terms_df is None or bp_terms_df.empty:
        print("Warning: No canonical BP top_30 terms available - cannot create MegaGO files")
        return None

    print("Creating MegaGO files with canonical BP top_30 terms for each module...")

    for module, module_data in bp_terms_df.groupby('Module', sort=True):
        cleaned_terms = []
        for term in module_data['Term'].tolist():
            if not isinstance(term, str):
                continue

            go_id = extract_go_id(term)
            if go_id:
                cleaned_terms.append(go_id)
                continue

            fallback = re.split(r' - |\(|;|,', term)[0].strip()
            if fallback:
                cleaned_terms.append(fallback)

        cleaned_terms = list(dict.fromkeys(cleaned_terms))
        if cleaned_terms:
            filepath = os.path.join(megago_dir, f"{module}_BP_terms.txt")
            with open(filepath, 'w') as handle:
                handle.write('GO_TERM\n')
                handle.write("\n".join(cleaned_terms))
            files_created += 1

    print(f"Created {files_created} GO BP term files for MegaGO clustering")
    
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

def megago_cluster_modules(module_pathways, bp_terms_df, n_clusters=5, output_dir='.'):
    """
    Cluster modules using megaGO on canonical BP top_30 terms when available,
    otherwise fall back to pathway similarity.

    Returns cluster labels as 'Cluster_N' strings for consistency with downstream code.

    Parameters:
    module_pathways (dict): Dictionary mapping module IDs to lists of pathway terms
    n_clusters (int): Number of clusters to create
    bp_terms_df (pd.DataFrame): Canonical BP top_30 terms used for MegaGO and rrvgo
    output_dir (str): Output directory for temporary files

    Returns:
    tuple: (module_clusters dict, similarity_matrix)
           module_clusters maps module_id -> 'Cluster_N' string
    """
    def _to_cluster_str(raw_clusters):
        """Convert integer cluster labels to 'Cluster_N' strings."""
        return {mid: f"Cluster_{int(lbl)}" for mid, lbl in raw_clusters.items()}

    if not module_pathways:
        print("No pathway data available for clustering")
        return {}, None

    modules = list(module_pathways.keys())

    # Skip clustering if requested or not enough clusters
    if n_clusters <= 1:
        print("Functional clustering disabled - assigning all modules to single cluster")
        return {mod: 'Cluster_1' for mod in modules}, None

    print("\nAttempting canonical MegaGO clustering...")
    megago_binary = shutil.which('megago')
    if megago_binary and bp_terms_df is not None and not bp_terms_df.empty:
        print(f"   Found canonical BP top_30 terms with {len(bp_terms_df)} rows")
        megago_dir = create_megago_files(bp_terms_df, output_dir)

        if megago_dir:
            print(f"   Created MegaGO files in: {megago_dir}")
            similarity_matrix, megago_module_ids = run_megago_clustering(megago_dir)

            if similarity_matrix is not None:
                print("Successfully obtained MegaGO similarity matrix")

                distance_matrix = 1 - similarity_matrix
                if len(megago_module_ids) >= 2:
                    condensed_dist = squareform(distance_matrix, checks=False)
                    linkage_matrix = linkage(condensed_dist, method='ward')

                    if n_clusters == 'auto':
                        n_clusters = min(5, max(2, len(megago_module_ids) // 3))

                    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

                    raw = {}
                    for i, module_id in enumerate(megago_module_ids):
                        raw[module_id] = cluster_labels[i]
                    for module in modules:
                        if module not in raw:
                            raw[module] = 0

                    module_clusters = _to_cluster_str(raw)
                    print(f"Created {n_clusters} module clusters using MegaGO semantic similarity")

                    cluster_counts = Counter(module_clusters.values())
                    for cid, count in sorted(cluster_counts.items()):
                        print(f"  {cid}: {count} modules")
                    return module_clusters, similarity_matrix

            print("MegaGO clustering failed - falling back to pathway similarity")
        else:
            print("Could not create MegaGO files - falling back to pathway similarity")
    elif bp_terms_df is None or bp_terms_df.empty:
        print("No canonical BP top_30 terms available - falling back to pathway similarity")
    else:
        print("MegaGO command-line tool not available - falling back to pathway similarity")

    # Fall back to enhanced pathway similarity clustering
    print("Using pathway similarity for clustering...")

    jaccard_sim = calculate_pathway_similarity_matrix(module_pathways)
    similarity_matrix = jaccard_sim.values
    valid_modules = list(module_pathways.keys())

    print("Using enhanced pathway similarity clustering (Jaccard)")

    distance_matrix = 1 - similarity_matrix

    if len(valid_modules) < 2:
        return {mod: 'Cluster_0' for mod in modules}, None

    condensed_dist = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_dist, method='ward')

    if n_clusters == 'auto':
        n_clusters = min(5, max(2, len(valid_modules) // 3))

    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

    raw = {}
    for i, module in enumerate(valid_modules):
        raw[module] = cluster_labels[i]
    for module in modules:
        if module not in raw:
            raw[module] = 0

    module_clusters = _to_cluster_str(raw)
    print(f"Created {n_clusters} module clusters based on functional similarity")

    cluster_counts = Counter(module_clusters.values())
    for cid, count in sorted(cluster_counts.items()):
        print(f"  {cid}: {count} modules")

    return module_clusters, similarity_matrix

def create_interactive_network_visualization(module_data, module_clusters, output_dir,
                                              enrichment_data=None, go_similarity_matrix=None,
                                              module_overview_df=None, enrichment_all_df=None,
                                              cluster_label_lookup=None, pkn_lookup=None,
                                              name_to_hmdb=None, module_genes_map=None,
                                              name_lookup=None, run_id=None):
    """
    Create interactive network visualization using Plotly with MegaGO cluster coloring,
    cluster labels, and PKN-based edge categorization.

    Generates two HTML files:
    - interactive_module_network.html (standard)
    - interactive_module_network_movable.html (draggable nodes)
    
    Parameters:
    module_data (list): List of module data dictionaries
    module_clusters (dict): module_id -> 'Cluster_N' string
    output_dir (str): Output directory
    enrichment_data (dict): Old per-database enrichment dict (fallback for hover)
    go_similarity_matrix: similarity matrix for layout
    module_overview_df (pd.DataFrame): Module overview for hover text
    enrichment_all_df (pd.DataFrame): Combined enrichment data with Database column
    cluster_label_lookup (dict or None): MegaGO cluster -> representative label
    pkn_lookup (dict or None): PKN edge lookup dict
    name_to_hmdb (dict or None): metabolite name -> HMDB mapping
    module_genes_map (dict or None): module_id -> gene list
    """
    print("Creating interactive network visualization with MegaGO cluster coloring...")
    
    # -- Build nodes and edges --
    nodes = []
    edges = []
    regulator_modules = {}
    
    for module_info in module_data:
        module_id = str(module_info['Module'])
        
        n_genes = len(module_info.get('Module_genes', '').split('|')) if module_info.get('Module_genes', '') != 'NA' else 0
        
        module_node = {
            'id': f"Module_{module_id}",
            'label': f"M{module_id}",
            'type': 'module',
            'n_genes': n_genes,
        }
        nodes.append(module_node)
        
        # Collect regulators dynamically for all regulator types
        regulator_columns = [col for col in module_info.keys() if col.endswith('_regulators')]
        for reg_col in regulator_columns:
            if module_info.get(reg_col, 'NA') != 'NA':
                reg_type = reg_col.replace('_regulators', '')
                regs = module_info[reg_col].split('|')
                
                if reg_type not in regulator_modules:
                    regulator_modules[reg_type] = {}
                
                for reg in regs:
                    reg = reg.strip()
                    if reg:
                        if reg not in regulator_modules[reg_type]:
                            regulator_modules[reg_type][reg] = []
                        regulator_modules[reg_type][reg].append(module_id)
    
    # Create regulator nodes and edges
    for reg_type, regulators in regulator_modules.items():
        for regulator, target_modules in regulators.items():
            modules_str = ', '.join(sorted(target_modules))
            display_name = name_lookup.get(regulator, regulator) if name_lookup else regulator
            regulator_node = {
                'id': regulator,
                'label': display_name[:15] if len(display_name) > 15 else display_name,
                'type': reg_type,
                'hover_info': f"<b>{display_name}</b><br>Type: {reg_type}<br>Targets ({len(target_modules)}): M{', M'.join(target_modules)}"
            }
            nodes.append(regulator_node)
            
            for module_id in target_modules:
                edge = {
                    'source': regulator,
                    'target': f"Module_{module_id}",
                    'type': f"{reg_type}_regulation"
                }
                edges.append(edge)
    
    # Annotate edges with PKN-based categories
    if pkn_lookup and module_genes_map:
        edges = annotate_edges_with_category(edges, pkn_lookup, name_to_hmdb or {}, module_genes_map)
    else:
        for edge in edges:
            edge['category'] = 'Other'
    
    # Build enriched hover info for modules
    for node in nodes:
        if node['type'] == 'module':
            module_id = node['id'].replace('Module_', '')
            node['hover_info'] = build_enriched_hover_text(
                module_id, module_overview_df, enrichment_all_df, edges, name_lookup=name_lookup
            )
    
    # -- Create layout positions --
    pos, cluster_positions = _create_cluster_layout(nodes, edges, module_clusters)
    
    # Apply overlap avoidance
    pos = adjust_positions_to_avoid_overlap(pos, min_dist=1.5)
    
    # -- Build color maps --
    all_clusters = sorted(set(module_clusters.values()))
    cluster_palette = ['#EE6677', '#4477AA', '#228833', '#AA3377', '#66CCEE',
                       '#CCBB44', '#EE99AA', '#44BB99', '#BBCC33', '#AAAA00']
    cluster_color_map = {cl: cluster_palette[i % len(cluster_palette)]
                         for i, cl in enumerate(all_clusters)}
    cluster_color_map['Unassigned'] = '#BBBBBB'

    regulator_color = '#FF8C00'

    # -- Edge styles by interaction category --
    edge_styles = {
        'Causal': {'color': 'rgba(80,80,80,0.9)', 'width': 3.0, 'label': 'Causal'},
        'Metabolic_pathway': {'color': 'rgba(80,80,80,0.85)', 'width': 2.5, 'label': 'Metabolic pathway'},
        'Other': {'color': 'rgba(80,80,80,0.6)', 'width': 1.5, 'label': 'Other'},
    }

    # -- Helper to build figures for both standard and movable variants --
    def _build_figure(movable=False):
        # Build edge traces per category
        temp_coords = {cat: {'x': [], 'y': []} for cat in edge_styles}
        for edge in edges:
            if edge['source'] in pos and edge['target'] in pos:
                x0, y0 = pos[edge['source']]
                x1, y1 = pos[edge['target']]
                cat = edge.get('category', 'Other')
                if cat not in temp_coords:
                    cat = 'Other'
                temp_coords[cat]['x'].extend([x0, x1, None])
                temp_coords[cat]['y'].extend([y0, y1, None])

        fig = go.Figure()

        # Add edge traces (Other first so thicker edges overlay)
        for cat in ['Other', 'Metabolic_pathway', 'Causal']:
            style = edge_styles[cat]
            fig.add_trace(go.Scatter(
                x=temp_coords[cat]['x'], y=temp_coords[cat]['y'],
                line=dict(width=style['width'], color=style['color']),
                hoverinfo='none', mode='lines',
                name=style['label'], showlegend=True
            ))

        # Add edge decorations (arrows and dots)
        for edge in edges:
            if edge['source'] in pos and edge['target'] in pos:
                x0, y0 = pos[edge['source']]
                x1, y1 = pos[edge['target']]
                cat = edge.get('category', 'Other')
                if cat == 'Causal':
                    fig.add_annotation(
                        x=x1, y=y1, ax=x0, ay=y0,
                        xref='x', yref='y', axref='x', ayref='y',
                        showarrow=True, arrowhead=3, arrowsize=1,
                        arrowwidth=edge_styles[cat]['width'],
                        arrowcolor=edge_styles[cat]['color'],
                        opacity=0.8
                    )
                elif cat == 'Metabolic_pathway':
                    fig.add_trace(go.Scatter(
                        x=[x1], y=[y1], mode='markers',
                        marker=dict(symbol='circle', size=6,
                                    color=edge_styles[cat]['color'], line=dict(width=0)),
                        hoverinfo='none', showlegend=False
                    ))

        if movable:
            # One trace per node for individual dragging
            for n in nodes:
                if n['id'] not in pos:
                    continue
                x, y = pos[n['id']]
                label = n.get('label', n['id'])
                txt = label
                symbol = 'circle' if n['type'] == 'module' else 'diamond-wide'
                if n['type'] == 'module':
                    m_id = n['id'].replace('Module_', '')
                    color = cluster_color_map.get(module_clusters.get(m_id, 'Unassigned'), '#BBBBBB')
                    size = 50
                    border = '#FFFFFF'
                    text_size = 12
                else:
                    color = regulator_color
                    size = min(80, max(25, len(txt) * 4))
                    border = 'darkorange'
                    text_size = 8

                fig.add_trace(go.Scatter(
                    x=[x], y=[y], mode='markers+text',
                    marker=dict(size=size, color=color,
                                line=dict(width=3, color=border),
                                symbol=symbol),
                    text=[txt], textposition='middle center',
                    textfont=dict(size=text_size, color='black', family='Arial Black'),
                    hovertext=[n.get('hover_info', txt)], hoverinfo='text',
                    name=n['id'], showlegend=False
                ))
        else:
            # -- Grouped module traces by cluster --
            for cluster in all_clusters:
                cat_nodes = [n for n in nodes if n['type'] == 'module' and
                             module_clusters.get(n['id'].replace('Module_', ''), 'Unassigned') == cluster
                             and n['id'] in pos]
                if not cat_nodes:
                    continue
                x_vals = [pos[n['id']][0] for n in cat_nodes]
                y_vals = [pos[n['id']][1] for n in cat_nodes]
                hover_texts = [n.get('hover_info', n['label']) for n in cat_nodes]
                labels = [n['id'].replace('Module_', '') for n in cat_nodes]

                line_colours = ['white'] * len(cat_nodes)

                cluster_label = (cluster_label_lookup or {}).get(cluster, '')
                trace_name = f"{cluster}: {cluster_label} ({len(cat_nodes)})" if cluster_label else f"{cluster} ({len(cat_nodes)})"

                trace = go.Scatter(
                    x=x_vals, y=y_vals, mode='markers+text',
                    marker=dict(size=50, color=cluster_color_map.get(cluster, '#CCCCCC'),
                                line=dict(width=3, color=line_colours), symbol='circle'),
                    text=labels, textposition='middle center',
                    textfont=dict(size=12, color='black', family='Arial Black'),
                    hovertext=hover_texts, hoverinfo='text',
                    name=trace_name, legendgroup=cluster
                )
                fig.add_trace(trace)

            # Unassigned modules
            unassigned_nodes = [n for n in nodes if n['type'] == 'module' and
                                module_clusters.get(n['id'].replace('Module_', ''), 'Unassigned') == 'Unassigned'
                                and n['id'] in pos]
            if unassigned_nodes:
                x_vals = [pos[n['id']][0] for n in unassigned_nodes]
                y_vals = [pos[n['id']][1] for n in unassigned_nodes]
                hover_texts = [n.get('hover_info', n['label']) for n in unassigned_nodes]
                labels = [n['id'].replace('Module_', '') for n in unassigned_nodes]
                fig.add_trace(go.Scatter(
                    x=x_vals, y=y_vals, mode='markers+text',
                    marker=dict(size=50, color='#BBBBBB',
                                line=dict(width=3, color=['white'] * len(x_vals)), symbol='circle'),
                    text=labels, textposition='middle center',
                    textfont=dict(size=12, color='black', family='Arial Black'),
                    hovertext=hover_texts, hoverinfo='text',
                    name=f'Unassigned ({len(unassigned_nodes)})', legendgroup='Unassigned'
                ))

            # -- Regulator traces grouped by type --
            for reg_type in sorted(regulator_modules.keys()):
                reg_nodes = [n for n in nodes if n['type'] == reg_type and n['id'] in pos]
                if not reg_nodes:
                    continue
                x_vals = [pos[n['id']][0] for n in reg_nodes]
                y_vals = [pos[n['id']][1] for n in reg_nodes]
                hover_texts = [n.get('hover_info', n['label']) for n in reg_nodes]
                labels = [n['label'] for n in reg_nodes]
                sizes = [min(80, max(25, len(label) * 4)) for label in labels]
                fig.add_trace(go.Scatter(
                    x=x_vals, y=y_vals, mode='markers+text',
                    marker=dict(size=sizes, color=regulator_color, symbol='diamond-wide',
                                line=dict(width=1, color='darkorange')),
                    text=labels, textposition='middle center',
                    textfont=dict(size=8, color='black'),
                    hovertext=hover_texts, hoverinfo='text',
                    name=f"{reg_type.capitalize()} ({len(reg_nodes)})", legendgroup=f'reg_{reg_type}'
                ))

        # -- Cluster label annotations --
        annotations = []
        for cl, (cx, cy) in cluster_positions.items():
            cluster_label = (cluster_label_lookup or {}).get(cl, '')
            display_label = cluster_label if len(cluster_label) <= 50 else cluster_label[:47] + '...'
            annotation_text = f'<b>{cl}</b>'
            if display_label:
                annotation_text += f'<br><span style="font-size:10px">{display_label}</span>'
            annotations.append(dict(
                x=cx * 1.3, y=cy * 1.3,
                text=annotation_text, showarrow=False,
                font=dict(size=14, color=cluster_color_map.get(cl, '#000000')),
                bgcolor='rgba(255,255,255,0.8)', borderpad=4
            ))

        fig.update_layout(
            title=dict(
                text='Module-Regulator Network<br><sub>Modules colored by canonical MegaGO clusters</sub>',
                x=0.5, font=dict(size=20)),
            showlegend=True, hovermode='closest',
            width=1400, height=1000,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            annotations=annotations,
            legend=dict(x=1.02, y=1, bgcolor='rgba(255,255,255,0.9)',
                        bordercolor='black', borderwidth=1, font=dict(size=10))
        )
        return fig

    # -- Generate Cytoscape.JS movable network --
    fig_standard = _build_figure(movable=False)
    output_movable = os.path.join(output_dir, 'interactive_module_network_movable.html')
    create_cytoscape_html_network(
        nodes, edges, module_clusters,
        config_name=run_id or 'Lemonite',
        variant_name='filtered',
        output_file=output_movable,
        module_overview=module_overview_df,
        module_enrichment_all=enrichment_all_df,
        name_lookup=name_lookup,
    )

    return fig_standard

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

def _create_cluster_layout(nodes, edges, module_clusters):
    """
    Position modules in a circle-per-cluster arrangement and place regulators
    near their connected modules with type-based angular offsets.

    Parameters:
        nodes (list): node dicts with 'id' and 'type'
        edges (list): edge dicts with 'source' and 'target'
        module_clusters (dict): module_id (str) -> 'Cluster_N' string

    Returns:
        (pos, cluster_positions):
            pos  – dict  node_id -> (x, y)
            cluster_positions – dict  cluster_name -> (cx, cy) centre
    """
    print("Creating cluster-based layout...")

    module_nodes = [n for n in nodes if n['type'] == 'module']

    # -- Determine unique clusters and their centres --
    unique_clusters = sorted(set(module_clusters.values()))
    n_clusters = len(unique_clusters)
    radius = 8

    if n_clusters <= 1:
        cluster_positions = {unique_clusters[0] if unique_clusters else 'Cluster_1': (0, 0)}
    elif n_clusters <= 6:
        cluster_positions = {
            c: (radius * np.cos(2 * np.pi * i / n_clusters),
                radius * np.sin(2 * np.pi * i / n_clusters))
            for i, c in enumerate(unique_clusters)
        }
    else:
        # Two concentric rings
        inner = n_clusters // 2
        outer = n_clusters - inner
        cluster_positions = {}
        for i, c in enumerate(unique_clusters[:inner]):
            angle = 2 * np.pi * i / inner
            cluster_positions[c] = (radius * 0.6 * np.cos(angle),
                                    radius * 0.6 * np.sin(angle))
        for i, c in enumerate(unique_clusters[inner:]):
            angle = 2 * np.pi * i / outer + np.pi / outer
            cluster_positions[c] = (radius * 1.2 * np.cos(angle),
                                    radius * 1.2 * np.sin(angle))

    np.random.seed(42)

    # -- Position module nodes around their cluster centres --
    pos = {}
    for node in module_nodes:
        mid = node['id'].replace('Module_', '')
        cluster = module_clusters.get(mid, unique_clusters[0] if unique_clusters else 'Cluster_1')
        cx, cy = cluster_positions.get(cluster, (0, 0))
        jx = np.random.normal(0, 1.2)
        jy = np.random.normal(0, 1.2)
        pos[node['id']] = (cx + jx, cy + jy)

    # -- Position regulators near connected modules with type offset --
    for node in nodes:
        if node['type'] == 'module' or node['id'] in pos:
            continue
        connected = [e['target'] for e in edges if e['source'] == node['id'] and e['target'].startswith('Module_')]
        connected += [e['source'] for e in edges if e['target'] == node['id'] and e['source'].startswith('Module_')]
        connected_pos = [pos[m] for m in connected if m in pos]

        if connected_pos:
            avg_x = np.mean([p[0] for p in connected_pos])
            avg_y = np.mean([p[1] for p in connected_pos])
            offset_dist = 2.5
            # Type-based angular offset
            if 'TF' in node['type']:
                base_angle = 0
            elif 'etabolite' in node['type'] or 'Metabolite' in node['type']:
                base_angle = 2 * np.pi / 3
            else:
                base_angle = 4 * np.pi / 3
            angle = base_angle + np.random.normal(0, 0.6)
            pos[node['id']] = (avg_x + offset_dist * np.cos(angle),
                                avg_y + offset_dist * np.sin(angle))
        else:
            pos[node['id']] = (np.random.normal(0, 4), np.random.normal(0, 4))

    print(f"Layout: {len(pos)} nodes positioned in {n_clusters} cluster(s)")
    return pos, cluster_positions

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


def _find_coherence_file(input_dir):
    """Search multiple possible locations for Module_coherence_scores.txt."""
    candidates = [
        os.path.join(input_dir, 'Module_coherence_scores.txt'),
        os.path.join(input_dir, 'Networks', 'Module_coherence_scores.txt'),
        'Module_coherence_scores.txt',  # current directory (Nextflow flat staging)
    ]
    # Also search recursively
    for root, dirs, files in os.walk(input_dir):
        if 'Module_coherence_scores.txt' in files:
            candidates.append(os.path.join(root, 'Module_coherence_scores.txt'))
    for candidate in candidates:
        if os.path.exists(candidate):
            print(f"Found coherence scores file: {candidate}")
            return candidate
    # Return default path (will trigger FileNotFoundError in caller)
    print(f"Warning: Module_coherence_scores.txt not found in any searched location")
    return os.path.join(input_dir, 'Module_coherence_scores.txt')


def main():
    parser = argparse.ArgumentParser(description='Generate interactive module overview with functional clustering')
    parser.add_argument('--input_dir', type=str, required=True,
                       help='Input directory containing module files')
    parser.add_argument('--output_dir', type=str, default='.',
                       help='Output directory for results')
    parser.add_argument('--regulator_files', type=str, required=False, default='',
                       help='Comma-separated list of regulator files (format: Type:Path,Type:Path)')
    parser.add_argument('--regulator_score_files', type=str, required=False, default='',
                       help='Comma-separated list of regulator score files (format: Type:Path,Type:Path)')
    parser.add_argument('--enrichment_method', type=str, default='auto',
                       choices=['EnrichR', 'GSEA', 'both', 'auto'],
                       help='Enrichment analysis method to use')
    parser.add_argument('--n_clusters', type=int, default=5,
                       help='Number of functional clusters to create')
    parser.add_argument('--organism', type=str, default='human',
                       help='Organism for canonical rrvgo labeling (human/mouse)')
    parser.add_argument('--pkn_file', type=str, default=None,
                       help='Path to PKN file (Lemonite_PKN.tsv) for edge categorization')
    parser.add_argument('--metabolite_mapping', type=str, default=None,
                       help='Path to metabolite name mapping file (name_map.csv) for HMDB resolution')
    parser.add_argument('--name_mapping', type=str, default=None,
                       help='Path to name_mapping.tsv (cleaned->original name mapping for restoring original feature names)')
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
    parser.add_argument('--run_id', type=str, default='Lemonite',
                       help='Run identifier used as network title')
    
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
    print(f"Organism: {args.organism}")
    print(f"Coherence threshold: {args.coherence_threshold}")
    print(f"Expression prioritization: {args.prioritize_by_expression}")
    print(f"Clustering workflow: canonical MegaGO {CANONICAL_CLUSTER_COLUMN} + rrvgo")
    if args.pkn_file:
        print(f"PKN file: {args.pkn_file}")
    if args.metabolite_mapping:
        print(f"Metabolite mapping: {args.metabolite_mapping}")
    if args.name_mapping:
        print(f"Name mapping: {args.name_mapping}")
    if args.prioritize_by_expression:
        print(f"Grouping column: {args.group_column}")
    
    # Load module data
    print("\nLoading module data...")
    
    # Load name mapping for restoring original feature names
    name_lookup = {}
    if args.name_mapping and os.path.exists(args.name_mapping):
        try:
            nm_df = pd.read_csv(args.name_mapping, sep='\t')
            if 'cleaned' in nm_df.columns and 'original' in nm_df.columns:
                name_lookup = dict(zip(nm_df['cleaned'], nm_df['original']))
                print(f"Loaded {len(name_lookup)} name mappings for original name restoration")
        except Exception as e:
            print(f"Warning: Could not load name mapping: {e}")
    
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
    
    # Parse and load regulator score files
    regulator_scores_dict = {}
    if args.regulator_score_files:
        for score_config in args.regulator_score_files.split(','):
            if ':' in score_config:
                reg_type, score_path = score_config.split(':', 1)
                print(f"Reading {reg_type} regulator scores from {score_path.strip()}")
                regulator_scores_dict[reg_type] = load_regulator_scores(score_path.strip())
    else:
        # Auto-discover score files from input_dir
        for reg_type in regulator_configs:
            score_file = os.path.join(args.input_dir, f'{reg_type}.selected_regulators_scores.txt')
            if os.path.exists(score_file):
                print(f"Auto-discovered {reg_type} regulator scores: {score_file}")
                regulator_scores_dict[reg_type] = load_regulator_scores(score_file)
    
    if regulator_scores_dict:
        for reg_type, df in regulator_scores_dict.items():
            print(f"  {reg_type}: {len(df)} regulator-module pairs loaded")
    
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
    
    enrichment_data, enrichment_metadata = load_enrichment_data(enrichment_dir, args.enrichment_method)

    # Add direction column to enrichment data for hover information
    for key in enrichment_data:
        if not enrichment_data[key].empty and '__direction__' not in enrichment_data[key].columns:
            # If no direction column, assume all are up-regulated for compatibility
            enrichment_data[key]['__direction__'] = 'Up'
    
    # Load PPI enrichment data if available
    print("\nLoading PPI enrichment data...")
    ppi_enrichment_data = {}
    moduleviewer_dir = os.path.join(args.input_dir, 'ModuleViewer_files')
    
    # Search multiple locations for PPI enrichment file (Nextflow may stage files flat)
    ppi_enrichment_candidates = [
        os.path.join(moduleviewer_dir, 'PPI_enrichment_results.csv'),
        os.path.join(args.input_dir, 'PPI_enrichment_results.csv'),  # flat staging
        'PPI_enrichment_results.csv',  # current directory
    ]
    # Also search recursively in input_dir
    for root, dirs, files in os.walk(args.input_dir):
        if 'PPI_enrichment_results.csv' in files:
            ppi_enrichment_candidates.append(os.path.join(root, 'PPI_enrichment_results.csv'))
    
    ppi_enrichment_file = None
    for candidate in ppi_enrichment_candidates:
        if os.path.exists(candidate):
            ppi_enrichment_file = candidate
            break
    
    if ppi_enrichment_file:
        try:
            ppi_df = pd.read_csv(ppi_enrichment_file)
            print(f"✓ Loaded PPI enrichment data from {ppi_enrichment_file}: {len(ppi_df)} modules")
            
            # Create dictionary for easy lookup - convert to int first to handle float module IDs
            for _, row in ppi_df.iterrows():
                module_id = str(int(row['Module']))
                ppi_enrichment_data[module_id] = {
                    'PPI_FDR': row['FDR'],
                    'PPI_fold_enrichment': row['Fold_enrichment']
                }
            
            # Summary statistics
            n_significant = (ppi_df['FDR'] < 0.05).sum()
            print(f"  - Modules significantly enriched for PPIs (FDR < 0.05): {n_significant}")
            print(f"  - Mean fold enrichment: {ppi_df['Fold_enrichment'].mean():.2f}")
        except Exception as e:
            print(f"Warning: Could not load PPI enrichment data: {e}")
            ppi_enrichment_data = {}
    else:
        print(f"PPI enrichment file not found in any searched location:")
        for c in ppi_enrichment_candidates[:3]:
            print(f"  - {c}")
        print("Continuing without PPI enrichment data...")
        ppi_enrichment_data = {}
    
    # Load metabolite-gene interaction enrichment data if available
    print("\nLoading metabolite-gene interaction enrichment data...")
    metgene_enrichment_data = {}
    
    metgene_enrichment_candidates = [
        os.path.join(moduleviewer_dir, 'Metabolite_Gene_enrichment_results.csv'),
        os.path.join(args.input_dir, 'Metabolite_Gene_enrichment_results.csv'),  # flat staging
        'Metabolite_Gene_enrichment_results.csv',  # current directory
    ]
    for root, dirs, files in os.walk(args.input_dir):
        if 'Metabolite_Gene_enrichment_results.csv' in files:
            metgene_enrichment_candidates.append(os.path.join(root, 'Metabolite_Gene_enrichment_results.csv'))
    
    metgene_enrichment_file = None
    for candidate in metgene_enrichment_candidates:
        if os.path.exists(candidate):
            metgene_enrichment_file = candidate
            break
    
    if metgene_enrichment_file:
        try:
            metgene_df = pd.read_csv(metgene_enrichment_file)
            print(f"✓ Loaded metabolite-gene enrichment data from {metgene_enrichment_file}: {len(metgene_df)} modules")
            
            for _, row in metgene_df.iterrows():
                module_id = str(int(row['Module']))
                metgene_enrichment_data[module_id] = {
                    'MetGene_FDR': row['FDR'],
                    'MetGene_fold_enrichment': row['Fold_enrichment'],
                    'MetGene_N_interactions': int(row['N_interactions_observed'])
                }
            
            n_significant = (metgene_df['FDR'] < 0.05).sum()
            print(f"  - Modules significantly enriched (FDR < 0.05): {n_significant}")
            print(f"  - Mean fold enrichment: {metgene_df['Fold_enrichment'].mean():.2f}")
        except Exception as e:
            print(f"Warning: Could not load metabolite-gene enrichment data: {e}")
            metgene_enrichment_data = {}
    else:
        print("Metabolite-gene enrichment file not found in any searched location")
        print("Continuing without metabolite-gene enrichment data...")
        metgene_enrichment_data = {}
    
    # Load coherence filtering - prioritize Networks/specific_modules.txt over coherence scores file
    
    print("\nLoading coherence-filtered modules...")    # First, try to load from Networks/specific_modules.txt (created by lemontree_to_network.py)
    filtered_modules_from_networks = load_filtered_modules_from_networks(args.input_dir)
    
    if filtered_modules_from_networks is not None:
        # Use the pre-filtered modules from Networks/specific_modules.txt
        modules_to_process = filtered_modules_from_networks.intersection(set(module_genes.keys()))
        print(f"Using coherence-filtered modules from Networks/specific_modules.txt")
        print(f"   Processing {len(modules_to_process)} modules that passed coherence filtering")
        
        # Load coherence scores file for informational purposes only
        # Search multiple possible locations for the coherence scores file
        coherence_file = _find_coherence_file(args.input_dir)
        _, coherence_df = load_coherence_filtered_modules(coherence_file, args.coherence_threshold)
        
    else:
        # Fallback to coherence scores file if Networks/specific_modules.txt is not available
        print("Networks/specific_modules.txt not found, falling back to coherence scores file")
        coherence_file = _find_coherence_file(args.input_dir)
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
        
        # Add regulator data dynamically, restoring original names via name_lookup
        for reg_type, regulators in regulators_dict.items():
            column_name = f'{reg_type}_regulators'
            raw_regs = regulators.get(module_id, [])
            if raw_regs:
                if name_lookup:
                    module_entry[column_name] = '|'.join(name_lookup.get(r, r) for r in raw_regs)
                else:
                    module_entry[column_name] = '|'.join(raw_regs)
            else:
                module_entry[column_name] = 'NA'
        
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
        
        # Add PPI enrichment data if available
        if module_id in ppi_enrichment_data:
            module_entry['PPI_FDR'] = ppi_enrichment_data[module_id]['PPI_FDR']
            module_entry['PPI_fold_enrichment'] = ppi_enrichment_data[module_id]['PPI_fold_enrichment']
        else:
            module_entry['PPI_FDR'] = 'NA'
            module_entry['PPI_fold_enrichment'] = 'NA'
        
        # Add metabolite-gene interaction enrichment data if available
        if module_id in metgene_enrichment_data:
            module_entry['MetGene_FDR'] = metgene_enrichment_data[module_id]['MetGene_FDR']
            module_entry['MetGene_fold_enrichment'] = metgene_enrichment_data[module_id]['MetGene_fold_enrichment']
            module_entry['MetGene_N_interactions'] = metgene_enrichment_data[module_id]['MetGene_N_interactions']
        else:
            module_entry['MetGene_FDR'] = 'NA'
            module_entry['MetGene_fold_enrichment'] = 'NA'
            module_entry['MetGene_N_interactions'] = 'NA'
        
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
    
    print(f"\nPerforming functional clustering with canonical MegaGO {CANONICAL_CLUSTER_COLUMN} workflow...")
    bp_terms_top30_df, bp_terms_top30_path, top30_dir = prepare_top30_bp_terms(
        enrichment_data,
        output_dir,
        enrichment_metadata.get('selected_source')
    )
    module_clusters, go_similarity_matrix = megago_cluster_modules(
        module_pathways,
        bp_terms_top30_df,
        args.n_clusters,
        output_dir=top30_dir
    )
    cluster_assignments_df, cluster_assignments_path = write_cluster_assignments(
        module_clusters,
        top30_dir,
        cluster_column=CANONICAL_CLUSTER_COLUMN
    )
    rrvgo_results = run_rrvgo_labeler(
        bp_terms_top30_path,
        cluster_assignments_path,
        top30_dir,
        args.organism,
        cluster_column=CANONICAL_CLUSTER_COLUMN
    )

    cluster_label_lookup = {}
    if not rrvgo_results['cluster_labels'].empty:
        cluster_labels_df = rrvgo_results['cluster_labels'].copy()
        cluster_label_lookup = {
            str(row['MegaGO_cluster']): str(row['MegaGO_label'])
            for _, row in cluster_labels_df.iterrows()
            if pd.notna(row.get('MegaGO_label')) and str(row.get('MegaGO_label')).strip()
        }

    module_label_lookup = {}
    if not rrvgo_results['module_labels'].empty:
        module_labels_df = rrvgo_results['module_labels'].copy()
        module_labels_df['Module'] = module_labels_df['Module'].astype(str)
        module_label_lookup = module_labels_df.set_index('Module').to_dict(orient='index')
    
    # -- Load PKN data for edge categorization --
    pkn_lookup = None
    name_to_hmdb = None
    if args.pkn_file and os.path.exists(args.pkn_file):
        mapping_path = args.metabolite_mapping if (args.metabolite_mapping and os.path.exists(args.metabolite_mapping)) else None
        pkn_lookup, name_to_hmdb = load_pkn_data(args.pkn_file, mapping_path)
    
    # Build module_genes_map for PKN edge matching
    module_genes_map = {}
    for mid, genes in module_genes.items():
        module_genes_map[mid] = genes
    
    # Add cluster information to module data
    for module_entry in module_data:
        module_id = str(module_entry['Module'])
        megago_cluster = module_clusters.get(module_id, 'Unassigned')
        label_info = module_label_lookup.get(module_id, {})

        megago_label = label_info.get('MegaGO_label', cluster_label_lookup.get(megago_cluster, 'NA'))
        if pd.isna(megago_label) or str(megago_label).strip() == '':
            megago_label = 'NA'

        representative_go_id = label_info.get('representative_go_id', 'NA')
        if pd.isna(representative_go_id) or str(representative_go_id).strip() == '':
            representative_go_id = 'NA'

        label_source = label_info.get('label_source', 'NA')
        if pd.isna(label_source) or str(label_source).strip() == '':
            label_source = 'NA'

        module_entry['Functional_Cluster'] = megago_cluster
        module_entry['MegaGO_Cluster'] = megago_cluster
        module_entry['MegaGO_Label'] = megago_label
        module_entry['MegaGO_Representative_GO_ID'] = representative_go_id
        module_entry['MegaGO_Label_Source'] = label_source
        
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

    # Build module overview DataFrame for enriched hover text
    module_overview_df = pd.DataFrame(module_data)

    # Build combined enrichment DataFrame with Database column
    enrichment_all_frames = []
    db_name_map = {
        'bp': 'GO_BP', 'mf': 'GO_MF', 'cc': 'GO_CC',
        'reactome': 'Reactome', 'kegg': 'KEGG'
    }
    for key, df in enrichment_data.items():
        if df is not None and not df.empty:
            df_copy = df.copy()
            df_copy['Database'] = db_name_map.get(key, key)
            enrichment_all_frames.append(df_copy)
    enrichment_all_df = pd.concat(enrichment_all_frames, ignore_index=True) if enrichment_all_frames else pd.DataFrame()

    # -- Create interactive visualizations --
    print(f"\nCreating interactive visualizations...")

    network_fig = create_interactive_network_visualization(
        module_data, module_clusters, output_dir,
        enrichment_data=enrichment_data,
        go_similarity_matrix=go_similarity_matrix,
        module_overview_df=module_overview_df,
        enrichment_all_df=enrichment_all_df,
        cluster_label_lookup=cluster_label_lookup,
        pkn_lookup=pkn_lookup,
        name_to_hmdb=name_to_hmdb,
        module_genes_map=module_genes_map,
        name_lookup=name_lookup,
        run_id=args.run_id,
    )

    # Export Cytoscape-compatible TSV files
    nodes_for_cytoscape = []
    edges_for_cytoscape = []
    for module_info in module_data:
        mid = str(module_info['Module'])
        nodes_for_cytoscape.append({'id': f"Module_{mid}", 'type': 'module'})
        for col in module_info:
            if col.endswith('_regulators') and module_info.get(col, 'NA') != 'NA':
                rtype = col.replace('_regulators', '')
                for reg in module_info[col].split('|'):
                    reg = reg.strip()
                    if reg:
                        nodes_for_cytoscape.append({'id': reg, 'type': rtype})
                        edges_for_cytoscape.append({'source': reg, 'target': f"Module_{mid}", 'type': f"{rtype}_regulation"})
    if pkn_lookup and module_genes_map:
        edges_for_cytoscape = annotate_edges_with_category(edges_for_cytoscape, pkn_lookup, name_to_hmdb or {}, module_genes_map)
    else:
        for e in edges_for_cytoscape:
            e['category'] = 'Other'

    export_cytoscape_files(nodes_for_cytoscape, edges_for_cytoscape, module_clusters,
                           module_overview_df, output_dir)

    cluster_stats = pd.DataFrame()  # Empty DataFrame to avoid errors
    
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
    xlsx_file = os.path.join(output_dir, 'Module_Overview.xlsx')
    module_df.to_excel(xlsx_file, index=False)
    print(f"Enhanced module overview saved to: {xlsx_file}")
    
    # Generate regulator ranking tables HTML
    regulator_tables_html = None
    if regulator_scores_dict:
        regulator_tables_path = generate_regulator_tables_html(regulator_scores_dict, output_dir, name_lookup=name_lookup)
        # Read the generated HTML for embedding in comprehensive report
        try:
            with open(regulator_tables_path, 'r') as f:
                regulator_tables_html = f.read()
        except Exception as e:
            print(f"Warning: Could not read regulator tables HTML: {e}")
    
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
    
    # Print PPI enrichment summary
    if ppi_enrichment_data:
        modules_with_ppi = sum(1 for entry in module_data if entry.get('PPI_FDR') != 'NA')
        significant_ppi = sum(1 for entry in module_data if entry.get('PPI_FDR') != 'NA' and entry.get('PPI_FDR') < 0.05)
        print(f"PPI enrichment analysis: {modules_with_ppi} modules analyzed")
        print(f"Significantly enriched modules (FDR < 0.05): {significant_ppi}")
    
    # Print metabolite-gene enrichment summary
    if metgene_enrichment_data:
        modules_with_metgene = sum(1 for entry in module_data if entry.get('MetGene_FDR') != 'NA')
        significant_metgene = sum(1 for entry in module_data if entry.get('MetGene_FDR') != 'NA' and entry.get('MetGene_FDR') < 0.05)
        print(f"Metabolite-gene enrichment analysis: {modules_with_metgene} modules analyzed")
        print(f"Significantly enriched modules (FDR < 0.05): {significant_metgene}")
    
    # Print cluster summary
    cluster_counts = Counter(module_clusters.values())
    print(f"\nCluster distribution:")
    for cluster_id in sorted(cluster_counts.keys()):
        count = cluster_counts[cluster_id]
        print(f"  {cluster_id}: {count} modules")
    
    print(f"\nFiles created:")
    print(f"  - {output_file}")
    
    movable_html_path = os.path.join(output_dir, 'interactive_module_network_movable.html')
    if os.path.exists(movable_html_path):
        print(f"  - {movable_html_path}")
    
    cytoscape_edges = os.path.join(output_dir, 'module_network_edges.txt')
    cytoscape_nodes = os.path.join(output_dir, 'module_network_node_attributes.txt')
    if os.path.exists(cytoscape_edges):
        print(f"  - {cytoscape_edges}")
        print(f"  - {cytoscape_nodes}")
    
    regulator_html_path = None
    if regulator_scores_dict:
        regulator_html_path = os.path.join(output_dir, 'regulator_rankings.html')
        print(f"  - {regulator_html_path}")

    
    if args.prioritize_by_expression and not expression_results_df.empty:
        print(f"  - {os.path.join(output_dir, 'module_expression_analysis.csv')}")
    
    print("="*80)
    print("                     SUCCESS! Analysis Complete")
    print("              Open HTML files in web browser for visualizations!")
    print("="*80)

if __name__ == "__main__":
    main()
