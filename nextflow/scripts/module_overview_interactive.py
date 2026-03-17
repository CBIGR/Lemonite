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
 

def generate_regulator_tables_html(regulator_scores_dict, output_dir):
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
            html_parts.append(
                f'<tr><td>{rank}</td>'
                f'<td>{row["Regulator"]}</td>'
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
            html_parts.append(
                f'<tr><td>{rank}</td>'
                f'<td>{row["Regulator"]}</td>'
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

        # Look for enrichment combined files (both EnrichR and GSEA patterns)
        up_candidates = [p for p in csv_files
                        if ('enrichr' in os.path.basename(p).lower() or 'gsea' in os.path.basename(p).lower())
                        and 'up' in os.path.basename(p).lower()
                        and 'per_module' in os.path.basename(p).lower()]
        down_candidates = [p for p in csv_files
                          if ('enrichr' in os.path.basename(p).lower() or 'gsea' in os.path.basename(p).lower())
                          and 'down' in os.path.basename(p).lower()
                          and 'per_module' in os.path.basename(p).lower()]

        # Also check exact filename fallbacks (both EnrichR and GSEA)
        for pattern in ['enrichr_top_10_enriched_pathways_up_per_module.csv',
                        'gsea_top_10_enriched_pathways_up_per_module.csv']:
            fallback = next((p for p in csv_files if os.path.basename(p).lower() == pattern), None)
            if fallback and fallback not in up_candidates:
                up_candidates.append(fallback)
        for pattern in ['enrichr_top_10_enriched_pathways_down_per_module.csv',
                        'gsea_top_10_enriched_pathways_down_per_module.csv']:
            fallback = next((p for p in csv_files if os.path.basename(p).lower() == pattern), None)
            if fallback and fallback not in down_candidates:
                down_candidates.append(fallback)

        # When multiple candidates exist per method+direction (e.g. from different runs),
        # try to select the best match per method using module count
        def _select_best_per_method(candidates, module_count=None):
            """Group candidates by method and select best file per method."""
            by_method = {}
            for c in candidates:
                bn = os.path.basename(c).lower()
                m = 'enrichr' if 'enrichr' in bn else ('gsea' if 'gsea' in bn else 'other')
                by_method.setdefault(m, []).append(c)

            selected = []
            for m, method_cands in by_method.items():
                if len(method_cands) == 1:
                    selected.append(method_cands[0])
                elif module_count:
                    matched = [c for c in method_cands if str(module_count) in os.path.basename(c)]
                    selected.append(matched[0] if matched else method_cands[0])
                else:
                    selected.append(method_cands[0])
            return selected

        # Determine module count for disambiguation if needed (>1 file per method)
        current_module_count = None
        needs_disambiguation = any(
            sum(1 for c in cands if kw in os.path.basename(c).lower()) > 1
            for cands in [up_candidates, down_candidates]
            for kw in ['enrichr', 'gsea']
        )
        if needs_disambiguation:
            for root, dirs, files in os.walk(os.path.dirname(enrichment_dir)):
                for f in files:
                    if 'module' in f.lower() and f.endswith('.txt'):
                        try:
                            with open(os.path.join(root, f), 'r') as file:
                                lines = file.readlines()
                                if lines:
                                    modules = set()
                                    for line in lines[1:]:
                                        parts_line = line.strip().split('\t')
                                        if len(parts_line) >= 2:
                                            modules.add(parts_line[1])
                                    if modules:
                                        current_module_count = len(modules)
                                        print(f"Detected {current_module_count} modules from {f}")
                                        break
                        except:
                            continue

            up_candidates = _select_best_per_method(up_candidates, current_module_count)
            down_candidates = _select_best_per_method(down_candidates, current_module_count)

        # Deduplicate and load all enrichment files (EnrichR + GSEA, up + down)
        all_enrichment_paths = list(dict.fromkeys(up_candidates + down_candidates))

        if all_enrichment_paths:
            parts = []
            for pth in all_enrichment_paths:
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
                    print(f"Loaded enrichment file: {os.path.basename(pth)} ({len(df)} rows)")
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

def categorize_modules_by_keywords(enrichment_data):
    """
    Naive heuristic to assign each module to a broad GO category based on keyword
    matching in its enriched Biological Process terms.
    
    Parameters:
    enrichment_data (dict): Enrichment data dictionary (key 'bp' has BP DataFrame)
    
    Returns:
    dict: Mapping module id (string) -> broad GO category string
    """
    go_categories = {
        'Immune': ['immune', 'inflammatory', 'cytokine', 'interferon', 'lymphocyte',
                   'leukocyte', 'antigen', 'defense', 'innate immunity', 'adaptive immunity',
                   'inflammation'],
        'Metabolic': ['metabolic', 'metabolism', 'biosynthetic', 'catabolic',
                      'glycolysis', 'oxidation', 'fatty acid', 'lipid metabolism',
                      'glucose', 'ATP'],
        'Cell Cycle': ['cell cycle', 'mitotic', 'division', 'proliferation',
                       'DNA replication', 'chromosome', 'cytokinesis', 'G1/S', 'G2/M'],
        'Signaling': ['signal transduction', 'signaling', 'receptor', 'kinase',
                      'phosphorylation', 'MAPK', 'cascade', 'pathway', 'GTPase'],
        'Development': ['development', 'differentiation', 'morphogenesis',
                        'embryonic', 'organogenesis', 'pattern specification'],
        'Apoptosis': ['apoptosis', 'cell death', 'programmed cell death',
                      'caspase'],
        'Adhesion_Migration': ['cell adhesion', 'migration', 'motility',
                               'locomotion', 'extracellular matrix', 'integrin'],
        'Transcription': ['transcription', 'RNA processing', 'gene expression',
                          'chromatin', 'histone', 'epigenetic'],
        'Transport': ['transport', 'localization', 'secretion', 'export',
                      'import', 'vesicle'],
        'Stress_Response': ['stress', 'response to stimulus', 'oxidative stress',
                            'DNA damage', 'hypoxia', 'heat shock']
    }

    module_categories = {}

    if 'bp' not in enrichment_data or enrichment_data['bp'].empty:
        print("No BP enrichment data available for keyword clustering")
        return module_categories

    bp_df = enrichment_data['bp']
    term_col = 'Term' if 'Term' in bp_df.columns else None
    mod_col = 'Module' if 'Module' in bp_df.columns else None
    if term_col is None or mod_col is None:
        return module_categories

    for module in bp_df[mod_col].unique():
        mod_df = bp_df[bp_df[mod_col] == module]
        bp_terms = mod_df[term_col].tolist()
        if not bp_terms:
            module_categories[str(module)] = 'Other'
            continue

        scores = {cat: 0 for cat in go_categories}
        for term in bp_terms:
            if not isinstance(term, str):
                continue
            tl = term.lower()
            for cat, kws in go_categories.items():
                for kw in kws:
                    if kw in tl:
                        scores[cat] += 1

        if max(scores.values()) > 0:
            module_categories[str(module)] = max(scores, key=scores.get)
        else:
            module_categories[str(module)] = 'Other'

    print(f"Keyword clustering assigned {len(module_categories)} modules to {len(set(module_categories.values()))} categories")
    return module_categories


def rand_index_from_labels(lbl1, lbl2):
    """Compute the Rand Index between two label vectors."""
    n = len(lbl1)
    assert n == len(lbl2)
    agreements = 0
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            same1 = lbl1[i] == lbl1[j]
            same2 = lbl2[i] == lbl2[j]
            if same1 == same2:
                agreements += 1
            total += 1
    return agreements / total if total else 0.0


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


def build_enriched_hover_text(module_id, module_overview_df, enrichment_all_df, edges):
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
            hover_text += f"<b>Genes:</b> {gene_count} genes<br><br>"

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
            hover_text += f"<b>{reg_type.capitalize()} ({len(regs)}):</b> "
            hover_text += ', '.join(regs[:10])
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


def export_cytoscape_files(nodes, edges, module_clusters, module_overview_df, naive_categories, output_dir):
    """
    Export Cytoscape-compatible edge list and node attribute TSV files.
    
    Parameters:
    nodes (list): Node dicts with 'id', 'type', etc.
    edges (list): Edge dicts with 'source', 'target', 'category'
    module_clusters (dict): module_id -> 'Cluster_N' string
    module_overview_df (pd.DataFrame): Module overview for attribute lookup
    naive_categories (dict or None): module_id -> naive keyword category
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
        header = "Node\tMegaGO_Cluster\tNaive_Category\tNode_Type\tExpression_significant\tPPI_significant\tModule_genes_count\n"
        f.write(header)
        for node in nodes:
            node_id = node['id']
            node_type = node['type']
            if node_type == 'module':
                module_num = node_id.replace('Module_', '')
                cluster = module_clusters.get(module_num, 'Unassigned')
                naive_cat = naive_categories.get(module_num, '') if naive_categories else ''

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
                            genes = mod_row.iloc[0].get('Module_genes', '')
                            if genes != 'NA' and pd.notna(genes):
                                gene_count = str(len(str(genes).split('|')))
                        except Exception:
                            pass
                f.write(f"{node_id}\t{cluster}\t{naive_cat}\t{node_type}\t{expr_flag}\t{ppi_flag}\t{gene_count}\n")
            else:
                naive_cat = ''
                f.write(f"{node_id}\t\t{naive_cat}\t{node_type}\t\t\t\n")
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
    Cluster modules using megaGO command-line tool if available, otherwise use pathway similarity.

    Returns cluster labels as 'Cluster_N' strings for consistency with downstream code.

    Parameters:
    module_pathways (dict): Dictionary mapping module IDs to lists of pathway terms
    n_clusters (int): Number of clusters to create
    use_megago (bool): Whether to use megago for clustering
    enrichment_data (dict): Enrichment data for creating megaGO files
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
                                              naive_categories=None, pkn_lookup=None,
                                              name_to_hmdb=None, module_genes_map=None):
    """
    Create interactive network visualization using Plotly with MegaGO cluster coloring,
    PKN-based edge categorization, and optional naive keyword category borders.

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
    naive_categories (dict or None): module_id -> naive keyword category
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
            regulator_node = {
                'id': regulator,
                'label': regulator[:15] if len(regulator) > 15 else regulator,
                'type': reg_type,
                'hover_info': f"<b>{regulator}</b><br>Type: {reg_type}<br>Targets ({len(target_modules)}): M{', M'.join(target_modules)}"
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
                module_id, module_overview_df, enrichment_all_df, edges
            )
            cluster = module_clusters.get(module_id, 'Unassigned')
            node['hover_info'] += f"<br><b>MegaGO Cluster:</b> {cluster}"
            if naive_categories:
                naive_cat = naive_categories.get(module_id, 'Other')
                node['hover_info'] += f"<br><b>Keyword Category:</b> {naive_cat}"
    
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

    naive_color_map = {}
    if naive_categories is not None:
        naives = sorted(set(naive_categories.values()))
        naive_palette = ['#FFD700', '#ADFF2F', '#00FA9A', '#FF69B4',
                         '#BA55D3', '#FFA500', '#87CEFA', '#40E0D0']
        naive_color_map = {cl: naive_palette[i % len(naive_palette)]
                           for i, cl in enumerate(naives)}
        naive_color_map['Other'] = '#DDDDDD'

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
                    if naive_categories is not None:
                        border = naive_color_map.get(naive_categories.get(m_id, 'Other'), '#FFFFFF')
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

                # Compute border colors per module
                line_colours = []
                for n in cat_nodes:
                    m_id = n['id'].replace('Module_', '')
                    if naive_categories is not None:
                        line_colours.append(naive_color_map.get(naive_categories.get(m_id, 'Other'), 'white'))
                    else:
                        line_colours.append('white')

                trace = go.Scatter(
                    x=x_vals, y=y_vals, mode='markers+text',
                    marker=dict(size=50, color=cluster_color_map.get(cluster, '#CCCCCC'),
                                line=dict(width=3, color=line_colours), symbol='circle'),
                    text=labels, textposition='middle center',
                    textfont=dict(size=12, color='black', family='Arial Black'),
                    hovertext=hover_texts, hoverinfo='text',
                    name=f"{cluster} ({len(cat_nodes)})", legendgroup=cluster
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
            annotations.append(dict(
                x=cx * 1.3, y=cy * 1.3,
                text=f'<b>{cl}</b>', showarrow=False,
                font=dict(size=14, color=cluster_color_map.get(cl, '#000000')),
                bgcolor='rgba(255,255,255,0.8)', borderpad=4
            ))

        # Add naive legend entries
        if naive_categories is not None and not movable:
            for cat, col in naive_color_map.items():
                fig.add_trace(go.Scatter(
                    x=[None], y=[None], mode='markers',
                    marker=dict(size=20, color='white', line=dict(width=3, color=col)),
                    name=f'Keyword: {cat}', showlegend=True
                ))

        fig.update_layout(
            title=dict(
                text='Module-Regulator Network<br><sub>Modules colored by MegaGO cluster</sub>',
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

    # -- Generate standard network --
    fig_standard = _build_figure(movable=False)
    output_file = os.path.join(output_dir, 'interactive_module_network.html')
    pyo.plot(fig_standard, filename=output_file, auto_open=False)
    print(f"Interactive network visualization saved to: {output_file}")

    # -- Generate movable network --
    fig_movable = _build_figure(movable=True)
    output_movable = os.path.join(output_dir, 'interactive_module_network_movable.html')
    fig_movable.write_html(output_movable, config={'editable': True})
    print(f"Movable network visualization saved to: {output_movable}")

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
                       choices=['EnrichR', 'GSEA', 'auto'],
                       help='Enrichment analysis method to use')
    parser.add_argument('--n_clusters', type=int, default=5,
                       help='Number of functional clusters to create')
    parser.add_argument('--clustering_method', type=str, default='megago',
                       choices=['megago', 'keyword', 'both'],
                       help='Clustering method: megago (default), keyword (naive GO keywords), or both')
    # Keep --use_megago/--no_megago for backward compatibility
    parser.add_argument('--use_megago', action='store_true', default=None,
                       help='(Deprecated) Use megago clustering. Prefer --clustering_method.')
    parser.add_argument('--no_megago', dest='use_megago', action='store_false',
                       help='(Deprecated) Disable megago clustering. Prefer --clustering_method keyword.')
    parser.add_argument('--pkn_file', type=str, default=None,
                       help='Path to PKN file (Lemonite_PKN.tsv) for edge categorization')
    parser.add_argument('--metabolite_mapping', type=str, default=None,
                       help='Path to metabolite name mapping file (name_map.csv) for HMDB resolution')
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
    
    # -- Backward compatibility: map deprecated --use_megago/--no_megago to --clustering_method --
    if args.use_megago is not None:
        if args.use_megago:
            args.clustering_method = 'megago'
        else:
            args.clustering_method = 'keyword'
        print(f"Note: --use_megago/--no_megago is deprecated. Mapped to --clustering_method {args.clustering_method}")
    
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
    print(f"Clustering method: {args.clustering_method}")
    if args.pkn_file:
        print(f"PKN file: {args.pkn_file}")
    if args.metabolite_mapping:
        print(f"Metabolite mapping: {args.metabolite_mapping}")
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
    
    enrichment_data = load_enrichment_data(enrichment_dir, args.enrichment_method)
    
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
    
    # -- Determine clustering method and run clustering --
    print(f"\nPerforming functional clustering (method={args.clustering_method})...")
    
    use_megago_flag = args.clustering_method in ('megago', 'both')
    module_clusters, go_similarity_matrix = megago_cluster_modules(
        module_pathways, args.n_clusters, use_megago_flag, enrichment_data, output_dir
    )
    
    # Naive keyword clustering (always run when method is 'keyword' or 'both')
    naive_categories = None
    if args.clustering_method in ('keyword', 'both'):
        print("Running naive keyword clustering on GO terms...")
        naive_categories = categorize_modules_by_keywords(enrichment_data)
        if args.clustering_method == 'keyword':
            # When keyword-only, use naive categories as the primary module_clusters
            module_clusters = {mid: cat for mid, cat in naive_categories.items()}
    
    # If 'both', compute Rand Index for comparison
    if args.clustering_method == 'both' and naive_categories:
        shared = sorted(set(module_clusters.keys()) & set(naive_categories.keys()))
        if len(shared) > 1:
            ri = rand_index_from_labels(
                [module_clusters[m] for m in shared],
                [naive_categories[m] for m in shared]
            )
            print(f"Rand Index (megaGO vs keyword): {ri:.3f}")
    
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
        naive_categories=naive_categories,
        pkn_lookup=pkn_lookup,
        name_to_hmdb=name_to_hmdb,
        module_genes_map=module_genes_map
    )
    
    # Export Cytoscape-compatible TSV files
    nodes_for_cytoscape = []
    edges_for_cytoscape = []
    # Re-build minimal nodes/edges for export (same logic as in visualization)
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
    # Annotate edges for export
    if pkn_lookup and module_genes_map:
        edges_for_cytoscape = annotate_edges_with_category(edges_for_cytoscape, pkn_lookup, name_to_hmdb or {}, module_genes_map)
    else:
        for e in edges_for_cytoscape:
            e['category'] = 'Other'
    
    export_cytoscape_files(nodes_for_cytoscape, edges_for_cytoscape, module_clusters,
                           module_overview_df, naive_categories, output_dir)
    
    cluster_stats = pd.DataFrame()  # Empty DataFrame to avoid errors
    
    # Add cluster information to module data
    for module_entry in module_data:
        module_id = str(module_entry['Module'])
        module_entry['Functional_Cluster'] = module_clusters.get(module_id, 'Cluster_0')
        
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
    
    # Generate regulator ranking tables HTML
    regulator_tables_html = None
    if regulator_scores_dict:
        regulator_tables_path = generate_regulator_tables_html(regulator_scores_dict, output_dir)
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
    
    network_html_path = os.path.join(output_dir, 'interactive_module_network.html')
    print(f"  - {network_html_path}")
    
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

    
    # Generate comprehensive HTML report combining network and regulator tables
    comprehensive_report = create_comprehensive_html_report(
        regulator_tables_html_path=regulator_html_path,
        network_html_path=network_html_path,
        output_dir=output_dir
    )
    print(f"  - {comprehensive_report}")
    
    if args.prioritize_by_expression and not expression_results_df.empty:
        print(f"  - {os.path.join(output_dir, 'module_expression_analysis.csv')}")
    
    print("="*80)
    print("                     SUCCESS! Analysis Complete")
    print("              Open HTML files in web browser for visualizations!")
    print("="*80)

if __name__ == "__main__":
    main()
