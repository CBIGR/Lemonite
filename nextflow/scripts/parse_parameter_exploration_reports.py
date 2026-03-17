#!/usr/bin/env python3

"""
Parse Lemonite Parameter Exploration HTML Reports

This script parses multiple HTML summary reports generated during parameter
exploration and creates an overview table with key metrics and visualizations.

Author: Boris Vandemoortele
"""

import os
import sys
import re
import glob
import argparse
import warnings
from pathlib import Path
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

warnings.filterwarnings('ignore')

# Try to import visualization libraries
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOT_AVAILABLE = True
except ImportError:
    PLOT_AVAILABLE = False
    print("Warning: matplotlib/seaborn not available. Visualizations will be disabled.")


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Parse Lemonite parameter exploration HTML reports'
    )
    parser.add_argument(
        '--input_dir',
        required=True,
        help='Directory containing HTML reports (will search recursively)'
    )
    parser.add_argument(
        '--output_prefix',
        default='parameter_exploration',
        help='Output file prefix (default: parameter_exploration)'
    )
    parser.add_argument(
        '--pattern',
        default='*Summary_Report.html',
        help='File pattern to match (default: *Summary_Report.html)'
    )
    parser.add_argument(
        '--no-plots',
        action='store_true',
        help='Disable plot generation'
    )
    return parser.parse_args()


def extract_stat_value(soup, section_id, stat_label):
    """Extract a statistic value from a specific section by label
    
    Args:
        soup: BeautifulSoup object
        section_id: ID of the section to search in
        stat_label: Label text to search for
        
    Returns:
        Extracted value as string or 'N/A'
    """
    section = soup.find('section', {'id': section_id})
    if not section:
        return 'N/A'
    
    # Find stat card with matching label
    for card in section.find_all('div', class_='stat-card'):
        label_div = card.find('div', class_='stat-label')
        if label_div and stat_label.lower() in label_div.get_text().lower():
            value_div = card.find('div', class_='stat-value')
            if value_div:
                return value_div.get_text().strip()
    
    return 'N/A'


def extract_regulator_counts(soup):
    """Extract regulator counts per type from the network section
    
    Returns:
        Dictionary with regulator type counts
    """
    regulators = {}
    
    # Find the regulators table in the network section
    network_section = soup.find('section', {'id': 'network-section'})
    if not network_section:
        return regulators
    
    # Look for the regulators table
    tables = network_section.find_all('table')
    for table in tables:
        headers = [th.get_text().strip() for th in table.find_all('th')]
        if 'Regulator Type' in headers:
            for row in table.find_all('tr')[1:]:  # Skip header
                cells = row.find_all('td')
                if len(cells) >= 2:
                    reg_type = cells[0].get_text().strip()
                    reg_count = cells[1].get_text().strip()
                    # Remove badge formatting if present
                    reg_type = re.sub(r'[^\w\s]', '', reg_type)
                    regulators[reg_type] = reg_count
    
    return regulators


def extract_run_id(soup):
    """Extract run ID from the header"""
    header = soup.find('div', class_='run-info')
    if header:
        text = header.get_text()
        match = re.search(r'Run ID:\s*([^\|]+)', text)
        if match:
            return match.group(1).strip()
    return 'Unknown'


def extract_folder_name(html_path):
    """Extract folder name from HTML path"""
    # Get the parent directory of the HTML file
    path = Path(html_path)
    # Usually it's 2 levels up (run_id/run_id/Lemonite_Summary_Report.html)
    folder_name = path.parent.parent.name
    return folder_name


def extract_parameters(soup):
    """Extract key parameters from the parameters section
    
    Returns:
        Dictionary of parameters
    """
    params = {}
    
    params_section = soup.find('section', {'id': 'parameters-section'})
    if not params_section:
        return params
    
    # Find all parameter items
    for param_item in params_section.find_all('div', class_='param-item'):
        name_span = param_item.find('span', class_='param-name')
        value_span = param_item.find('span', class_='param-value')
        
        if name_span and value_span:
            param_name = name_span.get_text().strip()
            param_value = value_span.get_text().strip()
            params[param_name] = param_value
    
    return params


def parse_html_report(html_path):
    """Parse a single HTML report and extract key statistics
    
    Args:
        html_path: Path to HTML report file
        
    Returns:
        Dictionary with extracted statistics
    """
    print(f"Parsing: {html_path}")
    
    with open(html_path, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f.read(), 'html.parser')
    
    stats = {}
    
    # Extract run ID and folder name
    stats['Run_ID'] = extract_run_id(soup)
    stats['Folder_Name'] = extract_folder_name(html_path)
    stats['Report_Path'] = str(html_path)
    
    # Extract key metrics from different sections
    
    # === Input/Summary Section ===
    stats['Total_Samples'] = extract_stat_value(soup, 'summary-section', 'Samples Analyzed')
    stats['Genes_Analyzed'] = extract_stat_value(soup, 'summary-section', 'Genes Analyzed')
    
    # === Network Section ===
    stats['Genes_in_Network'] = extract_stat_value(soup, 'network-section', 'Genes in Network')
    stats['TFs_in_Network'] = extract_stat_value(soup, 'network-section', 'TFs in Network')
    stats['Metabolites_in_Network'] = extract_stat_value(soup, 'network-section', 'Metabolites in Network')
    stats['Proteins_in_Network'] = extract_stat_value(soup, 'network-section', 'Proteins in Network')
    stats['Network_Edges'] = extract_stat_value(soup, 'network-section', 'Network Edges')
    
    # Extract regulator counts
    regulators = extract_regulator_counts(soup)
    for reg_type, count in regulators.items():
        stats[f'Regulators_{reg_type}'] = count
    
    # === Clustering Section ===
    stats['Consensus_Modules'] = extract_stat_value(soup, 'clustering-section', 'Consensus Modules')
    stats['Mean_Module_Size'] = extract_stat_value(soup, 'clustering-section', 'Mean Module Size')
    stats['Median_Module_Size'] = extract_stat_value(soup, 'clustering-section', 'Median Module Size')
    
    # === Network Section (Module counts) ===
    stats['Initial_Modules'] = extract_stat_value(soup, 'network-section', 'Initial Modules')
    stats['Final_Modules'] = extract_stat_value(soup, 'network-section', 'Final Modules')
    stats['DE_Modules'] = extract_stat_value(soup, 'network-section', 'DE Modules')
    stats['PPI_Enriched_Modules'] = extract_stat_value(soup, 'network-section', 'PPI Enriched')
    
    # === Extract key parameters ===
    params = extract_parameters(soup)
    
    # Add select parameters to stats
    key_params = [
        'top_n_genes', 'n_clusters', 'coherence_threshold',
        'regulator_selection_method', 'top_n_percent_regulators',
        'min_module_size', 'min_targets'
    ]
    
    for param in key_params:
        if param in params:
            stats[f'Param_{param}'] = params[param]
    
    return stats


def clean_numeric_column(series):
    """Clean and convert a pandas series to numeric values"""
    def clean_value(val):
        if pd.isna(val) or val == 'N/A':
            return np.nan
        # Remove any non-numeric characters except decimal point and minus
        val_str = str(val).strip()
        # Handle scientific notation
        val_str = re.sub(r'[^\d\.\-e]', '', val_str)
        try:
            return float(val_str)
        except:
            return np.nan
    
    return series.apply(clean_value)


def create_overview_table(html_files):
    """Parse all HTML files and create overview DataFrame
    
    Args:
        html_files: List of HTML file paths
        
    Returns:
        pandas DataFrame with all statistics
    """
    all_stats = []
    
    for html_file in html_files:
        try:
            stats = parse_html_report(html_file)
            all_stats.append(stats)
        except Exception as e:
            print(f"Error parsing {html_file}: {e}")
            continue
    
    if not all_stats:
        print("No reports successfully parsed!")
        return None
    
    df = pd.DataFrame(all_stats)
    
    # Convert numeric columns
    numeric_cols = [col for col in df.columns if any(x in col for x in [
        'Samples', 'Genes', 'TFs', 'Metabolites', 'Proteins', 'Modules',
        'Edges', 'Regulators', 'Size'
    ])]
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = clean_numeric_column(df[col])
    
    # Calculate percentages for quality metrics
    # DE modules as % of final modules
    if 'Final_Modules' in df.columns and 'DE_Modules' in df.columns:
        df['DE_Modules_Pct'] = (df['DE_Modules'] / df['Final_Modules'] * 100).fillna(0)
    
    # PPI enriched modules as % of final modules
    if 'Final_Modules' in df.columns and 'PPI_Enriched_Modules' in df.columns:
        df['PPI_Enriched_Modules_Pct'] = (df['PPI_Enriched_Modules'] / df['Final_Modules'] * 100).fillna(0)
    
    # Module retention after coherence filtering
    if 'Initial_Modules' in df.columns and 'Final_Modules' in df.columns:
        df['Module_Retention_Pct'] = (df['Final_Modules'] / df['Initial_Modules'] * 100).fillna(0)
    
    # Calculate composite quality scores that balance discovery and specificity
    # DE Module Discovery Score: (# DE modules) * (DE quality %)
    if 'DE_Modules' in df.columns and 'DE_Modules_Pct' in df.columns:
        df['DE_Discovery_Score'] = df['DE_Modules'] * (df['DE_Modules_Pct'] / 100)
    
    # PPI Module Discovery Score: (# PPI modules) * (PPI quality %)
    if 'PPI_Enriched_Modules' in df.columns and 'PPI_Enriched_Modules_Pct' in df.columns:
        df['PPI_Discovery_Score'] = df['PPI_Enriched_Modules'] * (df['PPI_Enriched_Modules_Pct'] / 100)
    
    # Overall Quality Score: (DE + PPI enriched modules) as count, weighted by their percentages
    if 'DE_Modules' in df.columns and 'PPI_Enriched_Modules' in df.columns:
        total_high_quality = df['DE_Modules'] + df['PPI_Enriched_Modules']
        avg_quality_pct = ((df['DE_Modules_Pct'].fillna(0) + df['PPI_Enriched_Modules_Pct'].fillna(0)) / 2)
        df['Overall_Quality_Score'] = total_high_quality * (avg_quality_pct / 100)
    
    return df


def create_visualizations(df, output_prefix):
    """Create comprehensive visualizations of parameter exploration results
    
    Args:
        df: DataFrame with statistics
        output_prefix: Prefix for output files
    """
    if not PLOT_AVAILABLE:
        print("Skipping visualizations (matplotlib/seaborn not available)")
        return
    
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.facecolor'] = 'white'
    
    # Create output directory for plots
    plot_dir = f"{output_prefix}_plots"
    os.makedirs(plot_dir, exist_ok=True)
    
    # Use folder names for x-axis labels
    x_labels = df['Folder_Name'].values if 'Folder_Name' in df.columns else df.index
    x_pos = np.arange(len(df))
    
    # === Plot 1: Module Quality Metrics (Percentages) ===
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Module Quality Metrics Across Parameter Settings', fontsize=16, fontweight='bold')
    
    # DE Modules as % of Final Modules
    ax1 = axes[0, 0]
    if 'DE_Modules_Pct' in df.columns:
        bars = ax1.bar(x_pos, df['DE_Modules_Pct'], alpha=0.8, color='#FDD835')
        ax1.set_ylabel('% of Final Modules')
        ax1.set_title('DE Modules (% of Final Modules)')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(x_labels, rotation=45, ha='right')
        ax1.grid(axis='y', alpha=0.3)
        # Add value labels on bars
        for i, v in enumerate(df['DE_Modules_Pct']):
            ax1.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # PPI Enriched Modules as % of Final Modules
    ax2 = axes[0, 1]
    if 'PPI_Enriched_Modules_Pct' in df.columns:
        bars = ax2.bar(x_pos, df['PPI_Enriched_Modules_Pct'], alpha=0.8, color='#FF9800')
        ax2.set_ylabel('% of Final Modules')
        ax2.set_title('PPI Enriched Modules (% of Final Modules)')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(x_labels, rotation=45, ha='right')
        ax2.grid(axis='y', alpha=0.3)
        # Add value labels on bars
        for i, v in enumerate(df['PPI_Enriched_Modules_Pct']):
            ax2.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # Module Retention Rate
    ax3 = axes[1, 0]
    if 'Module_Retention_Pct' in df.columns:
        bars = ax3.bar(x_pos, df['Module_Retention_Pct'], alpha=0.8, color='#2196F3')
        ax3.set_ylabel('Retention Rate (%)')
        ax3.set_title('Module Retention After Coherence Filtering')
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(x_labels, rotation=45, ha='right')
        ax3.axhline(y=50, color='r', linestyle='--', alpha=0.5, label='50% threshold')
        ax3.legend()
        ax3.grid(axis='y', alpha=0.3)
        # Add value labels on bars
        for i, v in enumerate(df['Module_Retention_Pct']):
            ax3.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # Overall Quality Score (balances discovery and specificity)
    ax4 = axes[1, 1]
    if 'Overall_Quality_Score' in df.columns:
        bars = ax4.bar(x_pos, df['Overall_Quality_Score'], alpha=0.8, color='#4CAF50')
        ax4.set_ylabel('Quality Score')
        ax4.set_title('Overall Quality Score\n(Balances Discovery & Specificity)')
        ax4.set_xticks(x_pos)
        ax4.set_xticklabels(x_labels, rotation=45, ha='right')
        ax4.grid(axis='y', alpha=0.3)
        # Add value labels on bars
        for i, v in enumerate(df['Overall_Quality_Score']):
            ax4.text(i, v + 0.1, f'{v:.2f}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/module_quality.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # === Plot 1b: Discovery Score - Balancing Absolute Counts with Quality ===
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Discovery Scores: Balancing Quantity & Quality', fontsize=16, fontweight='bold')
    
    # DE Discovery Score
    ax1 = axes[0]
    if 'DE_Discovery_Score' in df.columns:
        bars = ax1.bar(x_pos, df['DE_Discovery_Score'], alpha=0.8, color='#FDD835', edgecolor='#F57F17', linewidth=2)
        ax1.set_ylabel('Discovery Score\n(#DE modules × %DE specificity)')
        ax1.set_title('DE Module Discovery Score')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(x_labels, rotation=45, ha='right')
        ax1.grid(axis='y', alpha=0.3)
        # Add value labels
        for i, v in enumerate(df['DE_Discovery_Score']):
            ax1.text(i, v + 0.1, f'{v:.2f}', ha='center', va='bottom', fontweight='bold')
        # Add secondary info showing components
        for i in range(len(df)):
            comp_text = f"({df['DE_Modules'].iloc[i]:.0f} × {df['DE_Modules_Pct'].iloc[i]:.0f}%)"
            ax1.text(i, -0.5, comp_text, ha='center', va='top', fontsize=8, style='italic')
    
    # PPI Discovery Score
    ax2 = axes[1]
    if 'PPI_Discovery_Score' in df.columns:
        bars = ax2.bar(x_pos, df['PPI_Discovery_Score'], alpha=0.8, color='#FF9800', edgecolor='#E65100', linewidth=2)
        ax2.set_ylabel('Discovery Score\n(#PPI modules × %PPI specificity)')
        ax2.set_title('PPI Enrichment Discovery Score')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(x_labels, rotation=45, ha='right')
        ax2.grid(axis='y', alpha=0.3)
        # Add value labels
        for i, v in enumerate(df['PPI_Discovery_Score']):
            ax2.text(i, v + 0.05, f'{v:.2f}', ha='center', va='bottom', fontweight='bold')
        # Add secondary info showing components
        for i in range(len(df)):
            comp_text = f"({df['PPI_Enriched_Modules'].iloc[i]:.0f} × {df['PPI_Enriched_Modules_Pct'].iloc[i]:.0f}%)"
            ax2.text(i, -0.15, comp_text, ha='center', va='top', fontsize=8, style='italic')
    
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/discovery_scores.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # === Plot 2: Module Counts (Absolute and Percentages) ===
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Module Counts Across Parameter Settings', fontsize=16, fontweight='bold')
    
    # Absolute counts
    ax1 = axes[0]
    if 'Initial_Modules' in df.columns and 'Final_Modules' in df.columns:
        width = 0.35
        bars1 = ax1.bar(x_pos - width/2, df['Initial_Modules'], width, label='Initial Modules', alpha=0.8, color='#81C784')
        bars2 = ax1.bar(x_pos + width/2, df['Final_Modules'], width, label='Final Modules', alpha=0.8, color='#2E7D32')
        ax1.set_ylabel('Number of Modules')
        ax1.set_title('Initial vs Final Modules (Absolute Counts)')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(x_labels, rotation=45, ha='right')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
    
    # Absolute counts - DE and PPI enriched
    ax2 = axes[1]
    if 'DE_Modules' in df.columns and 'PPI_Enriched_Modules' in df.columns:
        width = 0.35
        bars1 = ax2.bar(x_pos - width/2, df['DE_Modules'], width, label='DE Modules', alpha=0.8, color='#FDD835')
        bars2 = ax2.bar(x_pos + width/2, df['PPI_Enriched_Modules'], width, label='PPI Enriched', alpha=0.8, color='#FF9800')
        ax2.set_ylabel('Number of Modules')
        ax2.set_title('DE and PPI Enriched Modules (Absolute Counts)')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(x_labels, rotation=45, ha='right')
        ax2.legend()
        ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/module_counts.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # === Plot 3: Module Size Distribution ===
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Module Size Metrics', fontsize=16, fontweight='bold')
    
    ax1 = axes[0]
    if 'Mean_Module_Size' in df.columns:
        ax1.plot(x_pos, df['Mean_Module_Size'], marker='o', linewidth=2, markersize=8, color='#2E7D32')
        ax1.set_xlabel('Parameter Settings')
        ax1.set_ylabel('Mean Module Size')
        ax1.set_title('Mean Module Size Across Runs')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(x_labels, rotation=45, ha='right')
        ax1.grid(alpha=0.3)
    
    ax2 = axes[1]
    if 'Median_Module_Size' in df.columns:
        ax2.plot(x_pos, df['Median_Module_Size'], marker='s', linewidth=2, markersize=8, color='#81C784')
        ax2.set_xlabel('Parameter Settings')
        ax2.set_ylabel('Median Module Size')
        ax2.set_title('Median Module Size Across Runs')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(x_labels, rotation=45, ha='right')
        ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/module_sizes.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # === Plot 4: Network Composition ===
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Network Composition', fontsize=16, fontweight='bold')
    
    # Genes and Regulators in Network
    ax1 = axes[0]
    network_cols = ['Genes_in_Network', 'TFs_in_Network', 'Metabolites_in_Network', 'Proteins_in_Network']
    network_data = df[[col for col in network_cols if col in df.columns]]
    if not network_data.empty:
        network_data.index = x_labels
        network_data.plot(kind='bar', ax=ax1, width=0.8)
        ax1.set_xlabel('Parameter Settings')
        ax1.set_ylabel('Count')
        ax1.set_title('Network Components Across Runs')
        ax1.legend(title='Component', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(axis='y', alpha=0.3)
        plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Regulator Counts
    ax2 = axes[1]
    regulator_cols = [col for col in df.columns if col.startswith('Regulators_')]
    if regulator_cols:
        regulator_data = df[regulator_cols].copy()
        regulator_data.index = x_labels
        # Rename columns for better labels
        regulator_data.columns = [col.replace('Regulators_', '') for col in regulator_data.columns]
        regulator_data.plot(kind='bar', ax=ax2, width=0.8)
        ax2.set_xlabel('Parameter Settings')
        ax2.set_ylabel('Number of Regulators')
        ax2.set_title('Regulators Assigned per Type')
        ax2.legend(title='Regulator Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.grid(axis='y', alpha=0.3)
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/network_composition.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # === Plot 5: Heatmap of All Metrics ===
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Select numeric columns for heatmap
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    # Exclude index and irrelevant columns
    exclude_patterns = ['Param_', 'Report_']
    heatmap_cols = [col for col in numeric_cols if not any(p in col for p in exclude_patterns)]
    
    if heatmap_cols:
        # Normalize data for better visualization
        heatmap_data = df[heatmap_cols].copy()
        heatmap_data.index = x_labels
        # Replace inf and -inf with NaN
        heatmap_data = heatmap_data.replace([np.inf, -np.inf], np.nan)
        # Normalize each column to 0-1 scale
        for col in heatmap_data.columns:
            col_min = heatmap_data[col].min()
            col_max = heatmap_data[col].max()
            if col_max > col_min:
                heatmap_data[col] = (heatmap_data[col] - col_min) / (col_max - col_min)
        
        sns.heatmap(heatmap_data.T, annot=False, cmap='RdYlGn', cbar_kws={'label': 'Normalized Value'},
                   ax=ax, linewidths=0.5, linecolor='white')
        ax.set_xlabel('Parameter Settings', fontsize=12)
        ax.set_ylabel('Metric', fontsize=12)
        ax.set_title('Normalized Metrics Heatmap (Green=High, Red=Low)', fontsize=14, fontweight='bold')
        plt.setp(ax.yaxis.get_majorticklabels(), rotation=0)
    
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/metrics_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # === Plot 6: Parameter Impact (if parameters are available) ===
    param_cols = [col for col in df.columns if col.startswith('Param_')]
    if param_cols and 'Final_Modules' in df.columns:
        # Try to create scatter plots showing parameter impact
        n_params = min(len(param_cols), 6)  # Limit to 6 parameters
        if n_params > 0:
            fig, axes = plt.subplots(2, 3, figsize=(16, 10))
            fig.suptitle('Parameter Impact on Module Quality (DE% + PPI%)', fontsize=16, fontweight='bold')
            axes = axes.flatten()
            
            quality_score = df['DE_Modules_Pct'].fillna(0) + df['PPI_Enriched_Modules_Pct'].fillna(0)
            
            for i, param_col in enumerate(param_cols[:n_params]):
                ax = axes[i]
                # Try to convert parameter to numeric
                param_values = clean_numeric_column(df[param_col])
                if not param_values.isna().all():
                    ax.scatter(param_values, quality_score, s=100, alpha=0.6, color='#2E7D32')
                    ax.set_xlabel(param_col.replace('Param_', ''), fontsize=10)
                    ax.set_ylabel('Quality Score (%)')
                    ax.grid(alpha=0.3)
                    
                    # Try to add trend line if enough points
                    valid_mask = ~(param_values.isna() | quality_score.isna())
                    if valid_mask.sum() > 2:
                        try:
                            z = np.polyfit(param_values[valid_mask], quality_score[valid_mask], 1)
                            p = np.poly1d(z)
                            x_line = np.linspace(param_values[valid_mask].min(), 
                                               param_values[valid_mask].max(), 100)
                            ax.plot(x_line, p(x_line), "r--", alpha=0.5, linewidth=2)
                        except:
                            pass
            
            # Hide unused subplots
            for i in range(n_params, 6):
                axes[i].set_visible(False)
            
            plt.tight_layout()
            plt.savefig(f"{plot_dir}/parameter_impact.png", dpi=300, bbox_inches='tight')
            plt.close()
    
    print(f"\n✓ Visualizations saved to: {plot_dir}/")
    print(f"  - module_quality.png (Quality metrics as percentages)")
    print(f"  - discovery_scores.png (NEW: Balances quantity & quality)")
    print(f"  - module_counts.png")
    print(f"  - module_sizes.png")
    print(f"  - network_composition.png")
    print(f"  - metrics_heatmap.png")
    if param_cols:
        print(f"  - parameter_impact.png")


def main():
    args = parse_arguments()
    
    print("="*80)
    print("LEMONITE PARAMETER EXPLORATION REPORT PARSER")
    print("="*80)
    print(f"Input directory: {args.input_dir}")
    print(f"Search pattern: {args.pattern}")
    print(f"Output prefix: {args.output_prefix}")
    print("="*80)
    
    # Find all HTML reports
    search_pattern = os.path.join(args.input_dir, '**', args.pattern)
    html_files = glob.glob(search_pattern, recursive=True)
    
    if not html_files:
        print(f"\nERROR: No HTML files found matching pattern: {search_pattern}")
        sys.exit(1)
    
    print(f"\nFound {len(html_files)} HTML report(s):")
    for f in html_files:
        print(f"  - {f}")
    
    # Parse all reports
    print("\nParsing reports...")
    df = create_overview_table(html_files)
    
    if df is None or df.empty:
        print("ERROR: No data extracted from reports!")
        sys.exit(1)
    
    print(f"\n✓ Successfully parsed {len(df)} report(s)")
    
    # Save overview table
    csv_output = f"{args.output_prefix}_overview.csv"
    df.to_csv(csv_output, index=False)
    print(f"\n✓ Overview table saved to: {csv_output}")
    
    # Also save as Excel for easier viewing (if openpyxl available)
    try:
        excel_output = f"{args.output_prefix}_overview.xlsx"
        df.to_excel(excel_output, index=False, engine='openpyxl')
        print(f"✓ Excel file saved to: {excel_output}")
    except ImportError:
        print("  (Install openpyxl to enable Excel export)")
    except Exception as e:
        print(f"  (Could not save Excel: {e})")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    # Key metrics summary
    key_metrics = [
        'Final_Modules', 'DE_Modules', 'PPI_Enriched_Modules',
        'DE_Modules_Pct', 'PPI_Enriched_Modules_Pct', 'Module_Retention_Pct',
        'DE_Discovery_Score', 'PPI_Discovery_Score', 'Overall_Quality_Score',
        'Mean_Module_Size', 'Median_Module_Size',
        'Genes_in_Network', 'TFs_in_Network', 'Metabolites_in_Network'
    ]
    
    for metric in key_metrics:
        if metric in df.columns:
            values = df[metric].dropna()
            if not values.empty:
                print(f"\n{metric}:")
                print(f"  Mean: {values.mean():.2f}")
                print(f"  Median: {values.median():.2f}")
                print(f"  Min: {values.min():.2f}")
                print(f"  Max: {values.max():.2f}")
                print(f"  Std: {values.std():.2f}")
    
    # Create visualizations
    if not args.no_plots:
        print("\n" + "="*80)
        print("GENERATING VISUALIZATIONS")
        print("="*80)
        create_visualizations(df, args.output_prefix)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)


if __name__ == '__main__':
    main()
