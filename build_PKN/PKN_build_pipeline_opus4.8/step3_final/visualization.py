"""
Step-3 visualization: figures summarising the assembled PKN.

Re-derives the reproducible figures from the notebooks (Collect_PKNdata_metabolites,
Collect_PKNdata_proteins, Build_final_PKN):

  - per-source edge-count barplots for each layer (metabolite-gene, PPI, hPTM);
  - UpSet plots of the database overlap for the metabolite-gene and PPI layers
    (which interactions are shared vs. unique across the source databases);
  - a combined node/edge-count summary barplot for the final PKN.

The notebook's MetalinksDB / MEBOCOST comparison and HMDB-superclass figures are
omitted: they require extra external database files that are not part of building or
validating the network (matching the note in analysis.py).

Every figure is written to config.FIGURES_DIR. This is best-effort: a plotting
failure (e.g. a missing optional dependency, or a headless display quirk) logs a
warning and is skipped — it must never fail a completed build. Uses the non-
interactive Agg backend so it runs unattended on HPC without a display.

The overlap/UpSet plots read the long-format per-source network files
(metabolite_gene_PKN.tsv, PPI_network.tsv, hPTM_network.tsv), where an edge may
appear once per supporting source — NOT the de-duplicated LemonIte_PKN.tsv, where
each edge carries only one Source.
"""

import logging
import os
import warnings

import pandas as pd

import matplotlib
matplotlib.use('Agg')  # headless: no display needed on HPC
import matplotlib.pyplot as plt  # noqa: E402

import config

logger = logging.getLogger('pkn.visualization')

# Sources shown in the metabolite-gene overlap plot (order = notebook order).
_MG_SOURCES = ['IntAct', 'chEMBL', 'UniProtKB', 'STITCH', 'BioGRID', 'LINCS',
               'L1000', 'Human1_GEM_dist1', 'Human1_GEM_dist2', 'MetalinksDB']
# Human-GEM distance-1/2 collapse to one label in the count barplot.
_MG_MERGE = {'Human1_GEM_dist1': 'Human1-GEM', 'Human1_GEM_dist2': 'Human1-GEM'}


def _fig_path(name):
    return os.path.join(config.FIGURES_DIR, name)


def _undirected_pairs(df, a, b):
    """Set of sorted (node, node) tuples — treats an edge as undirected."""
    return set(tuple(sorted((str(x), str(y)))) for x, y in zip(df[a], df[b]))


def _barplot_counts(counts, title, ylabel, path, log_scale=False):
    """Simple labelled bar chart of a name->count mapping."""
    counts = counts[counts > 0].sort_values(ascending=False)
    if counts.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(counts) + 2), 6))
    ax.bar(range(len(counts)), counts.values, color='#3b6ea5')
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=45, ha='right')
    ax.set_title(title, fontsize=14)
    ax.set_ylabel(ylabel + (' (log scale)' if log_scale else ''), fontsize=12)
    if log_scale:
        ax.set_yscale('log')
    else:
        for i, v in enumerate(counts.values):
            ax.text(i, v + max(counts.values) * 0.01, f'{int(v):,}',
                    ha='center', va='bottom', fontsize=9)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info("  wrote %s", path)


def _upset(source_sets, title, path):
    """UpSet plot of interaction-set overlap across sources (needs upsetplot)."""
    source_sets = {k: v for k, v in source_sets.items() if v}
    if len(source_sets) < 2:
        logger.info("  skipping %s: fewer than 2 non-empty sources", os.path.basename(path))
        return
    try:
        from upsetplot import from_memberships, UpSet
    except ImportError:
        logger.warning("  upsetplot not installed; skipping %s", os.path.basename(path))
        return

    all_ixns = set().union(*source_sets.values())
    memberships = [[src for src, s in source_sets.items() if ixn in s] for ixn in all_ixns]

    fig = plt.figure(figsize=(12, 6))
    with warnings.catch_warnings():
        # upsetplot 0.9 emits pandas-3.0 chained-assignment / downcasting
        # FutureWarnings internally (in both from_memberships and plot).
        warnings.simplefilter('ignore', FutureWarning)
        data = from_memberships(memberships)
        UpSet(data, subset_size='count', show_counts=True,
              sort_by='cardinality', min_subset_size=1).plot(fig)
    fig.suptitle(title, fontsize=14)
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info("  wrote %s", path)


def _metabolite_gene_figures():
    if not os.path.exists(config.METABOLITE_GENE_PKN):
        return
    mg = pd.read_csv(config.METABOLITE_GENE_PKN, sep='\t')
    if mg.empty:
        return

    # Per-source unique-interaction counts (Human-GEM dist1/2 merged).
    merged = mg.copy()
    merged['Source'] = merged['Source'].replace(_MG_MERGE)
    counts = merged.groupby('Source').apply(
        lambda x: len(_undirected_pairs(x, 'Metabolite', 'Gene'))).astype(int)
    _barplot_counts(counts, 'Metabolite-gene interactions per source',
                    'Unique interactions',
                    _fig_path('metabolite_gene_interactions_per_source.png'))

    # UpSet overlap across the individual source databases.
    present = [s for s in _MG_SOURCES if s in set(mg['Source'])]
    source_sets = {s: _undirected_pairs(mg[mg['Source'] == s], 'Metabolite', 'Gene')
                   for s in present}
    _upset(source_sets, 'Metabolite-gene interaction overlap across databases',
           _fig_path('metabolite_gene_overlap_upset.png'))


def _ppi_figures():
    if not os.path.exists(config.PPI_NETWORK_FILE):
        return
    ppi = pd.read_csv(config.PPI_NETWORK_FILE, sep='\t')
    if ppi.empty:
        return

    counts = ppi.groupby('Source').apply(
        lambda x: len(_undirected_pairs(x, 'GeneA', 'GeneB'))).astype(int)
    _barplot_counts(counts, 'Protein-protein interactions per source',
                    'Unique interactions', _fig_path('PPI_interactions_per_source.png'))

    source_sets = {s: _undirected_pairs(ppi[ppi['Source'] == s], 'GeneA', 'GeneB')
                   for s in sorted(set(ppi['Source']))}
    _upset(source_sets, 'Protein-protein interaction overlap across databases',
           _fig_path('PPI_overlap_upset.png'))


def _hptm_figures():
    if not os.path.exists(config.HPTM_NETWORK_FILE):
        return
    hptm = pd.read_csv(config.HPTM_NETWORK_FILE, sep='\t')
    if hptm.empty:
        return

    # Edges per source (QuickGO / UniProtKB_GO).
    src_counts = hptm.groupby('Source').apply(
        lambda x: len(_undirected_pairs(x, 'Enzyme', 'Mark'))).astype(int)
    _barplot_counts(src_counts, 'Enzyme-histone-mark interactions per source',
                    'Unique interactions', _fig_path('hPTM_interactions_per_source.png'))

    # UpSet: source overlap of enzyme->mark edges.
    source_sets = {s: _undirected_pairs(hptm[hptm['Source'] == s], 'Enzyme', 'Mark')
                   for s in sorted(set(hptm['Source']))}
    _upset(source_sets, 'Enzyme-histone-mark overlap across sources',
           _fig_path('hPTM_overlap_upset.png'))

    # Edges per activity (writer / eraser / reader), if the column is present.
    if 'Activity' in hptm.columns:
        act = hptm[['Enzyme', 'Mark', 'Activity']].drop_duplicates()
        act_counts = act['Activity'].value_counts()
        _barplot_counts(act_counts, 'Enzyme-histone-mark interactions by activity',
                        'Interactions', _fig_path('hPTM_interactions_by_activity.png'))


def _summary_figure():
    """Node and edge totals of the final combined PKN, per type (log scale)."""
    if not os.path.exists(config.FINAL_PKN_FILE):
        return
    pkn = pd.read_csv(config.FINAL_PKN_FILE, sep='\t')
    if pkn.empty:
        return

    mg = pkn[pkn['Type'] == 'metabolite-gene']
    ppi = pkn[pkn['Type'] == 'PPI']
    hptm = pkn[pkn['Type'] == 'histone-modification']
    metabolites = set(mg['Node1'])
    genes = set(mg['Node2']) | set(ppi['Node1']) | set(ppi['Node2']) | set(hptm['Node1'])
    marks = set(hptm['Node2'])

    counts = pd.Series({
        'Metabolites': len(metabolites),
        'Genes': len(genes),
        'Histone marks': len(marks),
        'Metabolite-gene edges': len(mg),
        'PPI edges': len(ppi),
        'Histone-mod. edges': len(hptm),
    })
    _barplot_counts(counts, 'LemonIte PKN composition',
                    'Count', _fig_path('PKN_composition.png'), log_scale=True)


def make_figures():
    """Generate all PKN summary figures. Best-effort; never raises."""
    try:
        os.makedirs(config.FIGURES_DIR, exist_ok=True)
        logger.info("Generating PKN figures -> %s", config.FIGURES_DIR)
        _metabolite_gene_figures()
        _ppi_figures()
        _hptm_figures()
        _summary_figure()
        logger.info("PKN figures complete")
    except Exception:  # noqa: BLE001
        logger.exception("Figure generation failed (build itself is unaffected)")
