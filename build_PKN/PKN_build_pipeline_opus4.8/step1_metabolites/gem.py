"""
Human-GEM metabolic-model metabolite-gene links (NetworkX graph traversal).

Builds a directed metabolic reaction graph from the Human-GEM model, excluding
ubiquitous currency metabolites. For each metabolite (matched by ChEBI id) it
finds genes reachable within ``distance`` reaction hops and records the actual
reaction *path* connecting metabolite -> gene, for provenance:

  * at distance 1 a single reaction R connects the metabolite to the gene, so
    only the first reaction column is filled;
  * at distance 2 the connection is metabolite -> intermediate metabolite ->
    gene and spans *two* reactions, so both reaction columns are filled — a
    single URL cannot represent a distance-2 edge, which is why the two
    reaction ids are tracked (Reaction_ID_1, Reaction_ID_2) with a Metabolic
    Atlas URL for each.

Re-derived from Collect_PKNdata_metabolites.ipynb cells 52-57.
"""

import logging
import re

import networkx as nx
import numpy as np
import pandas as pd

import config

logger = logging.getLogger('pkn.gem')


def _build_graph():
    """Build the metabolic reaction graph and the supporting lookup tables."""
    model = pd.read_csv(config.GEM_PATH, sep='\t')
    model_metabolites = pd.read_csv(config.GEM_METABOLITES_PATH, sep='\t')
    model_genes = pd.read_csv(config.GEM_GENES_PATH, sep='\t')

    metabolites_full = set(model_metabolites['mets'].tolist())

    exclude = {f'{m}{s}' for s in config.GEM_COMPARTMENT_SUFFIXES
               for m in config.GEM_UBIQUITOUS_METABOLITES}

    G = nx.DiGraph()
    for _, row in model.iterrows():
        formula = row['Formula']
        rxn = row['Rxn name']
        if '<=>' in formula:
            reactants, products = formula.split(' <=> ')
            reversible = True
        else:
            reactants, products = formula.split(' -> ')
            reversible = False
        reactants = reactants.strip().split(' + ')
        products = products.strip().split(' + ')
        for met in reactants + products:
            if met not in exclude and met in metabolites_full:
                G.add_node(met)
        for r in reactants:
            for p in products:
                if r in exclude or p in exclude:
                    continue
                G.add_edge(r, p, reaction=rxn, reversible=reversible)
                if reversible:
                    G.add_edge(p, r, reaction=rxn, reversible=True)

    # ChEBI -> metabolite base name
    model_metabolites['metChEBIID'] = model_metabolites['metChEBIID'].str.replace('CHEBI:', '', regex=False)
    chebi_to_name = dict(zip(model_metabolites['metChEBIID'], model_metabolites['metsNoComp']))

    # Ensembl gene -> symbol
    ensembl_to_symbol = dict(zip(model_genes['genes'], model_genes['geneSymbols']))

    # Reaction -> GPR gene string (cleaned of operators)
    model_clean = model.replace(
        to_replace=[r'\+', '->', '<-', '<=>', '=', '=>', '<=', 'and', 'or', r'\(', r'\)'],
        value='', regex=True)
    model_clean = model_clean.dropna(subset=['Gene-reaction association'])
    reaction_to_genes = dict(zip(model_clean['Rxn name'], model_clean['Gene-reaction association']))

    return G, metabolites_full, chebi_to_name, ensembl_to_symbol, reaction_to_genes


def _reaction_genes(rxn, ensembl_to_symbol, reaction_to_genes):
    """Gene symbols catalysing ``rxn`` (via its GPR association)."""
    gene_str = reaction_to_genes.get(rxn)
    if not gene_str:
        return []
    return [ensembl_to_symbol[g] for g in gene_str.split() if g in ensembl_to_symbol]


def _incident_reactions(G, node):
    """(neighbour, reaction) for every reaction incident to ``node``, treated as
    undirected — matching the original ego-graph traversal (undirected=True), which
    reaches reactions regardless of the stored edge direction. Reversible reactions
    already carry both directions; this also walks irreversible ones backward."""
    seen = []
    for _, nbr, d in G.out_edges(node, data=True):
        if 'reaction' in d:
            seen.append((nbr, d['reaction']))
    for nbr, _, d in G.in_edges(node, data=True):
        if 'reaction' in d:
            seen.append((nbr, d['reaction']))
    return seen


def _genes_within_distance(chebi, distance, G, metabolites_full, chebi_to_name,
                           ensembl_to_symbol, reaction_to_genes, store):
    """Genes reachable within ``distance`` hops; record the connecting reaction path.

    ``store`` maps (chebi, gene_symbol) -> set of reaction-path tuples. A path is
    ``(R1,)`` at distance 1 (one reaction connects the metabolite to the gene) or
    ``(R1, R2)`` at distance 2 (metabolite -R1-> intermediate -R2-> gene). Tracking
    the full path is what lets a distance-2 edge keep *both* reaction ids rather
    than a single, unrepresentative URL.
    """
    if pd.isna(chebi):
        return np.nan
    chebi = str(int(chebi))
    if chebi not in chebi_to_name:
        return np.nan
    base = chebi_to_name[chebi]
    pattern = re.compile(r'^' + re.escape(base) + r'[a-z]$')
    nodes = [m for m in metabolites_full if pattern.match(m)]
    if not nodes:
        return np.nan

    genes = set()

    def _record(symbol, path):
        genes.add(symbol)
        store.setdefault((chebi, symbol), set()).add(path)

    for start in nodes:
        if start not in G:
            continue

        # Distance-1: reactions directly incident to the metabolite. Any gene
        # catalysing such a reaction is one hop away.
        first_edges = _incident_reactions(G, start)
        for _, r1 in first_edges:
            for symbol in _reaction_genes(r1, ensembl_to_symbol, reaction_to_genes):
                _record(symbol, (r1,))

        if distance < 2:
            continue

        # Distance-2: metabolite -R1-> intermediate -R2-> gene. Walk one further
        # hop from each intermediate and record the (R1, R2) reaction pair.
        for intermediate, r1 in first_edges:
            if intermediate not in G:
                continue
            for _, r2 in _incident_reactions(G, intermediate):
                if r2 == r1:
                    continue
                for symbol in _reaction_genes(r2, ensembl_to_symbol, reaction_to_genes):
                    _record(symbol, (r1, r2))

    return '|'.join(sorted(genes)) if genes else np.nan


def run(df):
    """Build Human-GEM dist1 + dist2 interactions; write processed + link files."""
    G, metabolites_full, chebi_to_name, ensembl_to_symbol, reaction_to_genes = _build_graph()
    logger.info("GEM graph: %d nodes, %d edges", G.number_of_nodes(), G.number_of_edges())

    store_d1, store_d2 = {}, {}
    results = {}
    for dist, store, key in [(1, store_d1, 'Human1_GEM_dist1'),
                             (2, store_d2, 'Human1_GEM_dist2')]:
        col = df['ChEBI'].apply(
            lambda c: _genes_within_distance(c, dist, G, metabolites_full, chebi_to_name,
                                             ensembl_to_symbol, reaction_to_genes, store))
        df_proc = df[['HMDB', 'ChEBI']].copy()
        df_proc[key] = col
        df_proc.to_csv(config.DB_OUTPUT_FILES[key], index=False, sep='\t')
        results[key] = df_proc
        logger.info("GEM %s: %d metabolites with genes", key, col.notna().sum())

    # Reaction-level provenance links (dist1 + dist2 in one file)
    chebi_df = df.dropna(subset=['ChEBI']).copy()
    chebi_df['ChEBI_str'] = chebi_df['ChEBI'].astype(int).astype(str)
    chebi_to_hmdb = chebi_df.groupby('ChEBI_str')['HMDB'].apply(list).to_dict()

    # One row per (HMDB, gene, reaction-path, distance). A distance-1 edge is
    # connected by a single reaction (Reaction_ID_1); a distance-2 edge spans two
    # reactions, metabolite -R1-> intermediate -R2-> gene, so BOTH reaction ids are
    # kept (Reaction_ID_1, Reaction_ID_2), each with its own Metabolic Atlas URL.
    # A single URL cannot represent the two-reaction distance-2 connection, which is
    # the reason for the two-column schema. URL_dataset/URL_paper are retained for
    # general provenance.
    rxn_url = lambda r: config.URL_TEMPLATES['GEM_reaction'].format(reaction_id=r)
    link_rows = []
    for dist, store in [(1, store_d1), (2, store_d2)]:
        for (chebi, gene), paths in store.items():
            for hmdb in chebi_to_hmdb.get(chebi, []):
                for path in paths:
                    r1 = path[0]
                    r2 = path[1] if len(path) > 1 else None
                    link_rows.append({
                        'HMDB': hmdb, 'ChEBI': chebi, 'Gene': gene,
                        'Reaction_ID_1': r1,
                        'Reaction_ID_2': r2,
                        'Distance': dist,
                        'URL_1': rxn_url(r1),
                        'URL_2': rxn_url(r2) if r2 is not None else None,
                        'URL_dataset': config.URL_TEMPLATES['GEM_dataset'],
                        'URL_paper': config.URL_TEMPLATES['GEM_paper'],
                    })
    pd.DataFrame(link_rows).to_csv(config.URL_FILES['Human1_GEM_dist1'], index=False, sep='\t')
    logger.info("GEM: %d reaction-path links", len(link_rows))
    return results
