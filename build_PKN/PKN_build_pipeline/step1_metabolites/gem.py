"""
Human-GEM metabolic pathway retriever.

Builds NetworkX graph from Human-GEM model and finds genes within pathway distance.
"""

import pandas as pd
import numpy as np
import networkx as nx
import re
import logging
from typing import List, Dict
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pipeline import DatabaseRetriever
import config


class GEMRetriever(DatabaseRetriever):
    """
    Retrieve metabolite-gene interactions from Human-GEM metabolic model.
    
    Builds a NetworkX graph of metabolic reactions and finds genes within
    a specified pathway distance from each metabolite.
    """
    
    def __init__(self, distance=1, cache_file=None):
        """
        Initialize GEM retriever.
        
        Parameters:
        -----------
        distance : int
            Pathway distance (number of reactions away)
        """
        db_name = f'Human1_GEM_dist{distance}'
        
        if cache_file is None:
            cache_file = config.DB_OUTPUT_FILES.get(db_name)
        
        super().__init__(db_name=db_name, cache_file=cache_file)
        
        self.distance = distance
        self._initialized = False
        self.graph = None
        self.chebi_to_name = {}
        self.ensembl_mappings = {}
        self.reaction_to_genes = {}
        self.metabolites_full_name = []
    
    def _initialize_graph(self):
        """
        Initialize Human-GEM metabolic network.
        
        Builds NetworkX directed graph from model files:
        - Reactions with metabolites as nodes
        - Edges represent reactions
        - Gene-reaction associations
        """
        if self._initialized:
            return
        
        self.logger.info("="*80)
        self.logger.info(f"Initializing Human-GEM metabolic network (distance={self.distance})...")
        self.logger.info("="*80)
        
        # Load model files
        self.logger.info(f"Loading model from {config.GEM_PATH}")
        model = pd.read_csv(config.GEM_PATH, sep='\t')
        
        self.logger.info(f"Loading metabolites from {config.GEM_METABOLITES_PATH}")
        model_metabolites = pd.read_csv(config.GEM_METABOLITES_PATH, sep='\t')
        
        self.logger.info(f"Loading genes from {config.GEM_GENES_PATH}")
        model_genes = pd.read_csv(config.GEM_GENES_PATH, sep='\t')
        
        # Build ChEBI ID → metabolite name mapping
        model_metabolites['metChEBIID'] = model_metabolites['metChEBIID'].str.replace('CHEBI:', '')
        self.chebi_to_name = dict(zip(model_metabolites['metChEBIID'], model_metabolites['metsNoComp']))
        
        # Build Ensembl → gene symbol mapping
        self.ensembl_mappings = dict(zip(model_genes['genes'], model_genes['geneSymbols']))
        
        # Get list of all metabolites
        self.metabolites_full_name = list(set(model_metabolites['mets'].tolist()))
        self.logger.info(f"✓ Human-GEM contains {len(self.metabolites_full_name)} unique metabolites")
        
        # Define ubiquitous metabolites to exclude (cofactors, common currency)
        metabolites_ubiquitous = [
            'MAM02039', 'MAM02040', 'MAM01371', 'MAM02751', 'MAM01285',
            'MAM02555', 'MAM0254', 'MAM02630', 'MAM01597', 'MAM02552',
            'MAM02553', 'MAM02554', 'MAM02759', 'MAM02046', 'MAM01334',
            '2 MAM02039', '2 MAM02040'
        ]
        
        # Add compartment suffixes to exclusion list
        metabolites_to_exclude = []
        for suffix in ['c', 'g', 'l', 'm', 'n', 'x', 'r', 'e']:
            for metabolite in metabolites_ubiquitous:
                metabolites_to_exclude.append(f'{metabolite}{suffix}')
        
        self.logger.info(f"✓ Excluding {len(metabolites_to_exclude)} ubiquitous metabolites")
        
        # Build NetworkX graph
        self.logger.info("Building metabolic network graph...")
        self.graph = nx.DiGraph()
        
        for _, row in model.iterrows():
            rxn_name = row['Rxn name']
            formula = row['Formula']
            
            # Parse reaction formula
            if '<=>' in formula:
                reactants, products = formula.split(' <=> ')
                reversible = True
            else:
                reactants, products = formula.split(' -> ')
                reversible = False
            
            reactants = reactants.strip().split(' + ')
            products = products.strip().split(' + ')
            
            # Add nodes (exclude ubiquitous metabolites)
            for metabolite in reactants + products:
                if metabolite not in metabolites_to_exclude and metabolite in self.metabolites_full_name:
                    self.graph.add_node(metabolite)
            
            # Add edges
            for reactant in reactants:
                for product in products:
                    if reactant not in metabolites_to_exclude and product not in metabolites_to_exclude:
                        if reversible:
                            self.graph.add_edge(reactant, product, reaction=rxn_name, reversible=True)
                            self.graph.add_edge(product, reactant, reaction=rxn_name, reversible=True)
                        else:
                            self.graph.add_edge(reactant, product, reaction=rxn_name, reversible=False)
        
        self.logger.info(f"✓ Graph built: {self.graph.number_of_nodes()} nodes, {self.graph.number_of_edges()} edges")
        
        # Build reaction → genes mapping
        model_clean = model.replace(
            to_replace=['\+', '->', '<-', '<=>', '=', '=>', '<='],
            value='',
            regex=True
        )
        model_clean = model_clean.replace(
            to_replace=['and', 'or', '\(', '\)'],
            value='',
            regex=True
        )
        model_filtered = model_clean.dropna(subset=['Gene-reaction association'])
        
        self.reaction_to_genes = dict(
            zip(model_filtered['Rxn name'], model_filtered['Gene-reaction association'])
        )
        
        self.logger.info(f"✓ Mapped {len(self.reaction_to_genes)} reactions to genes")
        self.logger.info("="*80)
        self.logger.info("Human-GEM initialization complete!")
        self.logger.info("="*80)
        
        self._initialized = True
    
    def _find_metabolites_with_regex(self, lst: List[str], regex_pattern: str) -> List[str]:
        """Match metabolite names with compartment suffixes using regex."""
        pattern = re.compile(regex_pattern)
        return [element for element in lst if pattern.match(element)]
    
    def get_genes_for_chebi(self, chebi_id: str) -> List[str]:
        """
        Get genes within pathway distance from a metabolite.
        
        Parameters:
        -----------
        chebi_id : str
            ChEBI ID (numeric only, without 'CHEBI:' prefix)
        
        Returns:
        --------
        list of str
            Gene symbols within specified distance
        """
        # Ensure graph is initialized
        self._initialize_graph()
        
        if pd.isna(chebi_id):
            return []
        
        # Convert to string
        try:
            chebi_id = str(int(float(chebi_id)))
        except (ValueError, TypeError):
            return []
        
        # Check if ChEBI is in model
        if chebi_id not in self.chebi_to_name:
            return []
        
        metabolite = self.chebi_to_name[chebi_id]
        
        # Find metabolite variants with compartment suffixes
        metabolites = self._find_metabolites_with_regex(self.metabolites_full_name, metabolite)
        
        if not metabolites:
            return []
        
        # Collect genes from all compartment variants
        genes = set()
        
        for node in metabolites:
            if node not in self.graph:
                continue
            
            try:
                # Create ego graph (subgraph within distance)
                ego_graph = nx.ego_graph(self.graph, node, radius=self.distance, undirected=True)
                
                # Collect all reactions in ego graph
                reactions_within_distance = {
                    data['reaction']
                    for u, v, data in ego_graph.edges(data=True)
                    if 'reaction' in data
                }
                
                # Get genes for these reactions
                for reaction in reactions_within_distance:
                    if reaction not in self.reaction_to_genes:
                        continue
                    
                    gene_str = self.reaction_to_genes[reaction]
                    if not gene_str:
                        continue
                    
                    # Parse genes from reaction association
                    for gene in gene_str.split():
                        if gene in self.ensembl_mappings:
                            genes.add(self.ensembl_mappings[gene])
            
            except (nx.NetworkXError, nx.NodeNotFound):
                continue
            except Exception as e:
                self.logger.warning(f"Unexpected error for ChEBI {chebi_id}, node {node}: {e}")
                continue
        
        return list(genes)
    
    def get_interactions(self, metabolites: List[Dict]) -> pd.DataFrame:
        """
        Get GEM interactions for metabolites.
        
        Parameters:
        -----------
        metabolites : list of dict
            Metabolite records with 'HMDB_ID' and 'ChEBI'
        
        Returns:
        --------
        pd.DataFrame
            Columns: ['HMDB_ID', 'Gene', 'Source']
        """
        # Try cache first
        cached = self.load_cache()
        if cached is not None:
            return cached
        
        # Initialize graph
        self._initialize_graph()
        
        # Process metabolites
        self.logger.info(f"Processing {len(metabolites)} metabolites for GEM interactions...")
        
        interactions = []
        metabolites_with_genes = 0
        
        for metabolite in metabolites:
            hmdb_id = metabolite.get('HMDB_ID')
            chebi_id = metabolite.get('ChEBI')
            
            if pd.isna(chebi_id):
                continue
            
            genes = self.get_genes_for_chebi(chebi_id)
            
            if genes:
                metabolites_with_genes += 1
                for gene in genes:
                    interactions.append({
                        'HMDB_ID': hmdb_id,
                        'Gene': gene,
                        'Source': self.db_name
                    })
        
        df = pd.DataFrame(interactions)
        df = df.drop_duplicates(subset=['HMDB_ID', 'Gene'])
        
        self.logger.info(
            f"Found {len(df)} GEM interactions for "
            f"{metabolites_with_genes} metabolites "
            f"(distance={self.distance})"
        )
        
        # Save cache
        self.save_cache(df)
        
        return df


if __name__ == '__main__':
    # Test the retriever
    import logging
    logging.basicConfig(level=logging.INFO)
    
    from utils.file_io import load_hmdb_metabolites
    
    # Load metabolites
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    print(f"Loaded {len(metabolites)} metabolites")
    
    # Test with first 100 metabolites
    print("\nTesting GEM retriever with first 100 metabolites...")
    test_metabolites = metabolites[:100]
    
    # Test distance 1
    retriever1 = GEMRetriever(distance=1)
    interactions1 = retriever1.get_interactions(test_metabolites)
    print(f"\nGEM dist=1 interactions (test): {len(interactions1)}")
    print(interactions1.head())
    
    # Test distance 2
    retriever2 = GEMRetriever(distance=2)
    interactions2 = retriever2.get_interactions(test_metabolites)
    print(f"\nGEM dist=2 interactions (test): {len(interactions2)}")
    print(interactions2.head())
