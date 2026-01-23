# Methods: Module-Regulator Network Visualization

## Network Construction and Node Selection

Module-regulator networks were constructed and visualized for both the Lloyd-Price IBD dataset and the Wang GBM dataset using the `Module_Network_Visualization.ipynb` scripts. The networks display regulatory relationships between gene expression modules identified by LemonTree and their predicted regulatory molecules (transcription factors and metabolites).

### Node Types and Selection Criteria

Networks consisted of two node types:

1. **Module nodes**: All gene co-expression modules identified by LemonTree were included as nodes. For the Lloyd-Price IBD dataset, 87 modules were included, while the Wang GBM dataset contained all modules identified in that analysis.

2. **Regulator nodes**: Regulatory molecules (transcription factors and metabolites for Lloyd-Price; transcription factors, metabolites, and lipids for Wang) were included based on a two-step filtering strategy:
   
   **Step 1 - Connectivity-based filtering**: Regulators were first selected based on the number of modules they regulate. Only regulators that regulate at least 3 modules (min_targets=3) were retained. Within each regulator type (transcription factors, metabolites, or lipids), a maximum of 10 regulators (max_regulators_per_type=10) were kept, prioritizing those with the highest number of target modules and highest average LemonTree regulatory scores.
   
   **Step 2 - Score-based enrichment**: To ensure inclusion of the strongest individual regulatory relationships, the top 5 highest-scoring regulator-module pairs (top_pairs_by_score=5) were added for each regulator type, based on their LemonTree regulatory scores. This ensured that highly specific regulatory relationships were represented even if the regulator did not meet the connectivity threshold.

### Network Composition

For the **Lloyd-Price IBD dataset**, the final filtered network contained:
- 87 module nodes
- A subset of transcription factor nodes (filtered as described above)
- A subset of metabolite nodes (filtered as described above)
- Regulatory edges connecting regulators to their target modules

For the **Wang GBM dataset**, the final filtered network contained:
- All identified module nodes
- A subset of transcription factor, metabolite, and lipid regulator nodes (filtered as described above)
- Regulatory edges connecting regulators to their target modules

The exact number of nodes in the final networks was determined by the filtering parameters and the specific regulatory relationships identified by LemonTree in each dataset. Network statistics (total nodes and edges) were reported in the analysis output for each dataset.

### Network Visualization

Networks were visualized using both interactive (Plotly) and publication-quality (Cytoscape) approaches. Modules were spatially organized based on their Gene Ontology (GO) biological process enrichment, grouping functionally related modules into categories including Immune, Metabolic, Cell Cycle, Signaling, Development, Transport, Adhesion/Migration, and Other. Module nodes were displayed as circles colored by their GO functional category using a color-blind friendly palette (Paul Tol's bright scheme). Regulator nodes were displayed as rectangles with an orange color gradient reflecting their connectivity (darker orange indicates higher connectivity), positioned near their target modules. Edge thickness and node size reflected the strength and importance of regulatory relationships.

This visualization approach enabled identification of functional patterns in the regulatory architecture, highlighting which types of regulators control specific biological processes and revealing the multi-omics regulatory landscape underlying each disease condition.

---

## Protein-Protein Interaction (PPI) Enrichment Analysis

### Overview

To evaluate the biological coherence of gene co-expression modules identified by LemonTree, we performed protein-protein interaction (PPI) enrichment analysis. This analysis tests whether genes within each module have more physical interactions than expected by chance, indicating functional relatedness beyond co-expression.

### PPI Network Construction

PPI data were extracted from the LemonIte Prior Knowledge Network (PKN), which integrates interactions from multiple databases. The PKN was filtered to retain only protein-protein interactions (edge type = "PPI"), excluding other interaction types (e.g., transcription factor-gene, metabolite-gene). For each module, all pairwise gene combinations were examined to identify PPIs documented in the PKN. Since PPI networks are undirected, both orientations of each interaction (gene1-gene2 and gene2-gene1) were stored for efficient lookup.

### Statistical Testing: Hypergeometric Enrichment

PPI enrichment was assessed using the hypergeometric distribution, which models sampling without replacement from a finite population. This test evaluates whether the observed number of PPIs within a module is significantly higher than expected by random chance.

**Hypergeometric Test Parameters:**

For each module, the hypergeometric test was parameterized as follows:

- **Population (M)**: Total possible gene pairs across all modules  
  $M = \frac{N_{genes} \times (N_{genes} - 1)}{2}$  
  where $N_{genes}$ is the total number of genes across all modules

- **Success states in population (n)**: Total number of gene pairs with documented PPIs in the PKN among all module genes

- **Sample size (N)**: Number of possible gene pairs within the module  
  $N = \frac{k \times (k - 1)}{2}$  
  where $k$ is the number of genes in the module

- **Observed successes (x)**: Number of PPIs observed within the module

The hypergeometric probability mass function is:

$$P(X = x) = \frac{\binom{n}{x} \binom{M-n}{N-x}}{\binom{M}{N}}$$

To test for enrichment (more PPIs than expected), we calculated the upper-tail p-value using the survival function:

$$p\text{-value} = P(X \geq x) = \sum_{i=x}^{\min(n,N)} \frac{\binom{n}{i} \binom{M-n}{N-i}}{\binom{M}{N}}$$

This represents the probability of observing at least $x$ PPIs by random chance, given the background PPI density across all module genes.

### Enrichment Metrics

For each module, we calculated:

1. **Expected PPIs**: The expected number of PPIs under the null hypothesis of random gene assignment:  
   $E[X] = N \times \frac{n}{M}$

2. **Fold Enrichment**: The ratio of observed to expected PPIs:  
   $\text{Fold Enrichment} = \frac{x}{E[X]}$  
   Values > 1 indicate enrichment; values < 1 indicate depletion.

3. **PPI Density**: The proportion of possible gene pairs within the module that have documented interactions:  
   $\text{PPI Density} = \frac{x}{N}$

### Multiple Testing Correction

To account for multiple hypothesis testing across all modules, raw p-values were adjusted using the Benjamini-Hochberg false discovery rate (FDR) procedure. Modules with FDR < 0.05 were considered significantly enriched for PPIs.

### Interpretation

Significant PPI enrichment indicates that module genes are more functionally related than expected by chance, supporting the biological validity of the co-expression modules. Higher fold enrichment values suggest tighter functional coupling, while PPI density reflects the overall connectivity within each module. Modules with significant PPI enrichment are likely to represent coherent biological processes or protein complexes.
