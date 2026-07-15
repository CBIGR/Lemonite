"""
Step 2b: enzyme - histone-modification interactions.

Builds edges between a histone-modifying protein (writer / eraser / reader) and
the specific histone mark it acts on or binds, from Gene Ontology molecular-
function annotations retrieved via two REST APIs (QuickGO and UniProtKB). This is
enzyme-substrate specificity, distinct from the metabolite-gene (step 1) and
protein-protein (step 2) edge types. See hptm_integration.build_hptm_network.
"""
