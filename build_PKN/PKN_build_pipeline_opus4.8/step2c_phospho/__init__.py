"""
Step 2c: phospho-regulator prior-knowledge network (STANDALONE).

Builds a phosphosite/phosphoprotein -> target-gene network for the phospho-regulator
analysis. Two edge types:

  * kinase-substrate   (Node1 = kinase gene, Node2 = substrate gene)   -- OmniPath enzsub
  * phospho-TF-target  (Node1 = TF gene,     Node2 = target gene)      -- CollecTRI

This is deliberately kept SEPARATE from the main LemonIte_PKN.tsv (step 3 does NOT merge it):
it is an independent prior-knowledge layer consumed by the phospho-regulator evaluation
(Wang_GBM/Lemonite_phospho/evaluate_phospho_against_pkn.py). Output: phospho_pkn.tsv with
columns Node1, Node2, Source, Type, Residue, URL.

See phospho_integration.build_phospho_network.
"""
