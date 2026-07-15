# Phosphosite prior-knowledge evaluation & visualization — comprehensive report

**Case study:** Wang GBM (CPTAC), Lemonite phospho-regulator extension
**Branch:** `phospho-regulator-analysis` (nothing committed/pushed)
**Date:** 2026-07-04

---

## 1. Objective

Do for **phosphosite regulators** what the Lemonite pipeline already does for metabolites:
test whether a phosphosite assigned to a co-expression module has *known* prior-knowledge
connections to that module's genes — and visualize those known interactions on the module
heatmaps and subnetwork figures. This complements the earlier work that added
phosphoproteomics as a 4th regulator layer (see `README.md`); here we ask **which phospho
assignments are supported by external prior knowledge.**

The current shared PKN (`nextflow/PKN/Lemonite_PKN.tsv`) contains **no phosphosite edges**
(only PPI + metabolite-gene), so a phospho prior-knowledge network had to be built first.

## 2. What was built

### 2.1 Standalone phospho-PKN — `phospho_pkn.tsv` (65,717 edges)

Kept **standalone** (the shared PKN is not modified), same `Node1/Node2/Source/Type` schema:

| Edge type | Node1 → Node2 | Source | Edges |
|---|---|---|---|
| `kinase-substrate` | kinase → substrate gene | **OmniPath KSN** (`decoupleR::get_ksn_omnipath()`) | 20,319 |
| `phospho-TF-target` | TF → transcriptional target | **CollecTRI** (local) | 45,398 |

Plus edge type **C = PPI**, read at evaluation time from the main PKN (phosphoprotein nodes
only; 124,635–166,894 pairs depending on variant) — not copied into `phospho_pkn.tsv`.

**Why OmniPath KSN is the key add (resource research, see `PHOSPHO_PKN_RESOURCES.md`):**
its enzyme-substrate layer already aggregates ~18 curated + predicted databases —
PhosphoSitePlus, SIGNOR, DEPOD, dbPTM, phospho.ELM, HPRD, NetworKIN(MIMP), etc. — so one
license-clean pull (39,350 raw site-level edges → 20,319 gene-level after collapsing the
`GENE_residue` substrate ids) subsumes what would otherwise be many separate integrations.
The Johnson **Kinase Library** atlas is the only genuinely-additive external resource beyond
OmniPath but its bulk download is non-commercial-gated → left for a manual academic download.

### 2.2 Evaluation — `evaluate_phospho_against_pkn.py`

Adapts the pipeline's metabolite hypergeometric test (`evaluate_against_PKN.py:499`):
per module, count observed (phosphoprotein × module-gene) pairs that are known phospho-PKN
edges vs. expected under global density; hypergeometric SF + Benjamini-Hochberg FDR.
Matching is **gene-level** (phosphosite `GENE:RefSeq:Residue` → phosphoprotein gene).
`kinase-substrate` is scored direction-agnostically (our phosphoprotein may be the kinase
*or* the substrate); `phospho-TF-target` is directional (phosphoprotein must be the TF).

Run for **both variants**, scoring each edge type separately and combined.

## 3. Results

| | Top-2000 variable | All phosphosites |
|---|---|---|
| Modules tested | 46 | 46 |
| **PPI** — modules FDR<0.05 | **15/46** | **16/46** |
| kinase-substrate — modules FDR<0.05 | 0 | 3 |
| phospho-TF-target — modules FDR<0.05 | 0 | 2 |
| combined — modules FDR<0.05 | 15 | 15 |
| Supported same-module phospho→gene edges | 602 (596 PPI, 6 KS) | 731 (697 PPI, 20 TF, 14 KS) |
| Cross-module phospho-TF→target links | 82 | 261 |

**Interpretation.**
- **PPI is the dominant, strongly-significant evidence:** phosphoproteins assigned to a
  module physically interact with that module's genes far more than chance (e.g. top2000
  module 8: 38 observed, 6.9× fold, FDR 1e-18; module 39: 12.7× fold). 15–16 of 46 modules
  are significantly enriched.
- **Kinase-substrate and phospho-TF-target are sparser but higher-specificity.** They reach
  significance only in the larger `all_phosphosites` regulator pool (3 and 2 modules), which
  is expected — curated directional edges are rarer. Concrete supported examples:
  `EGFR:S1153 → GPNMB, DOK2` (kinase-substrate, module 13); `KCNQ2:S476 → CAMK2A` (module 0);
  `ADD2 → PRKCG/PRKCZ`, `PLCB1 → PRKCG` (module 6).
- **Cross-module phospho-TF links are the mechanistically interesting finding.** The
  phospho-TFs assigned as regulators (SOX2, SOX10, WWTR1/TAZ, PML, NFIA, MYT1 …; the
  all-phosphosites set adds RB1, SMARCA4, DNMT3A, SIRT1, STAT1, FOXO1, YY1, CHD4 — core GBM
  chromatin/TF machinery) have **almost no CollecTRI targets within their own assigned
  module, but 82 (top2000) / 261 (all) known targets in *other* network modules.** i.e. a
  phospho-event on SOX2 co-varies with module M, while SOX2's known transcriptional targets
  sit in modules M'. This is a genuine, testable hypothesis about phospho-mediated
  cross-module regulation, surfaced by the analysis rather than assumed.

Output tables per variant in `<variant>/run/PhosphoPKN_Evaluation/`:
`Phospho_PKN_enrichment_{kinase-substrate,phospho-TF-target,PPI,combined}.csv`,
`Phospho_PKN_supported_edges.csv`, `Phospho_TF_target_links.csv`.

## 4. Visualization (both variants, both figure types)

Known phospho interactions now render on the module figures using the **same mechanism** as
metabolite/PPI/HumanNet interactions (not a parallel path):

- **Module heatmaps** (`module_viewer.py`): the evaluation writes
  `phospho_LemoniteKG_interactions.mvf` (same `.mvf` format, `COLOR=ORANGE`), and a guarded,
  additive edit loads it and draws an **orange phospho-interaction panel** on the right of
  each module's expression heatmap — columns are the module's phosphosites, orange cells mark
  the module genes each phosphosite has a known phospho-PKN edge to. *Verified visually*
  (e.g. module 0: ANK2/GNG3/KCNQ2/… × CAMK2A/HPCA/RYR2/SPTB…). All 46 heatmaps regenerated
  per variant.
- **Subnetwork figures** (`create_subnetwork_graphs.py`): fed a **viz-only combined PKN**
  (`Lemonite_PKN_plus_phospho.tsv` = main PKN + phospho edges; the shared PKN is untouched).
  Added two edge categories — `Kinase_substrate` (orange solid) and `Phospho_TF_target`
  (dark-orange dashed) — with legend entries, so phospho regulator→target edges render
  distinctly alongside metabolite/PPI edges. Phospho nodes appear in orange. 38–39
  subnetworks regenerated per variant.

Guarded edits mean pipeline runs *without* a phospho `.mvf`/phospho-PKN are unaffected.

## 5. Files (all on branch `phospho-regulator-analysis`)

New in `Wang_GBM/Lemonite_phospho/`:
- `build_phospho_pkn.py` → `phospho_pkn.tsv` (+ `phospho_pkn_cache/omnipath_ksn.tsv`)
- `evaluate_phospho_against_pkn.py` → per-variant `PhosphoPKN_Evaluation/` + phospho `.mvf`
- `Lemonite_PKN_plus_phospho.tsv` (viz-only combined PKN)
- `PHOSPHO_PKN_RESOURCES.md` (resource review), `PHOSPHO_PKN_ANALYSIS_PLAN.md` (plan/resume),
  this report.

Edited (guarded/additive, isolated on branch):
- `nextflow/scripts/module_viewer.py` — phospho heatmap panel
- `nextflow/scripts/create_subnetwork_graphs.py` — phospho edge categories + styles
- `Wang_GBM/Lemonite_phospho/run_downstream.sh` — `PKN_VIZ` for subnetworks

Results per variant under
`…/results/LemonTree/noProteomics_percentile2_divide_by_sum_phospho/<variant>/run/`:
`PhosphoPKN_Evaluation/`, `heatmaps/` (46, with phospho panel),
`Networks/subnetworks_graphviz/` (with phospho edges),
`ModuleViewer_files/phospho_LemoniteKG_interactions.mvf`.

## 6. Caveats & next steps

- **PPI dominates** because it's undirected and dense; the directional kinase-substrate /
  TF-target layers are the biologically specific signal and are best read from the supported-
  edge and cross-module TF-link tables rather than the enrichment counts alone.
- **Namespace:** matching is on HGNC-style symbols; phosphosite RefSeq gene symbols aligned
  well with module symbols and OmniPath/CollecTRI (spot-checked), but a formal alias map
  (the PKN build ships `ensembl_mapping_*`) would tighten coverage.
- **Optional additions:** (a) the Johnson Kinase Library atlas (site-resolved predicted
  kinases for our exact phosphosites — manual academic download); (b) a shuffled-regulator
  null to put an empirical p on the PPI enrichment; (c) promote the builder to a first-class
  `step2c_phospho/` package in `build_PKN/PKN_build_pipeline_opus4.8/` (resumable, monthly).
