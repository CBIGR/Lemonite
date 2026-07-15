# Phospho-PKN resource review (online research, 2026-07-03)

Which prior-knowledge resources — beyond OmniPath and CollecTRI — are worth adding to the
knowledge graph for a **phosphosite / phospho-regulator** layer. Goal: known
phosphoprotein/kinase → target-gene (and phospho-TF → target) edges to validate phospho
regulator→module assignments, mirroring the metabolite-gene evaluation.

## TL;DR — recommended additions (in priority order)

1. **OmniPath enzyme–substrate (kinase–substrate) network** — *do this first; highest value,
   lowest effort, license-clean.* Available directly in the `.sif` via
   `decoupleR::get_ksn_omnipath()` / `OmnipathR::enzyme_substrate()`. It already
   **aggregates 18 resources** (115,215 enzyme-PTM interactions), including PhosphoSitePlus,
   SIGNOR, DEPOD, dbPTM, phospho.ELM, HPRD, KEA, PhosphoNetworks, NetworKIN(via MIMP),
   Reactome, NCI-PID. So one pull subsumes most of what we'd otherwise fetch individually.
2. **Johnson Kinase Library atlas** — Ser/Thr kinome (Nature 2023) + Tyr kinome (Nature 2024).
   Motif-based *predicted* kinase → phosphosite assignments covering ~84% of active human
   kinases and (in principle) *every* reported human phosphosite. **Site-resolved** and
   complementary to curated data (fills coverage gaps). Distributed via the PhosphoSitePlus
   "Kinase Library". License: academic/non-commercial OK (see licensing note). Adds
   predictive kinase→substrate edges keyed on the exact residue in our phosphosite IDs.
3. **SIGNOR 4.0 (2025 update, phosphorylation-focused)** — manually curated *causal*
   phosphorylation/dephosphorylation with residues + mechanism; clean **TSV API**
   (`https://signor.uniroma2.it`). Partly inside OmniPath already, but the direct pull is
   current and adds sign/mechanism. Good for a curated, high-confidence subset.
4. **DEPOD** — curated human **phosphatase**→substrate + dephosphorylation sites. The
   phosphatase counterpart to kinase-substrate (also partly in OmniPath). Optional; adds the
   "eraser" direction.

Already have locally (no fetch): **CollecTRI** (`nextflow/PKN/CollecTRI_network.txt`) for
phospho-TF → transcriptional-target edges, and the **PPI** layer of `Lemonite_PKN.tsv` for
phosphoprotein↔module-gene physical interactions.

## What each resource contributes

| Resource | Edge | Node1 → Node2 | Site-resolved? | Access | License |
|---|---|---|---|---|---|
| OmniPath enzsub | kinase/enzyme → substrate | enzyme gene → substrate gene (+residue) | yes (substrate side) | `.sif` OmnipathR / decoupleR, no manual DL | free (aggregates curated+predicted) |
| Kinase Library (Johnson) | predicted kinase → phosphosite | kinase → substrate gene @residue | yes | PhosphoSitePlus download | non-commercial/academic |
| SIGNOR 4.0 | curated causal (de)phosphorylation | enzyme → target (+residue, +sign) | yes | REST TSV API | free academic |
| DEPOD | phosphatase → substrate | phosphatase → substrate (+site) | yes | download | free |
| CollecTRI *(have)* | TF → transcriptional target | TF → gene (+mode) | n/a | local file | free |
| PPI in Lemonite_PKN *(have)* | physical interaction | protein → protein | n/a | local file | — |

## Coverage note (why OmniPath first, avoid duplication)

OmniPath's enzyme-substrate domain lists these 18 source resources:
`DEPOD, HPRD, HPRD_MIMP, KEA, Li2012, MIMP, PhosphoNetworks, PhosphoSite (=PhosphoSitePlus),
PhosphoSite_MIMP, ProtMapper (+BEL/NCI-PID/REACH/RLIMS-P/Reactome/Sparser variants),
SIGNOR, SIGNOR_ProtMapper, dbPTM, phosphoELM (+MIMP variant)`.
→ Pulling OmniPath gives us DEPOD, PhosphoSitePlus, SIGNOR, dbPTM, phospho.ELM, HPRD,
NetworKIN(MIMP) etc. in one shot. The **only genuinely additive** external fetch is the
**Kinase Library atlas** (motif-based genome-wide predictions not in OmniPath enzsub) and,
if we want the freshest curated causal layer, a **direct SIGNOR 4.0** pull.

## Licensing / feasibility for an unattended run

- **OmniPath / decoupleR** — free, programmatic, already in the `.sif`. ✅ safe to auto-run.
- **SIGNOR 4.0** — free academic; REST TSV API. ✅ auto-fetchable; be polite (rate-limit,
  cache).
- **DEPOD** — free; file download. ✅
- **PhosphoSitePlus bulk / Kinase Library** — *web use free; **bulk data download requires
  agreeing to non-commercial terms** (a click-through / registration), and commercial use is
  a paid license.* ⚠️ **Do NOT auto-scrape.** For the overnight run, obtain Kinase Library
  atlas edges via OmniPath where present; if the full atlas is wanted, flag it for the user
  to download manually under their academic agreement, then point the builder at the file.

## How this maps onto the build pipeline

Extend `build_PKN/PKN_build_pipeline_opus4.8/` with a new **`step2c_phospho/`** package
(mirroring `step2b_hPTM/`): fetch → write `phospho_pkn_network.tsv`
(`Node1, Node2, Source, Type[, Residue, url]`, `Type ∈ {kinase-substrate,
phospho-TF-target, phosphatase-substrate}`). Because the user chose a **standalone**
phospho-PKN, we do NOT merge it into `Lemonite_PKN.tsv`; we keep the output beside the
phospho analysis and load it only in `evaluate_phospho_against_pkn.py`. Use the pipeline's
existing `utils/cache.py` (resumable) + `utils/api_retry.py` (429-aware) so an overnight run
survives interruptions.

## Sources
- OmniPath 2025 (NAR): https://academic.oup.com/nar/article/54/D1/D652/8326458
- decoupleR `get_ksn_omnipath`: https://saezlab.github.io/decoupleR/reference/get_ksn_omnipath.html
- OmnipathR `get_enzsub_resources`: https://r.omnipathdb.org/reference/get_enzsub_resources.html
- Johnson et al. Ser/Thr kinome atlas (Nature 2023): https://www.nature.com/articles/s41586-022-05575-3
- Human tyrosine kinome specificity (Nature 2024): https://www.nature.com/articles/s41586-024-07407-y
- The Kinase Library: https://kinase-library.phosphosite.org/kinase-library/
- SIGNOR 4.0 (NAR 2025): https://academic.oup.com/nar/article/54/D1/D682/8324960
- DEPOD: https://depod.bioss.uni-freiburg.de/
- PhosphoSitePlus licensing: https://www.phosphosite.org/staticLicensing.action
