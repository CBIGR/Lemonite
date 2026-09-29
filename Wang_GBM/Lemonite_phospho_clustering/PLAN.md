# Phospho-as-clustering-input Lemonite variant — plan

**Goal (user, 2026-07-06).** A *second* phospho extension, distinct from
`Wang_GBM/Lemonite_phospho/` (where phospho was a regulator). Here:

- **Phosphoproteomics is the CLUSTERING INPUT** (primary matrix), not a regulator.
- **Regulators assigned to the phospho-driven modules:** TFs, Metabolites, Proteins, Lipids,
  **+ Kinase activity** (decoupleR + KSN).
- **Activity inference:** keep TFA (from expression/primary matrix) AND add **kinase
  activity via decoupleR** with an OmniPath kinase-substrate network (KSN).
- **Two clustering variants** (like the earlier extension):
  1. `all_sites` — all phosphosites (no feature selection)
  2. `top5000_hvg` — top-5000 most-variable phosphosites
- First analysis + **biological interpretation** of results when done.
- **Do NOT change nextflow pipeline code.** If the pipeline can't handle these data types,
  **list the issues** instead of patching it.

Everything separate: scripts in `Wang_GBM/Lemonite_phospho_clustering/`; results in a new
results dir (mirroring the earlier `…_phospho/` layout).

## Inputs (all verified present, GBM data dir)
- Phospho: `mmc3.xlsx` sheet `phosphoproteome_normalized` (70,330 sites × 109 samples,
  already log2/median-polished/ComBat; ids `GENE:RefSeq:Residue`).
- Proteome: `Proteome_normalized.csv` (10,999 × ~109, tab-sep, `symbol,refseq_prot_id,hgnc_id,<samples>`).
- Metabolome: `metabolome.csv` (134 × samples). Lipids: `lipidome_pos.csv`/`neg` (334/248).
- TFs: `lovering_TF_list.txt`. KSN: OmniPath via decoupleR in the `.sif`.

## Key pipeline facts (verified — reuse, don't re-derive)
- `preprocessing_type=proteomics` (`preprocessing_prescaled_proteomics.R`) accepts a
  **pre-scaled primary matrix** (no DESeq2) and writes `LemonPreprocessed_expression.txt`
  (the clustering input) + `LemonPreprocessed_complete.txt`. → route phospho through this.
- TFA: `Preprocessing_TFA_Proteomics.R` uses `decoupleR::decouple(mat, net, .source, .target)`
  + `run_consensus()`, net=CollecTRI. **Kinase activity = identical call with net=KSN**
  (`source=kinase, target=GENE_residue`). Mirror lines 208–277.
- Clustering = LemonTree Gibbs on `LemonPreprocessed_expression.txt`; regulator assignment =
  `-task regulators` per layer against `tight_clusters.txt` (see earlier extension's
  `run_phospho_regulators.sh`).
- Phosphosite → KSN target key: `GENE:RefSeq:Res[;Res2]` → `GENE_Res` (expand multi-residue).

## Steps
1. ☐ Build phospho clustering matrix (2 variants) in LemonTree format
   (`LemonPreprocessed_expression.txt`), NA→0, z-score per site. Reuse
   `../Lemonite_phospho/prepare_phospho_regulator.py` logic (already does exactly this;
   `--top-var 5000` for variant 2, none for variant 1).
2. ☐ Build regulator abundance files aligned to the phospho sample set:
   Proteins (Proteome_normalized), Metabolites, Lipids, TFs (lovering list) — z-scored,
   rows appended to `LemonPreprocessed_complete.txt` (LemonTree reg lookup).
3. ☐ Kinase activity: `kinase_activity.R` — decoupleR decouple/consensus on the phospho
   matrix (rows re-keyed `GENE_residue`) with OmniPath KSN → per-sample kinase-activity
   matrix; also run TFA on the primary matrix. Add kinase activity as a regulator layer
   (abundance = activity score per sample) + standalone heatmap/table.
4. ☐ Clustering + regulator assignment per variant (LemonTree in the `.sif`). If nextflow
   `main.nf` can't ingest phospho-as-primary, run the equivalent LemonTree steps directly
   (as the earlier extension did) and **record the pipeline issues** in `PIPELINE_ISSUES.md`.
5. ☐ First analysis + biological interpretation (module sizes, which kinases/TFs/mets
   regulate which phospho-modules, GBM relevance).

## Open decisions — RESOLVED (user 2026-07-06)
- Cluster input: **both** all-sites and top-5000 HVG (2 variants).
- Regulators: **all** (TFs, Metabolites, Proteins, Lipids) **+ kinase activity via decoupleR**.
- Kinase activity: **decoupleR on phospho + KSN**.

## Progress log
- 2026-07-06 — Plan created; inputs + pipeline mechanism verified (proteomics-prescaled path,
  TFA/decoupleR pattern, KSN site format `GENE_residue`).
- 2026-07-06 — Built both clustering matrices; kinase activity (run_ulm + KSN) + TFA working
  (GBM-relevant kinases: CDK1/2/4/5, CSNK2, EGF/AKT, GSK3B, MAPK14…). Built all 5 regulator
  layers + complete matrix. Documented pipeline issues (main.nf can't ingest phospho-primary).
  **Two bugs found+fixed:** (1) phospho sheet mixes CPTAC discovery (C3L/C3N, 99) with 10
  pediatric CBTTC PT-* samples that lack multi-omics — now filtered to 99 CPTAC via
  `--sample-prefix C3L,C3N`; (2) LemonTree `regulators` crashed with
  ArrayIndexOutOfBounds because the complete matrix had DUPLICATE row names (a gene appears as
  Protein + TF + Kinase) — fixed by deduping rows and suffixing kinase-activity features
  `_kact`. (3) **LemonTree splits data_file rows on whitespace**, so metabolite/lipid feature
  names containing spaces (e.g. `2-hydroxybutyric acid`) shifted columns → the SAME
  ArrayIndexOutOfBounds at "Reading expression data". Fixed by sanitizing all regulator
  feature names (spaces/tabs → `_`) in build_regulators.py. Confirmed `-task regulators` now
  parses. Next: run all 5 regulator layers on existing tight_clusters; then all_sites; then
  first analysis/interpretation.

## ▶ RESUME TOMORROW (state as of 2026-07-06 evening)
**Everything expensive is on disk and safe to stop.** To resume: `bash resume.sh` (skips
finished steps, re-launches only what's missing).

- **top5000_hvg: COMPLETE** — clustering (`Lemon_out/tight_clusters.txt`, 178 modules) + all
  5 regulator layers (`{TFs,Metabolites,Proteins,Lipids,KinaseActivity}.topreg.txt`).
  (Lipids/KinaseActivity were run with a `_par` suffix during parallelization — resume.sh /
  analysis should treat `Lipids_par`≈`Lipids`, `KinaseActivity_par`≈`KinaseActivity`, or just
  rename them.)
- **all_sites: prep COMPLETE, clustering IN PROGRESS** — `Preprocessing/` (46,060-row complete
  matrix) + `Activity/` (kinase+TF) are saved. Clustering on 34,090 sites is the long pole
  (~1–3 h) and is NOT resumable from partial (LemonTree ganesh restarts fresh) — but its
  inputs are saved, so `resume.sh` just relaunches clustering + regulators. No prep is lost.

Nothing to rebuild by hand tomorrow — all matrices/activity are persisted; only LemonTree
compute needs to finish.

### STOP STATE 2026-07-06 ~22:00 (user stopped for the day)
- **top5000_hvg:** clustering done (179 modules, 12,705 sites). Regulator layers on disk:
  **TFs, Metabolites, KinaseActivity DONE**; **Proteins + Lipids were on the final module
  (178/178)** when stopped — likely finished or ~1 min short. Tomorrow: check for
  `Proteins.topreg.txt` / `Lipids_par.topreg.txt`; if missing, `resume.sh` re-runs just those.
  NOTE: Lipids output is named `Lipids_par.*` (parallel-launch suffix) — rename to `Lipids.*`.
- **all_sites:** prep + activity SAVED; clustering 0/10 (barely started — nothing lost).
  Tomorrow: `bash resume.sh` launches all_sites clustering (~1–3 h) then its 5 regulator
  layers, and finishes any missing top5000 layers.
- Sanity already good: TF regs (PURG/MEOX2/REL/KLF16), kinase-activity regs
  (LYN/UHMK1/CAMK2B/BUB1/CSNK1D) assigned to phospho modules — GBM-plausible.
- **Still TODO next time:** finish all_sites; canonicalize `_par` names; write
  `FIRST_ANALYSIS.md` (biological interpretation of both variants).
- Running java procs may still be finishing; safe to let them complete or kill — all inputs
  are persisted, nothing to lose.
