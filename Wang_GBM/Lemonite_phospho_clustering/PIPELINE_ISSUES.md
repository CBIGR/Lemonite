# nextflow pipeline issues for phospho-as-clustering-input

Earlier version of this doc concluded `main.nf` couldn't ingest this variant and ran LemonTree
directly instead (see git history / `run_lemontree.sh`). That conclusion was **too pessimistic
on two of the four points** and has since been corrected: the analysis now runs through the
real `nextflow run main.nf`, with one small, authorized adaptation (below). Kept the earlier
issues here since #3 and #4 are still real and still not attempted.

## What actually blocked `main.nf`, and what was fixed

1. **`Preprocessing_TFA_Proteomics.R` had no generic `--regulator_types` loader.**
   `Preprocessing_TFA_RNA.R` has a fully generic `parse_regulator_config()` reading arbitrary
   `Prefix:File[:c|d]` pairs (`nextflow/scripts/Preprocessing_TFA_RNA.R:48-70`) — exactly what
   the CLAUDE.md "config-driven, no code changes" description promises. The proteomics script
   only had two fixed hooks (`--hptm_file`, `--metabolomics_file`), so 4 of our 5 regulator
   layers (TFs-activity, Proteins, Lipids, KinaseActivity) had no way in. **Fixed** (user
   authorized adapting this specific script): ported the same generic loader into
   `Preprocessing_TFA_Proteomics.R`, modeled on the RNA script's pattern but z-scoring like the
   existing hPTM path (no log/pareto — this script's convention is "pre-scaled data in,
   minimal reprocessing"), plus de-duplication of feature names both within a layer and across
   layers (LemonTree corrupts its array indexing on duplicate row names — hit this immediately
   with `Proteome_normalized.csv`, which has repeat gene symbols across protein isoforms).
   `modules/preprocessing.nf`'s proteomics branch was updated to pass `--regulator_types`
   through (previously only derived `--hptm_file`/`--metabolomics_file` from it).

2. **Original "re-applies its own scaling" claim was wrong.** `Preprocessing_TFA_Proteomics.R`
   does *not* rescale the primary matrix — it only does `--top_n_genes` variable-feature
   selection on whatever is handed in. Feeding it an already z-scored, full (untruncated)
   phospho matrix and letting `--top_n_genes` do the selection (5000 for the HVG variant, a
   number ≥ site count for `all_sites`) works cleanly and matches how the pipeline treats any
   other proteomics-type primary matrix.

3. **`--expression_file`/`--metadata_file` already support explicit overrides** — the "no way
   to hand it a pre-built matrix" claim was also wrong; auto-detection is only the fallback
   path when these aren't set.

4. **Container couldn't see edited scripts without `-profile dev`.** Unrelated latent bug:
   the proteomics branch's script-selection guard (`modules/preprocessing.nf`) only checks
   `${projectDir}/scripts/Preprocessing_TFA_Proteomics.R`, with no `/app/scripts/` fallback
   (the RNA branch has both). Since `projectDir` isn't bound into the singularity container by
   default, that check silently failed and execution fell through to the `else` branch, which
   *does* have an `/app/scripts/` fallback — i.e. it silently ran `Preprocessing_TFA_RNA.R`
   instead, with no error. `-profile singularity,dev` (bind-mounts `scripts/` and `PKN/` into
   the container, meant for script iteration) fixes this by making `projectDir` visible. Not
   patched further since `-profile dev` is a documented, existing escape hatch — but worth
   flagging: **any** `preprocessing_type=proteomics` run without `dev` silently mis-runs as RNA
   preprocessing instead of erroring. That asymmetry is a real pipeline bug, left unfixed here
   (out of scope of "adapt Preprocessing_TFA_Proteomics.R").

## Bug found after `top5000_hvg` completed: metabolomics scaled on the wrong axis

`pareto_scale()` in `Preprocessing_TFA_Proteomics.R` (used for the `Metabolites` regulator
layer via `--metabolomics_file`) is a row-wise (per-feature) function — same design as the
reference implementation in `Preprocessing_TFA_RNA.R:1041-1052`, which is called directly on
a `features × samples` matrix. But the metabolomics call site wrapped it in the
`t(fn(t(x)))` idiom that's only correct for column-wise functions like R's built-in `scale()`
(copied from the neighboring hPTM line, `t(scale(t(hPTM)))`). Feeding `t(metabolomics)`
(`samples × metabolites`) into a row-wise function computes each **sample's** mean/sd across
all metabolites and centers/scales per sample instead of per metabolite — flattening real
cross-sample variation for every metabolite (confirmed with a toy R example: buggy output was
~9.5/9.6/9.5/9.5 across 4 samples for one feature; correct output was -0.92/1.31/0.20/-0.60).
This is what caused the "values look highly similar across samples" symptom in the
metabolomics heatmaps. **Fixed**: call `pareto_scale(metabolomics)` directly, no transpose
wrapping (`Preprocessing_TFA_Proteomics.R:365-370`).

Only `Metabolites` was affected (uses the native `--metabolomics_file` hook → `pareto_scale`).
`Proteins`/`Lipids`/`TFs`/`KinaseActivity` go through the generic loader's `t(scale(t(layer)))`
path, which is the correct column-wise idiom, and hPTM uses the same correct pattern — neither
was affected. The already-completed `top5000_hvg` run's Metabolites layer used the buggy
scaling; it needs to be re-preprocessed (and CLUSTERING re-run, since the combined matrix baked
the wrong values in) to get correct results. `preprocessing_prescaled_proteomics.R` has the
identical bug but is dead code — not referenced by any `modules/*.nf` file — left as-is.

## Still real, still not attempted

5. **No kinase-activity step in the pipeline.** TF activity and kinase activity for this
   variant both come from `kinase_and_tf_activity.R` (decoupleR `run_ulm`, CollecTRI + OmniPath
   KSN), run standalone against each variant's real (post `--top_n_genes` selection) primary
   matrix, then fed back in as two more `--regulator_types` layers. No `main.nf` hook for this
   kind of inference exists or was added.

6. **Regulator layers assume gene-symbol primary rows.** Downstream network/enrichment steps
   (`lemontree_to_network.py`, enrichment, MegaGO) assume module members are genes; here they're
   phosphosites. Clustering + regulator assignment are agnostic to row identity and work fine,
   but network/enrichment/overview would need phosphosite→gene mapping to be meaningful — out
   of scope for this first analysis, so runs use `--stop_after_network`.

## What actually runs now

Real `nextflow run main.nf -profile singularity,dev,local`, `preprocessing_type=proteomics`,
`--expression_file`/`--metadata_file` pointing at a full z-scored phospho matrix + metadata,
`--regulator_types` listing all 5 layers, `--perform_tfa false` (CollecTRI targets are gene
symbols — meaningless matched against phosphosite ids), `--stop_after_network`, pipeline
defaults otherwise (**n_clusters=100, node_clustering=true** — the two defaults an earlier ad
hoc run had wrong, at 10 and `false`). External `-c nf_6cores.config` caps the local executor
at 6 concurrent cores; no pipeline file is touched by that cap. Per-variant TFs/KinaseActivity
raw files are precomputed once (`kinase_and_tf_activity.R` on each variant's real
`--top_n_genes`-selected primary matrix) since main.nf regenerates the primary matrix inside a
single preprocessing run and there's no way to feed it stage-by-stage in one pass.
