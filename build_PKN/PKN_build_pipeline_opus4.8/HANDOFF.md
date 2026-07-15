# Handoff: running the PKN pipeline on a new machine / HPC

This note carries the context needed to continue the PKN build pipeline in a fresh
session (for example on an HPC cluster, where the original chat history is not
available). It records what was built, what was fixed and why, the current
verification status, and exactly how to launch the full run.

Last updated: 2026-06-23.


## Current status

- The pipeline lives in `build_PKN/PKN_build_pipeline_opus4.8/` and is committed to
  git on `main`. It is a faithful, modular re-implementation of the three reference
  notebooks in `build_PKN/` (`Collect_PKNdata_metabolites.ipynb`,
  `Collect_PKNdata_proteins.ipynb`, `Build_final_PKN.ipynb`).
- A 500-metabolite test build was validated against the original `PKN/` output. The
  network is correct: 7 of 9 metabolite-gene sources matched the reference exactly;
  the 2 that differed (UniProtKB, ChEMBL) were real upstream issues that have since
  been fixed and re-verified (details below).
- The full all-metabolites run has NOT been done yet. That is the next step, intended
  to run on HPC.


## What was fixed this session (do not regress these)

1. STRING endpoint: default switched from the flaky `version-11-5` mirror to the
   current stable `https://string-db.org/api` (config.py, `STRING_API_URL`). Same
   request format; override with `PKN_STRING_API_URL` if needed.

2. Provenance URLs now point to evidence for each specific interaction, not landing
   pages. Notably: ChEMBL uses `explore/compound|target` (not the dead report-card
   URLs); LINCS falls back to the ChEMBL compound page (its data is ChEMBL-derived);
   MetalinksDB points to the metabolite's HMDB page; Human-GEM points to the specific
   Metabolic Atlas reaction page, distinct for distance 1 and distance 2. A build-time
   audit (`utils/url_audit.py`) checks these and writes `url_audit_report.csv`.

3. UniProtKB retriever rewritten (`step1_metabolites/uniprot.py`). The UniProt
   `inchikey:` search field is dead (returns nothing, even for ATP), which made the
   old query produce zero edges. It now queries the `chebi:` field, expanded to the
   metabolite's conjugate-acid/base/tautomer ChEBI species (UniProt annotates charged
   forms, e.g. biotin(1-) CHEBI:57586 not neutral CHEBI:15956). Expansion is resolved
   from the ChEBI ontology via OLS and cached. A per-metabolite cap
   (`MAX_PROTEINS_PER_METABOLITE = 200`) excludes promiscuous cofactors. See
   DOCUMENTATION.md section 3.2.1. Verified: recovers 92-100% of the reference genes
   on a 10-metabolite spot check.

4. ChEMBL activity filter broadened (`step1_metabolites/chembl.py`). The old exact
   `activity_comment in {active, substrate}` rule dropped real interactions. It now
   keeps activities whose comment contains an interaction keyword (active, substrate,
   inhibitor, agonist, antagonist, binder, inducer) and is not an explicit negative,
   OR that have a measured `pchembl_value >= 5.0` (about 10 micromolar). See
   DOCUMENTATION.md section 3.2.2.

5. Preprocessing hardened (`step1_metabolites/preprocessing.py`): it no longer
   re-merges the ChEMBL mapping when `ChEMBL_id` is already present, and strips any
   `_x`/`_y` merge-suffix columns before persisting the annotated table (a prior bug
   had polluted `HMDB_metabolites_annotated.csv` with `ChEMBL_id_x`/`_y`).


## Environment notes (important for HPC)

- Python deps are in `requirements.txt`. Locally the working interpreter was a conda
  env named `LemonIte` (pandas, rdkit, mygene, chembl-webresource-client, networkx).
  On HPC, either `bash setup.sh` to build a venv, or use the pre-existing
  `venv_pkn_doduo` referenced by `submit_pkn_hpc.sh`.
- IMPORTANT: do NOT launch via `nohup conda run ...` or `nohup` with a bare `conda`
  command. `conda` is a shell function and is not on PATH in a detached subshell, so
  it fails immediately. Activate the environment first (`source .../bin/activate`, as
  `schedule/run_monthly.sh` and `submit_pkn_hpc.sh` do) or call the interpreter by
  absolute path.
- Configuration is via environment variables (no file edits needed):
  `PKN_WORKDIR`, `PKN_OUTPUT_DIR_NAME`, `PKN_DB_DIR`, `PKN_GEM_DIR`, `PKN_CONFIG`,
  `PKN_STRING_API_URL`. `submit_pkn_hpc.sh` already sets these for VSC Ghent; adjust
  the paths in it to the actual HPC locations of the database files and Human-GEM
  model.


## Required input data on HPC

The local-file sources must be present under `PKN_DB_DIR` (and Human-GEM under
`PKN_GEM_DIR`). On the local machine these lived on an external drive at
`/run/media/borisvdm/5250-D000/resources/databases`. The expected files are listed in
DOCUMENTATION.md section 2 and include: HMDB `hmdb_metabolites.xml`, BioGRID chemicals
and ALL tab3, STITCH `9606.protein_chemical.links.v5.0.tsv`, `metalinks.csv`, the four
LINCS `lsp_*.csv` files, `HuRI.tsv`, `HumanNet/HS-PI.symbol.tsv`, an Ensembl symbol
mapping, and the Human-GEM model files. Confirm these exist on HPC before launching.


## How to run the full build on HPC

1. Ensure the database files are present and the env variables in `submit_pkn_hpc.sh`
   point at them.
2. Activate the environment.
3. Submit:  `qsub submit_pkn_hpc.sh`   (runs `python main.py --all --resume`).

The run is checkpointed and resumable: if the job hits the wall-time, re-submit and it
continues from the last completed database without re-querying. A 500-metabolite build
took about 13 minutes locally; the full HMDB set (~200k metabolites) is dominated by
the rate-limited API stages (UniProtKB, IntAct, ChEMBL) and the GEM graph traversal,
so expect a long, resumable run.

For a quick HPC smoke test first:
`python main.py --step 1 --max-metabolites 500 --output-dir PKN_test_500`


## How to re-verify correctness on HPC (optional)

If a reference `PKN/` directory is available on HPC, the comparison method used locally
was: restrict both `metabolite_gene_PKN.tsv` files to the test metabolites (by HMDB id
parsed from the `Metabolite` column), then compare `(HMDB, Gene, Source)` edge sets per
source. Expect the local-file sources to match closely; UniProtKB and ChEMBL will
differ from an older reference because of the upstream fixes above and ongoing database
updates (the pipeline is now more current and, for ChEMBL, more selective about weak
activities).


## Open / possible follow-ups

- Full all-metabolites run (the immediate next step).
- The UniProtKB 200-protein cap truncates alphabetically for ~5 percent of metabolites
  (highly-connected cofactors). If a relevance-ranked cap is wanted, that is a future
  refinement.
- IntAct partner labels are sometimes complexes or protein-chain identifiers rather
  than single gene symbols; by decision these are retained as-is.
