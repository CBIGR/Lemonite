# LemonIte Prior Knowledge Network (PKN) Build Pipeline

A reproducible pipeline that constructs the LemonIte Prior Knowledge Network: a
heterogeneous interaction network combining metabolite-gene, protein-protein, and
enzyme-histone-modification interactions assembled from curated databases and public
web APIs. The network serves as a causal scaffold for mechanistic multi-omics
modelling in the LemonIte project.

The assembled network can be queried and explored interactively at
**https://www.lemonite.ugent.be/Knowledge_Graph_Exploration**.

For data sources, output schema, provenance-URL design and methodological detail,
see [DOCUMENTATION.md](DOCUMENTATION.md).


## Overview

The pipeline builds a network with three edge types:

| Edge type              | Meaning                                                         |
| ---------------------- | -------------------------------------------------------------- |
| `metabolite-gene`      | A metabolite regulates, binds, or is metabolised by a gene product |
| `PPI`                  | A protein-protein interaction between two gene products        |
| `histone-modification` | A histone-modifying enzyme or reader acts on or binds a specific histone mark |

Every edge carries a `Source` label (the database it was derived from) and a
provenance URL that links back to the record supporting the interaction.

It runs in the following steps:

1. **Metabolite-gene interactions** from nine sources: HMDB (annotation), BioGRID,
   STITCH, UniProtKB, IntAct, ChEMBL, LINCS, Human-GEM (reaction distances 1 and 2),
   and MetalinksDB.
2. **Protein-protein interactions** from STRING, BioGRID, HuRI, and HumanNet, seeded
   by the genes found in step 1.
2b. **Enzyme-histone-modification interactions** from Gene Ontology molecular-function
   annotations, queried from QuickGO and UniProtKB — which protein writes, erases, or
   reads a specific histone mark (for example, EZH2 writes H3K27me3). This step is
   independent of steps 1 and 2.
2c. **Phospho-regulator network (standalone)** — kinase-substrate edges from the
   OmniPath enzyme-substrate service (which aggregates PhosphoSitePlus, SIGNOR, DEPOD,
   dbPTM, phospho.ELM, HPRD, NetworKIN, … ~18 resources) plus phospho-TF-target edges from
   CollecTRI. Written to `phospho_pkn.tsv`. **This layer is deliberately kept separate: step
   3 does NOT merge it into `LemonIte_PKN.tsv`.** It is consumed by the phospho-regulator
   analysis (`Wang_GBM/Lemonite_phospho/`). Independent of steps 1, 2, 2b.
3. **Integration**: the networks are merged and de-duplicated, provenance URLs are
   attached, and a URL audit is performed.

The pipeline is designed for unattended monthly re-runs to keep the network current,
on a workstation or on an HPC cluster.


## Installation

```bash
cd PKN_build_pipeline_opus4.8
bash setup.sh                 # create a virtual environment and install requirements
source venv/bin/activate
```

Requirements: Python 3.9 or newer. Dependencies are listed in `requirements.txt`
(pandas, numpy, requests, networkx, rdkit, mygene, chembl-webresource-client).


## Usage

```bash
# Full build
python main.py --all

# Individual steps
python main.py --step 1       # metabolite-gene interactions
python main.py --step 2       # protein-protein interactions
python main.py --step 2b      # enzyme-histone-modification interactions
python main.py --step 2c      # phospho-regulator network (standalone phospho_pkn.tsv)
python main.py --step 3       # integration, annotation, URL audit

# A fast end-to-end test on a metabolite subset, in a separate output directory
python main.py --step 1 --max-metabolites 500 --output-dir PKN_test_500
```

Selected options:

| Option                   | Effect                                                       |
| ------------------------ | ----------------------------------------------------------- |
| `--databases a,b,c`      | Run only the named step-1 databases                         |
| `--resume` / `--no-resume` | Resume from checkpoints (default on) or force fresh queries |
| `--max-metabolites N`    | Limit step 1 to the first N metabolites (test runs)         |
| `--output-dir DIR`       | Write outputs to an alternate directory                     |
| `--workers N`            | Override the per-database thread count                       |
| `--audit-urls`           | Run only the provenance-URL audit                           |


## Configuration

All paths, thresholds and API settings are defined in `config.py`. The settings most
likely to change between environments are read from environment variables, so HPC and
container runs do not require editing the file:

| Variable              | Default                       | Description                          |
| --------------------- | ----------------------------- | ------------------------------------ |
| `PKN_WORKDIR`         | project working directory     | Root working directory               |
| `PKN_OUTPUT_DIR_NAME` | `PKN`                         | Output sub-directory name            |
| `PKN_DB_DIR`          | local database root           | Root directory of local database files |
| `PKN_GEM_DIR`         | Human-GEM model directory     | Human-GEM model files                |
| `PKN_CONFIG`          | `config`                      | Config module (`config_hpc` on HPC)  |
| `PKN_STRING_API_URL`  | `https://string-db.org/api`   | STRING REST endpoint                 |


## Outputs

Written to `$PKN_WORKDIR/$PKN_OUTPUT_DIR_NAME/`:

| File                          | Columns                          | Description                                  |
| ----------------------------- | -------------------------------- | -------------------------------------------- |
| `LemonIte_PKN.tsv`            | `Node1, Node2, Source, Type`     | The combined, de-duplicated network          |
| `LemonIte_PKN_with_URLs.tsv`  | `HMDB, Gene, Source, URL`        | Metabolite-gene edges with provenance URLs   |
| `metabolite_gene_PKN.tsv`     | `Metabolite, Gene, Source, url`  | Step-1 long-format metabolite-gene network   |
| `PPI_network.tsv`             | `GeneA, GeneB, combined_score, Source` | Step-2 protein-protein network        |
| `hPTM_network.tsv`            | `Enzyme, Mark, Activity, Source, GO_ID, URL` | Step-2b enzyme-histone-mark network with provenance |
| `phospho_pkn.tsv`             | `Node1, Node2, Source, Type, Residue, URL` | Step-2c standalone phospho network (kinase-substrate + phospho-TF-target); NOT merged into LemonIte_PKN.tsv |
| `OmniPath_kinase_substrate.tsv` | raw OmniPath enzsub TSV | Step-2c resumable cache of the OmniPath kinase-substrate pull |
| `url_audit_report.csv`        | audit results                    | Per-source reachability of sampled URLs      |
| `PKN_run_history.jsonl`       | run ledger                       | One record per run (see Run ledger below)    |
| `changes/<run_id>_diff.tsv.gz`| edge diff                        | Edges added/removed vs the previous run      |
| `figures/*.png`               | summary figures                  | Per-source counts, database-overlap UpSet plots, composition (see Figures below) |

Per-database intermediate caches (`HMDB_metabolites_<DB>_processed.csv`) and
per-database link files (`HMDB_metabolites_annotated_<DB>_links.csv`) are also
written, and are reused by `--resume`.

## Figures

At the end of step 3 the pipeline writes summary figures to `figures/`, re-deriving
the reproducible plots from the notebooks:

| Figure                                          | Shows                                                     |
| ----------------------------------------------- | --------------------------------------------------------- |
| `metabolite_gene_interactions_per_source.png`   | Unique metabolite-gene interactions per source database   |
| `metabolite_gene_overlap_upset.png`             | UpSet plot: overlap of metabolite-gene edges across sources |
| `PPI_interactions_per_source.png`               | Unique PPI interactions per source database                |
| `PPI_overlap_upset.png`                         | UpSet plot: overlap of PPI edges across sources            |
| `hPTM_interactions_per_source.png`              | Enzyme-histone-mark interactions per source (step 2b)      |
| `hPTM_overlap_upset.png`                        | UpSet plot: enzyme-mark overlap across sources (step 2b)   |
| `hPTM_interactions_by_activity.png`             | Enzyme-histone-mark interactions by writer/eraser/reader   |
| `PKN_composition.png`                           | Node and edge totals of the final PKN, per type            |

The overlap (UpSet) plots read the long-format per-source network files (where an
edge may appear once per supporting source), not the de-duplicated `LemonIte_PKN.tsv`.
Figure generation is best-effort: a plotting failure is logged and skipped, and never
fails a completed build. The notebook's MetalinksDB/MEBOCOST comparison and
HMDB-superclass figures are omitted here because they require external database files
not needed to build or validate the network.

## Run ledger and change tracking

Every build records itself automatically at the end of step 3: a run entry (id,
timestamp, git commit, parameters, per-source edge counts, checksum) is appended to
`PKN_run_history.jsonl`, the network is diffed against the previous run (full
added/removed edge list written to `changes/`), and a human-readable summary is
appended to the repo-versioned `PKN_CHANGELOG.md`. This gives a permanent, auditable
history of every run and of how the PKN changes over time. See DOCUMENTATION.md
section 9 for detail.


## Robustness and rate limiting

The pipeline queries several public APIs that impose rate limits and occasionally
time out. Reliability is handled by three deterministic layers (no external services
or language models are involved):

1. **Connection level** (`utils/http.py`): a shared HTTP session with automatic
   retry on transient connection failures and on the `429`, `500`, `502`, `503` and
   `504` status codes, honouring the `Retry-After` header.
2. **Application level** (`utils/api_retry.py`): a per-database retry policy with
   exponential backoff and jitter, explicit handling of HTTP 429, and logging to
   `api_errors.log`. After exhausting retries a request returns a sentinel value
   rather than aborting the run.
3. **Persistence level** (`utils/cache.py`): each database is processed in
   checkpointed chunks written to disk. A run interrupted by a failure, a timeout, or
   an exhausted cluster wall-time resumes from the last completed chunk with
   `--resume`, without re-querying completed databases.


## Provenance URLs

Each edge links to a record that provides evidence for that specific interaction,
rather than to a database landing page. URL templates are centralised in
`config.URL_TEMPLATES` and were validated against the live databases. At the end of
step 3, `utils/url_audit.py` samples URLs per source, checks their reachability, and
writes `url_audit_report.csv`; a future run that hits a changed URL scheme is flagged
automatically. See [DOCUMENTATION.md](DOCUMENTATION.md) for the per-source URL design.


## Scheduled (monthly) re-runs

The build is intended to be re-run periodically. A wrapper script and an example
cron configuration are provided in `schedule/`:

```bash
# 1. Set the PKN_* paths in schedule/run_monthly.sh (or export them in the environment)
# 2. Install a cron entry (see schedule/crontab.example), for example on the 1st of
#    each month at 03:00:
#       0 3 1 * *  /path/to/schedule/run_monthly.sh >> /path/to/pkn_cron.log 2>&1
```

`run_monthly.sh` runs `python main.py --all --resume` and is idempotent: a re-run
after a successful build is a fast no-op, and a re-run after an interrupted build
completes it.


## HPC deployment

```bash
qsub submit_pkn_hpc.sh
```

The PBS job sets the `PKN_*` environment variables and `PKN_CONFIG=config_hpc`, then
runs the full build. If the job reaches its wall-time, re-submitting it resumes from
the last checkpoint.


## Repository layout

```
PKN_build_pipeline_opus4.8/
  main.py                 Command-line entry point
  config.py               Paths, thresholds, API settings, URL templates
  config_hpc.py           HPC overrides (PKN_CONFIG=config_hpc)
  requirements.txt
  setup.sh
  README.md               This file
  DOCUMENTATION.md        Data sources, schema, provenance and methods
  utils/
    http.py               Shared HTTP session with connection-level retry
    api_retry.py          Per-database retry decorator
    cache.py              Checkpointed, resumable processing
    file_io.py            HMDB XML parser, SMILES canonicalisation
    url_audit.py          Provenance-URL reachability audit
  step1_metabolites/      Preprocessing, nine source modules, integration
  step2_proteins/         STRING, BioGRID, HuRI, HumanNet, integration
  step2b_hPTM/            QuickGO + UniProtKB enzyme-histone-mark retrieval, integration
  step2c_phospho/         OmniPath enzsub (kinase-substrate) + CollecTRI (phospho-TF-target) -> standalone phospho_pkn.tsv
  step3_final/            Combiner, annotator, analysis, visualization, provenance
  schedule/               run_monthly.sh, crontab.example
  submit_pkn_hpc.sh       PBS/qsub job script
```


## Citation and access

The LemonIte Prior Knowledge Network can be explored interactively at
https://www.lemonite.ugent.be/Knowledge_Graph_Exploration. When using the network or
this pipeline, please cite the LemonIte project and the underlying source databases
listed in [DOCUMENTATION.md](DOCUMENTATION.md).
