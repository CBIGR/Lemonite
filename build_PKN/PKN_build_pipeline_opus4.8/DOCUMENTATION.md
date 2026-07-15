# LemonIte Prior Knowledge Network: Data Sources, Schema, and Methods

This document describes the data sources, construction methods, output schema, and
provenance-URL design of the LemonIte Prior Knowledge Network (PKN). For installation
and usage, see [README.md](README.md). The assembled network can be explored at
https://www.lemonite.ugent.be/Knowledge_Graph_Exploration.


## 1. Network definition

The PKN is a heterogeneous network with three node types (metabolites, identified by
HMDB accession; genes, identified by HGNC symbol; and histone marks, identified by a
controlled vocabulary such as `H3K27me3`) and three edge types:

- `metabolite-gene`: a metabolite regulates, binds, transports, or is metabolised by
  a gene product.
- `PPI`: a protein-protein interaction between two gene products.
- `histone-modification`: a histone-modifying enzyme (writer or eraser) or reader
  protein acts on or binds a specific histone mark (Node1 = enzyme gene, Node2 =
  histone mark).

Each edge is annotated with the source database from which it was derived and a
provenance URL pointing to the supporting record.


## 2. Data sources

### 2.1 Metabolite-gene interactions (step 1)

| Source             | Access            | Metabolite key | Description                                                                 |
| ------------------ | ----------------- | -------------- | --------------------------------------------------------------------------- |
| HMDB               | Local XML         | -              | Metabolite annotation (names, identifiers, structures, taxonomy)            |
| BioGRID (chemicals)| Local file        | InChIKey       | Curated chemical-protein interactions                                       |
| STITCH             | Local file + API  | PubChem CID    | Chemical-protein associations; Ensembl protein IDs mapped to symbols via MyGene.info |
| UniProtKB          | REST API          | ChEBI ID       | Proteins annotated with the metabolite as a cofactor, ligand, or substrate (see 3.2.1) |
| IntAct             | REST API          | ChEBI ID       | Curated molecular interactions involving the metabolite                     |
| ChEMBL             | REST API          | ChEMBL ID      | Measured bioactivities of the metabolite against human targets (see 3.2.2)   |
| LINCS              | Local files       | ChEMBL ID      | Biochemical binding affinities (filtered by an IC50 threshold)              |
| Human-GEM (dist 1) | Local model       | ChEBI ID       | Genes catalysing reactions involving the metabolite (one reaction hop)      |
| Human-GEM (dist 2) | Local model       | ChEBI ID       | Genes within two reaction hops of the metabolite                            |
| MetalinksDB        | Local file        | HMDB ID        | Metabolite-receptor interactions aggregated from CellPhoneDB, NicheNet, scConnect, and others |

### 2.2 Protein-protein interactions (step 2)

| Source   | Access     | Description                                                       |
| -------- | ---------- | ---------------------------------------------------------------- |
| STRING   | REST API   | Functional associations (combined score), human, queried per gene set |
| BioGRID  | Local file | Curated human physical interactions                              |
| HuRI     | Local file | Human Reference Interactome; Ensembl IDs mapped to HGNC symbols  |
| HumanNet | Local file | Functional gene network (gene symbols)                           |

The PPI network is seeded by the set of genes appearing in the step-1 metabolite-gene
network, retaining interactions in which at least one partner is a seed gene.

### 2.3 Enzyme-histone-modification interactions (step 2b)

| Source       | Access   | Query key            | Description                                                              |
| ------------ | -------- | -------------------- | ------------------------------------------------------------------------ |
| QuickGO (EBI)| REST API | GO molecular function| Gene products annotated with a residue-specific histone writer/eraser/reader GO term |
| UniProtKB    | REST API | GO molecular function| Reviewed human proteins carrying the same GO terms (independent cross-source) |

These edges capture enzyme-substrate specificity — which protein writes, erases, or
reads a given histone mark (for example, EZH2 writes H3K27me3, KDM6A erases H3K27me3,
BRD4 reads H3K27ac). They are derived from Gene Ontology molecular-function
annotations, not from ChIP-seq genomic-proximity data. The two sources are queried
independently and both are retained, so an edge supported by both carries two source
labels and two provenance URLs (see 3.5). This step is independent of steps 1 and 2:
it queries the GO term set directly and does not require the metabolite-gene seed.


## 3. Construction method

### 3.1 Preprocessing

The HMDB metabolite XML is stream-parsed into an annotation table. SMILES strings are
canonicalised with RDKit, and canonical SMILES are mapped to ChEMBL molecule
identifiers through the ChEMBL web client. The ChEMBL mapping is cached and reused by
both the ChEMBL and LINCS modules.

### 3.2 Metabolite-gene retrieval

Each source is queried independently using the metabolite identifier appropriate to
that source (see the table in section 2.1). API-backed sources (UniProtKB, IntAct,
ChEMBL) are processed in checkpointed chunks with a per-database thread pool; local
and model-based sources are processed in a single pass. Human-GEM interactions are
derived by building a directed reaction graph (excluding ubiquitous currency
metabolites) and, for each metabolite, collecting the genes that catalyse reactions
within one or two reaction hops, recording the connecting reaction *path*: one reaction
at distance 1, and the two-reaction (metabolite -> intermediate -> gene) path at
distance 2.

#### 3.2.1 UniProtKB ChEBI query and charge-state expansion

UniProtKB associations are retrieved by querying the UniProtKB `chebi:` cross-reference
field for human entries whose binding-site, cofactor, or catalytic-activity annotations
reference the metabolite's ChEBI identifier. (The reference notebooks queried the
UniProtKB `inchikey:` field; that field no longer returns results from the UniProt REST
API and produced empty results for every metabolite, so the pipeline now queries by
ChEBI.)

UniProt annotates ligands under specific charge and tautomeric states rather than the
neutral species. For example, UniProt references biotin(1-) (CHEBI:57586) rather than
neutral biotin (CHEBI:15956), and the L-phenylalanine zwitterion (CHEBI:58095) rather
than CHEBI:17295. Querying the neutral HMDB ChEBI alone therefore misses many
associations, particularly for amino acids and central metabolites. To address this,
each metabolite's ChEBI identifier is expanded to include its conjugate-acid,
conjugate-base, protonated, deprotonated, and tautomeric species, resolved from the
ChEBI ontology through the EBI Ontology Lookup Service and cached per ChEBI. The
UniProtKB query then covers the metabolite's ChEBI and all of these related species.

To prevent a small number of highly promiscuous cofactors (ATP, NAD+, magnesium, and
similar), which are bound by thousands of proteins, from dominating the network, the
number of proteins accepted per metabolite is capped (default 200), in line with the
reference network's per-metabolite distribution.

#### 3.2.2 ChEMBL activity filter

ChEMBL associations are kept when an activity record indicates a genuine interaction.
The reference notebooks kept only activities whose `activity_comment` was exactly
`active` or `substrate`. In current ChEMBL this exact-match rule discards many real
interactions: substrate and inhibitor records whose comment carries an embedded
constant (for example `substrate [Km=330uM]` or `inhibitor [487uM]`), and the large
fraction of activities that have no comment but a measured potency. The filter is
therefore broadened to keep an activity when either (a) its `activity_comment` contains
an interaction keyword (active, substrate, inhibitor, agonist, antagonist, binder,
inducer) and does not contain an explicit negative (not active, inactive, not
determined, inconclusive), or (b) it has a measured `pchembl_value` of at least 5.0
(a potency of about 10 micromolar or better, consistent with the LINCS IC50 threshold).
Weak or explicitly negative activities (for example "Inhibition < 50% at 10 uM") are
excluded.

Each source writes two files: a processed table (one row per metabolite, with
pipe-separated gene symbols) and a link table (one row per interaction, with the
supporting provenance URL). These are merged into the long-format
`metabolite_gene_PKN.tsv`, with provenance URLs attached per edge.

### 3.3 Protein-protein retrieval

STRING is queried in chunks for the seed gene set and cached. BioGRID, HuRI, and
HumanNet are loaded from local files and pruned to interactions involving the seed
genes. The four sources are concatenated into `PPI_network.tsv`.

### 3.4 Enzyme-histone-modification retrieval

A curated set of residue-specific histone Gene Ontology molecular-function terms
(`config.HISTONE_GO_TERMS`) maps each GO term to a histone mark and an activity class
(writer, eraser, or reader). Examples: `GO:0046976` (histone H3K27 methyltransferase
activity) -> (`H3K27`, writer); `GO:0071558` (histone H3K27me2/me3 demethylase
activity) -> (`H3K27`, eraser); `GO:0140119` (histone H3K27ac reader activity) ->
(`H3K27ac`, reader). The GO identifiers were verified against QuickGO's ontology
search.

For each term, QuickGO's annotation-search endpoint is queried for human protein
annotations (`goUsage=descendants`, so enzymes annotated only to a more specific child
term are still captured under the queried parent's mark and activity), and UniProtKB's
REST API is queried for reviewed human proteins carrying the same term (its `go:`
filter matches the term and its descendants over the GO hierarchy). Each annotation
becomes one row: enzyme gene, histone mark, activity, source, GO identifier, and a
provenance URL. The two sources are concatenated into `hPTM_network.tsv` without
merging, so cross-source support is preserved.

### 3.5 Integration

The metabolite-gene, protein-protein, and enzyme-histone-modification networks are
reconciled to a common `Node1, Node2, Source, Type` schema, concatenated, and
de-duplicated on the node pair. The metabolite-gene edges with their provenance URLs
are exported separately; the full per-source enzyme-histone-mark table (with both
QuickGO and UniProtKB provenance URLs) is retained in `hPTM_network.tsv`. A URL audit
samples and checks the provenance URLs across all sources, including the two hPTM
sources.

### 3.6 Figures

At the end of step 3, `step3_final/visualization.py` writes summary figures to
`figures/`, re-deriving the reproducible plots from the notebooks
(`Collect_PKNdata_metabolites`, `Collect_PKNdata_proteins`, `Build_final_PKN`):

- per-source unique-interaction barplots for each layer (metabolite-gene, PPI, and
  the step-2b enzyme-histone-mark layer, the last also broken down by
  writer/eraser/reader activity);
- UpSet plots of the database overlap for each layer — which interactions are shared
  vs. unique across the source databases;
- a composition barplot of the final PKN's node and edge totals per type.

The overlap (UpSet) plots are computed from the long-format per-source network files
(`metabolite_gene_PKN.tsv`, `PPI_network.tsv`, `hPTM_network.tsv`), where an edge may
appear once per supporting source — not from the de-duplicated `LemonIte_PKN.tsv`.
Edges are treated as undirected (node pairs are sorted) so that A-B and B-A count as
one interaction. Figure generation uses the headless Agg backend and is best-effort:
a plotting failure (e.g. a missing optional dependency) is logged and skipped and
never fails a completed build. The notebook's MetalinksDB/MEBOCOST comparison and
HMDB-superclass figures are omitted, as they require external database files that are
not part of building or validating the network.


## 4. Output schema

| File                          | Columns                                | Description                              |
| ----------------------------- | -------------------------------------- | ---------------------------------------- |
| `LemonIte_PKN.tsv`            | `Node1, Node2, Source, Type`           | Combined, de-duplicated network          |
| `LemonIte_PKN_with_URLs.tsv`  | `HMDB, Gene, Source, URL`              | Metabolite-gene edges with provenance    |
| `metabolite_gene_PKN.tsv`     | `Metabolite, Gene, Source, url`        | Long-format step-1 network               |
| `PPI_network.tsv`             | `GeneA, GeneB, combined_score, Source` | Step-2 protein-protein network           |
| `hPTM_network.tsv`            | `Enzyme, Mark, Activity, Source, GO_ID, URL` | Step-2b enzyme-histone-mark network with provenance |
| `url_audit_report.csv`        | `Source, URL, status_code, reachable, spa_host` | URL audit results               |
| `figures/*.png`               | —                                      | Per-source count, database-overlap (UpSet), and composition figures (see §3.6) |

In `LemonIte_PKN.tsv`, metabolite nodes are written as `MetaboliteName_HMDBID`, gene
nodes as HGNC symbols, and histone-mark nodes as their controlled-vocabulary label
(for example `H3K27me3`). In `LemonIte_PKN_with_URLs.tsv`, the metabolite is reduced to
its HMDB accession. In `hPTM_network.tsv`, an enzyme-mark pair supported by both
QuickGO and UniProtKB appears as two rows (one per source); the combined
`LemonIte_PKN.tsv` de-duplicates it to a single edge.


## 5. Provenance URLs

Each edge links to a record that provides evidence for that specific interaction, not
to a database landing page. The templates are centralised in `config.URL_TEMPLATES`.

| Source             | URL target                | Evidence the user sees                                            |
| ------------------ | ------------------------- | ----------------------------------------------------------------- |
| BioGRID            | Chemical page             | The chemical's curated gene interactions                          |
| STITCH             | Chemical network page     | The chemical's protein interaction network                        |
| UniProtKB          | Protein entry             | The protein matched to the metabolite                             |
| IntAct             | Interaction record        | The specific binary interaction, by EBI accession                 |
| ChEMBL             | Compound report           | The compound's measured bioactivities and their targets          |
| Human-GEM (d1)     | Reaction page             | The single Metabolic Atlas reaction connecting metabolite to gene |
| Human-GEM (d2)     | Two reaction pages        | Both Metabolic Atlas reactions on the metabolite -> intermediate -> gene path (joined by ` -> `) |
| LINCS              | ChEMBL compound report    | The measured binding (LINCS biochemical data is ChEMBL-derived)   |
| MetalinksDB        | HMDB metabolite page      | The metabolite's protein associations                             |
| QuickGO            | Annotation viewer         | The gene product's annotations for that exact GO term (evidence codes, references), or the PubMed reference directly when the annotation carries a PMID |
| UniProtKB (GO)     | Protein entry, function   | The protein's GO molecular-function annotations for the histone activity |

Notes:

- For Human-GEM, a distance-1 edge is connected by a single reaction and links to that
  one reaction page. A distance-2 edge spans two reactions
  (metabolite -> intermediate metabolite -> gene); a single URL cannot represent it, so
  both reaction ids are tracked (`Reaction_ID_1`, `Reaction_ID_2`) and both reaction
  pages are kept (`URL_1`, `URL_2` in the link file, joined by ` -> ` in the final
  table's URL column). The URL audit checks each component independently.
- For IntAct, the molecule label reported by the API is used as the partner node;
  this is occasionally a protein complex or chain identifier rather than a single gene
  symbol. These records are retained: the interaction and its URL remain valid.
- For MetalinksDB, the originating source database name and homepage are retained in a
  secondary `Source_DB_URL` column in the link file.
- A `403` response (for example, hmdb.ca behind bot protection) is treated as
  reachable, as the page resolves in a browser; only `404` and `5xx` are treated as
  failures by the audit.


## 6. Reference network statistics

Indicative figures from a 500-metabolite test build (`PKN_test_500`). A full build
across all HMDB metabolites produces a substantially larger network.

Combined network: 1,077,337 edges (985,335 PPI; 92,002 metabolite-gene).

Metabolite-gene edges by source:

| Source            | Edges  |
| ----------------- | ------ |
| Human1_GEM_dist2  | 81,723 |
| Human1_GEM_dist1  | 5,174  |
| MetalinksDB       | 2,984  |
| BioGRID           | 758    |
| STITCH            | 628    |
| LINCS             | 276    |
| ChEMBL            | 230    |
| IntAct            | 229    |

Protein-protein edges by source:

| Source         | Edges   |
| -------------- | ------- |
| STRING         | 496,300 |
| BioGRID_genes  | 346,376 |
| HumanNet       | 141,256 |
| HuRI           | 1,403   |

Provenance-URL coverage of the metabolite-gene edges: 100 percent.

Enzyme-histone-modification edges (step 2b) are queried from GO annotations and are
independent of the metabolite-gene seed. An indicative run yields on the order of a
thousand annotations (roughly 500 QuickGO and 400 UniProtKB) covering about 500
distinct enzymes across some 30 histone marks and the writer, eraser, and reader
activity classes; after de-duplication on the node pair these collapse to the unique
enzyme-mark edges in `LemonIte_PKN.tsv`. Exact counts vary between runs as the GO
annotations are updated. Provenance-URL coverage of these edges is 100 percent.


## 7. Key parameters

| Parameter                    | Default | Meaning                                                  |
| ---------------------------- | ------- | -------------------------------------------------------- |
| `PATHWAY_DISTANCE`           | 2       | Maximum reaction hops for Human-GEM traversal            |
| `LINCS_IC50_THRESHOLD`       | 10000   | Maximum IC50 (nM) for a LINCS binding edge (10 micromolar) |
| `CHUNK_SIZE`                 | 800     | Metabolites per checkpointed processing chunk            |
| STRING confidence (implicit) | -       | STRING combined score is retained per edge for filtering |
| `HISTONE_GO_TERMS`           | 48 terms| Curated residue-specific histone GO terms queried in step 2b, each mapped to (mark, activity) |
| `QUICKGO_TAXON`              | 9606    | Taxon for the step-2b GO annotation queries (human)      |
| `QUICKGO_PAGE_LIMIT`         | 200     | Annotations per QuickGO page (the API maximum)           |

`GEM_UBIQUITOUS_METABOLITES` lists the currency metabolites (ATP, NAD+, water, and
similar) excluded from the Human-GEM reaction graph so that traversal reflects
specific biochemistry rather than shared cofactors.


## 8. Reproducibility and updates

The pipeline is deterministic given fixed input databases and is intended to be
re-run periodically (for example, monthly) to incorporate database updates. Local
database files are versioned by filename; API-derived data reflects the state of each
service at run time. Each run records a full log (`pipeline_progress.log`), an API
error log (`api_errors.log`), and the URL audit report (`url_audit_report.csv`).


## 9. Run ledger and change tracking

Every build automatically records itself and how the network changed, at the end of
step 3 (implemented in `step3_final/provenance.py`). This provides a permanent,
auditable history without any manual step. Four artefacts are produced:

| Artefact | Location | Content |
| -------- | -------- | ------- |
| `PKN_run_history.jsonl` | output dir | Append-only ledger, one JSON record per run: run id, UTC timestamp, git commit of the pipeline code, run parameters, node/edge totals, per-source edge counts, a content checksum, and the change summary versus the previous run. |
| `changes/<run_id>_diff.tsv.gz` | output dir | The full list of added and removed edges (Node1, Node2, Source) relative to the previous run. Written only when the network changed. |
| `PKN_edge_snapshot.tsv.gz` | output dir | The current run's edge set, used as the baseline for the next run's diff. |
| `PKN_CHANGELOG.md` | pipeline (repo) dir | A concise, human-readable entry per run, summarising totals and the per-source added/removed counts, with a pointer to the full diff. Being in the pipeline directory, it is git-versioned across runs. |

The run id is the UTC timestamp (`YYYYMMDDThhmmssZ`). The content checksum is an
order-independent SHA-256 of the `(Node1, Node2, Source, Type)` edge set, so two runs
that produce the same network share a checksum regardless of row order. Change
detection is on the `(Node1, Node2, Source)` key. The first recorded run establishes a
baseline (no diff). Provenance recording is best-effort: if it fails, the completed
build is unaffected.

To inspect history: read `PKN_CHANGELOG.md` for a summary, `PKN_run_history.jsonl` for
machine-readable per-run records, and the `changes/` diffs for exact edge-level changes.


## 10. Source database references

When using the network, please cite the underlying databases:

- HMDB: Wishart et al., Human Metabolome Database.
- BioGRID: Oughtred et al., Biological General Repository for Interaction Datasets.
- STITCH: Szklarczyk et al., Search Tool for Interactions of Chemicals.
- UniProt: The UniProt Consortium.
- IntAct: Del Toro et al., IntAct molecular interaction database.
- ChEMBL: Zdrazil et al., ChEMBL bioactivity database.
- LINCS / Laboratory of Systems Pharmacology biochemical data.
- Human-GEM: Robinson et al., Science Signaling 13 (2020), doi:10.1126/scisignal.aaz1482.
- MetalinksDB / OmniPath: Türei et al., and the contributing resources (CellPhoneDB,
  NicheNet, scConnect, and others).
- STRING: Szklarczyk et al., Search Tool for the Retrieval of Interacting Genes.
- HuRI: Luck et al., Human Reference Interactome.
- HumanNet: Kim et al., HumanNet functional gene network.
- Gene Ontology / QuickGO: The Gene Ontology Consortium; Binns et al., QuickGO.
- UniProt (GO annotations): The UniProt Consortium.
