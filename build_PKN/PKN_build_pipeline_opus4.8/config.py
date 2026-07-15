"""
Configuration for the LemonIte PKN (Prior Knowledge Network) build pipeline.

All paths, thresholds, API-retry settings and provenance-URL templates live here.
The most common settings are overridable via environment variables so that HPC /
container runs never need to edit this file:

    PKN_WORKDIR          root working directory (default below)
    PKN_OUTPUT_DIR_NAME  output sub-directory name (default 'PKN')
    PKN_DB_DIR           root directory holding the local database files
    PKN_GEM_DIR          directory holding the Human-GEM model files
    PKN_CONFIG           python module name to use as config (e.g. 'config_hpc')

This module is the single source of truth re-derived from the reference notebooks
(Collect_PKNdata_metabolites/proteins.ipynb, Build_final_PKN.ipynb). The provenance
URL templates in URL_TEMPLATES were audited against the live databases on 2026-06-23;
the broken ChEMBL report-card and LINCS lsp-cheminformatics URLs from the notebooks
have been corrected here (see URL_TEMPLATES below).
"""

import os
from pathlib import Path

# ===== Root configuration (env-overridable) =============================
WORKDIR = os.environ.get('PKN_WORKDIR', '/home/borisvdm/Documents/PhD/Lemonite')
OUTPUT_DIR_NAME = os.environ.get('PKN_OUTPUT_DIR_NAME', 'PKN')
DB_DIR = os.environ.get('PKN_DB_DIR', '/home/borisvdm/Documents/PhD/resources/databases')
GEM_DIR = os.environ.get(
    'PKN_GEM_DIR',
    '/home/borisvdm/Documents/PhD/resources/models/Human1-GEM/model',
)
# ========================================================================

OUTPUT_DIR = os.path.join(WORKDIR, OUTPUT_DIR_NAME)
FIGURES_DIR = os.path.join(OUTPUT_DIR, 'figures')
CACHE_DIR = os.path.join(OUTPUT_DIR, 'cache')

# =========================================================================
# STEP 1: METABOLITE-GENE INTERACTION SOURCES
# =========================================================================

HMDB_METABOLITES_XML = os.path.join(DB_DIR, 'HMDB/hmdb_metabolites.xml')
BIOGRID_LOCATION = os.path.join(DB_DIR, 'BioGRID/BIOGRID-CHEMICALS-4.4.238.chemtab.txt')
STITCH_LOCATION = os.path.join(DB_DIR, 'STITCHdb/9606.protein_chemical.links.v5.0.tsv')
METALINKS_PATH = os.path.join(DB_DIR, 'metalinks/metalinks.csv')

# Human-GEM metabolic model files
GEM_PATH = os.path.join(GEM_DIR, 'Human-GEM.txt')
GEM_METABOLITES_PATH = os.path.join(GEM_DIR, 'metabolites.tsv')
GEM_REACTIONS_PATH = os.path.join(GEM_DIR, 'reactions.tsv')
GEM_GENES_PATH = os.path.join(GEM_DIR, 'genes.tsv')

# LINCS / Laboratory of Systems Pharmacology biochemical binding data
LINCS_COMPOUND_MAPPING = os.path.join(DB_DIR, 'LINCS/lsp_compound_mapping.csv')
LINCS_TARGET_MAPPING = os.path.join(DB_DIR, 'LINCS/lsp_target_mapping.csv')
LINCS_TARGET_DICTIONARY = os.path.join(DB_DIR, 'LINCS/lsp_target_dictionary.csv')
LINCS_BIOCHEM_AGG = os.path.join(DB_DIR, 'LINCS/lsp_biochem_agg.csv')
LINCS_REFERENCES = os.path.join(DB_DIR, 'LINCS/lsp_references.csv')

# Thresholds (re-derived from notebook defaults)
PATHWAY_DISTANCE = 2          # GEM: max reaction hops metabolite -> gene
LINCS_IC50_THRESHOLD = 10000  # nM (= 10 uM) maximum IC50 to accept a LINCS edge
CHUNK_SIZE = 800              # metabolites per processing/checkpoint chunk
PAUSE_BETWEEN_CHUNKS = 8      # seconds to pause between checkpoint chunks

# Ubiquitous / currency metabolites excluded from the GEM graph traversal
# (Human-GEM metabolite IDs without compartment suffix)
GEM_UBIQUITOUS_METABOLITES = [
    'MAM02039', 'MAM02040', 'MAM01371', 'MAM02751', 'MAM01285',
    'MAM02555', 'MAM0254', 'MAM02630', 'MAM01597', 'MAM02552',
    'MAM02553', 'MAM02554', 'MAM02759', 'MAM02046', 'MAM01334',
    '2 MAM02039', '2 MAM02040',
]
GEM_COMPARTMENT_SUFFIXES = ['c', 'g', 'l', 'm', 'n', 'x', 'r', 'e']

# =========================================================================
# STEP 2: PROTEIN-PROTEIN INTERACTION SOURCES
# =========================================================================

BIOGRID_PPI_LOCATION = os.path.join(DB_DIR, 'BioGRID/BIOGRID-ALL-4.4.238.tab3.txt')
HURI_LOCATION = os.path.join(DB_DIR, 'HuRi/HuRI.tsv')
HUMANNET_LOCATION = os.path.join(DB_DIR, 'HumanNet/HS-PI.symbol.tsv')
ENSEMBL_MAPPING_FILE = os.environ.get(
    'PKN_ENSEMBL_MAPPING',
    os.path.join(DB_DIR, 'Ensembl/ensembl_mapping_jan2024.txt'),
)

# STRING REST API. Default to the current stable endpoint (string-db.org/api,
# v12) — the version-11-5 mirror the reference notebook pinned is slow/flaky and
# frequently times out. Same request format & TSV layout. Override via env if a
# specific STRING version is required.
STRING_API_URL = os.environ.get('PKN_STRING_API_URL', 'https://string-db.org/api')
STRING_SPECIES = 9606  # Homo sapiens
STRING_CALLER_IDENTITY = 'LemonIte'
STRING_CHUNK_SIZE = 1000
STRING_MAX_WORKERS = 10

# =========================================================================
# STEP 2b: ENZYME - HISTONE-MODIFICATION INTERACTIONS
# =========================================================================
# Edges between a histone-modifying protein (writer / eraser / reader) and the
# specific histone mark it acts on or binds — e.g. "EZH2 writes H3K27me3",
# "KDM6A erases H3K27me3", "BRD4 reads H3K27ac". This is enzyme-substrate
# specificity, NOT ChIP-seq genomic proximity.
#
# Retrieved from Gene Ontology molecular-function annotations via two REST APIs:
#   - QuickGO   (EBI)      : ebi.ac.uk/QuickGO/services/annotation/search
#   - UniProtKB (rest API) : rest.uniprot.org/uniprotkb/search  (go: filter)
# Both are queried per residue-specific GO term below; QuickGO auto-expands each
# term to its descendants, so this is a curated set of parent terms, not every
# leaf. taxon = human (9606). Each annotation carries a GO evidence code (ECO)
# and a reference (GO_REF or PMID).
#
# HISTONE_GO_TERMS maps each GO id -> (mark, activity) where `mark` is the
# controlled-vocabulary histone-mark node used as Node2 in the final network
# (e.g. 'H3K27me3') and `activity` is one of 'writer' | 'eraser' | 'reader'.
# The GO ids were verified live against QuickGO's ontology search on 2026-07-03.

HISTONE_GO_TERMS = {
    # ---- Writers: methyltransferases ----
    'GO:0042800': ('H3K4',   'writer'),   # histone H3K4 methyltransferase activity
    'GO:0046974': ('H3K9',   'writer'),   # histone H3K9 methyltransferase activity
    'GO:0140947': ('H3K9me2','writer'),   # histone H3K9me2 methyltransferase activity
    'GO:0046976': ('H3K27',  'writer'),   # histone H3K27 methyltransferase activity
    'GO:0046975': ('H3K36',  'writer'),   # histone H3K36 methyltransferase activity
    'GO:0031151': ('H3K79',  'writer'),   # histone H3K79 methyltransferase activity
    'GO:0140759': ('H3K56',  'writer'),   # histone H3K56 methyltransferase activity
    'GO:0042799': ('H4K20',  'writer'),   # histone H4K20 methyltransferase activity
    'GO:0140941': ('H4K20me','writer'),   # histone H4K20me methyltransferase activity
    # ---- Writers: acetyltransferases ----
    'GO:0004402': ('histone_ac', 'writer'),  # histone acetyltransferase activity (general)
    'GO:0043992': ('H3K9ac',  'writer'),  # histone H3K9 acetyltransferase activity
    'GO:0036408': ('H3K14ac', 'writer'),  # histone H3K14 acetyltransferase activity
    'GO:0044017': ('H3K27ac', 'writer'),  # histone H3K27 acetyltransferase activity
    'GO:0044018': ('H3K36ac', 'writer'),  # histone H3K36 acetyltransferase activity
    'GO:0032931': ('H3K56ac', 'writer'),  # histone H3K56 acetyltransferase activity
    'GO:0046972': ('H4K16ac', 'writer'),  # histone H4K16 acetyltransferase activity
    # ---- Writers: ubiquitin ligases ----
    'GO:0140862': ('H2AK119ub', 'writer'),  # histone H2AK119 ubiquitin ligase activity
    'GO:0141054': ('H2Bub',     'writer'),  # histone H2B ubiquitin ligase activity
    # ---- Erasers: demethylases ----
    'GO:0032453': ('H3K4',   'eraser'),   # histone H3K4 demethylase activity
    'GO:0034647': ('H3K4me', 'eraser'),   # histone H3K4me/me2/me3 demethylase activity
    'GO:0032454': ('H3K9',   'eraser'),   # histone H3K9 demethylase activity
    'GO:0071558': ('H3K27',  'eraser'),   # histone H3K27me2/me3 demethylase activity
    'GO:0051864': ('H3K36',  'eraser'),   # histone H3K36 demethylase activity
    'GO:0035575': ('H4K20',  'eraser'),   # histone H4K20 demethylase activity
    'GO:0140760': ('H3K56',  'eraser'),   # histone H3K56me2/me3 demethylase activity
    # ---- Erasers: deacetylases ----
    'GO:0004407': ('histone_ac', 'eraser'),  # histone deacetylase activity (general)
    'GO:0141050': ('H3Kac',   'eraser'),  # histone H3K deacetylase activity
    'GO:0141051': ('H4Kac',   'eraser'),  # histone H4K deacetylase activity
    'GO:0046969': ('H3K9ac',  'eraser'),  # histone H3K9 deacetylase activity, NAD-dependent
    'GO:0032041': ('H3K14ac', 'eraser'),  # histone H3K14 deacetylase activity, NAD-dependent
    'GO:0097372': ('H3K18ac', 'eraser'),  # histone H3K18 deacetylase activity, NAD-dependent
    'GO:0140765': ('H3K56ac', 'eraser'),  # histone H3K56 deacetylase activity, NAD-dependent
    'GO:0046970': ('H4K16ac', 'eraser'),  # histone H4K16 deacetylase activity, NAD-dependent
    # ---- Readers: methyl / acetyl / ubiquitin binding ----
    'GO:0042393': ('histone', 'reader'),  # histone binding (general)
    'GO:0140566': ('histone', 'reader'),  # histone reader activity (general)
    'GO:0140002': ('H3K4me3',  'reader'), # histone H3K4me3 reader activity
    'GO:0140109': ('H3K4me1',  'reader'), # histone H3K4me1 reader activity
    'GO:0062072': ('H3K9me2/3','reader'), # histone H3K9me2/3 reader activity
    'GO:0061628': ('H3K27me3', 'reader'), # histone H3K27me3 reader activity
    'GO:0140003': ('H3K36me3', 'reader'), # histone H3K36me3 reader activity
    'GO:0140072': ('H3K9ac',   'reader'), # histone H3K9ac reader activity
    'GO:0140119': ('H3K27ac',  'reader'), # histone H3K27ac reader activity
    'GO:0140015': ('H3K14ac',  'reader'), # histone H3K14ac reader activity
    'GO:0140046': ('H4K16ac',  'reader'), # histone H4K16ac reader activity
    'GO:0140117': ('H4K20me1', 'reader'), # histone H4K20me1 reader activity
    'GO:0140005': ('H4K20me2', 'reader'), # histone H4K20me2 reader activity
    'GO:1990889': ('H4K20me3', 'reader'), # histone H4K20me3 reader activity
    'GO:0061649': ('histone_ub','reader'),# ubiquitin-modified histone reader activity
}

# QuickGO annotation query: batch several GO ids per request (comma-separated),
# human only, protein gene products, one evidence line per row.
QUICKGO_TAXON = 9606
QUICKGO_BATCH_SIZE = 10       # GO ids per annotation-search request
QUICKGO_PAGE_LIMIT = 200      # annotations per page (QuickGO max is 200)

# =========================================================================
# STEP 3: FINAL PKN INTEGRATION (comparison databases, optional)
# =========================================================================

METALINKS_DB = METALINKS_PATH
MEBOCOST_DB = os.path.join(
    DB_DIR, 'MEBOCOST/metabolite_associated_gene_reaction_HMDB_summary.tsv'
)

# =========================================================================
# OUTPUT FILES
# =========================================================================

# Parsed HMDB annotation table (also holds Canonical_smiles + ChEMBL_id after preprocessing)
HMDB_ANNOTATED_FILE = os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated.csv')
CHEMBL_MAPPING_FILE = os.path.join(OUTPUT_DIR, 'HMDB_metabolites_ChEMBL_mapping.csv')

# Step-1 per-database processed caches (one row per metabolite, pipe-joined genes)
DB_OUTPUT_FILES = {
    'IntAct': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_IntAct_processed.csv'),
    'UniProtKB': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_UniProtKB_processed.csv'),
    'chEMBL': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_chEMBL_processed.csv'),
    'STITCH': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_STITCH_processed.csv'),
    'BioGRID': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_BioGRID_processed.csv'),
    'LINCS': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_LINCS_processed.csv'),
    'Human1_GEM_dist1': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_Human1_GEM_dist1_processed.csv'),
    'Human1_GEM_dist2': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_Human1_GEM_dist2_processed.csv'),
    'MetalinksDB': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_MetalinksDB_processed.csv'),
}

# Step-1 per-database link files (one row per interaction, with provenance URL)
URL_FILES = {
    'UniProtKB': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_UniProtKB_links.csv'),
    'IntAct': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_IntAct_links.csv'),
    'chEMBL': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_chEMBL_links.csv'),
    'STITCH': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_STITCH_links.csv'),
    'BioGRID': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_BioGRID_links.csv'),
    'LINCS': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_LINCS_links.csv'),
    'MetalinksDB': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_MetalinksDB_links.csv'),
    'Human1_GEM_dist1': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_Human1GEM_links.csv'),
    'Human1_GEM_dist2': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_Human1GEM_links.csv'),
}

# Step-1 merged outputs
OUTPUT_FILE_FINAL = os.path.join(OUTPUT_DIR, 'HMDB_metabolites_gene_interactions.csv')
METABOLITE_GENE_PKN = os.path.join(OUTPUT_DIR, 'metabolite_gene_PKN.tsv')

# Step-2 outputs
STRING_INTERACTIONS_FILE = os.path.join(OUTPUT_DIR, 'STRING_interactions.csv')
BIOGRID_GENES_FILE = os.path.join(OUTPUT_DIR, 'BioGRID_genes_interactions.csv')
HUMANNET_GENES_FILE = os.path.join(OUTPUT_DIR, 'HumanNet_interactions.csv')
HURI_GENES_FILE = os.path.join(OUTPUT_DIR, 'HuRI_pruned_interactions.csv')
PPI_NETWORK_FILE = os.path.join(OUTPUT_DIR, 'PPI_network.tsv')

# Step-2b outputs (enzyme - histone-modification edges). Per-source caches are
# whole-file (like STRING): one query pass over the GO term set, cached to disk.
# The merged long-format network carries the provenance URL directly.
QUICKGO_HPTM_FILE = os.path.join(OUTPUT_DIR, 'QuickGO_hPTM_interactions.csv')
UNIPROT_HPTM_FILE = os.path.join(OUTPUT_DIR, 'UniProt_hPTM_interactions.csv')
HPTM_NETWORK_FILE = os.path.join(OUTPUT_DIR, 'hPTM_network.tsv')

# Step-2c phospho-regulator network (STANDALONE — deliberately NOT merged into the main
# LemonIte_PKN.tsv; it is a separate prior-knowledge layer consumed by the phospho-regulator
# analysis). Edge types: kinase-substrate (OmniPath enzsub) + phospho-TF-target (CollecTRI).
OMNIPATH_ENZSUB_URL = os.environ.get(
    'PKN_OMNIPATH_ENZSUB_URL',
    'https://omnipathdb.org/enzsub')
# CollecTRI TF->target network shipped with the repo (source, target, mor).
COLLECTRI_FILE = os.environ.get(
    'PKN_COLLECTRI_FILE',
    os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                 'nextflow', 'PKN', 'CollecTRI_network.txt'))
OMNIPATH_KSN_FILE = os.path.join(OUTPUT_DIR, 'OmniPath_kinase_substrate.tsv')
PHOSPHO_PKN_FILE = os.path.join(OUTPUT_DIR, 'phospho_pkn.tsv')

# Step-3 final outputs
FINAL_PKN_FILE = os.path.join(OUTPUT_DIR, 'LemonIte_PKN.tsv')
FINAL_PKN_WITH_LINKS_FILE = os.path.join(OUTPUT_DIR, 'LemonIte_PKN_with_URLs.tsv')
URL_AUDIT_REPORT = os.path.join(OUTPUT_DIR, 'url_audit_report.csv')

# Run ledger + PKN change tracking (written automatically at the end of step 3).
# RUN_HISTORY is an append-only JSON-lines ledger (one record per run). CHANGES_DIR
# holds the per-run full edge diffs. PKN_SNAPSHOT is the previous run's edge set,
# used to compute the next diff. CHANGELOG is a human-readable summary that is meant
# to be committed to the repo (see also the pipeline-dir copy in reconfigure).
RUN_HISTORY_FILE = os.path.join(OUTPUT_DIR, 'PKN_run_history.jsonl')
CHANGES_DIR = os.path.join(OUTPUT_DIR, 'changes')
PKN_SNAPSHOT_FILE = os.path.join(OUTPUT_DIR, 'PKN_edge_snapshot.tsv.gz')
# Human-readable changelog. Lives in the pipeline (repo) directory so it is
# git-versioned across runs, independent of the (gitignored) output directory.
PKN_CHANGELOG_FILE = os.path.join(os.path.dirname(__file__), 'PKN_CHANGELOG.md')

# =========================================================================
# PROVENANCE URL TEMPLATES  (audited 2026-06-23 against live databases)
# =========================================================================
# These are the single source of truth for every clickable source URL.  Each is
# a python str.format template.  url_audit.py samples emitted URLs and checks them.
#
# FIXES applied relative to the reference notebooks:
#   - ChEMBL: compound_report_card / target_report_card serve only an empty SPA
#     shell -> replaced with the current canonical explore/* browse URLs.
#   - LINCS: labsyspharm.github.io/.../compounds/{id} is 404 (toolkit docs, no
#     per-compound pages) -> PubMed reference is preferred; the fallback is the
#     HMS-LINCS small-molecule browser root.
#   - HMDB: the malformed 'http://www.hmdb.ca}metabolite' template -> proper
#     per-metabolite URL.
URL_TEMPLATES = {
    'UniProtKB':        'https://www.uniprot.org/uniprotkb/{accession}/entry',
    'IntAct':           'https://www.ebi.ac.uk/intact/details/interaction/{ebi_id}',
    'chEMBL_compound':  'https://www.ebi.ac.uk/chembl/explore/compound/{chembl_id}',
    'chEMBL_target':    'https://www.ebi.ac.uk/chembl/explore/target/{target_id}',
    'STITCH':           'http://stitch.embl.de/cgi/show_network_section.pl?identifier=CID{cid}&species=9606',
    'BioGRID':          'https://thebiogrid.org/chemical/{chemical_id}/{chemical_name}.html',
    'PubMed':           'https://pubmed.ncbi.nlm.nih.gov/{pmid}',
    'LINCS_fallback':   'https://lincs.hms.harvard.edu/db/sm/',
    'HMDB':             'https://hmdb.ca/metabolites/{hmdb_id}',
    'GEM_dataset':      'https://github.com/SysBioChalmers/Human-GEM',
    'GEM_paper':        'https://doi.org/10.1126/scisignal.aaz1482',
    # Per-reaction Metabolic Atlas page — traces a metabolite->gene edge back to
    # the specific Human-GEM reaction(s) that connect them.
    'GEM_reaction':     'https://metabolicatlas.org/explore/Human-GEM/gem-browser/reaction/{reaction_id}',
    # Enzyme - histone-mark edges (step 2b). The QuickGO annotation viewer,
    # filtered to this gene product + GO term, shows the user exactly the
    # evidence (evidence codes + references) supporting that specific activity.
    'QuickGO_annotation': 'https://www.ebi.ac.uk/QuickGO/annotations?geneProductId={accession}&goId={go_id}',
    # UniProt entry's function/GO section for the same enzyme.
    'UniProt_entry':      'https://www.uniprot.org/uniprotkb/{accession}/entry#function',
}

# MetalinksDB per-source provenance URL roots
METALINKS_SOURCE_URLS = {
    'CellPhoneDB': 'https://www.cellphonedb.org',
    'NicheNet':    'https://github.com/saeyslab/nichenetr',
    'HMDB':        'https://hmdb.ca',
    'MetalinksDB': 'https://metalinks.omnipathdb.org',
}
METALINKS_DEFAULT_URL = 'https://metalinks.omnipathdb.org'

# Aliases kept for backward-compatibility with the GEM module
GEM_GITHUB_URL = URL_TEMPLATES['GEM_dataset']
GEM_PAPER_URL = URL_TEMPLATES['GEM_paper']

# =========================================================================
# API RETRY CONFIGURATION
# =========================================================================
# Per-database retry/throttle profiles. backoff_factor controls exponential
# backoff (wait = backoff_factor ** attempt, with jitter); pause_after/duration
# inject a periodic cooldown so threaded request bursts stay under rate limits.
API_RETRY_CONFIG = {
    'UniProtKB':      {'max_retries': 10, 'backoff_factor': 3, 'timeout': 30, 'max_workers': 3, 'pause_after': 50,  'pause_duration': 15},
    'IntAct':         {'max_retries': 10, 'backoff_factor': 2, 'timeout': 25, 'max_workers': 4, 'pause_after': 75,  'pause_duration': 10},
    'chEMBL':         {'max_retries': 10, 'backoff_factor': 2, 'timeout': 25, 'max_workers': 5, 'pause_after': 100, 'pause_duration': 8},
    'ChEMBL_Mapping': {'max_retries': 10, 'backoff_factor': 2, 'timeout': 25, 'max_workers': 5, 'pause_after': 100, 'pause_duration': 8},
    'STRING':         {'max_retries': 5,  'backoff_factor': 2, 'timeout': 30, 'max_workers': 10,'pause_after': 10,  'pause_duration': 5},
    'MyGene':         {'max_retries': 5,  'backoff_factor': 2, 'timeout': 30, 'max_workers': 1, 'pause_after': 5,   'pause_duration': 3},
    'QuickGO':        {'max_retries': 8,  'backoff_factor': 2, 'timeout': 30, 'max_workers': 4, 'pause_after': 40,  'pause_duration': 8},
    'UniProt_hPTM':   {'max_retries': 8,  'backoff_factor': 3, 'timeout': 30, 'max_workers': 3, 'pause_after': 40,  'pause_duration': 10},
}
MAX_WORKERS_DEFAULT = 4

# =========================================================================
# LOGGING
# =========================================================================
LOG_FILE_API_ERRORS = os.path.join(OUTPUT_DIR, 'api_errors.log')
LOG_FILE_PIPELINE = os.path.join(OUTPUT_DIR, 'pipeline_progress.log')
LOG_LEVEL = 'INFO'
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'


def ensure_dirs():
    """Create the output, figures and cache directories if missing."""
    for d in (OUTPUT_DIR, FIGURES_DIR, CACHE_DIR):
        os.makedirs(d, exist_ok=True)


def reconfigure_output_dir(output_dir_name: str):
    """
    Repoint OUTPUT_DIR and every derived output path to a new directory.

    Accepts a name relative to WORKDIR or an absolute path. Used by the
    ``--output-dir`` CLI flag for test runs without editing this module.
    """
    import sys
    mod = sys.modules[__name__]

    new_dir = output_dir_name if os.path.isabs(output_dir_name) \
        else os.path.join(mod.WORKDIR, output_dir_name)

    mod.OUTPUT_DIR_NAME = output_dir_name
    mod.OUTPUT_DIR = new_dir
    mod.FIGURES_DIR = os.path.join(new_dir, 'figures')
    mod.CACHE_DIR = os.path.join(new_dir, 'cache')

    mod.HMDB_ANNOTATED_FILE = os.path.join(new_dir, 'HMDB_metabolites_annotated.csv')
    mod.CHEMBL_MAPPING_FILE = os.path.join(new_dir, 'HMDB_metabolites_ChEMBL_mapping.csv')
    mod.DB_OUTPUT_FILES = {k: os.path.join(new_dir, os.path.basename(v))
                           for k, v in mod.DB_OUTPUT_FILES.items()}
    mod.URL_FILES = {k: os.path.join(new_dir, os.path.basename(v))
                     for k, v in mod.URL_FILES.items()}
    mod.OUTPUT_FILE_FINAL = os.path.join(new_dir, 'HMDB_metabolites_gene_interactions.csv')
    mod.METABOLITE_GENE_PKN = os.path.join(new_dir, 'metabolite_gene_PKN.tsv')

    mod.STRING_INTERACTIONS_FILE = os.path.join(new_dir, 'STRING_interactions.csv')
    mod.BIOGRID_GENES_FILE = os.path.join(new_dir, 'BioGRID_genes_interactions.csv')
    mod.HUMANNET_GENES_FILE = os.path.join(new_dir, 'HumanNet_interactions.csv')
    mod.HURI_GENES_FILE = os.path.join(new_dir, 'HuRI_pruned_interactions.csv')
    mod.PPI_NETWORK_FILE = os.path.join(new_dir, 'PPI_network.tsv')

    mod.QUICKGO_HPTM_FILE = os.path.join(new_dir, 'QuickGO_hPTM_interactions.csv')
    mod.UNIPROT_HPTM_FILE = os.path.join(new_dir, 'UniProt_hPTM_interactions.csv')
    mod.HPTM_NETWORK_FILE = os.path.join(new_dir, 'hPTM_network.tsv')

    mod.OMNIPATH_KSN_FILE = os.path.join(new_dir, 'OmniPath_kinase_substrate.tsv')
    mod.PHOSPHO_PKN_FILE = os.path.join(new_dir, 'phospho_pkn.tsv')

    mod.FINAL_PKN_FILE = os.path.join(new_dir, 'LemonIte_PKN.tsv')
    mod.FINAL_PKN_WITH_LINKS_FILE = os.path.join(new_dir, 'LemonIte_PKN_with_URLs.tsv')
    mod.URL_AUDIT_REPORT = os.path.join(new_dir, 'url_audit_report.csv')

    mod.RUN_HISTORY_FILE = os.path.join(new_dir, 'PKN_run_history.jsonl')
    mod.CHANGES_DIR = os.path.join(new_dir, 'changes')
    mod.PKN_SNAPSHOT_FILE = os.path.join(new_dir, 'PKN_edge_snapshot.tsv.gz')
    # PKN_CHANGELOG_FILE intentionally NOT repointed: the changelog is a single
    # repo-versioned file shared across output directories.

    mod.LOG_FILE_API_ERRORS = os.path.join(new_dir, 'api_errors.log')
    mod.LOG_FILE_PIPELINE = os.path.join(new_dir, 'pipeline_progress.log')

    # Create the dirs from `mod` directly. Do NOT call ensure_dirs() here: when an
    # alternate config (PKN_CONFIG=config_hpc) is active, `mod` is that aliased
    # module, but ensure_dirs()'s globals are this base module — so it would create
    # the stale OUTPUT_DIR instead of new_dir.
    for d in (mod.OUTPUT_DIR, mod.FIGURES_DIR, mod.CACHE_DIR):
        os.makedirs(d, exist_ok=True)
