"""
HPC configuration for the PKN Build Pipeline (Ghent University HPC / doduo cluster).

Inherits all defaults from config.py and overrides paths that differ on HPC.
Activated by setting the environment variable:
    export PKN_CONFIG=config_hpc

The submit scripts (submit_pkn_hpc.sh / submit_pkn_hpc_test.sh) set this
automatically, so no manual action is needed when submitting jobs.

HPC-specific env vars (set in submit scripts):
    PKN_WORKDIR        – root Lemonite directory  (default below)
    PKN_OUTPUT_DIR_NAME– output sub-folder name
    PKN_DB_DIR         – databases root directory (default below)
    PKN_GEM_DIR        – Human-GEM model directory (default below)
"""

import os

# ── bring in all shared defaults ──────────────────────────────────────────────
from config import *   # noqa: F401, F403  (intentional star-import for config)

# ── HPC root paths ────────────────────────────────────────────────────────────
_LEMONITE = '/user/gent/435/vsc43501/25BVM_Lemonite'

WORKDIR    = os.environ.get('PKN_WORKDIR',        _LEMONITE)
DB_DIR     = os.environ.get('PKN_DB_DIR',         os.path.join(_LEMONITE, 'databases'))
_GEM_DIR   = os.environ.get('PKN_GEM_DIR',        os.path.join(_LEMONITE, 'models/Human1-GEM/model'))

OUTPUT_DIR_NAME = os.environ.get('PKN_OUTPUT_DIR_NAME', 'PKN')
OUTPUT_DIR      = os.path.join(WORKDIR, OUTPUT_DIR_NAME)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Step 1 database paths (DB_DIR-relative, same structure as local) ──────────
HMDB_METABOLITES_XML = os.path.join(DB_DIR, 'HMDB/hmdb_metabolites.xml')
BIOGRID_LOCATION     = os.path.join(DB_DIR, 'BioGRID/BIOGRID-CHEMICALS-4.4.238.chemtab.txt')
STITCH_LOCATION      = os.path.join(DB_DIR, 'STITCHdb/9606.protein_chemical.links.v5.0.tsv')
UNIPROT_LOCATION     = os.path.join(DB_DIR, 'UniProt/uniprot_sprot_human.xml')
METALINKS_PATH       = os.path.join(DB_DIR, 'metalinks/metalinks.csv')

LINCS_COMPOUND_MAPPING   = os.path.join(DB_DIR, 'LINCS/lsp_compound_mapping.csv')
LINCS_TARGET_MAPPING     = os.path.join(DB_DIR, 'LINCS/lsp_target_mapping.csv')
LINCS_TARGET_DICTIONARY  = os.path.join(DB_DIR, 'LINCS/lsp_target_dictionary.csv')
LINCS_BIOCHEM_AGG        = os.path.join(DB_DIR, 'LINCS/lsp_biochem_agg.csv')
LINCS_REFERENCES         = os.path.join(DB_DIR, 'LINCS/lsp_references.csv')

# ── Human-GEM paths ───────────────────────────────────────────────────────────
GEM_PATH            = os.path.join(_GEM_DIR, 'Human-GEM.txt')
GEM_METABOLITES_PATH = os.path.join(_GEM_DIR, 'metabolites.tsv')
GEM_REACTIONS_PATH  = os.path.join(_GEM_DIR, 'reactions.tsv')
GEM_GENES_PATH      = os.path.join(_GEM_DIR, 'genes.tsv')

# ── Step 2 database paths ─────────────────────────────────────────────────────
BIOGRID_PPI_LOCATION  = os.path.join(DB_DIR, 'BioGRID/BIOGRID-ALL-4.4.238.tab3.txt')
HURI_LOCATION         = os.path.join(DB_DIR, 'HuRI/HuRI_Tong_2021.tsv')
HUMANNET_LOCATION     = os.path.join(DB_DIR, 'HumanNet/HS-PI.symbol.tsv')
ENSEMBL_MAPPING_FILE  = os.path.join(DB_DIR, 'Ensembl/ensembl_to_hgnc.tsv')

# ── Output files (re-derived from the HPC OUTPUT_DIR) ────────────────────────
DB_OUTPUT_FILES = {
    'IntAct':          os.path.join(OUTPUT_DIR, 'HMDB_metabolites_IntAct_processed.csv'),
    'UniProtKB':       os.path.join(OUTPUT_DIR, 'HMDB_metabolites_UniProtKB_processed.csv'),
    'chEMBL':          os.path.join(OUTPUT_DIR, 'HMDB_metabolites_chEMBL_processed.csv'),
    'STITCH':          os.path.join(OUTPUT_DIR, 'HMDB_metabolites_STITCH_processed.csv'),
    'BioGRID':         os.path.join(OUTPUT_DIR, 'HMDB_metabolites_BioGRID_processed.csv'),
    'LINCS':           os.path.join(OUTPUT_DIR, 'HMDB_metabolites_LINCS_processed.csv'),
    'Human1_GEM_dist1': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_Human1_GEM_dist1_processed.csv'),
    'Human1_GEM_dist2': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_Human1_GEM_dist2_processed.csv'),
    'MetalinksDB':     os.path.join(OUTPUT_DIR, 'HMDB_metabolites_MetalinksDB_processed.csv'),
    'HumanNet':        os.path.join(OUTPUT_DIR, 'HMDB_metabolites_HumanNet_processed.csv'),
}

URL_FILES = {
    'UniProtKB': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_UniProtKB_links.csv'),
    'IntAct':    os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_IntAct_links.csv'),
    'chEMBL':    os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_chEMBL_links.csv'),
    'STITCH':    os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_STITCH_links.csv'),
    'BioGRID':   os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_BioGRID_links.csv'),
    'LINCS':     os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_LINCS_links.csv'),
}

OUTPUT_FILE_FINAL = os.path.join(OUTPUT_DIR, 'HMDB_metabolites_gene_interactions.csv')

# Step 2 / 3 outputs
PPI_OUTPUT_FILE    = os.path.join(OUTPUT_DIR, 'PPI_network.tsv')
PPI_STRING_CACHE   = os.path.join(OUTPUT_DIR, 'cache/string_ppi_cache.csv')
PPI_BIOGRID_CACHE  = os.path.join(OUTPUT_DIR, 'cache/biogrid_ppi_cache.csv')
PPI_HUMANNET_CACHE = os.path.join(OUTPUT_DIR, 'cache/humannet_ppi_cache.csv')
PPI_HURI_CACHE     = os.path.join(OUTPUT_DIR, 'cache/huri_ppi_cache.csv')

METABOLITE_GENE_PKN      = os.path.join(OUTPUT_DIR, 'metabolite_gene_PKN.tsv')
FINAL_PKN_OUTPUT         = os.path.join(OUTPUT_DIR, 'LemonIte_PKN.tsv')
FINAL_PKN_WITH_URLS      = os.path.join(OUTPUT_DIR, 'LemonIte_PKN_with_URLs.tsv')

# Aliases
BIOGRID_PPI_FILE          = BIOGRID_PPI_LOCATION
HURI_FILE                 = HURI_LOCATION
PPI_NETWORK               = PPI_OUTPUT_FILE
FINAL_PKN_FILE            = FINAL_PKN_OUTPUT
FINAL_PKN_WITH_LINKS_FILE = FINAL_PKN_WITH_URLS

# Logging
LOG_FILE_API_ERRORS = os.path.join(OUTPUT_DIR, 'api_errors.log')
LOG_FILE_PIPELINE   = os.path.join(OUTPUT_DIR, 'pipeline_progress.log')
