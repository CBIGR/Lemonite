"""
Configuration file for the PKN (Prior Knowledge Network) pipeline.

This module contains all paths, settings, and parameters used across
the three-step pipeline:
  - Step 1: Metabolite-gene interactions
  - Step 2: Protein-protein interactions
  - Step 3: Final PKN integration
"""

import os
from pathlib import Path

# ===== CONFIGURATION - Change output directory here =====
# These can be overridden via environment variables for container use:
#   PKN_WORKDIR, PKN_OUTPUT_DIR_NAME, PKN_DB_DIR
WORKDIR = os.environ.get('PKN_WORKDIR', '/home/borisvdm/Documents/PhD/Lemonite')
OUTPUT_DIR_NAME = os.environ.get('PKN_OUTPUT_DIR_NAME', 'PKN')
# ========================================================

# Set up paths
OUTPUT_DIR = os.path.join(WORKDIR, OUTPUT_DIR_NAME)
DB_DIR = os.environ.get('PKN_DB_DIR', '/run/media/borisvdm/5250-D000/resources/databases')

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================================================================
# STEP 1: METABOLITE-GENE INTERACTION SOURCES
# =========================================================================

# HMDB metabolites XML database
HMDB_METABOLITES_XML = os.path.join(DB_DIR, 'HMDB/hmdb_metabolites.xml')

# BioGRID chemical-protein interactions
BIOGRID_LOCATION = os.path.join(DB_DIR, 'BioGRID/BIOGRID-CHEMICALS-4.4.238.chemtab.txt')

# STITCH protein-chemical links
STITCH_LOCATION = os.path.join(DB_DIR, 'STITCHdb/9606.protein_chemical.links.v5.0.tsv')

# UniProt human proteome XML
UNIPROT_LOCATION = os.path.join(DB_DIR, 'UniProt/uniprot_sprot_human.xml')

# Human-GEM metabolic model (https://github.com/SysBioChalmers/Human-GEM)
GEM_PATH = '/home/borisvdm/Documents/PhD/resources/models/Human1-GEM/model/Human-GEM.txt'
GEM_METABOLITES_PATH = '/home/borisvdm/Documents/PhD/resources/models/Human1-GEM/model/metabolites.tsv'
GEM_REACTIONS_PATH = '/home/borisvdm/Documents/PhD/resources/models/Human1-GEM/model/reactions.tsv'
GEM_GENES_PATH = '/home/borisvdm/Documents/PhD/resources/models/Human1-GEM/model/genes.tsv'

# GEM provenance URLs
GEM_GITHUB_URL = 'https://github.com/SysBioChalmers/Human-GEM'
GEM_PAPER_URL = 'https://doi.org/10.1126/scisignal.aaz1482'

# MetalinksDB (from liana+ package)
METALINKS_PATH = os.path.join(DB_DIR, 'metalinks/metalinks.csv')

# LINCS biochemical binding data
LINCS_COMPOUND_MAPPING = os.path.join(DB_DIR, 'LINCS/lsp_compound_mapping.csv')
LINCS_TARGET_MAPPING = os.path.join(DB_DIR, 'LINCS/lsp_target_mapping.csv')
LINCS_TARGET_DICTIONARY = os.path.join(DB_DIR, 'LINCS/lsp_target_dictionary.csv')
LINCS_BIOCHEM_AGG = os.path.join(DB_DIR, 'LINCS/lsp_biochem_agg.csv')
LINCS_REFERENCES = os.path.join(DB_DIR, 'LINCS/lsp_references.csv')

# Pathway distance for GEM (how many reactions away to include)
PATHWAY_DISTANCE = 2

# Ubiquitous metabolites to exclude from GEM graph (cofactors, common currency)
# These are Human-GEM metabolite IDs (without compartment suffix)
GEM_UBIQUITOUS_METABOLITES = [
    'MAM02039', 'MAM02040', 'MAM01371', 'MAM02751', 'MAM01285',
    'MAM02555', 'MAM0254', 'MAM02630', 'MAM01597', 'MAM02552',
    'MAM02553', 'MAM02554', 'MAM02759', 'MAM02046', 'MAM01334',
    '2 MAM02039', '2 MAM02040'
]

# Compartment suffixes used in Human-GEM
GEM_COMPARTMENT_SUFFIXES = ['c', 'g', 'l', 'm', 'n', 'x', 'r', 'e']

# =========================================================================
# STEP 2: PROTEIN-PROTEIN INTERACTION SOURCES
# =========================================================================

# BioGRID PPI database
BIOGRID_PPI_LOCATION = os.path.join(DB_DIR, 'BioGRID/BIOGRID-ALL-4.4.238.tab3.txt')

# HuRI human protein-protein interactome
HURI_LOCATION = os.path.join(DB_DIR, 'HuRI/HuRI_Tong_2021.tsv')

# HumanNet PPI interactions
HUMANNET_LOCATION = os.path.join(DB_DIR, 'HumanNet/HS-PI.symbol.tsv')

# Ensembl ID to gene symbol mapping (for HuRI)
ENSEMBL_MAPPING_FILE = os.path.join(DB_DIR, 'Ensembl/ensembl_to_hgnc.tsv')

# STRING API endpoint
STRING_API_URL = "https://string-db.org/api"
STRING_VERSION = "12.0"
STRING_SPECIES = "9606"  # Human

# =========================================================================
# STEP 3: FINAL PKN INTEGRATION
# =========================================================================

# Comparative databases for validation
METADB_PATH = os.path.join(DB_DIR, 'MetalinksDB/metadb_interactions.csv')
MEBOCOST_PATH = os.path.join(DB_DIR, 'MEBOCOST/mebocost_interactions.csv')

# =========================================================================
# OUTPUT FILES (STEP 1: METABOLITES)
# =========================================================================

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
    'HumanNet': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_HumanNet_processed.csv')
}

URL_FILES = {
    'UniProtKB': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_UniProtKB_links.csv'),
    'IntAct': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_IntAct_links.csv'),
    'chEMBL': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_chEMBL_links.csv'),
    'STITCH': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_STITCH_links.csv'),
    'BioGRID': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_BioGRID_links.csv'),
    'LINCS': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_annotated_LINCS_links.csv')
}

OUTPUT_FILE_FINAL = os.path.join(OUTPUT_DIR, 'HMDB_metabolites_gene_interactions.csv')

# =========================================================================
# OUTPUT FILES (STEP 2: PPI)
# =========================================================================

PPI_OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'PPI_network.tsv')
PPI_STRING_CACHE = os.path.join(OUTPUT_DIR, 'cache/string_ppi_cache.csv')
PPI_BIOGRID_CACHE = os.path.join(OUTPUT_DIR, 'cache/biogrid_ppi_cache.csv')
PPI_HUMANNET_CACHE = os.path.join(OUTPUT_DIR, 'cache/humannet_ppi_cache.csv')
PPI_HURI_CACHE = os.path.join(OUTPUT_DIR, 'cache/huri_ppi_cache.csv')

# =========================================================================
# OUTPUT FILES (STEP 3: FINAL PKN)
# =========================================================================

METABOLITE_GENE_PKN = os.path.join(OUTPUT_DIR, 'metabolite_gene_PKN.tsv')
FINAL_PKN_OUTPUT = os.path.join(OUTPUT_DIR, 'LemonIte_PKN.tsv')
FINAL_PKN_WITH_URLS = os.path.join(OUTPUT_DIR, 'LemonIte_PKN_with_URLs.tsv')

# Aliases used by various modules
BIOGRID_PPI_FILE = BIOGRID_PPI_LOCATION
HURI_FILE = HURI_LOCATION
PPI_NETWORK = PPI_OUTPUT_FILE
FINAL_PKN_FILE = FINAL_PKN_OUTPUT
FINAL_PKN_WITH_LINKS_FILE = FINAL_PKN_WITH_URLS

# =========================================================================
# API RETRY CONFIGURATION
# =========================================================================

API_RETRY_CONFIG = {
    'UniProtKB': {
        'max_retries': 10,
        'backoff_factor': 3,  # Wait time multiplier (exponential backoff)
        'timeout': 30,  # Request timeout in seconds
        'max_workers': 3,  # Concurrent threads (lower = gentler on API)
        'pause_after': 50,  # Pause every N requests
        'pause_duration': 15  # Pause duration in seconds
    },
    'IntAct': {
        'max_retries': 10,
        'backoff_factor': 2,
        'timeout': 25,
        'max_workers': 4,
        'pause_after': 75,
        'pause_duration': 10
    },
    'chEMBL': {
        'max_retries': 10,
        'backoff_factor': 2,
        'timeout': 25,
        'max_workers': 5,
        'pause_after': 100,
        'pause_duration': 8
    },
    'ChEMBL_Mapping': {
        'max_retries': 10,
        'backoff_factor': 2,
        'timeout': 25,
        'max_workers': 5,
        'pause_after': 100,
        'pause_duration': 8
    },
    'STRING': {
        'max_retries': 5,
        'backoff_factor': 2,
        'timeout': 30,
        'max_workers': 3,
        'pause_after': 10,
        'pause_duration': 5
    },
    'MyGene': {
        'max_retries': 5,
        'backoff_factor': 2,
        'timeout': 30,
        'max_workers': 1,
        'pause_after': 5,
        'pause_duration': 3
    }
}

# =========================================================================
# PROCESSING PARAMETERS
# =========================================================================

CHUNK_SIZE = 800  # For batching metabolite processing
MAX_WORKERS_DEFAULT = 4  # Default thread pool size

# =========================================================================
# LOGGING CONFIGURATION
# =========================================================================

LOG_FILE_API_ERRORS = os.path.join(OUTPUT_DIR, 'api_errors.log')
LOG_FILE_PIPELINE = os.path.join(OUTPUT_DIR, 'pipeline_progress.log')
LOG_LEVEL = 'INFO'
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'


def reconfigure_output_dir(output_dir_name: str):
    """
    Update OUTPUT_DIR and all derived output paths to use a new directory name.

    Parameters
    ----------
    output_dir_name : str
        Directory name (relative to WORKDIR) or an absolute path.
    """
    import sys
    mod = sys.modules[__name__]

    if os.path.isabs(output_dir_name):
        new_dir = output_dir_name
    else:
        new_dir = os.path.join(mod.WORKDIR, output_dir_name)

    os.makedirs(new_dir, exist_ok=True)
    mod.OUTPUT_DIR = new_dir
    mod.OUTPUT_DIR_NAME = output_dir_name

    # Step 1 output files
    mod.DB_OUTPUT_FILES = {k: os.path.join(new_dir, os.path.basename(v))
                           for k, v in mod.DB_OUTPUT_FILES.items()}
    mod.URL_FILES = {k: os.path.join(new_dir, os.path.basename(v))
                     for k, v in mod.URL_FILES.items()}
    mod.OUTPUT_FILE_FINAL = os.path.join(new_dir, 'HMDB_metabolites_gene_interactions.csv')

    # Step 2 output files
    mod.PPI_OUTPUT_FILE = os.path.join(new_dir, 'PPI_network.tsv')
    mod.PPI_STRING_CACHE = os.path.join(new_dir, 'cache/string_ppi_cache.csv')
    mod.PPI_BIOGRID_CACHE = os.path.join(new_dir, 'cache/biogrid_ppi_cache.csv')
    mod.PPI_HURI_CACHE = os.path.join(new_dir, 'cache/huri_ppi_cache.csv')

    # Step 3 output files
    mod.METABOLITE_GENE_PKN = os.path.join(new_dir, 'metabolite_gene_PKN.tsv')
    mod.FINAL_PKN_OUTPUT = os.path.join(new_dir, 'LemonIte_PKN.tsv')
    mod.FINAL_PKN_WITH_URLS = os.path.join(new_dir, 'LemonIte_PKN_with_URLs.tsv')

    # Aliases
    mod.PPI_NETWORK = mod.PPI_OUTPUT_FILE
    mod.FINAL_PKN_FILE = mod.FINAL_PKN_OUTPUT
    mod.FINAL_PKN_WITH_LINKS_FILE = mod.FINAL_PKN_WITH_URLS

    # Logging
    mod.LOG_FILE_API_ERRORS = os.path.join(new_dir, 'api_errors.log')
    mod.LOG_FILE_PIPELINE = os.path.join(new_dir, 'pipeline_progress.log')
