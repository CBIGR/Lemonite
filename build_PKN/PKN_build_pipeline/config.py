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
WORKDIR = '/home/borisvdm/Documents/PhD/Lemonite'
OUTPUT_DIR_NAME = 'PKN'  # Change this to use a different output folder
# ========================================================

# Set up paths
OUTPUT_DIR = os.path.join(WORKDIR, OUTPUT_DIR_NAME)
DB_DIR = '/home/borisvdm/Documents/PhD/resources/databases'

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================================================================
# STEP 1: METABOLITE-GENE INTERACTION SOURCES
# =========================================================================

# HMDB metabolites XML database
HMDB_METABOLITES_XML = '/home/borisvdm/Documents/PhD/resources/databases/HMDB/hmdb_metabolites.xml'

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

# MetalinksDB (from liana+ package)
METALINKS_PATH = os.path.join(DB_DIR, 'metalinks/metalinks.csv')

# LINCS biochemical binding data
LINCS_COMPOUND_MAPPING = os.path.join(DB_DIR, 'LINCS/lsp_compound_mapping.csv')
LINCS_TARGET_MAPPING = os.path.join(DB_DIR, 'LINCS/lsp_target_mapping.csv')
LINCS_TARGET_DICTIONARY = os.path.join(DB_DIR, 'LINCS/lsp_target_dictionary.csv')
LINCS_BIOCHEM_AGG = os.path.join(DB_DIR, 'LINCS/lsp_biochem_agg.csv')
LINCS_REFERENCES = os.path.join(DB_DIR, 'LINCS/lsp_references.csv')

# L1000 gene expression signatures (2.1GB GMT file)
L1000_GMT_LOCATION = os.path.join(DB_DIR, 'LINCS/L1000/l1000_cp.gmt')
L1000_COMPOUNDS_LOCATION = os.path.join(DB_DIR, 'LINCS/L1000/LINCS_small_molecules.tsv')

# Pathway distance for GEM (how many reactions away to include)
PATHWAY_DISTANCE = 2

# =========================================================================
# STEP 2: PROTEIN-PROTEIN INTERACTION SOURCES
# =========================================================================

# BioGRID PPI database
BIOGRID_PPI_LOCATION = os.path.join(DB_DIR, 'BioGRID/BIOGRID-ALL-4.4.238.tab3.txt')

# HuRI human protein-protein interactome
HURI_LOCATION = os.path.join(DB_DIR, 'HuRI/HuRI_Tong_2021.tsv')

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
    'L1000': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_L1000_processed.csv'),
    'Human1_GEM_dist1': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_Human1_GEM_dist1_processed.csv'),
    'Human1_GEM_dist2': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_Human1_GEM_dist2_processed.csv'),
    'MetalinksDB': os.path.join(OUTPUT_DIR, 'HMDB_metabolites_MetalinksDB_processed.csv')
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
PPI_HURI_CACHE = os.path.join(OUTPUT_DIR, 'cache/huri_ppi_cache.csv')

# =========================================================================
# OUTPUT FILES (STEP 3: FINAL PKN)
# =========================================================================

METABOLITE_GENE_PKN = os.path.join(OUTPUT_DIR, 'metabolite_gene_PKN.tsv')
FINAL_PKN_OUTPUT = os.path.join(OUTPUT_DIR, 'LemonIte_PKN.tsv')
FINAL_PKN_WITH_URLS = os.path.join(OUTPUT_DIR, 'LemonIte_PKN_with_URLs.tsv')

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
    }
}

# =========================================================================
# PROCESSING PARAMETERS
# =========================================================================

CHUNK_SIZE = 800  # For batching metabolite processing
MAX_WORKERS_DEFAULT = 4  # Default thread pool size
RESUME_SAVE_INTERVAL = 1000  # Save progress every N metabolites (for L1000)

# =========================================================================
# LOGGING CONFIGURATION
# =========================================================================

LOG_FILE_API_ERRORS = os.path.join(OUTPUT_DIR, 'api_errors.log')
LOG_FILE_PIPELINE = os.path.join(OUTPUT_DIR, 'pipeline_progress.log')
LOG_LEVEL = 'INFO'
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
