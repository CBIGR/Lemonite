"""
HPC config overrides (e.g. VSC Ghent / Doduo).

Activated by setting ``PKN_CONFIG=config_hpc`` before running main.py. It imports
everything from ``config`` and then overrides only what differs on the cluster.
In practice all paths are already env-driven (PKN_WORKDIR / PKN_DB_DIR /
PKN_GEM_DIR), so this module mostly exists to document the HPC contract and to be
a place to pin cluster-specific values (e.g. a STRING mirror) without editing the
base config.
"""

import os
from config import *  # noqa: F401,F403  (re-export the full base config)

# All real paths come from environment variables set by submit_pkn_hpc.sh:
#   PKN_WORKDIR, PKN_OUTPUT_DIR_NAME, PKN_DB_DIR, PKN_GEM_DIR
# config.py already reads those, so nothing path-related needs repeating here.

# Example cluster-specific override (uncomment / adjust if a mirror is required):
# STRING_API_URL = os.environ.get('PKN_STRING_API_URL', 'https://version-11-5.string-db.org/api')
