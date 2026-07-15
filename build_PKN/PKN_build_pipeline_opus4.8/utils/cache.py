"""
Persistence / resume layer for the rate-limit-safe pipeline.

Two responsibilities:
  * ``apply_function_with_multithreading`` / ``process_single_database`` —
    chunked, checkpointed processing of an API-backed retriever. Each chunk is
    saved to ``DB_OUTPUT_FILES[db_name]`` immediately, so a job that dies or is
    throttled for hours resumes from the last saved chunk via ``--resume`` instead
    of re-querying everything. Metabolites that error are recorded so a later pass
    retries only them.
  * ``cache_valid`` / ``write_meta`` — a generic ``.meta`` JSON sidecar (timestamp,
    row count) written next to local-file caches, used by retrievers to decide
    whether to recompute.

Re-derived from Collect_PKNdata_metabolites.ipynb cell 3.
"""

import json
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import config
from utils.api_retry import rate_limit_pause

logger = logging.getLogger('pkn.cache')

# Per-database accumulation of interaction-level URL rows, flushed to URL_FILES.
url_data_storage: dict = {}


def is_processed(value):
    """A metabolite is processed if its cell holds anything (incl. 'none'/'error')."""
    return pd.notna(value) and value != '' and value != ' '


def save_url_data(db_name):
    """Flush accumulated interaction-level URL rows for ``db_name`` to its link file."""
    rows = url_data_storage.get(db_name)
    if rows and db_name in config.URL_FILES:
        pd.DataFrame(rows).to_csv(config.URL_FILES[db_name], index=False, sep='\t')
        logger.info("Saved %d URL records for %s", len(rows), db_name)


def apply_function_with_multithreading(df_subset, query_column, db_name,
                                       interaction_function, max_workers=None):
    """
    Run ``interaction_function`` over ``df_subset[query_column]`` with a thread pool.

    Uses the per-database worker count and periodic pause from API_RETRY_CONFIG so
    concurrent requests stay under the API rate limit. Returns {row_index: result}.
    """
    conf = config.API_RETRY_CONFIG.get(db_name, {})
    if max_workers is None:
        max_workers = conf.get('max_workers', config.MAX_WORKERS_DEFAULT)

    results, error_ids, count = {}, [], 0
    logger.info("[%s] %d workers", db_name, max_workers)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(interaction_function, val): idx
                   for idx, val in df_subset[query_column].items()}
        for future in as_completed(futures):
            idx = futures[future]
            try:
                results[idx] = future.result()
            except Exception as e:  # noqa: BLE001
                logger.warning("[%s] error at index %s: %s", db_name, idx, e)
                error_ids.append(idx)
                results[idx] = 'error_occurred'
            count += 1
            rate_limit_pause(count, db_name)

    logger.info("[%s] %d ok, %d errors of %d",
                db_name, len(results) - len(error_ids), len(error_ids), len(results))
    return results


def process_single_database(df_input, db_name, column, function, resume=True):
    """
    Checkpointed processing of one API database, with resume.

    Saves progress to DB_OUTPUT_FILES[db_name] after every chunk and after URL
    accumulation, so an interrupted run continues where it stopped. Returns the
    working dataframe with a ``db_name`` column of pipe-joined gene strings.
    """
    logger.info("=" * 70)
    logger.info("Processing %s (column=%s, resume=%s)", db_name, column, resume)

    url_data_storage[db_name] = []
    output_file = config.DB_OUTPUT_FILES[db_name]

    if resume and os.path.exists(output_file):
        df_existing = pd.read_csv(output_file, sep='\t', low_memory=False,
                                  na_values=['', ' ', 'NA'])
        df_work = df_input[['HMDB', column]].copy()
        if db_name in df_existing.columns:
            df_work = df_work.merge(df_existing[['HMDB', db_name]], on='HMDB', how='left')
        else:
            df_work[db_name] = np.nan
        done = df_work[db_name].apply(is_processed).sum()
        logger.info("Resuming %s: %d done, %d remaining", db_name, done, len(df_work) - done)
    else:
        df_work = df_input[['HMDB', column]].copy()
        df_work[db_name] = np.nan
        logger.info("Starting %s fresh: %d rows", db_name, len(df_work))

    # The result column holds gene strings; use object dtype so per-cell string
    # assignment never collides with an inferred float64 (all-NaN) column.
    df_work[db_name] = df_work[db_name].astype(object)

    total_chunks = (len(df_work) - 1) // config.CHUNK_SIZE + 1

    for chunk_num in range(total_chunks):
        start = chunk_num * config.CHUNK_SIZE
        end = min(start + config.CHUNK_SIZE, len(df_work))
        chunk = df_work.iloc[start:end]

        unprocessed = ~chunk[db_name].apply(is_processed)
        has_id = chunk[column].notna()

        # Metabolites without an identifier can never be queried -> mark 'none'
        for idx in chunk[unprocessed & ~has_id].index:
            df_work.loc[idx, db_name] = 'none'

        to_do = chunk[unprocessed & has_id]
        if to_do.empty:
            continue

        logger.info("[%s] chunk %d/%d: %d metabolites",
                    db_name, chunk_num + 1, total_chunks, len(to_do))
        results = apply_function_with_multithreading(to_do, column, db_name, function)
        for idx, res in results.items():
            df_work.loc[idx, db_name] = res

        save_url_data(db_name)
        df_work.to_csv(output_file, index=False, sep='\t')
        time.sleep(config.PAUSE_BETWEEN_CHUNKS)

    save_url_data(db_name)
    df_work.to_csv(output_file, index=False, sep='\t')
    logger.info("[%s] complete -> %s", db_name, output_file)
    return df_work


# --------------------------------------------------------------------------
# Generic .meta sidecar for local-file caches
# --------------------------------------------------------------------------

def _meta_path(path):
    return path + '.meta'


def write_meta(path, **extra):
    """Write a .meta JSON sidecar (timestamp + row count) next to ``path``."""
    meta = {'created': datetime.now(timezone.utc).isoformat(), **extra}
    try:
        meta.setdefault('rows', sum(1 for _ in open(path)) - 1)
    except OSError:
        pass
    with open(_meta_path(path), 'w') as fh:
        json.dump(meta, fh, indent=2)


def cache_valid(path):
    """True if a cache file and its .meta sidecar both exist."""
    return os.path.exists(path) and os.path.exists(_meta_path(path))
