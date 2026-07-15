"""
Run provenance and PKN change tracking.

Called at the end of every step-3 build. It:
  1. records one JSON line per run in an append-only ledger (PKN_run_history.jsonl):
     timestamp, run id, git commit, parameters, per-source edge counts, node/edge
     totals, and a content checksum of the network;
  2. diffs the current network against the previous run's snapshot and writes the
     full set of added / removed edges to changes/<run_id>_diff.tsv.gz, with per-source
     counts summarised in the ledger record;
  3. saves the current edge set as the snapshot for the next run;
  4. appends a human-readable entry to PKN_CHANGELOG.md (repo-versioned).

Together these give a permanent, auditable record of every run and of how the PKN
changed between runs.
"""

import hashlib
import json
import logging
import os
import subprocess
from datetime import datetime, timezone

import pandas as pd

import config

logger = logging.getLogger('pkn.provenance')

_EDGE_COLS = ['Node1', 'Node2', 'Source', 'Type']


def _git_commit():
    """Best-effort short git commit of the pipeline code (or 'unknown')."""
    try:
        out = subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=os.path.dirname(__file__), capture_output=True, text=True, timeout=10)
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:  # noqa: BLE001
        pass
    return 'unknown'


def _checksum(df):
    """Order-independent content checksum of the edge set."""
    keyed = df[['Node1', 'Node2', 'Source', 'Type']].astype(str)
    joined = keyed['Node1'] + '\t' + keyed['Node2'] + '\t' + keyed['Source'] + '\t' + keyed['Type']
    h = hashlib.sha256()
    for line in sorted(joined.tolist()):
        h.update(line.encode('utf-8'))
    return h.hexdigest()


def _source_counts(df):
    """Per-source edge counts as a plain dict (metabolite-gene + PPI sources)."""
    return {str(k): int(v) for k, v in df['Source'].value_counts().items()}


def _edge_keys(df):
    """Set of (Node1, Node2, Source) tuples for diffing."""
    return set(map(tuple, df[['Node1', 'Node2', 'Source']].itertuples(index=False)))


def record_run(params=None):
    """Record this run in the ledger and diff the PKN against the previous run.

    Returns the run record dict. Never raises: provenance is best-effort and must
    not fail a completed build.
    """
    try:
        return _record_run(params or {})
    except Exception:  # noqa: BLE001
        logger.exception("Provenance recording failed (build itself is unaffected)")
        return None


def _record_run(params):
    os.makedirs(config.CHANGES_DIR, exist_ok=True)

    pkn = pd.read_csv(config.FINAL_PKN_FILE, sep='\t')
    now = datetime.now(timezone.utc)
    run_id = now.strftime('%Y%m%dT%H%M%SZ')
    # Guard against a collision if two runs land in the same second.
    if os.path.exists(os.path.join(config.CHANGES_DIR, f'{run_id}_diff.tsv.gz')):
        run_id = f"{run_id}-{now.strftime('%f')[:3]}"

    mg = pkn[pkn['Type'] == 'metabolite-gene']
    ppi = pkn[pkn['Type'] == 'PPI']
    hptm = pkn[pkn['Type'] == 'histone-modification']
    metabolites = set(mg['Node1'].unique())
    # hPTM Node1 is the enzyme gene; Node2 is the histone mark (a distinct node type).
    genes = set(mg['Node2']) | set(ppi['Node1']) | set(ppi['Node2']) | set(hptm['Node1'])
    marks = set(hptm['Node2'].unique())

    record = {
        'run_id': run_id,
        'timestamp': now.isoformat(),
        'git_commit': _git_commit(),
        'output_dir': config.OUTPUT_DIR,
        'params': params,
        'checksum_sha256': _checksum(pkn),
        'totals': {
            'unique_metabolites': len(metabolites),
            'unique_genes': len(genes),
            'unique_histone_marks': len(marks),
            'total_nodes': len(metabolites) + len(genes) + len(marks),
            'metabolite_gene_edges': int(len(mg)),
            'ppi_edges': int(len(ppi)),
            'histone_modification_edges': int(len(hptm)),
            'total_edges': int(len(pkn)),
        },
        'source_counts': _source_counts(pkn),
    }

    # Diff against the previous snapshot, if any.
    change = _diff_against_snapshot(pkn, run_id)
    record['change'] = change

    # Append to the JSON-lines ledger.
    with open(config.RUN_HISTORY_FILE, 'a') as fh:
        fh.write(json.dumps(record) + '\n')

    # Save the current edge set as the snapshot for the next run.
    pkn[_EDGE_COLS].to_csv(config.PKN_SNAPSHOT_FILE, sep='\t', index=False,
                           compression='gzip')

    _append_changelog(record)

    logger.info("Run %s recorded -> %s", run_id, config.RUN_HISTORY_FILE)
    if change['is_first_run']:
        logger.info("  first recorded run: %d edges (baseline)", record['totals']['total_edges'])
    else:
        logger.info("  changes vs previous run: +%d / -%d edges (net %+d)",
                    change['added'], change['removed'], change['net'])
    return record


def _diff_against_snapshot(pkn, run_id):
    """Compare the current PKN to the previous snapshot; write the full edge diff."""
    if not os.path.exists(config.PKN_SNAPSHOT_FILE):
        return {'is_first_run': True, 'added': 0, 'removed': 0, 'net': 0,
                'added_by_source': {}, 'removed_by_source': {}, 'diff_file': None}

    prev = pd.read_csv(config.PKN_SNAPSHOT_FILE, sep='\t', compression='gzip')
    cur_keys = _edge_keys(pkn)
    prev_keys = _edge_keys(prev)

    added = cur_keys - prev_keys
    removed = prev_keys - cur_keys

    diff_rows = ([{'change': 'added', 'Node1': a, 'Node2': b, 'Source': s} for a, b, s in added]
                 + [{'change': 'removed', 'Node1': a, 'Node2': b, 'Source': s} for a, b, s in removed])
    diff_file = os.path.join(config.CHANGES_DIR, f'{run_id}_diff.tsv.gz')
    if diff_rows:
        pd.DataFrame(diff_rows).to_csv(diff_file, sep='\t', index=False, compression='gzip')
    else:
        diff_file = None  # no change -> no diff file

    def _by_source(keys):
        counts = {}
        for _, _, s in keys:
            counts[str(s)] = counts.get(str(s), 0) + 1
        return counts

    return {
        'is_first_run': False,
        'added': len(added),
        'removed': len(removed),
        'net': len(added) - len(removed),
        'added_by_source': _by_source(added),
        'removed_by_source': _by_source(removed),
        'diff_file': os.path.basename(diff_file) if diff_file else None,
    }


def _append_changelog(record):
    """Append a concise, human-readable entry to the repo-versioned changelog."""
    t = record['totals']
    c = record['change']
    new = not os.path.exists(config.PKN_CHANGELOG_FILE)
    with open(config.PKN_CHANGELOG_FILE, 'a') as fh:
        if new:
            fh.write("# LemonIte PKN change log\n\n"
                     "Automatically appended at the end of every build (step 3). "
                     "Each entry records the run and how the network changed versus "
                     "the previous run. Full per-run detail is in the output "
                     "directory's `PKN_run_history.jsonl` and `changes/` diffs.\n\n")
        fh.write(f"## {record['run_id']}  (git {record['git_commit']})\n\n")
        fh.write(f"- Timestamp: {record['timestamp']}\n")
        fh.write(f"- Output dir: {record['output_dir']}\n")
        if record['params']:
            fh.write(f"- Params: {record['params']}\n")
        fh.write(f"- Totals: {t['total_edges']:,} edges "
                 f"({t['metabolite_gene_edges']:,} metabolite-gene, {t['ppi_edges']:,} PPI, "
                 f"{t.get('histone_modification_edges', 0):,} histone-modification); "
                 f"{t['unique_metabolites']:,} metabolites, {t['unique_genes']:,} genes, "
                 f"{t.get('unique_histone_marks', 0):,} histone marks\n")
        if c['is_first_run']:
            fh.write("- Change: first recorded run (baseline)\n")
        else:
            fh.write(f"- Change vs previous run: +{c['added']:,} / -{c['removed']:,} edges "
                     f"(net {c['net']:+,})\n")
            if c['added_by_source'] or c['removed_by_source']:
                srcs = sorted(set(c['added_by_source']) | set(c['removed_by_source']))
                parts = [f"{s} +{c['added_by_source'].get(s, 0)}/-{c['removed_by_source'].get(s, 0)}"
                         for s in srcs]
                fh.write(f"  - By source: {', '.join(parts)}\n")
            if c['diff_file']:
                fh.write(f"  - Full edge diff: changes/{c['diff_file']}\n")
        fh.write(f"- Checksum (sha256): {record['checksum_sha256']}\n\n")
