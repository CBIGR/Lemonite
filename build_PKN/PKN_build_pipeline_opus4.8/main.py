"""
LemonIte PKN build pipeline — CLI entry point.

Orchestrates the pipeline steps:
  1.  metabolite-gene interactions (9 sources)
  2.  protein-protein interactions (STRING, BioGRID, HuRI, HumanNet)
  2b. enzyme-histone-modification interactions (QuickGO, UniProtKB GO annotations)
  3.  final assembly + URL annotation + URL audit + summary figures

Designed to be re-run monthly (see schedule/run_monthly.sh). All API work is
rate-limit-safe and checkpointed, so ``--resume`` continues an interrupted or
throttled run without re-querying completed databases.

Usage
-----
    python main.py --all                       # full build
    python main.py --step 1                     # one step
    python main.py --step 2b                    # enzyme-histone-mark edges only
    python main.py --step 2c                    # phospho-regulator network (standalone) only
    python main.py --step 1 --databases biogrid,stitch,lincs
    python main.py --step 1 --max-metabolites 50 --output-dir PKN_smoke
    python main.py --all --resume               # resume after interruption
    python main.py --audit-urls                 # re-run only the URL audit

Environment overrides (no file edits needed for HPC/containers):
    PKN_WORKDIR PKN_OUTPUT_DIR_NAME PKN_DB_DIR PKN_GEM_DIR PKN_CONFIG
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

# Allow an alternate config module (e.g. config_hpc) via PKN_CONFIG, resolved
# before any submodule imports `config`.
_config_name = os.environ.get('PKN_CONFIG', 'config')
if _config_name != 'config':
    import importlib
    sys.modules['config'] = importlib.import_module(_config_name)

import config  # noqa: E402

STEP1_DATABASES = ['biogrid', 'stitch', 'uniprot', 'intact', 'chembl',
                   'lincs', 'gem', 'metalinks']


def setup_logging():
    config.ensure_dirs()
    root = logging.getLogger()
    root.handlers.clear()
    logging.basicConfig(
        level=getattr(logging, config.LOG_LEVEL, logging.INFO),
        format=config.LOG_FORMAT,
        handlers=[logging.FileHandler(config.LOG_FILE_PIPELINE),
                  logging.StreamHandler(sys.stdout)],
    )
    # Keep noisy third-party clients quiet
    logging.getLogger('chembl_webresource_client').setLevel(logging.WARNING)
    logging.getLogger('urllib3').setLevel(logging.WARNING)


def run_step1(databases=None, resume=True, max_metabolites=None):
    from step1_metabolites import (preprocessing, integration, biogrid, stitch,
                                    uniprot, intact, chembl, lincs, gem, metalinks)
    log = logging.getLogger('main')
    log.info("=" * 70)
    log.info("STEP 1: metabolite-gene interactions")

    df = preprocessing.load_or_parse_hmdb(max_metabolites=max_metabolites)
    df = preprocessing.preprocess_metabolites(df, resume=resume)

    selected = databases or STEP1_DATABASES
    for name in selected:
        log.info("--- %s ---", name)
        try:
            if name == 'biogrid':
                biogrid.run(df)
            elif name == 'stitch':
                stitch.run(df)
            elif name == 'uniprot':
                uniprot.run(df, resume=resume)
            elif name == 'intact':
                intact.run(df, resume=resume)
            elif name == 'chembl':
                chembl.run(df, resume=resume)
            elif name == 'lincs':
                lincs.run(df)
            elif name == 'gem':
                gem.run(df)
            elif name == 'metalinks':
                metalinks.run(df)
            else:
                log.error("Unknown database: %s", name)
        except Exception:  # noqa: BLE001
            log.exception("%s failed; continuing with remaining databases", name)

    integration.integrate(df)
    log.info("STEP 1 complete -> %s", config.METABOLITE_GENE_PKN)


def run_step2(resume=True):
    from step2_proteins import ppi_integration
    log = logging.getLogger('main')
    log.info("=" * 70)
    log.info("STEP 2: protein-protein interactions")
    if not os.path.exists(config.METABOLITE_GENE_PKN):
        log.error("Missing %s — run step 1 first", config.METABOLITE_GENE_PKN)
        return
    ppi_integration.build_ppi_network(resume=resume)
    log.info("STEP 2 complete -> %s", config.PPI_NETWORK_FILE)


def run_step2b(resume=True):
    from step2b_hPTM import hptm_integration
    log = logging.getLogger('main')
    log.info("=" * 70)
    log.info("STEP 2b: enzyme - histone-modification interactions")
    # Independent of steps 1/2: queries GO terms directly, no gene seed needed.
    hptm_integration.build_hptm_network(resume=resume)
    log.info("STEP 2b complete -> %s", config.HPTM_NETWORK_FILE)


def run_step2c(resume=True):
    from step2c_phospho import phospho_integration
    log = logging.getLogger('main')
    log.info("=" * 70)
    log.info("STEP 2c: phospho-regulator network (kinase-substrate + phospho-TF-target)")
    # Independent standalone layer: OmniPath enzsub (HTTP) + CollecTRI (local). NOT merged
    # into the main PKN by step 3; consumed by the phospho-regulator analysis.
    phospho_integration.build_phospho_network(resume=resume)
    log.info("STEP 2c complete -> %s", config.PHOSPHO_PKN_FILE)


def run_step3(audit=True, run_params=None):
    from step3_final import combiner, annotator, analysis, visualization, provenance
    log = logging.getLogger('main')
    log.info("=" * 70)
    log.info("STEP 3: final assembly + annotation")
    missing = [p for p in (config.METABOLITE_GENE_PKN, config.PPI_NETWORK_FILE)
               if not os.path.exists(p)]
    if missing:
        log.error("Missing inputs: %s — run steps 1 and 2 first", missing)
        return
    combiner.combine_networks()
    annotator.annotate_pkn()
    analysis.analyze()
    visualization.make_figures()
    if audit:
        from utils import url_audit
        url_audit.audit()
    # Record this run in the ledger and diff the PKN against the previous run.
    provenance.record_run(params=run_params or {})
    log.info("STEP 3 complete -> %s", config.FINAL_PKN_WITH_LINKS_FILE)


def main():
    parser = argparse.ArgumentParser(
        description='LemonIte PKN build pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--all', action='store_true', help='Run all three steps')
    parser.add_argument('--step', type=str, choices=['1', '2', '2b', '2c', '3'],
                        help='Run a single step (1, 2, 2b, 2c, or 3)')
    parser.add_argument('--databases', type=str,
                        help='Comma-separated step-1 databases (default: all)')
    parser.add_argument('--resume', action='store_true', default=True,
                        help='Resume from checkpoints (default: on)')
    parser.add_argument('--no-resume', dest='resume', action='store_false',
                        help='Force fresh queries, ignoring caches')
    parser.add_argument('--max-metabolites', type=int, default=None,
                        help='Limit step 1 to the first N metabolites (test runs)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory name (relative to WORKDIR) or absolute path')
    parser.add_argument('--workers', type=int, default=None,
                        help='Override worker threads for every API database')
    parser.add_argument('--audit-urls', action='store_true',
                        help='Run only the provenance-URL audit and exit')
    args = parser.parse_args()

    if args.output_dir:
        config.reconfigure_output_dir(args.output_dir)
    if args.workers is not None:
        config.MAX_WORKERS_DEFAULT = args.workers
        for db in config.API_RETRY_CONFIG:
            config.API_RETRY_CONFIG[db]['max_workers'] = args.workers

    setup_logging()
    log = logging.getLogger('main')
    log.info("PKN pipeline started — output dir: %s", config.OUTPUT_DIR)
    start = time.time()

    databases = [d.strip().lower() for d in args.databases.split(',')] if args.databases else None

    # Parameters recorded in the run ledger for provenance.
    run_params = {
        'mode': 'all' if args.all else (f'step{args.step}' if args.step else 'other'),
        'databases': databases,
        'resume': args.resume,
        'max_metabolites': args.max_metabolites,
        'workers': args.workers,
    }

    if args.audit_urls:
        from utils import url_audit
        url_audit.audit()
    elif args.all:
        run_step1(databases, args.resume, args.max_metabolites)
        run_step2(args.resume)
        run_step2b(args.resume)
        run_step2c(args.resume)
        run_step3(run_params=run_params)
    elif args.step == '1':
        run_step1(databases, args.resume, args.max_metabolites)
    elif args.step == '2':
        run_step2(args.resume)
    elif args.step == '2b':
        run_step2b(args.resume)
    elif args.step == '2c':
        run_step2c(args.resume)
    elif args.step == '3':
        run_step3(run_params=run_params)
    else:
        parser.print_help()
        return

    log.info("PKN pipeline finished in %.1fs", time.time() - start)


if __name__ == '__main__':
    main()
