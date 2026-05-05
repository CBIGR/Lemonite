"""
Main entry point for the PKN (Prior Knowledge Network) pipeline.

This script orchestrates the three-step process:
1. Collect metabolite-gene interactions from multiple databases
2. Build protein-protein interaction network
3. Integrate and analyze final PKN

Usage:
------
# Run all steps
python main.py --all

# Run individual steps
python main.py --step 1
python main.py --step 2
python main.py --step 3

# Run with specific databases only
python main.py --step 1 --databases biogrid,stitch,lincs

# Resume from checkpoint
python main.py --step 1 --resume
"""

import argparse
import logging
import sys
import os
import time
from pathlib import Path

# Add current directory to Python path
sys.path.append(str(Path(__file__).parent))

import config
from utils.file_io import load_hmdb_metabolites


def setup_logging():
    """Configure logging for the pipeline."""
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)
    
    # Remove existing handlers to avoid duplicates on re-run
    root = logging.getLogger()
    root.handlers.clear()
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(config.LOG_FILE_PIPELINE),
            logging.StreamHandler(sys.stdout)
        ]
    )


def run_step1_metabolites(databases=None, resume=False, max_metabolites=None):
    """
    Run Step 1: Collect metabolite-gene interactions.
    
    Parameters:
    -----------
    databases : list, optional
        List of database names to run (default: all)
    resume : bool
        Whether to resume from checkpoint
    max_metabolites : int, optional
        Cap the number of metabolites processed (for test runs)
    """
    from step1_metabolites import preprocessing, integration
    from step1_metabolites.biogrid import BioGRIDRetriever
    from step1_metabolites.stitch import STITCHRetriever
    from step1_metabolites.uniprot import UniProtRetriever
    from step1_metabolites.intact import IntActRetriever
    from step1_metabolites.chembl import ChEMBLRetriever
    from step1_metabolites.lincs import LINCSRetriever
    from step1_metabolites.gem import GEMRetriever
    from step1_metabolites.metalinks import MetalinksRetriever
    
    logger = logging.getLogger('main')
    logger.info("="*80)
    logger.info("STEP 1: COLLECTING METABOLITE-GENE INTERACTIONS")
    logger.info("="*80)
    
    step_start = time.time()
    
    # Load metabolites from HMDB
    logger.info("Loading HMDB metabolites...")
    if not os.path.exists(config.HMDB_METABOLITES_XML):
        logger.error(f"HMDB metabolites file not found: {config.HMDB_METABOLITES_XML}")
        return None
    metabolites = load_hmdb_metabolites(config.HMDB_METABOLITES_XML)
    logger.info(f"Loaded {len(metabolites)} metabolites")

    if max_metabolites is not None:
        metabolites = metabolites[:max_metabolites]
        logger.info(f"Limiting to first {max_metabolites} metabolites for this run")

    # Preprocess metabolites (add ChEMBL IDs, etc.)
    logger.info("\nPreprocessing metabolites...")
    metabolites_df = preprocessing.preprocess_metabolites(metabolites)
    
    # Define available databases and their retriever classes
    retriever_classes = {
        'biogrid': BioGRIDRetriever,
        'stitch': STITCHRetriever,
        'uniprot': UniProtRetriever,
        'intact': IntActRetriever,
        'chembl': ChEMBLRetriever,
        'lincs': LINCSRetriever,
        'gem_dist1': lambda: GEMRetriever(distance=1),
        'gem_dist2': lambda: GEMRetriever(distance=2),
        'metalinks': MetalinksRetriever
    }
    
    # Use specified databases or all
    databases_to_run = databases if databases else list(retriever_classes.keys())
    
    logger.info(f"\nRunning databases: {', '.join(databases_to_run)}")
    
    # Run retrievers
    results = {}
    
    for db_name in databases_to_run:
        try:
            logger.info(f"\n{'='*80}")
            logger.info(f"Processing {db_name.upper()}")
            logger.info(f"{'='*80}")
            
            # Initialize retriever
            retriever_factory = retriever_classes.get(db_name)
            if retriever_factory is None:
                logger.error(f"Unknown database: {db_name}")
                continue
            
            # Create retriever instance
            if callable(retriever_factory) and not isinstance(retriever_factory, type):
                retriever = retriever_factory()  # Factory function (for GEM)
            else:
                retriever = retriever_factory()  # Regular class
            
            # Get interactions
            interactions = retriever.get_interactions(metabolites_df.to_dict('records'))
            results[db_name] = interactions
            
            if interactions is not None and len(interactions) > 0:
                logger.info(f"✓ {db_name}: {len(interactions)} interactions")
            else:
                logger.warning(f"✗ {db_name}: No interactions found")
        
        except Exception as e:
            logger.error(f"✗ {db_name}: Failed with error: {e}", exc_info=True)
            results[db_name] = None
    
    # Integrate results
    logger.info(f"\n{'='*80}")
    logger.info("INTEGRATING RESULTS")
    logger.info(f"{'='*80}")
    
    final_network = integration.integrate_databases(results)
    integration.create_visualizations(results)
    
    logger.info(f"\n{'='*80}")
    logger.info("STEP 1 COMPLETE!")
    logger.info(f"{'='*80}")
    logger.info(f"Final network: {len(final_network)} interactions")
    logger.info(f"Saved to: {config.METABOLITE_GENE_PKN}")
    logger.info(f"Step 1 elapsed time: {time.time() - step_start:.1f}s")
    
    return final_network


def run_step2_proteins():
    """Run Step 2: Build protein-protein interaction network."""
    from step2_proteins import ppi_integration
    
    logger = logging.getLogger('main')
    logger.info("\n"+"="*80)
    logger.info("STEP 2: BUILDING PROTEIN-PROTEIN INTERACTION NETWORK")
    logger.info("="*80)
    
    # Validate that step 1 output exists
    if not os.path.exists(config.METABOLITE_GENE_PKN):
        logger.error(
            f"Metabolite-gene PKN not found: {config.METABOLITE_GENE_PKN}\n"
            f"Run step 1 first: python main.py --step 1"
        )
        return None
    
    step_start = time.time()
    
    # Run PPI collection
    ppi_network = ppi_integration.build_ppi_network()
    
    logger.info(f"\nPPI network:")
    logger.info(f"  Total interactions: {len(ppi_network)}")
    logger.info(f"  Unique proteins: {ppi_network[['Node1', 'Node2']].stack().nunique()}")
    logger.info(f"  Output: {config.PPI_OUTPUT_FILE}")
    logger.info(f"Step 2 elapsed time: {time.time() - step_start:.1f}s")
    
    logger.info("\nSTEP 2 COMPLETE")


def run_step3_final():
    """Run Step 3: Integrate and analyze final PKN."""
    from step3_final import combiner, annotator, analysis, visualization
    
    logger = logging.getLogger('main')
    logger.info("\n"+"="*80)
    logger.info("STEP 3: FINAL PKN INTEGRATION AND ANALYSIS")
    logger.info("="*80)
    
    # Validate that step 1 and step 2 outputs exist
    missing = []
    if not os.path.exists(config.METABOLITE_GENE_PKN):
        missing.append(f"Metabolite-gene PKN: {config.METABOLITE_GENE_PKN}")
    if not os.path.exists(config.PPI_NETWORK):
        missing.append(f"PPI network: {config.PPI_NETWORK}")
    if missing:
        logger.error(
            f"Required input files not found:\n  " + "\n  ".join(missing) +
            "\nRun steps 1 and 2 first."
        )
        return None
    
    step_start = time.time()
    
    # Combine networks
    final_pkn = combiner.combine_networks()
    
    # Add URL annotations
    logger.info("\n" + "="*80)
    logger.info("ADDING ANNOTATIONS")
    logger.info("="*80)
    annotated_pkn = annotator.annotate_pkn()
    
    # Analyze coverage
    logger.info("\n" + "="*80)
    logger.info("ANALYZING COVERAGE")
    logger.info("="*80)
    stats = analysis.analyze_coverage()
    
    # Generate visualizations
    logger.info("\n" + "="*80)
    logger.info("GENERATING VISUALIZATIONS")
    logger.info("="*80)
    visualization.create_all_visualizations()
    
    logger.info(f"\n{'='*80}")
    logger.info("STEP 3 COMPLETE!")
    logger.info(f"{'='*80}")
    logger.info(f"Final PKN: {len(final_pkn)} interactions")
    logger.info(f"Saved to: {config.FINAL_PKN_FILE}")
    logger.info(f"With URLs: {config.FINAL_PKN_WITH_LINKS_FILE}")
    logger.info(f"Step 3 elapsed time: {time.time() - step_start:.1f}s")
    
    return final_pkn


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='PKN Pipeline: Build Prior Knowledge Network from metabolites',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--step',
        type=int,
        choices=[1, 2, 3],
        help='Run specific step (1=metabolites, 2=PPI, 3=final)'
    )
    
    parser.add_argument(
        '--all',
        action='store_true',
        help='Run all three steps sequentially'
    )
    
    parser.add_argument(
        '--databases',
        type=str,
        help='Comma-separated list of databases for step 1 (e.g., biogrid,stitch,lincs)'
    )
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from checkpoint (for long-running API tasks)'
    )

    parser.add_argument(
        '--workers',
        type=int,
        default=None,
        metavar='N',
        help='Number of concurrent worker threads (overrides per-database defaults)'
    )

    parser.add_argument(
        '--max-metabolites',
        type=int,
        default=None,
        metavar='N',
        help='Limit step 1 to the first N metabolites (useful for test runs)'
    )

    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        metavar='DIR',
        help='Output directory name (relative to WORKDIR) or absolute path'
    )
    
    args = parser.parse_args()

    # Apply output directory override BEFORE setup_logging so logs land in the right place
    if args.output_dir:
        config.reconfigure_output_dir(args.output_dir)

    # Apply workers override to every database entry and the global default
    if args.workers is not None:
        config.MAX_WORKERS_DEFAULT = args.workers
        for db in config.API_RETRY_CONFIG:
            config.API_RETRY_CONFIG[db]['max_workers'] = args.workers

    # Setup logging
    setup_logging()
    logger = logging.getLogger('main')
    
    logger.info("="*80)
    logger.info("PKN PIPELINE STARTED")
    logger.info("="*80)
    logger.info(f"Output directory: {config.OUTPUT_DIR}")
    
    # Parse databases list
    databases = None
    if args.databases:
        databases = [db.strip().lower() for db in args.databases.split(',')]
    
    # Run requested steps
    if args.all:
        run_step1_metabolites(databases=databases, resume=args.resume, max_metabolites=args.max_metabolites)
        run_step2_proteins()
        run_step3_final()
    elif args.step == 1:
        run_step1_metabolites(databases=databases, resume=args.resume, max_metabolites=args.max_metabolites)
    elif args.step == 2:
        run_step2_proteins()
    elif args.step == 3:
        run_step3_final()
    else:
        parser.print_help()
        return
    
    logger.info("\n"+"="*80)
    logger.info("PKN PIPELINE COMPLETE")
    logger.info("="*80)


if __name__ == '__main__':
    main()
