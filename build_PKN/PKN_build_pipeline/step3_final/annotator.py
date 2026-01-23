"""
PKN annotator module for adding database URLs to interactions.

Handles:
- Adding URLs for metabolite-gene interactions from various databases
- Creating annotated PKN with clickable references
"""

import pandas as pd
import logging
import sys
import os

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config


def generate_intact_url(chebi_id: str, gene: str) -> str:
    """Generate IntAct URL for ChEBI ID and gene."""
    if pd.isna(chebi_id):
        return None
    return f"https://www.ebi.ac.uk/intact/query/{gene}+AND+{chebi_id}"


def generate_uniprot_url(inchikey: str) -> str:
    """Generate UniProt URL for InChIKey search."""
    if pd.isna(inchikey):
        return None
    return f"https://www.uniprot.org/uniprotkb?query={inchikey}"


def generate_chembl_url(chembl_id: str) -> str:
    """Generate chEMBL URL for compound."""
    if pd.isna(chembl_id):
        return None
    return f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}"


def generate_biogrid_url(hmdb_id: str) -> str:
    """Generate BioGRID URL for chemical interactions."""
    if pd.isna(hmdb_id):
        return None
    return f"https://thebiogrid.org/search.php?search={hmdb_id}&organism=all"


def generate_stitch_url(pubchem_id: str) -> str:
    """Generate STITCH URL for compound."""
    if pd.isna(pubchem_id):
        return None
    return f"http://stitch.embl.de/cgi/show_network_section.pl?identifier={pubchem_id}"


def annotate_pkn() -> pd.DataFrame:
    """
    Add database URLs to final PKN.
    
    Returns:
    --------
    pd.DataFrame
        Annotated PKN with URL column
    """
    logger = logging.getLogger('annotator')
    logger.info("="*80)
    logger.info("ANNOTATING PKN WITH DATABASE URLS")
    logger.info("="*80)
    
    # Load final PKN
    logger.info(f"\nLoading final PKN from: {config.FINAL_PKN_FILE}")
    pkn = pd.read_csv(config.FINAL_PKN_FILE, sep='\t')
    
    logger.info(f"  Total interactions: {len(pkn)}")
    
    # Load HMDB metabolites for mapping
    logger.info(f"\nLoading HMDB annotations...")
    
    # For metabolite-gene interactions, add URLs based on source
    pkn['URL'] = None
    
    # Filter for metabolite-gene interactions
    met_gene_mask = pkn['Type'] == 'metabolite-gene'
    
    logger.info(f"\nGenerating URLs for {met_gene_mask.sum()} metabolite-gene interactions...")
    
    # For each source, generate appropriate URL
    # This is a simplified version - full implementation would load all metadata
    
    for idx, row in pkn[met_gene_mask].iterrows():
        hmdb_id = row['Node1']
        gene = row['Node2']
        source = row['Source']
        
        # Generate URL based on source
        if 'IntAct' in source:
            pkn.at[idx, 'URL'] = f"https://www.ebi.ac.uk/intact/query/{gene}"
        elif 'UniProtKB' in source:
            pkn.at[idx, 'URL'] = f"https://www.uniprot.org/uniprotkb?query={gene}+AND+{hmdb_id}"
        elif 'chEMBL' in source:
            pkn.at[idx, 'URL'] = f"https://www.ebi.ac.uk/chembl/g/#browse/activities/filter/target_organism%3A%22Homo%20sapiens%22"
        elif 'BioGRID' in source:
            pkn.at[idx, 'URL'] = generate_biogrid_url(hmdb_id)
        elif 'STITCH' in source:
            pkn.at[idx, 'URL'] = f"http://stitch.embl.de/cgi/network.pl?identifiers={hmdb_id}"
        elif 'MetalinksDB' in source:
            pkn.at[idx, 'URL'] = f"http://metalinks.csb.pitt.edu/"
        elif 'L1000' in source:
            pkn.at[idx, 'URL'] = f"https://clue.io/"
        elif 'GEM' in source:
            pkn.at[idx, 'URL'] = f"https://metabolicatlas.org/"
    
    urls_added = pkn['URL'].notna().sum()
    logger.info(f"  Added {urls_added} URLs")
    
    # Save annotated PKN
    output_file = config.FINAL_PKN_WITH_LINKS_FILE
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    pkn.to_csv(output_file, sep='\t', index=False)
    logger.info(f"  Saved to: {output_file}")
    
    return pkn


if __name__ == '__main__':
    # Test annotator
    import logging
    logging.basicConfig(level=logging.INFO)
    
    annotated = annotate_pkn()
    
    print(f"\nAnnotated PKN:")
    print(f"  Total interactions: {len(annotated)}")
    print(f"  With URLs: {annotated['URL'].notna().sum()}")
    print(annotated[annotated['URL'].notna()].head(10))
