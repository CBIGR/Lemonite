"""
Step-3 annotator: produce the final URL-annotated metabolite-gene table.

The long-format metabolite_gene_PKN.tsv already carries a per-edge ``url`` (attached
in step-1 integration from the per-database link files). This module reshapes it to
the published LemonIte_PKN_with_URLs.tsv (HMDB, Gene, Source, URL) and logs per-source
URL coverage. No URL is regenerated here — all URLs originate in the step-1 retrievers
(via config.URL_TEMPLATES), so the audited/fixed templates are the single source.

Re-derived from Build_final_PKN.ipynb cells 9-23.
"""

import logging

import pandas as pd

import config

logger = logging.getLogger('pkn.annotator')


def annotate_pkn():
    """Write LemonIte_PKN_with_URLs.tsv (HMDB, Gene, Source, URL) and log coverage."""
    pkn = pd.read_csv(config.METABOLITE_GENE_PKN, sep='\t')

    # Metabolite is "Name_HMDB"; extract the HMDB accession (last underscore field).
    out = pkn.copy()
    out['HMDB'] = out['Metabolite'].str.rsplit('_', n=1).str[-1]
    out = out.rename(columns={'url': 'URL'})
    out = out[['HMDB', 'Gene', 'Source', 'URL']]
    out.to_csv(config.FINAL_PKN_WITH_LINKS_FILE, sep='\t', index=False)

    logger.info("Final URL table: %d metabolite-gene edges -> %s",
                len(out), config.FINAL_PKN_WITH_LINKS_FILE)
    logger.info("Per-source URL coverage:")
    for source, sub in out.groupby('Source'):
        n_url = sub['URL'].notna().sum()
        logger.info("  %-18s %7d/%-7d (%.1f%%)",
                    source, n_url, len(sub), 100 * n_url / max(len(sub), 1))
    return out
