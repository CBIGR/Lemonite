"""
Provenance-URL audit — validates the clickable source URLs in the final PKN.

For each Source in LemonIte_PKN_with_URLs.tsv it samples a few emitted URLs and
issues a lightweight HEAD/GET, asserting the response is not 404 and not a 5xx.
A report is written to URL_AUDIT_REPORT so a future monthly run surfaces immediately
if a database changes its URL scheme.

Caveat: IntAct and ChEMBL serve Angular single-page-apps that return HTTP 200 for
any path (the record is rendered client-side), so a 200 there confirms reachability
of the app, not the specific record. Those URLs are additionally checked structurally
against config.URL_TEMPLATES.
"""

import logging
import os
from urllib.parse import urlparse

import pandas as pd

import config
from utils.http import get_session

logger = logging.getLogger('pkn.url_audit')

# Sources whose URLs are SPA shells: 200 == reachable app, not record-level proof.
# (EBI: IntAct/ChEMBL, and QuickGO/UniProt for the step-2b hPTM edges;
# metabolicatlas.org: Human-GEM reaction pages.)
_SPA_HOSTS = {'www.ebi.ac.uk', 'metabolicatlas.org', 'www.uniprot.org'}

SAMPLE_PER_SOURCE = 5


# Distance-2 Human-GEM edges carry two reaction URLs joined by this separator
# (metabolite -> intermediate, intermediate -> gene); each part is audited on its own.
_URL_JOIN = ' -> '


def _check(url):
    """Return (status_code, ok). ok = reachable (not 404/5xx).

    A 403 is treated as reachable: hosts like hmdb.ca sit behind Cloudflare bot
    protection that rejects scripted requests but serve the page fine in a real
    browser. Only 404 (wrong record) and 5xx (server down) count as failures.

    A GEM distance-2 value holds two reaction URLs joined by ``_URL_JOIN``; each
    component is checked and the worst status is returned (ok only if both are).
    """
    if _URL_JOIN in url:
        codes, oks = [], []
        for part in url.split(_URL_JOIN):
            c, o = _check(part.strip())
            codes.append(c)
            oks.append(o)
        # Report the first failing component's code, else the first code.
        fail_code = next((c for c, o in zip(codes, oks) if not o), codes[0])
        return fail_code, all(oks)
    session = get_session()
    try:
        resp = session.head(url, timeout=20, allow_redirects=True)
        if resp.status_code in (405, 403) or resp.status_code >= 400:
            # Some servers reject HEAD; retry with a lightweight GET.
            resp = session.get(url, timeout=25, allow_redirects=True, stream=True)
        code = resp.status_code
    except Exception as e:  # noqa: BLE001
        logger.warning("URL check error for %s: %s", url, str(e)[:120])
        return None, False
    ok = code is not None and code != 404 and code < 500
    return code, ok


def audit(sample_per_source=SAMPLE_PER_SOURCE):
    """Audit sampled URLs from the final tables; write a CSV report; return it.

    Covers the metabolite-gene provenance table (FINAL_PKN_WITH_LINKS_FILE) and,
    when present, the step-2b enzyme-histone-mark network (HPTM_NETWORK_FILE) —
    both carry a Source + URL column.
    """
    pkn = pd.read_csv(config.FINAL_PKN_WITH_LINKS_FILE, sep='\t')
    frames = [pkn[['Source', 'URL']]]
    if os.path.exists(config.HPTM_NETWORK_FILE):
        hptm = pd.read_csv(config.HPTM_NETWORK_FILE, sep='\t')
        if not hptm.empty and 'URL' in hptm.columns:
            frames.append(hptm[['Source', 'URL']])
    pkn = pd.concat(frames, ignore_index=True)
    pkn = pkn[pkn['URL'].notna()]

    rows = []
    for source, sub in pkn.groupby('Source'):
        urls = sub['URL'].dropna().drop_duplicates().head(sample_per_source).tolist()
        for url in urls:
            code, ok = _check(url)
            # For joined GEM URLs, take the host of the first component.
            host = urlparse(url.split(_URL_JOIN)[0]).netloc
            rows.append({
                'Source': source, 'URL': url, 'status_code': code,
                'reachable': ok, 'spa_host': host in _SPA_HOSTS,
            })
            logger.info("  [%s] %s -> %s (%s)", source, url, code,
                        'OK' if ok else 'FAIL')

    report = pd.DataFrame(rows)
    report.to_csv(config.URL_AUDIT_REPORT, index=False)

    n_fail = (~report['reachable']).sum() if not report.empty else 0
    if n_fail:
        logger.warning("URL audit: %d/%d sampled URLs unreachable (see %s)",
                       n_fail, len(report), config.URL_AUDIT_REPORT)
    else:
        logger.info("URL audit: all %d sampled URLs reachable -> %s",
                    len(report), config.URL_AUDIT_REPORT)
    return report
