"""
UniProtKB metabolite-protein interactions (REST API, ChEBI ligand query).

For each metabolite, queries the UniProtKB search API for human entries that
reference the metabolite's ChEBI identifier (as a bound ligand, cofactor, or
substrate), accumulating per-interaction link rows with the UniProt entry URL.
Runs through the checkpointed process_single_database so it resumes safely after
rate limiting.

NOTE ON THE QUERY FIELD: the reference notebook queried the UniProtKB ``inchikey:``
field. That field has since stopped returning results from the UniProt REST API
(even ATP yields zero hits), so every InChIKey query produced an empty result. The
current, working mechanism is the ``chebi:`` cross-reference field, which matches the
metabolite to UniProt entries whose binding-site / cofactor annotations reference that
ChEBI id. This was validated to recover the expected genes (e.g. ChEBI:1189 ->
COMT, UGT1A8, matching the reference output).

A per-metabolite cap is applied: a few ubiquitous cofactors (ATP, NAD+, Mg2+, ...)
are referenced by thousands of proteins, which are low-specificity associations that
would flood the network. The cap keeps the retrieval in line with the reference
distribution (whose maximum was ~100 proteins per metabolite).
"""

import logging

import pandas as pd
from requests.exceptions import RequestException

import config
from utils import cache
from utils.api_retry import retry_api_call
from utils.http import get_session

logger = logging.getLogger('pkn.uniprot')

# Maximum proteins to accept per metabolite. Cofactors such as ATP/NAD+ are bound by
# thousands of proteins; those are non-specific and would dominate the network. The
# reference network's maximum was ~101 proteins per metabolite, so 200 is comfortably
# above legitimate cases while excluding the runaway cofactors.
MAX_PROTEINS_PER_METABOLITE = 200
_PAGE_SIZE = 100

# ChEBI ontology relationships that connect a metabolite's (neutral) ChEBI id to the
# charged / tautomeric species that UniProt actually annotates on its binding sites.
_CHEBI_CONJUGATE_RELS = {
    'is tautomer of', 'is protonated form of', 'is deprotonated form of',
    'is conjugate acid of', 'is conjugate base of',
}
_OLS_GRAPH = ('https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms/'
              'http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCHEBI_{chebi}/graph')

# Process-wide cache: ChEBI id -> sorted list of [self + related charged species].
_chebi_expansion_cache: dict = {}


def _expand_chebi(chebi_id):
    """Return [chebi_id] plus its conjugate-acid/base/tautomer ChEBI ids.

    UniProt annotates ligands under specific charge/tautomer states (e.g. biotin(1-)
    CHEBI:57586, not neutral biotin CHEBI:15956), so a metabolite's neutral HMDB
    ChEBI alone misses many associations. The related species are resolved once per
    ChEBI from the ChEBI ontology (via OLS) and cached.
    """
    if chebi_id in _chebi_expansion_cache:
        return _chebi_expansion_cache[chebi_id]

    import re
    related = {chebi_id}
    try:
        resp = get_session().get(_OLS_GRAPH.format(chebi=chebi_id), timeout=20)
        if resp.status_code == 200:
            for edge in resp.json().get('edges', []):
                if edge.get('label') in _CHEBI_CONJUGATE_RELS:
                    for node in (edge.get('source'), edge.get('target')):
                        m = re.search(r'CHEBI_(\d+)', str(node))
                        if m:
                            related.add(int(m.group(1)))
    except Exception:  # noqa: BLE001 - expansion is best-effort; fall back to self
        pass

    result = sorted(related)
    _chebi_expansion_cache[chebi_id] = result
    return result


@retry_api_call(db_name='UniProtKB')
def get_uniprot_interactions(chebi_id, inchikey=None, reviewed=False):
    """Return pipe-joined gene symbols for a metabolite ChEBI id.

    Queries UniProtKB for the metabolite's ChEBI id and its related charged/tautomeric
    ChEBI species (see _expand_chebi). Stores per-interaction link rows in the shared
    URL accumulator, keyed by InChIKey so step-3 integration can remap them to HMDB
    ids (the link-file schema is unchanged from the InChIKey-based version).
    """
    if pd.isna(chebi_id):
        return 'none'
    try:
        chebi_id = int(float(chebi_id))
    except (ValueError, TypeError):
        return 'none'

    timeout = config.API_RETRY_CONFIG['UniProtKB']['timeout']
    chebi_ids = _expand_chebi(chebi_id)
    chebi_clause = ' OR '.join(f'chebi:{c}' for c in chebi_ids)
    query = f"({chebi_clause}) AND (organism_id:9606)"
    if reviewed:
        query += " AND (reviewed:true)"

    store = cache.url_data_storage.setdefault('UniProtKB', [])
    interactions = []
    seen_genes = set()
    session = get_session()

    # Paginate via the cursor in the Link header until the cap is reached.
    url = (f"https://rest.uniprot.org/uniprotkb/search"
           f"?query={query}&format=tsv&fields=accession,gene_names&size={_PAGE_SIZE}")
    while url and len(interactions) < MAX_PROTEINS_PER_METABOLITE:
        response = session.get(url, timeout=timeout)
        if response.status_code != 200:
            raise RequestException(f"HTTP {response.status_code}")

        lines = response.text.strip().split('\n')
        if len(lines) >= 2:
            for line in lines[1:]:
                cols = line.split('\t')
                if len(cols) < 2:
                    continue
                uni_id = cols[0].strip()
                genes = cols[1].strip().split(' ') if cols[1].strip() else []
                for gene in (g.strip() for g in genes):
                    if not gene or gene in seen_genes:
                        continue
                    seen_genes.add(gene)
                    interactions.append(gene)
                    store.append({
                        'InChIKey': inchikey, 'ChEBI_ID': chebi_id, 'Gene': gene,
                        'UniProt_ID': uni_id,
                        'URL': config.URL_TEMPLATES['UniProtKB'].format(accession=uni_id),
                    })
                    if len(interactions) >= MAX_PROTEINS_PER_METABOLITE:
                        break
                if len(interactions) >= MAX_PROTEINS_PER_METABOLITE:
                    break

        # Follow the rel="next" cursor link for the next page, if any.
        url = response.links.get('next', {}).get('url')

    return '|'.join(interactions) if interactions else 'none'


def run(df, resume=True):
    """Process UniProtKB for all metabolites (keyed by ChEBI) with checkpointing.

    The query column is ChEBI, but the InChIKey is passed through so the link rows
    stay keyed by InChIKey, preserving the link-file schema and the step-3
    InChIKey -> HMDB remapping.
    """
    # ChEBI (int) -> InChIKey, precomputed once for fast per-call lookup.
    chebi_to_inchikey = {}
    for _, r in df.dropna(subset=['ChEBI']).iterrows():
        try:
            chebi_to_inchikey.setdefault(int(float(r['ChEBI'])), r['InChIKey'])
        except (ValueError, TypeError):
            continue

    def _query(chebi_id):
        try:
            key = chebi_to_inchikey.get(int(float(chebi_id)))
        except (ValueError, TypeError):
            key = None
        return get_uniprot_interactions(chebi_id, inchikey=key)

    return cache.process_single_database(
        df, 'UniProtKB', 'ChEBI', _query, resume=resume)
