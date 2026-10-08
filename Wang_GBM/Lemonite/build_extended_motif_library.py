#!/usr/bin/env python3
"""
Build an EXTENDED HOMER motif library for the TF-motif validation of LemonTree regulators.

Why: HOMER's known-motif library has no motif for ~75% of the predicted module-TF pairs
(HOMER_motif_enrichment.R / HOMER_enhancer_enrichment.R report them as "no_motif_in_library").
Most of those TFs DO have a motif in JASPAR 2024 CORE and/or HOCOMOCO v12 (H12CORE). This script
converts those matrices into HOMER motif format and appends them to HOMER's known.motifs, so both
enrichment scripts can be re-run unchanged in logic (same backgrounds, same tiers, same q cut-off)
via `Rscript HOMER_*_enrichment.R extended`.

TRACEABILITY - every library motif <-> predicted-TF link is written to custom_motif_map.tsv:
  Relation = direct          the motif belongs to the predicted TF itself
  Relation = paralog_proxy   the predicted TF has NO motif anywhere; the motif of a same-family
                             paralog stands in for it (BATF2->BATF/BATF3 etc., PARALOG_PROXIES)
The enrichment scripts carry Relation / Motif_TF / Motif_DB / Matrix_ID into their output tables
and label paralog-derived .mvf entries "TF (via PARALOG)", so a paralog-derived call can never be
mistaken for a direct one. A proxy is only ever used when the TF has no direct motif.

Detection thresholds: HOMER scores a site as sum(ln(p_i/0.25)) (verified against `homer2 find
-mscore`). HOMER's own 472 known motifs have a median threshold p-value of 10^-3.98 under a uniform
background, so imported matrices get the score threshold at p = 1e-4 (exact DP, same background),
i.e. the same stringency as HOMER's library - with a floor of 0.44 x max score (HOMER's 5th
percentile) and a cap of 0.77 x max (95th percentile), so long near-deterministic and very short
matrices are not scanned at cutoffs outside anything in HOMER's own library.

Sources are fetched once and cached in the library directory.
"""

import csv, json, math, re, sys, time, urllib.request, collections
import concurrent.futures as cf
import numpy as np

ENRICH = '/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/transcriptomics_clustering/Enrichment'
LIBDIR = f'{ENRICH}/motif_library'
HOMER_KNOWN = '/home/borisvdm/Bioinformatics_software/homer/data/knownTFs/vertebrates/known.motifs'
PREV_RESULTS = [f'{ENRICH}/HOMER_motif/HOMER_motif_enrichment_46_modules.csv',
                f'{ENRICH}/HOMER_enhancer/HOMER_enhancer_enrichment_modules.csv']

P_THRESHOLD = 1e-4          # motif detection p-value (uniform background), matches HOMER's library
MIN_FRAC_OF_MAX = 0.44      # ...clamped to HOMER's own 5th-95th percentile of threshold / max score
MAX_FRAC_OF_MAX = 0.77      # (median 0.63). The floor binds for long near-deterministic matrices (15-22
                            # bp, 20-28 bits), where p=1e-4 alone would allow several mismatches and even
                            # give a negative threshold. The cap binds for 6-7 bp matrices (mostly
                            # homeodomain TAATTA-type), where p=1e-4 is unreachable (4^6 sequences) and
                            # the threshold would sit at the max score, i.e. exact-consensus-only.
PSEUDOCOUNT = 0.8           # total pseudocount spread over the 4 bases (JASPAR count matrices)
MIN_PROB = 0.001            # HOMER floors probabilities at 0.001
HOCOMOCO_QUALITY = set('ABC')   # grade D motifs are not used
JASPAR_SPECIES = ('Homo sapiens', 'Mus musculus')   # Drosophila/Ciona matrices are dropped
JASPAR_API = 'https://jaspar.elixir.no/api/v1'
HOCOMOCO_URL = 'https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12CORE/H12CORE_annotation.jsonl'

# Predicted TFs that have no motif in HOMER, JASPAR CORE or HOCOMOCO H12CORE, but whose gene family
# does. The paralogs' motifs stand in for them. Family-level motifs of these families are nearly
# identical (BATF/JUN-like TGA(C/G)TCA, NeuroD E-box, AP-2 GCCNNNGGC, CNC-bZIP MARE), so this is a
# family-level test, NOT evidence for the individual TF.
PARALOG_PROXIES = {
    'BATF2':  (['BATF', 'BATF3'],                'BATF family (bZIP)'),
    'NEUROD6': (['NEUROD1', 'NEUROD2'],          'NeuroD family (bHLH)'),
    'TFAP2D': (['TFAP2A', 'TFAP2B', 'TFAP2C', 'TFAP2E'], 'AP-2 family'),
    'NFE2L3': (['NFE2', 'NFE2L1', 'NFE2L2'],     'CNC-bZIP family'),
}

# ------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------

def fetch_json(url, cache=None):
    if cache:
        try:
            return json.load(open(cache))
        except FileNotFoundError:
            pass
    for i in range(4):
        try:
            d = json.load(urllib.request.urlopen(
                urllib.request.Request(url, headers={'Accept': 'application/json'}), timeout=90))
            if cache:
                json.dump(d, open(cache, 'w'))
            return d
        except Exception as e:
            err = e
            time.sleep(2 * (i + 1))
    raise err

def score_dist(P, bg=0.25, res=0.001):
    """Exact distribution of sum(ln(p/bg)) for a uniformly random sequence, on a `res` grid."""
    L = np.log(np.asarray(P) / bg)
    n = len(L)
    lo = int(np.floor(L.min(1).sum() / res)) - n
    hi = int(np.ceil(L.max(1).sum() / res)) + n
    dist = np.zeros(hi - lo + 1)
    dist[-lo] = 1.0
    for row in L:
        new = np.zeros_like(dist)
        for b in range(4):
            sh = int(round(row[b] / res))
            if sh >= 0:
                new[sh:] += dist[:len(dist) - sh] * bg
            else:
                new[:sh] += dist[-sh:] * bg
        dist = new
    return dist, lo, res

def threshold_for_p(P, p):
    dist, lo, res = score_dist(P)
    tail = np.cumsum(dist[::-1])[::-1]
    return (int(np.argmax(tail <= p)) + lo) * res

def counts_to_prob(counts):
    c = np.asarray(counts, float)                          # rows = positions, cols = A,C,G,T
    P = (c + PSEUDOCOUNT / 4) / (c.sum(1, keepdims=True) + PSEUDOCOUNT)
    P = np.maximum(P, MIN_PROB)
    return P / P.sum(1, keepdims=True)

def consensus(P):
    out = []
    for r in P:
        b = int(np.argmax(r))
        out.append('ACGT'[b] if r[b] >= 0.5 else 'N')
    return ''.join(out)

floored, capped = [], []    # motifs where the fraction-of-max clamp, not the p-value, set the threshold

def homer_block(P, name):
    max_score = float(np.log(P.max(1) / 0.25).sum())
    thr = threshold_for_p(P, P_THRESHOLD)
    if thr < MIN_FRAC_OF_MAX * max_score:
        floored.append((name, thr, MIN_FRAC_OF_MAX * max_score))
        thr = MIN_FRAC_OF_MAX * max_score
    elif thr > MAX_FRAC_OF_MAX * max_score:
        capped.append((name, thr, MAX_FRAC_OF_MAX * max_score))
        thr = MAX_FRAC_OF_MAX * max_score
    hdr = f'>{consensus(P)}\t{name}\t{thr:.6f}\t-1000\t0\t0\n'
    return hdr + ''.join('\t'.join(f'{x:.4f}' for x in row) + '\n' for row in P)

# ------------------------------------------------------------------------------------------------
# 1. which TFs need a custom motif
# ------------------------------------------------------------------------------------------------

import os
os.makedirs(LIBDIR, exist_ok=True)

predicted, no_motif, tested = set(), set(), set()
for f in PREV_RESULTS:
    for r in csv.DictReader(open(f)):
        predicted.add(r['TF'])
        (no_motif if r['Status'] == 'no_motif_in_library' else tested).add(r['TF'])
no_motif -= tested
print(f'{len(predicted)} predicted TFs | {len(tested)} already testable in HOMER | '
      f'{len(no_motif)} without a HOMER motif')

want_direct = set(no_motif)                                  # TFs whose own motif we look for
want_paralog = {p for _, (ps, _) in PARALOG_PROXIES.items() for p in ps}
wanted = want_direct | want_paralog

# ------------------------------------------------------------------------------------------------
# 2. JASPAR 2024 CORE (latest versions, monomers only) - by exact gene symbol
# ------------------------------------------------------------------------------------------------

listing, url = [], f'{JASPAR_API}/matrix/?collection=CORE&version=latest&page_size=500&format=json'
cache = f'{LIBDIR}/jaspar_CORE_listing.json'
if os.path.exists(cache):
    listing = json.load(open(cache))
else:
    while url:
        d = fetch_json(url)
        listing += d['results']
        url = d.get('next')
        url = url.replace('http://', 'https://') if url else url
    json.dump(listing, open(cache, 'w'))

jaspar_ids = collections.defaultdict(list)                   # TF -> [matrix ids]
for m in listing:
    nm = re.sub(r'\(var\.\d+\)', '', m['name']).strip().upper()
    if '::' not in nm and nm in wanted:
        jaspar_ids[nm].append(m['matrix_id'])

def jaspar_detail(mid):
    return mid, fetch_json(f'{JASPAR_API}/matrix/{mid}/?format=json', f'{LIBDIR}/jaspar_{mid}.json')

all_ids = sorted({i for v in jaspar_ids.values() for i in v})
with cf.ThreadPoolExecutor(8) as ex:
    jdetail = dict(ex.map(jaspar_detail, all_ids))

# ------------------------------------------------------------------------------------------------
# 3. HOCOMOCO v12 H12CORE (human gene symbol from the masterlist)
# ------------------------------------------------------------------------------------------------

hoc_file = f'{LIBDIR}/H12CORE_annotation.jsonl'
if not os.path.exists(hoc_file):
    urllib.request.urlretrieve(HOCOMOCO_URL, hoc_file)
hoc = collections.defaultdict(list)                          # TF -> [records]
for line in open(hoc_file):
    d = json.loads(line)
    h = d.get('masterlist_info', {}).get('species', {}).get('HUMAN')
    if h and h['gene_symbol'].upper() in wanted and d.get('quality') in HOCOMOCO_QUALITY:
        hoc[h['gene_symbol'].upper()].append(d)

# ------------------------------------------------------------------------------------------------
# 4. assemble library + traceability map
# ------------------------------------------------------------------------------------------------

entries = []      # (name, motif_tf, db, matrix_id, species, quality, P)
skipped = []
for tf, ids in sorted(jaspar_ids.items()):
    for mid in ids:
        d = jdetail[mid]
        sp = [s['name'] for s in d.get('species', [])]
        if not any(s in JASPAR_SPECIES for s in sp):
            skipped.append((tf, mid, 'JASPAR', ','.join(sp)))
            continue
        pfm = d['pfm']
        P = counts_to_prob(np.array([pfm[b] for b in 'ACGT']).T)
        entries.append((f'{tf}/JASPAR:{mid}', tf, 'JASPAR', mid, ','.join(sp), 'CORE', P))
for tf, recs in sorted(hoc.items()):
    for d in recs:
        P = counts_to_prob(d['pcm'])
        entries.append((f"{tf}/HOCOMOCO:{d['name'].split('.', 1)[1]}", tf, 'HOCOMOCO', d['name'],
                        'HUMAN', d['quality'], P))

map_rows = []
for name, tf, db, mid, sp, q, P in entries:
    if tf in predicted:
        map_rows.append((name, tf, db, mid, sp, q, tf, 'direct', ''))
    for target, (paralogs, basis) in PARALOG_PROXIES.items():
        if tf in paralogs:
            map_rows.append((name, tf, db, mid, sp, q, target, 'paralog_proxy', basis))

with open(f'{LIBDIR}/custom_motifs.motif', 'w') as f:
    for name, tf, db, mid, sp, q, P in entries:
        f.write(homer_block(P, name))
with open(f'{LIBDIR}/custom_motif_map.tsv', 'w') as f:
    f.write('Motif_Name\tMotif_TF\tMotif_DB\tMatrix_ID\tSpecies\tQuality\tPredicted_TF\tRelation\tParalog_basis\n')
    for r in map_rows:
        f.write('\t'.join(r) + '\n')
homer_txt = open(HOMER_KNOWN).read()
with open(f'{LIBDIR}/extended_known.motifs', 'w') as f:
    f.write(homer_txt if homer_txt.endswith('\n') else homer_txt + '\n')
    f.write(open(f'{LIBDIR}/custom_motifs.motif').read())

# ------------------------------------------------------------------------------------------------
# 5. report
# ------------------------------------------------------------------------------------------------

covered_direct = {r[6] for r in map_rows if r[7] == 'direct'}
covered_proxy = {r[6] for r in map_rows if r[7] == 'paralog_proxy'}
print(f'\nlibrary: {len(entries)} custom motifs '
      f'({sum(e[2] == "JASPAR" for e in entries)} JASPAR, {sum(e[2] == "HOCOMOCO" for e in entries)} HOCOMOCO)'
      f' appended to {homer_txt.count(chr(10) + ">") + 1} HOMER known motifs')
print(f'no-HOMER-motif TFs now covered directly: {len(no_motif & covered_direct)} of {len(no_motif)}')
print(f'via paralog proxy only: {sorted(covered_proxy - covered_direct)}')
print(f'still without any motif: {sorted(no_motif - covered_direct - covered_proxy)}')
if skipped:
    print('dropped (non vertebrate-model species):', skipped)
for target, (paralogs, basis) in PARALOG_PROXIES.items():
    got = sorted({r[1] for r in map_rows if r[6] == target and r[7] == 'paralog_proxy'})
    print(f'  proxy {target} <- {got} ({basis})')
thr = [float(l.split('\t')[2]) for l in open(f'{LIBDIR}/custom_motifs.motif') if l.startswith('>')]
print(f'threshold clamped to HOMER range of max score: floor {MIN_FRAC_OF_MAX} on {len(floored)}, '
      f'cap {MAX_FRAC_OF_MAX} on {len(capped)} of {len(entries)} motifs')
print(f'custom thresholds: min {min(thr):.2f} median {np.median(thr):.2f} max {max(thr):.2f}')
print(f'wrote: {LIBDIR}/{{extended_known.motifs,custom_motifs.motif,custom_motif_map.tsv}}')
