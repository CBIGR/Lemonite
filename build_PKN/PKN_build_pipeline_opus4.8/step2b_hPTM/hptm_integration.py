"""
Step-2b orchestration: build the enzyme - histone-modification network from the
QuickGO and UniProtKB GO-annotation sources.

Output: hPTM_network.tsv with columns
    Enzyme, Mark, Activity, Source, GO_ID, URL
where Enzyme is the modifying/binding gene, Mark is the controlled-vocabulary
histone-mark node (e.g. 'H3K27me3'), and Activity is 'writer' | 'eraser' |
'reader'. This is the input step 3's combiner turns into 'histone-modification'
edges (Node1=Enzyme, Node2=Mark).

An edge supported by both QuickGO and UniProtKB appears once per source (two rows,
two provenance URLs) so the cross-validation is visible in the final network.
"""

import logging

import pandas as pd

import config
from step2b_hPTM import quickgo, uniprot_go

logger = logging.getLogger('pkn.hptm')

_COLS = ['Enzyme', 'Mark', 'Activity', 'Source', 'GO_ID', 'URL']


def build_hptm_network(resume=True):
    """Collect QuickGO + UniProtKB annotations; write hPTM_network.tsv."""
    quickgo_df = quickgo.run(resume=resume)
    uniprot_df = uniprot_go.run(resume=resume)

    frames = []
    for df in (quickgo_df, uniprot_df):
        if df is None or df.empty:
            continue
        sub = df.rename(columns={'Gene': 'Enzyme'})[
            ['Enzyme', 'Mark', 'Activity', 'Source', 'GO_ID', 'URL']]
        frames.append(sub)

    if not frames:
        hptm = pd.DataFrame(columns=_COLS)
    else:
        hptm = pd.concat(frames, ignore_index=True)
        hptm['Enzyme'] = hptm['Enzyme'].astype(str)
        hptm = hptm[hptm['Enzyme'].str.strip().ne('') & hptm['Enzyme'].ne('nan')]
        hptm = hptm.drop_duplicates(subset=['Enzyme', 'Mark', 'Activity', 'Source', 'GO_ID'])

    hptm.to_csv(config.HPTM_NETWORK_FILE, sep='\t', index=False)
    if not hptm.empty:
        logger.info("hPTM network: %d edges, %d enzymes, %d marks (%d writer / %d eraser / %d reader) -> %s",
                    len(hptm), hptm['Enzyme'].nunique(), hptm['Mark'].nunique(),
                    (hptm['Activity'] == 'writer').sum(),
                    (hptm['Activity'] == 'eraser').sum(),
                    (hptm['Activity'] == 'reader').sum(),
                    config.HPTM_NETWORK_FILE)
    else:
        logger.warning("hPTM network: 0 edges -> %s", config.HPTM_NETWORK_FILE)
    return hptm
