"""
HMDB XML parsing and small data helpers.

Re-derived from Collect_PKNdata_metabolites.ipynb cells 5 & 13.
"""

import logging
import os
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

logger = logging.getLogger('pkn.io')

_HMDB_NS = {'hmdb': 'http://www.hmdb.ca'}


def _text(elem, path):
    found = elem.find(path, _HMDB_NS)
    return found.text if found is not None else None


def parse_hmdb(xml_path, max_metabolites=None):
    """
    Stream-parse the HMDB metabolites XML into a dataframe.

    Uses iterparse + elem.clear() so the multi-GB file is processed with bounded
    memory. Returns columns: Name, HMDB, ChEBI, KEGG, PubChem, IUPAC_Name, SMILES,
    InChIKey, PDB_ID, Kingdom, Super_Class, Sub_Class.
    """
    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"HMDB metabolites XML not found: {xml_path}")

    data = []
    context = ET.iterparse(xml_path, events=('end',))
    for _, elem in context:
        if elem.tag != '{http://www.hmdb.ca}metabolite':
            continue
        data.append({
            'Name': _text(elem, 'hmdb:name'),
            'HMDB': _text(elem, 'hmdb:accession'),
            'ChEBI': _text(elem, 'hmdb:chebi_id'),
            'KEGG': _text(elem, 'hmdb:kegg_id'),
            'PubChem': _text(elem, 'hmdb:pubchem_compound_id'),
            'IUPAC_Name': _text(elem, 'hmdb:iupac_name'),
            'SMILES': _text(elem, 'hmdb:smiles'),
            'InChIKey': _text(elem, 'hmdb:inchikey'),
            'PDB_ID': _text(elem, 'hmdb:pdb_id'),
            'Kingdom': _text(elem, 'hmdb:taxonomy/hmdb:kingdom'),
            'Super_Class': _text(elem, 'hmdb:taxonomy/hmdb:super_class'),
            'Sub_Class': _text(elem, 'hmdb:taxonomy/hmdb:sub_class'),
        })
        elem.clear()
        if max_metabolites is not None and len(data) >= max_metabolites:
            break

    df = pd.DataFrame(data)
    logger.info("Parsed %d metabolites from HMDB", len(df))
    return df


def smiles_to_canonical(smiles):
    """Canonicalize a SMILES string with RDKit; return np.nan on failure."""
    from rdkit import Chem  # imported lazily so the module loads without rdkit
    try:
        if smiles and smiles != 'NA':
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:  # noqa: BLE001
        pass
    return np.nan
