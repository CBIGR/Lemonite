#!/usr/bin/env python3
"""
Minimal test script for the PKN_build_pipeline.

Tests core pipeline logic using small mock datasets so that no real
database files or network access are required.  Real-file retrievers
are verified only when their input files exist on disk.

Usage:
    python test_minimal.py
    python test_minimal.py -v          # verbose (show full tracebacks)
    python test_minimal.py --real      # also run local-file retrievers
"""

import sys
import os
import logging
import argparse
import tempfile
import textwrap
import traceback
from pathlib import Path
from unittest.mock import patch

import pandas as pd

# ---------------------------------------------------------------------------
# Path setup – add the pipeline package to sys.path
# ---------------------------------------------------------------------------
PIPELINE_DIR = Path(__file__).parent / "PKN_build_pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.WARNING,          # suppress noisy pipeline INFO messages
    format="%(levelname)s | %(name)s | %(message)s",
)
root_logger = logging.getLogger()

# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------
RESULTS: list[dict] = []
VERBOSE = False


def run_test(name: str, fn):
    """Execute *fn* and record pass / fail."""
    try:
        fn()
        RESULTS.append({"name": name, "status": "PASS", "msg": ""})
        print(f"  [PASS] {name}")
    except Exception as exc:
        msg = traceback.format_exc() if VERBOSE else str(exc)
        RESULTS.append({"name": name, "status": "FAIL", "msg": msg})
        print(f"  [FAIL] {name}")
        if VERBOSE:
            print(textwrap.indent(traceback.format_exc(), "        "))
        else:
            print(f"        {exc}")


def skip_test(name: str, reason: str):
    RESULTS.append({"name": name, "status": "SKIP", "msg": reason})
    print(f"  [SKIP] {name}  ({reason})")


# ---------------------------------------------------------------------------
# Mock metabolite data  (5 well-known metabolites)
# ---------------------------------------------------------------------------
MOCK_METABOLITES = [
    {"HMDB_ID": "HMDB0000122", "name": "Glucose",    "InChIKey": "WQZGKKKJIJFFOK-GASJEMHNSA-N", "SMILES": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
    {"HMDB_ID": "HMDB0000161", "name": "Alanine",    "InChIKey": "QNAYBMKLOCPYGJ-REOHCLBHSA-N", "SMILES": "C[C@@H](N)C(O)=O"},
    {"HMDB_ID": "HMDB0000067", "name": "Citric acid","InChIKey": "KRKNYBCHXYNGOX-UHFFFAOYSA-N", "SMILES": "OC(CC(O)(C(O)=O)CC(O)=O)(C(O)=O)"},
    {"HMDB_ID": "HMDB0001294", "name": "AMP",        "InChIKey": "UDMBCSSLTHHNCD-KQYNXXCUSA-N", "SMILES": "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O"},
    {"HMDB_ID": "HMDB0000517", "name": "L-Carnitine","InChIKey": "PHIQHXFUZVPYII-ZCFIWIBFSA-N", "SMILES": "C[N+](C)(C)C[C@@H](O)CC([O-])=O"},
]

MOCK_METABOLITES_DF = pd.DataFrame(MOCK_METABOLITES)

# Minimal mock interactions – shape expected by integration functions
MOCK_MG_INTERACTIONS = pd.DataFrame({
    "HMDB_ID": ["HMDB0000122", "HMDB0000122", "HMDB0000161", "HMDB0000067"],
    "Gene":    ["HK1",         "GCK",          "BCAT1",        "ACO2"],
    "Source":  ["mock_db",     "mock_db",       "mock_db",      "mock_db"],
})

MOCK_PPI_INTERACTIONS = pd.DataFrame({
    "GeneA":  ["HK1",  "GCK",  "ACO2"],
    "GeneB":  ["GCK",  "ACO2", "IDH2"],
    "Source": ["mock", "mock", "mock"],
})


# ===========================================================================
# TEST DEFINITIONS
# ===========================================================================

def test_imports_utils():
    from utils.file_io import load_hmdb_metabolites, save_interactions  # noqa
    from utils.pipeline import DatabaseRetriever, LocalFileRetriever     # noqa
    from utils.api_retry import retry_api_call                           # noqa


def test_imports_step1():
    from step1_metabolites import preprocessing, integration             # noqa
    from step1_metabolites.biogrid import BioGRIDRetriever               # noqa
    from step1_metabolites.stitch import STITCHRetriever                 # noqa
    from step1_metabolites.metalinks import MetalinksRetriever           # noqa
    from step1_metabolites.lincs import LINCSRetriever                   # noqa
    from step1_metabolites.chembl import ChEMBLRetriever                 # noqa
    from step1_metabolites.intact import IntActRetriever                 # noqa
    from step1_metabolites.gem import GEMRetriever                       # noqa
    from step1_metabolites.l1000 import L1000Retriever                   # noqa
    from step1_metabolites.uniprot import UniProtRetriever               # noqa


def test_imports_step2():
    from step2_proteins import ppi_integration                           # noqa
    from step2_proteins.biogrid_ppi import BioGRIDPPIRetriever           # noqa
    from step2_proteins.huri import HuRIRetriever                        # noqa
    from step2_proteins.string_api import STRINGRetriever                # noqa


def test_imports_step3():
    from step3_final import combiner, annotator, analysis, visualization # noqa


def test_preprocessing_mock():
    """preprocess_metabolites should return a DataFrame with required columns."""
    from step1_metabolites.preprocessing import preprocess_metabolites  # noqa

    # preprocess_metabolites may call ChEMBL API; patch it to return 'none'
    with patch(
        "step1_metabolites.preprocessing.get_chembl_id_from_smiles",
        return_value="none",
    ):
        result = preprocess_metabolites(MOCK_METABOLITES)

    assert isinstance(result, pd.DataFrame), "Expected a DataFrame"
    assert "HMDB_ID" in result.columns, "Missing HMDB_ID column"
    assert len(result) == len(MOCK_METABOLITES), "Row count mismatch"


def test_step1_integration_mock():
    """integrate_databases should combine mock results correctly."""
    from step1_metabolites.integration import integrate_databases

    db_results = {
        "db_A": pd.DataFrame({
            "HMDB_ID": ["HMDB0000122", "HMDB0000161"],
            "Gene":    ["HK1",          "BCAT1"],
            "Source":  ["db_A",         "db_A"],
        }),
        "db_B": pd.DataFrame({
            "HMDB_ID": ["HMDB0000122", "HMDB0000067"],
            "Gene":    ["HK1",          "ACO2"],
            "Source":  ["db_B",         "db_B"],
        }),
        "db_C": None,   # simulate a retriever that returned nothing
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_metabolite_gene_pkn = os.path.join(tmpdir, "mg_pkn.tsv")
        with patch("step1_metabolites.integration.config") as mock_cfg:
            mock_cfg.METABOLITE_GENE_PKN = tmp_metabolite_gene_pkn
            result = integrate_databases(db_results)

    assert isinstance(result, pd.DataFrame), "Expected a DataFrame"
    assert set(result.columns) >= {"HMDB_ID", "Gene", "Source"}, "Missing columns"
    # HMDB0000122-HK1 appears in both db_A and db_B; sources should be merged
    hk1_row = result[(result["HMDB_ID"] == "HMDB0000122") & (result["Gene"] == "HK1")]
    assert len(hk1_row) == 1, "Duplicate pair should be deduplicated"
    assert "db_A" in hk1_row["Source"].values[0], "db_A source missing"
    assert "db_B" in hk1_row["Source"].values[0], "db_B source missing"
    assert len(result) == 3, f"Expected 3 unique pairs, got {len(result)}"


def test_step2_integration_mock():
    """integrate_ppi_databases should combine mock PPI results."""
    from step2_proteins.ppi_integration import integrate_ppi_databases

    ppi_results = {
        "STRING": pd.DataFrame({
            "GeneA":         ["HK1",  "GCK"],
            "GeneB":         ["GCK",  "ACO2"],
            "Source":        ["STRING", "STRING"],
            "combined_score": [700, 650],
        }),
        "HuRI": pd.DataFrame({
            "GeneA": ["ACO2"],
            "GeneB": ["IDH2"],
            "Source": ["HuRI"],
        }),
        "BioGRID": None,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_ppi_out = os.path.join(tmpdir, "ppi.tsv")
        with patch("step2_proteins.ppi_integration.config") as mock_cfg:
            mock_cfg.PPI_OUTPUT_FILE = tmp_ppi_out
            result = integrate_ppi_databases(ppi_results)

    assert isinstance(result, pd.DataFrame), "Expected a DataFrame"
    assert {"GeneA", "GeneB", "Source"} <= set(result.columns), "Missing PPI columns"
    assert len(result) >= 3, "Expected at least 3 PPI rows"


def test_step3_combiner_mock():
    """combine_networks should vertically join metabolite-gene and PPI data."""
    from step3_final.combiner import combine_networks

    with tempfile.TemporaryDirectory() as tmpdir:
        mg_path  = os.path.join(tmpdir, "mg_pkn.tsv")
        ppi_path = os.path.join(tmpdir, "ppi.tsv")
        out_path = os.path.join(tmpdir, "final_pkn.tsv")

        MOCK_MG_INTERACTIONS.to_csv(mg_path, sep="\t", index=False)

        ppi_df = MOCK_PPI_INTERACTIONS.rename(columns={"GeneA": "Node1", "GeneB": "Node2"})
        ppi_df["Type"] = "PPI"
        ppi_df.to_csv(ppi_path, sep="\t", index=False)

        with patch("step3_final.combiner.config") as mock_cfg:
            mock_cfg.METABOLITE_GENE_PKN = mg_path
            mock_cfg.PPI_NETWORK         = ppi_path
            mock_cfg.FINAL_PKN_FILE      = out_path
            result = combine_networks()

    assert isinstance(result, pd.DataFrame), "Expected a DataFrame"
    assert {"Node1", "Node2", "Type", "Source"} <= set(result.columns), "Missing columns"
    types = set(result["Type"].unique())
    assert "metabolite-gene" in types, "'metabolite-gene' edge type missing"
    assert "PPI" in types, "'PPI' edge type missing"
    assert len(result) == len(MOCK_MG_INTERACTIONS) + len(MOCK_PPI_INTERACTIONS), \
        "Row count mismatch in combined network"


def test_file_io_save_load():
    """save_interactions / load round-trip."""
    from utils.file_io import save_interactions

    df = MOCK_MG_INTERACTIONS.copy()
    with tempfile.TemporaryDirectory() as tmpdir:
        out = os.path.join(tmpdir, "test_interactions.csv")
        save_interactions(df, out, add_source=False)
        assert os.path.exists(out), "Output file not created"
        loaded = pd.read_csv(out)
        assert len(loaded) == len(df), "Saved/loaded row count mismatch"
        assert set(loaded.columns) >= {"HMDB_ID", "Gene"}, "Missing columns in saved file"


# ---------------------------------------------------------------------------
# Optional: real-file retriever smoke-tests (run only when files exist)
# ---------------------------------------------------------------------------

def _maybe_test_metalinks(real: bool):
    import config as cfg
    path = cfg.METALINKS_PATH
    if not real or not os.path.exists(path):
        skip_test("local_retriever:MetalinksDB", f"file not found: {path}")
        return

    def _test():
        from step1_metabolites.metalinks import MetalinksRetriever
        with tempfile.TemporaryDirectory() as tmpdir:
            r = MetalinksRetriever(cache_file=os.path.join(tmpdir, "cache.csv"))
            result = r.get_interactions(MOCK_METABOLITES)
        assert isinstance(result, pd.DataFrame)
        assert {"HMDB_ID", "Gene", "Source"} <= set(result.columns)

    run_test("local_retriever:MetalinksDB", _test)


def _maybe_test_biogrid(real: bool):
    import config as cfg
    path = cfg.BIOGRID_LOCATION
    if not real or not os.path.exists(path):
        skip_test("local_retriever:BioGRID_metabolites", f"file not found: {path}")
        return

    def _test():
        from step1_metabolites.biogrid import BioGRIDRetriever
        with tempfile.TemporaryDirectory() as tmpdir:
            r = BioGRIDRetriever(cache_file=os.path.join(tmpdir, "cache.csv"))
            result = r.get_interactions(MOCK_METABOLITES)
        assert isinstance(result, pd.DataFrame)

    run_test("local_retriever:BioGRID_metabolites", _test)


def _maybe_test_lincs(real: bool):
    import config as cfg
    files = [
        cfg.LINCS_COMPOUND_MAPPING,
        cfg.LINCS_TARGET_MAPPING,
        cfg.LINCS_BIOCHEM_AGG,
    ]
    if not real or not all(os.path.exists(f) for f in files):
        skip_test("local_retriever:LINCS", "one or more LINCS files not found")
        return

    def _test():
        from step1_metabolites.lincs import LINCSRetriever
        with tempfile.TemporaryDirectory() as tmpdir:
            r = LINCSRetriever(cache_file=os.path.join(tmpdir, "cache.csv"))
            result = r.get_interactions(MOCK_METABOLITES)
        assert isinstance(result, pd.DataFrame)

    run_test("local_retriever:LINCS", _test)


def _maybe_test_biogrid_ppi(real: bool):
    import config as cfg
    path = cfg.BIOGRID_PPI_LOCATION
    if not real or not os.path.exists(path):
        skip_test("local_retriever:BioGRID_PPI", f"file not found: {path}")
        return

    def _test():
        from step2_proteins.biogrid_ppi import BioGRIDPPIRetriever
        genes = ["HK1", "GCK", "ACO2", "IDH2", "BCAT1"]
        with tempfile.TemporaryDirectory() as tmpdir:
            r = BioGRIDPPIRetriever(cache_file=os.path.join(tmpdir, "cache.csv"))
            result = r.get_interactions(genes)
        assert isinstance(result, pd.DataFrame)

    run_test("local_retriever:BioGRID_PPI", _test)


def _maybe_test_huri(real: bool):
    import config as cfg
    path = cfg.HURI_LOCATION
    if not real or not os.path.exists(path):
        skip_test("local_retriever:HuRI", f"file not found: {path}")
        return

    def _test():
        from step2_proteins.huri import HuRIRetriever
        genes = ["HK1", "GCK", "ACO2", "IDH2", "BCAT1"]
        with tempfile.TemporaryDirectory() as tmpdir:
            r = HuRIRetriever(cache_file=os.path.join(tmpdir, "cache.csv"))
            result = r.get_interactions(genes)
        assert isinstance(result, pd.DataFrame)

    run_test("local_retriever:HuRI", _test)


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    global VERBOSE

    parser = argparse.ArgumentParser(description="Minimal PKN pipeline tests")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Show full tracebacks on failure")
    parser.add_argument("--real", action="store_true",
                        help="Also run local-file retriever tests (requires DB files)")
    args = parser.parse_args()
    VERBOSE = args.verbose

    print("=" * 64)
    print("PKN_build_pipeline  –  minimal test suite")
    print("=" * 64)

    # ------------------------------------------------------------------
    # Core tests (always run, no real files needed)
    # ------------------------------------------------------------------
    print("\n--- Module imports ---")
    run_test("imports:utils",  test_imports_utils)
    run_test("imports:step1",  test_imports_step1)
    run_test("imports:step2",  test_imports_step2)
    run_test("imports:step3",  test_imports_step3)

    print("\n--- Logic tests with mock data ---")
    run_test("preprocessing:mock_metabolites",     test_preprocessing_mock)
    run_test("step1:integrate_databases_mock",     test_step1_integration_mock)
    run_test("step2:integrate_ppi_databases_mock", test_step2_integration_mock)
    run_test("step3:combine_networks_mock",        test_step3_combiner_mock)
    run_test("utils:file_io_save_load",            test_file_io_save_load)

    # ------------------------------------------------------------------
    # Optional real-file tests
    # ------------------------------------------------------------------
    print("\n--- Local-file retriever tests ---")
    _maybe_test_metalinks(args.real)
    _maybe_test_biogrid(args.real)
    _maybe_test_lincs(args.real)
    _maybe_test_biogrid_ppi(args.real)
    _maybe_test_huri(args.real)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_pass = sum(1 for r in RESULTS if r["status"] == "PASS")
    n_fail = sum(1 for r in RESULTS if r["status"] == "FAIL")
    n_skip = sum(1 for r in RESULTS if r["status"] == "SKIP")

    print("\n" + "=" * 64)
    print(f"Results:  {n_pass} passed  |  {n_fail} failed  |  {n_skip} skipped")
    print("=" * 64)

    if n_fail:
        print("\nFailed tests:")
        for r in RESULTS:
            if r["status"] == "FAIL":
                print(f"  • {r['name']}")
                if r["msg"] and not VERBOSE:
                    first_line = r["msg"].split("\n")[0]
                    print(f"    {first_line}")
        sys.exit(1)


if __name__ == "__main__":
    main()
