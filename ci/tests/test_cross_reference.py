"""Tests for cross-reference checks."""

import pytest
from checks import CheckResult
from checks.cross_reference import check_cross_references, parse_provenance
from pathlib import Path


TEST_ENTRY = Path(__file__).parent.parent.parent / "entries" / "test-entry"


def test_cross_reference_valid():
    result = CheckResult()
    prov = parse_provenance(TEST_ENTRY / "provenance.tsv", result)
    assert result.passed
    assert prov is not None

    xref_result = CheckResult()
    protein_ids = ["Mos1_Dmaur", "SynMar1_Dpse", "SynMar2_Agam"]
    dna_ids = ["Mos1_Dmaur", "SynMar1_Dpse", "SynMar2_Agam"]
    check_cross_references(protein_ids, dna_ids, prov, xref_result)
    assert xref_result.passed


def test_cross_reference_missing_dna_id():
    result = CheckResult()
    check_cross_references(
        protein_ids=["A", "B", "C"],
        dna_ids=["A", "B"],  # missing C
        provenance_rows=[{"id": "A"}, {"id": "B"}, {"id": "C"}],
        result=result,
    )
    assert not result.passed
    assert any("missing from dna.sto" in m for m in result.messages)


def test_cross_reference_missing_provenance_id():
    result = CheckResult()
    check_cross_references(
        protein_ids=["A", "B", "C"],
        dna_ids=["A", "B", "C"],
        provenance_rows=[{"id": "A"}, {"id": "B"}],  # missing C
        result=result,
    )
    assert not result.passed


def test_cross_reference_too_few():
    result = CheckResult()
    check_cross_references(
        protein_ids=["A", "B"],
        dna_ids=["A", "B"],
        provenance_rows=[{"id": "A"}, {"id": "B"}],
        result=result,
    )
    assert not result.passed
    assert any("minimum is 3" in m for m in result.messages)


def test_provenance_missing_columns(tmp_path):
    bad = tmp_path / "provenance.tsv"
    bad.write_text("id\tfamily\n")
    result = CheckResult()
    rows = parse_provenance(bad, result)
    assert not result.passed
