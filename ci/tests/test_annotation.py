"""Tests for annotation consistency checks."""

import pytest
from pathlib import Path
from checks import CheckResult
from checks.stockholm import parse_protein_sto, iter_dna_blocks, StoBlock
from checks.annotation import (
    check_catalytic_triad,
    check_element_structure,
)


TEST_ENTRY = Path(__file__).parent.parent.parent / "entries" / "test-entry"


def test_catalytic_triad_valid():
    result = CheckResult()
    block = parse_protein_sto(TEST_ENTRY / "protein.sto", result)
    assert result.passed

    triad_result = CheckResult()
    columns = check_catalytic_triad(block, triad_result)
    assert triad_result.passed
    assert len(columns) == 3


def test_catalytic_triad_wrong_residue():
    """Triad column doesn't contain D in sequences."""
    block = StoBlock(
        sequences={"A": "MKAE", "B": "MKAE", "C": "MKAE"},
        gc={"catalytic_triad": ".D.d"},
    )
    result = CheckResult()
    columns = check_catalytic_triad(block, result)
    # D at column 1 -> K (wrong), d at column 3 -> E (wrong)
    assert not result.passed


def test_catalytic_triad_too_many_marks():
    block = StoBlock(
        sequences={"A": "DDDDD"},
        gc={"catalytic_triad": "DdDDE"},
    )
    result = CheckResult()
    check_catalytic_triad(block, result)
    assert not result.passed


def test_element_structure_valid():
    """Full integration: element_structure from test entry should pass."""
    fmt_result = CheckResult()
    protein_block = parse_protein_sto(TEST_ENTRY / "protein.sto", fmt_result)
    assert fmt_result.passed

    dna_result = CheckResult()
    blocks = list(iter_dna_blocks(TEST_ENTRY / "dna.sto", dna_result))
    assert dna_result.passed

    annot_result = CheckResult()
    for dna_block in blocks:
        sid = dna_block.seq_ids[0]
        prot_seq = protein_block.sequences.get(sid, "")
        issues = []
        check_element_structure(dna_block, prot_seq, annot_result, issues)

    assert annot_result.passed


def test_element_structure_bad_tsd():
    """TSD not spelling TA should fail."""
    block = StoBlock(
        sequences={"A": "AAAA" + "GG" + "AAAA" + "ATGATGTGA" + "AAAA" + "GG" + "AAAA"},
        gc={"element_structure": "5555" + "AA" + "<<<<" + "012012012" + ">>>>" + "BB" + "3333"},
    )
    result = CheckResult()
    issues = []
    check_element_structure(block, None, result, issues)
    assert not result.passed
    assert any("TSD" in m for m in result.messages)


def test_element_structure_invalid_chars():
    block = StoBlock(
        sequences={"A": "ACGTACGT"},
        gc={"element_structure": "5555XXXX"},
    )
    result = CheckResult()
    issues = []
    check_element_structure(block, None, result, issues)
    assert not result.passed
