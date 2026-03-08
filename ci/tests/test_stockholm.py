"""Tests for Stockholm format parsing."""

import pytest
from pathlib import Path
from checks import CheckResult
from checks.stockholm import parse_protein_sto, iter_dna_blocks


TEST_ENTRY = Path(__file__).parent.parent.parent / "entries" / "test-entry"


def test_parse_valid_protein_sto():
    result = CheckResult()
    block = parse_protein_sto(TEST_ENTRY / "protein.sto", result)
    assert result.passed
    assert block is not None
    assert len(block.sequences) == 3
    assert "catalytic_triad" in block.gc
    assert "Mos1_Dmaur" in block.sequences


def test_parse_protein_sto_missing_file(tmp_path):
    result = CheckResult()
    block = parse_protein_sto(tmp_path / "nonexistent.sto", result)
    assert not result.passed
    assert block is None


def test_parse_protein_sto_no_header(tmp_path):
    bad = tmp_path / "protein.sto"
    bad.write_text("not a stockholm file\n//\n")
    result = CheckResult()
    block = parse_protein_sto(bad, result)
    assert not result.passed


def test_parse_protein_sto_no_triad(tmp_path):
    bad = tmp_path / "protein.sto"
    bad.write_text("# STOCKHOLM 1.0\nSeqA  MKDDD\n//\n")
    result = CheckResult()
    block = parse_protein_sto(bad, result)
    assert not result.passed


def test_parse_protein_sto_no_sequences(tmp_path):
    bad = tmp_path / "protein.sto"
    bad.write_text("# STOCKHOLM 1.0\n#=GC catalytic_triad ...\n//\n")
    result = CheckResult()
    block = parse_protein_sto(bad, result)
    assert not result.passed


def test_iter_dna_blocks_valid():
    result = CheckResult()
    blocks = list(iter_dna_blocks(TEST_ENTRY / "dna.sto", result))
    assert result.passed
    assert len(blocks) == 3
    for block in blocks:
        assert len(block.sequences) == 1
        assert "element_structure" in block.gc


def test_iter_dna_blocks_missing_file(tmp_path):
    result = CheckResult()
    blocks = list(iter_dna_blocks(tmp_path / "missing.sto", result))
    assert not result.passed
    assert len(blocks) == 0


def test_iter_dna_blocks_bad_block(tmp_path):
    bad = tmp_path / "dna.sto"
    # Block with two sequences (invalid)
    bad.write_text(
        "# STOCKHOLM 1.0\n"
        "SeqA  ACGT\n"
        "SeqB  ACGT\n"
        "#=GC element_structure 5555\n"
        "//\n"
    )
    result = CheckResult()
    blocks = list(iter_dna_blocks(bad, result))
    assert not result.passed


def test_iter_dna_blocks_missing_annotation(tmp_path):
    bad = tmp_path / "dna.sto"
    bad.write_text("# STOCKHOLM 1.0\nSeqA  ACGT\n//\n")
    result = CheckResult()
    blocks = list(iter_dna_blocks(bad, result))
    assert not result.passed
