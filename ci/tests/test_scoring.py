"""Tests for scoring module."""

import pytest
from pathlib import Path
from checks import CheckResult
from checks.stockholm import parse_protein_sto
from checks.scoring import (
    score_diversity,
    score_annotation,
    score_functionality,
    score_size,
    compute_total_score,
)


TEST_ENTRY = Path(__file__).parent.parent.parent / "entries" / "test-entry"


def test_score_size_minimum():
    assert score_size(3) > 0
    assert score_size(2) == 0.0
    assert score_size(301) == 0.0


def test_score_size_logarithmic():
    s16 = score_size(16)
    s32 = score_size(32)
    s64 = score_size(64)
    assert s16 < s32 < s64
    assert abs(s64 - 1.0) < 0.01  # 64 is the target


def test_score_functionality_empty():
    assert score_functionality([]) == 0.0


def test_score_functionality_with_paralogs():
    rows = [
        {"max_paralog_identity": "0.995", "reference": "10.1234/test"},
        {"max_paralog_identity": "0.80"},
        {},
    ]
    score = score_functionality(rows)
    assert score > 0.0


def test_score_annotation_all_pass():
    checks = {
        "format": CheckResult(),
        "annotation": CheckResult(),
        "protein": CheckResult(),
        "dna": CheckResult(),
        "provenance": CheckResult(),
    }
    score = score_annotation(checks, [])
    assert score == 1.0


def test_score_annotation_with_warnings():
    checks = {
        "format": CheckResult(),
        "annotation": CheckResult(),
        "protein": CheckResult(),
        "dna": CheckResult(),
        "provenance": CheckResult(),
    }
    checks["protein"].warn("something minor")
    score = score_annotation(checks, [])
    assert 0.5 < score < 1.0


def test_score_diversity_valid():
    result = CheckResult()
    block = parse_protein_sto(TEST_ENTRY / "protein.sto", result)
    assert result.passed
    score = score_diversity(block, {"DD34D_mariner"})
    assert 0.0 <= score <= 1.0


def test_compute_total_score_hard_fail():
    """Hard fail in format or annotation should give total 0."""
    result = CheckResult()
    block = parse_protein_sto(TEST_ENTRY / "protein.sto", result)

    checks = {
        "format": CheckResult(),
        "annotation": CheckResult(),
        "protein": CheckResult(),
        "dna": CheckResult(),
        "provenance": CheckResult(),
    }
    checks["annotation"].fail("something bad")

    scores = compute_total_score(block, {"DD34D_mariner"}, checks, [], 3)
    assert scores["total"] == 0.0
