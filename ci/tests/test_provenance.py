"""Tests for provenance checks."""

from checks import CheckResult
from checks.provenance import check_provenance


def test_provenance_valid():
    rows = [
        {"id": "A", "host_taxid": "7225", "assembly": "GCF_004382145.1",
         "start": "100", "end": "200", "strand": "+"},
    ]
    result = CheckResult()
    issues = {}
    check_provenance(rows, result, issues)
    assert result.passed


def test_provenance_bad_coordinates():
    rows = [
        {"id": "A", "host_taxid": "7225", "assembly": "GCF_004382145.1",
         "start": "500", "end": "100", "strand": "+"},
    ]
    result = CheckResult()
    issues = {}
    check_provenance(rows, result, issues)
    assert result.status == "warn"
    assert any("start" in m for m in result.messages)


def test_provenance_bad_strand():
    rows = [
        {"id": "A", "host_taxid": "7225", "assembly": "GCF_004382145.1",
         "start": "100", "end": "200", "strand": "forward"},
    ]
    result = CheckResult()
    issues = {}
    check_provenance(rows, result, issues)
    assert result.status == "warn"


def test_provenance_bad_taxid():
    rows = [
        {"id": "A", "host_taxid": "abc", "assembly": "GCF_004382145.1",
         "start": "100", "end": "200", "strand": "+"},
    ]
    result = CheckResult()
    issues = {}
    check_provenance(rows, result, issues)
    assert result.status == "warn"


def test_provenance_bad_assembly():
    rows = [
        {"id": "A", "host_taxid": "7225", "assembly": "bad_format",
         "start": "100", "end": "200", "strand": "+"},
    ]
    result = CheckResult()
    issues = {}
    check_provenance(rows, result, issues)
    assert result.status == "warn"
