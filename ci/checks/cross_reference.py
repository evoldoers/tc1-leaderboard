"""Cross-reference checks: ID consistency across protein.sto, dna.sto, provenance.tsv."""

from __future__ import annotations

import csv
from pathlib import Path

from . import CheckResult

REQUIRED_COLUMNS = [
    "id", "family", "host_species", "host_taxid", "assembly",
    "chrom", "start", "end", "strand", "source", "reference",
]


def parse_provenance(path: Path, result: CheckResult) -> list[dict] | None:
    """Parse provenance.tsv and validate required columns.

    Returns list of row dicts, or None on hard failure.
    """
    try:
        text = path.read_text()
    except (OSError, UnicodeDecodeError) as e:
        result.fail(f"Cannot read {path.name}: {e}")
        return None

    lines = text.strip().splitlines()
    if not lines:
        result.fail(f"{path.name}: file is empty")
        return None

    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames is None:
        result.fail(f"{path.name}: cannot parse header")
        return None

    missing_cols = set(REQUIRED_COLUMNS) - set(reader.fieldnames)
    if missing_cols:
        result.fail(f"{path.name}: missing required columns: {', '.join(sorted(missing_cols))}")
        return None

    rows = list(reader)
    if not rows:
        result.fail(f"{path.name}: no data rows")
        return None

    return rows


def check_cross_references(
    protein_ids: list[str],
    dna_ids: list[str],
    provenance_rows: list[dict],
    result: CheckResult,
) -> None:
    """Check that IDs match across all three files."""
    prot_set = set(protein_ids)
    dna_set = set(dna_ids)
    prov_set = {row["id"] for row in provenance_rows}

    # protein.sto IDs not in dna.sto
    missing_dna = prot_set - dna_set
    if missing_dna:
        result.fail(
            f"IDs in protein.sto but missing from dna.sto: {', '.join(sorted(missing_dna))}"
        )

    # dna.sto IDs not in protein.sto
    extra_dna = dna_set - prot_set
    if extra_dna:
        result.fail(
            f"IDs in dna.sto but missing from protein.sto: {', '.join(sorted(extra_dna))}"
        )

    # protein.sto IDs not in provenance.tsv
    missing_prov = prot_set - prov_set
    if missing_prov:
        result.fail(
            f"IDs in protein.sto but missing from provenance.tsv: {', '.join(sorted(missing_prov))}"
        )

    # provenance.tsv IDs not in protein.sto
    extra_prov = prov_set - prot_set
    if extra_prov:
        result.fail(
            f"IDs in provenance.tsv but missing from protein.sto: {', '.join(sorted(extra_prov))}"
        )

    # Minimum entry count
    n = len(prot_set)
    if n < 3:
        result.fail(f"Only {n} sequences submitted; minimum is 3")
    if n > 300:
        result.fail(f"{n} sequences submitted; maximum is 300")
