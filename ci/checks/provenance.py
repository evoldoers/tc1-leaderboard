"""Provenance metadata checks."""

from __future__ import annotations

import re

from . import CheckResult

_ASSEMBLY_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")
_VALID_STRANDS = {"+", "-"}


def check_provenance(
    rows: list[dict],
    result: CheckResult,
    per_seq_issues: dict[str, list[str]],
) -> None:
    """Validate provenance metadata (scored checks, not hard fail).

    Args:
        rows: Parsed provenance.tsv rows.
        result: CheckResult to accumulate.
        per_seq_issues: Dict of seq_id -> issue list.
    """
    for row in rows:
        sid = row.get("id", "???")
        issues = per_seq_issues.setdefault(sid, [])

        # Coordinate consistency
        try:
            start = int(row.get("start", ""))
            end = int(row.get("end", ""))
            if start >= end:
                msg = f"start ({start}) >= end ({end})"
                result.warn(f"{sid}: {msg}")
                issues.append(msg)
        except ValueError:
            msg = "start/end coordinates are not valid integers"
            result.warn(f"{sid}: {msg}")
            issues.append(msg)

        # Strand
        strand = row.get("strand", "")
        if strand not in _VALID_STRANDS:
            msg = f"strand '{strand}' is not '+' or '-'"
            result.warn(f"{sid}: {msg}")
            issues.append(msg)

        # Assembly accession format
        assembly = row.get("assembly", "")
        if assembly and not _ASSEMBLY_RE.match(assembly):
            msg = f"assembly '{assembly}' does not match expected format GC[AF]_NNNNNNNNN.N"
            result.warn(f"{sid}: {msg}")
            issues.append(msg)

        # Taxonomy ID: should be a positive integer
        taxid = row.get("host_taxid", "")
        try:
            tid = int(taxid)
            if tid <= 0:
                raise ValueError
        except ValueError:
            msg = f"host_taxid '{taxid}' is not a valid positive integer"
            result.warn(f"{sid}: {msg}")
            issues.append(msg)
