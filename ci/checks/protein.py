"""Protein-level checks: length, triad spacing, domain predictions."""

from __future__ import annotations

from . import CheckResult
from .stockholm import StoBlock

# Expected transposase length ranges by family
FAMILY_LENGTH_RANGES: dict[str, tuple[int, int]] = {
    "DD34D_mariner": (250, 450),
    "DD34E_Tc1": (250, 450),
    "DDxD_pogo": (250, 500),
    "IS630": (180, 350),
}
DEFAULT_LENGTH_RANGE = (180, 500)

# Expected triad spacing (residues between 2nd and 3rd catalytic residue)
FAMILY_TRIAD_SPACING: dict[str, tuple[int, int]] = {
    "DD34D_mariner": (30, 38),
    "DD34E_Tc1": (30, 38),
    "DDxD_pogo": (25, 50),
    "IS630": (25, 50),
}
DEFAULT_TRIAD_SPACING = (25, 50)


def check_protein(
    protein_block: StoBlock,
    triad_columns: list[int],
    families: dict[str, str],
    result: CheckResult,
    per_seq_issues: dict[str, list[str]],
) -> None:
    """Run protein-level checks (scored, not hard fail).

    Args:
        protein_block: Parsed protein Stockholm block.
        triad_columns: 0-indexed column positions of catalytic triad.
        families: Mapping of seq_id -> family from provenance.
        result: CheckResult to accumulate.
        per_seq_issues: Dict of seq_id -> issue list.
    """
    for sid, seq in protein_block.sequences.items():
        ungapped = seq.replace("-", "").replace(".", "")
        aa_len = len(ungapped)
        family = families.get(sid, "")
        issues = per_seq_issues.setdefault(sid, [])

        # Length check
        lo, hi = FAMILY_LENGTH_RANGES.get(family, DEFAULT_LENGTH_RANGE)
        if aa_len < lo or aa_len > hi:
            msg = (
                f"protein length {aa_len} aa outside expected range "
                f"[{lo}, {hi}] for {family or 'unknown family'}"
            )
            result.warn(f"{sid}: {msg}")
            issues.append(msg)

    # Triad spacing check (uses alignment columns, so ungap to get spacing)
    if len(triad_columns) >= 3:
        for sid, seq in protein_block.sequences.items():
            family = families.get(sid, "")
            issues = per_seq_issues.setdefault(sid, [])

            # Count non-gap residues between 2nd and 3rd triad columns
            col2, col3 = triad_columns[1], triad_columns[2]
            between = seq[col2 + 1:col3]
            spacing = sum(1 for ch in between if ch not in ("-", "."))

            lo, hi = FAMILY_TRIAD_SPACING.get(family, DEFAULT_TRIAD_SPACING)
            if spacing < lo or spacing > hi:
                msg = (
                    f"triad spacing {spacing} residues between 2nd and 3rd "
                    f"catalytic residue (expected [{lo}, {hi}] for "
                    f"{family or 'unknown family'})"
                )
                result.warn(f"{sid}: {msg}")
                issues.append(msg)
