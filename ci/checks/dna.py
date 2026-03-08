"""DNA-level checks: element length, TIR length, GC content, ambiguous bases."""

from __future__ import annotations

import re

from . import CheckResult
from .stockholm import StoBlock

# Expected element length ranges (excluding flanking), by family
FAMILY_ELEMENT_RANGES: dict[str, tuple[int, int]] = {
    "DD34D_mariner": (1100, 2500),
    "DD34E_Tc1": (1200, 2500),
    "DDxD_pogo": (1200, 3000),
    "IS630": (800, 2000),
}
DEFAULT_ELEMENT_RANGE = (800, 3000)

# Expected TIR length ranges
FAMILY_TIR_RANGES: dict[str, tuple[int, int]] = {
    "DD34D_mariner": (15, 40),
    "DD34E_Tc1": (15, 200),
    "DDxD_pogo": (10, 200),
    "IS630": (10, 100),
}
DEFAULT_TIR_RANGE = (10, 200)


def check_dna(
    dna_block: StoBlock,
    family: str,
    result: CheckResult,
    seq_issues: list[str],
) -> None:
    """Run DNA-level checks (scored, not hard fail).

    Args:
        dna_block: Parsed single-sequence DNA Stockholm block.
        family: Family classification from provenance.
        result: CheckResult to accumulate.
        seq_issues: Issue list for this sequence.
    """
    sid = dna_block.seq_ids[0]
    seq = dna_block.sequences[sid].upper()
    annot = dna_block.gc["element_structure"]

    # Compute element region (everything between flanking, exclusive of TSDs)
    # Element = TIR_left + interior + TIR_right
    elem_chars = [c for c in annot if c in "<>012tn"]
    elem_len = len(elem_chars)

    lo, hi = FAMILY_ELEMENT_RANGES.get(family, DEFAULT_ELEMENT_RANGE)
    if elem_len < lo or elem_len > hi:
        msg = (
            f"element length {elem_len} bp outside expected range "
            f"[{lo}, {hi}] for {family or 'unknown family'}"
        )
        result.warn(f"{sid}: {msg}")
        seq_issues.append(msg)

    # TIR length
    tir_left_len = annot.count("<")
    tir_right_len = annot.count(">")

    tir_lo, tir_hi = FAMILY_TIR_RANGES.get(family, DEFAULT_TIR_RANGE)
    for label, tir_len in [("left", tir_left_len), ("right", tir_right_len)]:
        if tir_len < tir_lo or tir_len > tir_hi:
            msg = (
                f"{label} TIR length {tir_len} bp outside expected range "
                f"[{tir_lo}, {tir_hi}] for {family or 'unknown family'}"
            )
            result.warn(f"{sid}: {msg}")
            seq_issues.append(msg)

    # GC content of ORF
    orf_bases = [seq[i] for i, c in enumerate(annot) if c in "012"]
    if orf_bases:
        gc_count = sum(1 for b in orf_bases if b in "GC")
        gc_pct = gc_count / len(orf_bases)
        if gc_pct < 0.15 or gc_pct > 0.75:
            msg = f"ORF GC content {gc_pct:.1%} outside plausible range [15%, 75%]"
            result.warn(f"{sid}: {msg}")
            seq_issues.append(msg)

    # Ambiguous bases in ORF
    n_count = sum(1 for b in orf_bases if b == "N")
    if n_count > 0:
        msg = f"{n_count} ambiguous base(s) (N) in ORF region"
        result.warn(f"{sid}: {msg}")
        seq_issues.append(msg)
