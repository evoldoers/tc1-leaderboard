"""Annotation consistency checks for catalytic_triad and element_structure."""

from __future__ import annotations

import re

from . import CheckResult
from .stockholm import StoBlock

# Valid element_structure characters
_VALID_STRUCT_CHARS = set("53<>ABnt012.")

# Expected structure pattern: 5+AA<+[t012n.]+>+BB3+
# (with possible gaps/flexibility)
_STRUCT_REGION_RE = re.compile(
    r"^(5+)(AA)(<+)((?:[t012n.])+)(>+)(BB)(3+)$"
)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _reverse_complement(seq: str) -> str:
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def _translate(dna: str) -> str:
    """Translate DNA to protein (stops at first stop codon)."""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def _pairwise_identity(s1: str, s2: str) -> float:
    """Fraction of identical residues between two equal-length strings."""
    if not s1 or not s2:
        return 0.0
    matches = sum(a == b for a, b in zip(s1, s2))
    return matches / max(len(s1), len(s2))


def check_catalytic_triad(
    protein_block: StoBlock,
    result: CheckResult,
) -> list[int]:
    """Validate the catalytic_triad annotation in protein.sto.

    Returns 0-indexed column positions of the triad residues (for scoring use).
    """
    triad_str = protein_block.gc["catalytic_triad"]
    # Find marked positions
    marked = []
    for i, ch in enumerate(triad_str):
        if ch in "DdE":
            marked.append((i, ch))
        elif ch != ".":
            result.fail(
                f"catalytic_triad: unexpected character '{ch}' at column {i+1} "
                "(expected D, d, E, or .)"
            )

    if len(marked) < 2 or len(marked) > 3:
        result.fail(
            f"catalytic_triad: expected 2-3 marked positions, found {len(marked)}"
        )
        return []

    # Validate the pattern: first should be D, second should be d, third D or E
    expected_patterns = [
        ["D", "d", "E"],  # DD?E families (e.g. DD34E Tc1)
        ["D", "d", "D"],  # DDD families (e.g. DD34D mariner)
        ["D", "d"],       # 2-residue triad annotation (rare)
    ]
    actual_pattern = [ch for _, ch in marked]
    if actual_pattern not in expected_patterns:
        result.fail(
            f"catalytic_triad: unexpected marker pattern "
            f"{''.join(actual_pattern)} (expected D,d,E or D,d,D)"
        )
        return []

    # Check that marked positions are D/D/E (or D/D/D) in majority of sequences
    columns = [pos for pos, _ in marked]
    expected_aa = {"D": "D", "d": "D", "E": "E"}

    for col, marker in marked:
        exp_aa = expected_aa[marker]
        residues = []
        for sid, seq in protein_block.sequences.items():
            if col < len(seq):
                aa = seq[col]
                if aa not in (".", "-"):
                    residues.append((sid, aa))

        if not residues:
            result.fail(f"catalytic_triad: column {col+1} has no residues (all gaps)")
            continue

        matching = sum(1 for _, aa in residues if aa == exp_aa)
        total = len(residues)
        if matching < total:
            mismatches = [(sid, aa) for sid, aa in residues if aa != exp_aa]
            if matching == 0:
                result.fail(
                    f"catalytic_triad: column {col+1} expected {exp_aa} but none match: "
                    + ", ".join(f"{sid}={aa}" for sid, aa in mismatches)
                )
            elif matching / total < 0.5:
                result.fail(
                    f"catalytic_triad: column {col+1} expected {exp_aa} but "
                    f"only {matching}/{total} match: "
                    + ", ".join(f"{sid}={aa}" for sid, aa in mismatches)
                )
            else:
                result.warn(
                    f"catalytic_triad: column {col+1} {matching}/{total} are {exp_aa}: "
                    + ", ".join(f"{sid}={aa}" for sid, aa in mismatches)
                )

    return columns


def check_element_structure(
    dna_block: StoBlock,
    protein_seq: str | None,
    result: CheckResult,
    seq_issues: list[str],
) -> None:
    """Validate the element_structure annotation for a single dna.sto block.

    Args:
        dna_block: Parsed Stockholm block with one sequence.
        protein_seq: The corresponding protein sequence from protein.sto (ungapped).
        result: CheckResult to accumulate errors/warnings.
        seq_issues: Per-sequence issue list for this element.
    """
    sid = dna_block.seq_ids[0]
    seq = dna_block.sequences[sid].upper()
    annot = dna_block.gc["element_structure"]

    # Check valid characters
    invalid = set(annot) - _VALID_STRUCT_CHARS
    if invalid:
        result.fail(
            f"{sid}: element_structure contains invalid characters: "
            f"{''.join(sorted(invalid))}"
        )
        return

    # Check overall structure pattern
    m = _STRUCT_REGION_RE.match(annot)
    if m is None:
        result.fail(
            f"{sid}: element_structure does not match expected pattern "
            "5+AA<+[interior]+>+BB3+"
        )
        return

    flank5, tsd5, tir_left, interior, tir_right, tsd3, flank3 = (
        m.group(1), m.group(2), m.group(3), m.group(4),
        m.group(5), m.group(6), m.group(7),
    )

    # Check TSD sequences spell "TA"
    tsd5_start = len(flank5)
    tsd5_seq = seq[tsd5_start:tsd5_start + 2]
    if tsd5_seq != "TA":
        result.fail(f"{sid}: 5' TSD (A annotation) is '{tsd5_seq}', expected 'TA'")

    tsd3_start = len(annot) - len(flank3) - 2
    tsd3_seq = seq[tsd3_start:tsd3_start + 2]
    if tsd3_seq != "TA":
        result.fail(f"{sid}: 3' TSD (B annotation) is '{tsd3_seq}', expected 'TA'")

    # Extract TIR sequences and check reverse complement
    tir_l_start = len(flank5) + 2
    tir_l_end = tir_l_start + len(tir_left)
    tir_l_seq = seq[tir_l_start:tir_l_end]

    tir_r_end = len(annot) - len(flank3) - 2
    tir_r_start = tir_r_end - len(tir_right)
    tir_r_seq = seq[tir_r_start:tir_r_end]

    rc_right = _reverse_complement(tir_r_seq)
    tir_len = min(len(tir_l_seq), len(rc_right))
    if tir_len > 0:
        hamming = sum(a != b for a, b in zip(tir_l_seq[:tir_len], rc_right[:tir_len]))
        match_pct = 1.0 - hamming / tir_len
        if match_pct < 0.70:
            result.fail(
                f"{sid}: TIR reverse-complement match is {match_pct:.0%} "
                f"(below 70% threshold)"
            )
        elif match_pct < 0.85:
            msg = f"TIR reverse-complement match: {match_pct:.0%} (below 85% warning threshold)"
            result.warn(f"{sid}: {msg}")
            seq_issues.append(msg)

    # Extract ORF (positions annotated 0, 1, 2) and introns (n)
    orf_positions = [i for i, c in enumerate(annot) if c in "012"]
    if not orf_positions:
        result.fail(f"{sid}: no ORF annotation (012) found in element_structure")
        return

    # Check codon position annotation pattern
    expected_phase = 0
    for pos in orf_positions:
        actual_phase = int(annot[pos])
        if actual_phase != expected_phase:
            result.fail(
                f"{sid}: ORF codon phase error at position {pos+1}: "
                f"expected {expected_phase}, got {actual_phase}"
            )
            return
        expected_phase = (expected_phase + 1) % 3

    orf_seq = "".join(seq[i] for i in orf_positions)

    # ORF length must be divisible by 3
    if len(orf_seq) % 3 != 0:
        result.fail(
            f"{sid}: ORF length {len(orf_seq)} is not divisible by 3"
        )
        return

    # Translate and compare to protein
    translated = _translate(orf_seq)
    if not translated:
        result.fail(f"{sid}: ORF translation is empty")
        return

    if protein_seq is not None:
        # Remove gaps from protein sequence
        prot_clean = protein_seq.replace("-", "").replace(".", "")
        identity = _pairwise_identity(translated, prot_clean)
        if identity < 0.90:
            result.fail(
                f"{sid}: ORF translation vs protein.sto identity is {identity:.1%} "
                f"(below 90% threshold). "
                f"Translated length: {len(translated)}, protein length: {len(prot_clean)}"
            )
        elif identity < 0.95:
            msg = f"ORF translation vs protein identity: {identity:.1%}"
            result.warn(f"{sid}: {msg}")
            seq_issues.append(msg)
