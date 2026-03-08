"""Score computation for ITm seed leaderboard entries.

Scoring weights and parameters are defined here as module-level constants
for easy tuning. The total score is normalized to 0-100.
"""

from __future__ import annotations

import math

import numpy as np

from . import CheckResult
from .stockholm import StoBlock

# ── Weights ──────────────────────────────────────────────────────────────────
W_DIVERSITY = 0.40
W_ANNOTATION = 0.25
W_FUNCTIONALITY = 0.20
W_SIZE = 0.15

# ── Diversity parameters ─────────────────────────────────────────────────────
REDUNDANCY_THRESHOLD = 0.95  # protein identity above this triggers penalty
REDUNDANCY_PENALTY_PER_PAIR = 0.02  # subtracted from diversity score per pair

# ── Size parameters ──────────────────────────────────────────────────────────
SIZE_TARGET = 64  # target number of sequences (log2(64) = 6)

# ── Functionality parameters ─────────────────────────────────────────────────
PARALOG_IDENTITY_STRONG = 0.99  # near-identical paralog threshold
EXPERIMENTAL_PDB_BONUS = 0.10
PREDICTED_PDB_BONUS = 0.03
LITERATURE_BONUS = 0.05
TIER_SCORES = {"1": 1.0, "2": 0.7, "3": 0.4, "4": 0.2, "MITE": 0.1}


def compute_pairwise_identity_matrix(protein_block: StoBlock) -> np.ndarray:
    """Compute all-vs-all pairwise protein identity from the alignment.

    Returns an NxN matrix of identity fractions (0-1).
    """
    seqs = list(protein_block.sequences.values())
    n = len(seqs)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            si, sj = seqs[i], seqs[j]
            length = max(len(si), len(sj))
            if length == 0:
                mat[i, j] = mat[j, i] = 0.0
                continue
            matches = sum(
                a == b and a not in ("-", ".")
                for a, b in zip(si, sj)
            )
            # Denominator: aligned non-gap positions
            aligned = sum(
                a not in ("-", ".") or b not in ("-", ".")
                for a, b in zip(si, sj)
            )
            identity = matches / aligned if aligned > 0 else 0.0
            mat[i, j] = identity
            mat[j, i] = identity
    return mat


def score_diversity(
    protein_block: StoBlock,
    families: set[str],
) -> float:
    """Compute sequence diversity score (0-1).

    Components:
    - Mean pairwise distance (1 - identity)
    - Multi-family bonus
    - Redundancy penalty
    """
    identity_mat = compute_pairwise_identity_matrix(protein_block)
    n = identity_mat.shape[0]

    if n < 2:
        return 0.0

    # Mean pairwise distance
    upper_tri = identity_mat[np.triu_indices(n, k=1)]
    mean_distance = 1.0 - np.mean(upper_tri)

    # Multi-family bonus: each additional family adds 10% (up to 50% bonus)
    n_families = len(families)
    family_bonus = min(0.5, (n_families - 1) * 0.10) if n_families > 1 else 0.0

    # Redundancy penalty
    redundant_pairs = np.sum(upper_tri > REDUNDANCY_THRESHOLD)
    penalty = redundant_pairs * REDUNDANCY_PENALTY_PER_PAIR

    score = mean_distance * (1.0 + family_bonus) - penalty
    return max(0.0, min(1.0, score))


def score_annotation(check_results: dict[str, CheckResult], provenance_rows: list[dict]) -> float:
    """Compute annotation completeness score (0-1).

    Full marks if all checks pass; partial credit for warnings.
    Bonuses for optional metadata.
    """
    # Base score from check pass/warn/fail
    check_categories = ["format", "annotation", "protein", "dna", "provenance"]
    total_checks = len(check_categories)
    passed = 0
    warned = 0
    for cat in check_categories:
        cr = check_results.get(cat)
        if cr is None:
            continue
        if cr.status == "pass":
            passed += 1
        elif cr.status == "warn":
            warned += 1

    base = (passed + 0.5 * warned) / total_checks if total_checks > 0 else 0.0

    # Optional metadata bonuses
    bonus = 0.0
    for row in provenance_rows:
        if (row.get("paralog_hits") or "").strip():
            bonus += 0.01
        if (row.get("confidence_tier") or "").strip():
            bonus += 0.005
        if (row.get("pdb_file") or "").strip():
            bonus += 0.01

    return min(1.0, base + bonus)


def score_functionality(provenance_rows: list[dict]) -> float:
    """Compute evidence of functionality score (0-1).

    Rewards near-identical paralogs, experimental structures, literature references.
    """
    n = len(provenance_rows)
    if n == 0:
        return 0.0

    total = 0.0
    for row in provenance_rows:
        seq_score = 0.0

        # Near-identical paralogs
        try:
            max_id = float(row.get("max_paralog_identity") or "0")
            if max_id >= PARALOG_IDENTITY_STRONG:
                seq_score += 0.40
            elif max_id >= 0.95:
                seq_score += 0.15
        except (ValueError, TypeError):
            pass

        # PDB structures
        pdb_source = (row.get("pdb_source") or "").strip().lower()
        if pdb_source == "experimental":
            seq_score += EXPERIMENTAL_PDB_BONUS
        elif pdb_source in ("alphafold", "esmfold", "colabfold"):
            seq_score += PREDICTED_PDB_BONUS

        # Literature reference
        ref = (row.get("reference") or "").strip()
        if ref and ("10." in ref or "PMC" in ref or "PMID" in ref):
            seq_score += LITERATURE_BONUS

        # Confidence tier
        tier = (row.get("confidence_tier") or "").strip()
        seq_score += TIER_SCORES.get(tier, 0.0) * 0.10

        total += min(1.0, seq_score)

    return min(1.0, total / n)


def score_size(n_sequences: int) -> float:
    """Compute collection size score (0-1).

    Logarithmic scaling: score = min(1, log2(N) / log2(target)).
    """
    if n_sequences < 3 or n_sequences > 300:
        return 0.0
    if n_sequences <= 1:
        return 0.0
    return min(1.0, math.log2(n_sequences) / math.log2(SIZE_TARGET))


def compute_total_score(
    protein_block: StoBlock,
    families: set[str],
    check_results: dict[str, CheckResult],
    provenance_rows: list[dict],
    n_sequences: int,
) -> dict[str, float]:
    """Compute all score components and total.

    Returns dict with keys: diversity, annotation, functionality, size, total.
    """
    # If any hard-fail check, total is 0
    hard_fail_categories = ["format", "annotation"]
    for cat in hard_fail_categories:
        cr = check_results.get(cat)
        if cr and cr.status == "fail":
            return {
                "diversity": 0.0,
                "annotation": 0.0,
                "functionality": 0.0,
                "size": 0.0,
                "total": 0.0,
            }

    div = score_diversity(protein_block, families)
    ann = score_annotation(check_results, provenance_rows)
    func = score_functionality(provenance_rows)
    size = score_size(n_sequences)

    total = (
        W_DIVERSITY * div
        + W_ANNOTATION * ann
        + W_FUNCTIONALITY * func
        + W_SIZE * size
    ) * 100.0

    return {
        "diversity": round(div, 4),
        "annotation": round(ann, 4),
        "functionality": round(func, 4),
        "size": round(size, 4),
        "total": round(total, 1),
    }
