"""Stockholm format parsing for protein.sto and dna.sto files."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from . import CheckResult


@dataclass
class StoBlock:
    """A single Stockholm alignment block."""
    gf: dict[str, str] = field(default_factory=dict)  # #=GF lines
    sequences: dict[str, str] = field(default_factory=dict)  # seqid -> sequence
    gc: dict[str, str] = field(default_factory=dict)  # #=GC feature -> annotation

    @property
    def seq_ids(self) -> list[str]:
        return list(self.sequences.keys())

    def strip_id(self, raw_id: str) -> str:
        """Strip coordinate range suffix (e.g. '/1-345') from a Stockholm ID."""
        return raw_id.split("/")[0]


_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,80}$")


def _parse_block_lines(lines: list[str]) -> StoBlock:
    """Parse lines (excluding header and //) into a StoBlock."""
    block = StoBlock()
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith("# "):
            continue
        if line.startswith("#=GF"):
            parts = line.split(None, 2)
            if len(parts) >= 3:
                block.gf[parts[1]] = parts[2]
        elif line.startswith("#=GC"):
            parts = line.split(None, 2)
            if len(parts) >= 3:
                block.gc[parts[1]] = parts[2].strip()
        elif line.startswith("#"):
            continue  # other comment lines
        else:
            parts = line.split()
            if len(parts) >= 2:
                raw_id = parts[0]
                seq = parts[1]
                stripped = block.strip_id(raw_id)
                if stripped in block.sequences:
                    # Interleaved: append
                    block.sequences[stripped] += seq
                else:
                    block.sequences[stripped] = seq
    return block


def parse_protein_sto(path: Path, result: CheckResult) -> StoBlock | None:
    """Parse a single-block Stockholm protein alignment.

    Returns the parsed block, or None on hard failure.
    """
    try:
        text = path.read_text()
    except (OSError, UnicodeDecodeError) as e:
        result.fail(f"Cannot read {path.name}: {e}")
        return None

    lines = text.splitlines()
    if not lines or not lines[0].strip().startswith("# STOCKHOLM 1.0"):
        result.fail(f"{path.name}: missing '# STOCKHOLM 1.0' header")
        return None

    # Find the // terminator
    term_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "//":
            term_idx = i
            break
    if term_idx is None:
        result.fail(f"{path.name}: missing '//' terminator")
        return None

    block = _parse_block_lines(lines[1:term_idx])

    if not block.sequences:
        result.fail(f"{path.name}: no sequences found")
        return None

    if "catalytic_triad" not in block.gc:
        result.fail(f"{path.name}: missing #=GC catalytic_triad annotation line")
        return None

    # Validate sequence IDs
    for sid in block.seq_ids:
        if not _ID_RE.match(sid):
            result.fail(
                f"{path.name}: invalid sequence ID '{sid}' "
                "(must be alphanumeric, underscores, hyphens; max 80 chars)"
            )

    # Validate alignment: all sequences and GC annotations must be same length
    lengths = {sid: len(seq) for sid, seq in block.sequences.items()}
    gc_lengths = {feat: len(ann) for feat, ann in block.gc.items()}
    all_lengths = {**lengths, **gc_lengths}
    unique_lens = set(all_lengths.values())
    if len(unique_lens) > 1:
        result.fail(
            f"{path.name}: alignment columns have inconsistent lengths: "
            + ", ".join(f"{k}={v}" for k, v in all_lengths.items())
        )
        return None

    return block


def iter_dna_blocks(path: Path, result: CheckResult) -> Iterator[StoBlock]:
    """Parse a multi-Stockholm dna.sto file, yielding one StoBlock per element.

    Adds errors/warnings to result for malformed blocks.
    """
    try:
        text = path.read_text()
    except (OSError, UnicodeDecodeError) as e:
        result.fail(f"Cannot read {path.name}: {e}")
        return

    # Split on '//' lines
    raw_blocks = re.split(r"^//\s*$", text, flags=re.MULTILINE)

    block_count = 0
    for raw in raw_blocks:
        lines = raw.strip().splitlines()
        if not lines:
            continue

        # Check header
        if not lines[0].strip().startswith("# STOCKHOLM 1.0"):
            # Could be trailing whitespace after last //
            if any(line.strip() for line in lines):
                result.fail(f"{path.name}: block {block_count+1} missing '# STOCKHOLM 1.0' header")
            continue

        block = _parse_block_lines(lines[1:])
        block_count += 1

        if len(block.sequences) != 1:
            result.fail(
                f"{path.name} block {block_count}: expected exactly 1 sequence, "
                f"found {len(block.sequences)}"
            )
            continue

        if "element_structure" not in block.gc:
            result.fail(
                f"{path.name} block {block_count}: missing #=GC element_structure annotation"
            )
            continue

        sid = block.seq_ids[0]
        if not _ID_RE.match(sid):
            result.fail(
                f"{path.name} block {block_count}: invalid sequence ID '{sid}'"
            )

        # Check length consistency
        seq_len = len(block.sequences[sid])
        annot_len = len(block.gc["element_structure"])
        if seq_len != annot_len:
            result.fail(
                f"{path.name} block {block_count} ({sid}): sequence length ({seq_len}) "
                f"!= element_structure length ({annot_len})"
            )
            continue

        yield block

    if block_count == 0:
        result.fail(f"{path.name}: no valid Stockholm blocks found")
