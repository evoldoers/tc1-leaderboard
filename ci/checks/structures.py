"""PDB structure file validation."""

from __future__ import annotations

from pathlib import Path

from . import CheckResult


def check_structures(
    entry_dir: Path,
    provenance_rows: list[dict],
    protein_lengths: dict[str, int],
    result: CheckResult,
) -> None:
    """Validate PDB files referenced in provenance.tsv.

    Args:
        entry_dir: Path to the entry directory.
        provenance_rows: Parsed provenance rows.
        protein_lengths: Dict of seq_id -> ungapped protein length.
        result: CheckResult to accumulate.
    """
    pdb_count = 0
    for row in provenance_rows:
        pdb_file = row.get("pdb_file", "").strip()
        if not pdb_file:
            continue

        pdb_path = entry_dir / pdb_file
        if not pdb_path.exists():
            result.warn(f"PDB file '{pdb_file}' referenced in provenance.tsv does not exist")
            continue

        pdb_count += 1
        sid = row.get("id", "???")

        try:
            text = pdb_path.read_text()
        except (OSError, UnicodeDecodeError) as e:
            result.warn(f"Cannot read PDB file '{pdb_file}': {e}")
            continue

        # Basic PDB format check: should have ATOM or HETATM records
        has_atom = any(
            line.startswith(("ATOM", "HETATM")) for line in text.splitlines()
        )
        # Also accept mmCIF: _atom_site.
        has_mmcif = "_atom_site." in text

        if not has_atom and not has_mmcif:
            result.warn(f"'{pdb_file}': not parseable as PDB or mmCIF (no ATOM records found)")
            continue

        # Chain length check (rough): count unique residue numbers from ATOM records
        if has_atom:
            residue_ids = set()
            for line in text.splitlines():
                if line.startswith("ATOM"):
                    # PDB format: residue number at columns 22-26
                    try:
                        resnum = int(line[22:26].strip())
                        residue_ids.add(resnum)
                    except (ValueError, IndexError):
                        pass

            if residue_ids:
                chain_len = len(residue_ids)
                expected_len = protein_lengths.get(sid, 0)
                if expected_len > 0:
                    ratio = chain_len / expected_len
                    if ratio < 0.5 or ratio > 2.0:
                        result.warn(
                            f"'{pdb_file}': chain length {chain_len} residues vs "
                            f"protein length {expected_len} aa (ratio {ratio:.2f})"
                        )

        # Check pLDDT for predicted structures
        pdb_source = row.get("pdb_source", "").strip().lower()
        if pdb_source in ("alphafold", "esmfold", "colabfold") and has_atom:
            # B-factor column should contain pLDDT values (0-100)
            bfactors = []
            for line in text.splitlines():
                if line.startswith("ATOM"):
                    try:
                        bfactor = float(line[60:66].strip())
                        bfactors.append(bfactor)
                    except (ValueError, IndexError):
                        pass
            if bfactors:
                mean_bfactor = sum(bfactors) / len(bfactors)
                if mean_bfactor < 1.0 or mean_bfactor > 100.0:
                    result.warn(
                        f"'{pdb_file}': source is {pdb_source} but B-factor range "
                        f"({min(bfactors):.1f}-{max(bfactors):.1f}) doesn't look like pLDDT"
                    )

    if pdb_count > 0:
        result.info(f"{pdb_count} PDB file(s) validated")
