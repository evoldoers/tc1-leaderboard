# Scoring Policy — ITm Seed Dataset Leaderboard

This document describes the current scoring scheme and validation rules. It is the authoritative reference for what the CI pipeline checks and how scores are computed.

## Hard-fail checks (score = 0 if any fail)

### Format checks
- `protein.sto` must be valid Stockholm 1.0 with a `#=GC catalytic_triad` line.
- `dna.sto` must be valid multi-Stockholm with exactly one sequence per block and a `#=GC element_structure` line per block.
- `provenance.tsv` must be valid TSV with all required columns: `id`, `family`, `host_species`, `host_taxid`, `assembly`, `chrom`, `start`, `end`, `strand`, `source`, `reference`.
- All IDs must be consistent across `protein.sto`, `dna.sto`, and `provenance.tsv`.
- Minimum 3 sequences, maximum 300.

### Annotation consistency checks
- `catalytic_triad` must mark exactly 2–3 positions with D, d, and E (or D for DDD families). Marked columns must contain the expected residue in the majority of sequences.
- `element_structure` must:
  - Use only valid characters: `5 3 < > A B 0 1 2 n t .`
  - Follow the pattern: `5+AA<+[interior]+>+BB3+`
  - Have `A` and `B` positions (TSD) spelling `TA` in the sequence.
  - Have TIRs (`<` and `>`) that are approximate reverse complements (>70% match).
  - Have ORF annotation (`012`) with correct codon phasing and length divisible by 3.
  - Have ORF translation matching the protein in `protein.sto` at >=90% identity.

## Scored checks (reduce score but do not reject)

### Protein-level
- Transposase length within expected range for declared family.
- Catalytic triad spacing matches declared family (e.g. DD34D = 34 residues between 2nd and 3rd).

### DNA-level
- Element length (excluding flanking) within expected range.
- TIR length within expected range.
- ORF GC content between 15% and 75%.
- No ambiguous bases (N) in the ORF.

### Provenance
- `host_taxid` is a valid positive integer.
- `assembly` matches GC[AF]_NNNNNNNNN.N format.
- Coordinates are consistent (`start` < `end`, `strand` is `+` or `-`).

## Scoring weights

| Component | Weight | What it rewards |
|-----------|--------|-----------------|
| Sequence diversity | 40% | Phylogenetic breadth, multi-family coverage |
| Annotation completeness | 25% | Passing all checks, optional metadata |
| Evidence of functionality | 20% | Near-identical paralogs, experimental data, literature |
| Collection size | 15% | Enough sequences to seed a model (target: 64) |

## Scoring details

### Diversity (40%)
- Mean pairwise protein distance (1 − identity) across all sequence pairs.
- Multi-family bonus: +10% per additional family (up to +50%).
- Redundancy penalty: −2% per pair with >95% protein identity.

### Annotation completeness (25%)
- Base score from check pass/warn/fail rates.
- Bonuses for optional metadata: `paralog_hits` (+1% per entry), `confidence_tier` (+0.5% per entry), `pdb_file` (+1% per entry).

### Functionality (20%)
- Per-sequence scoring, averaged:
  - `max_paralog_identity` >= 0.99: +40%
  - `max_paralog_identity` >= 0.95: +15%
  - Experimental PDB: +10%
  - Predicted PDB: +3%
  - Literature reference with DOI: +5%
  - Confidence tier bonus (tier 1 = +10%, tier 2 = +7%, tier 3 = +4%, tier 4 = +2%)

### Size (15%)
- `min(1, log2(N) / log2(64))`
- Saturates at 64 sequences. Hard reject if N < 3 or N > 300.

## Tie-breaking

Entries with equal total scores are ranked by annotation completeness, then by diversity.

## Policy changes

This policy may be updated. The leaderboard always reflects the current policy. When the policy changes, all entries are re-scored on the next CI run.
