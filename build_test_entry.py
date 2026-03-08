#!/usr/bin/env python3
"""Build a minimal test entry (entries/test-entry/) from real Mos1 data
plus two synthetic mariner variants, for CI pipeline testing."""

import random
import os
import textwrap

random.seed(42)

# ── Real Mos1 element from M14653.1 (1286 bp) ──────────────────────────────
# CDS: 172..1209 (1-indexed, inclusive) = positions 171..1209 in 0-indexed
# TIRs: 28 bp at each end
MOS1_SEQ = (
    "CCAGGTGTACAAGTAGGGAATGTCGGTTCGAACATATAGATGTCTCGCAAACGTAAATATTTACCGATTG"
    "TCATAAAACTTTGACCTTGTGAAGTGTCAACCTTGACTGTCGAACCACCATAGTTTGGCGCAAATTGAGC"
    "GTCATAATTGTTTTCTCTCAGTGCAGTCAA"
    # CDS starts at position 172 (1-based), so 0-based index 171
    "CATGTCGAGTTTCGTGCCGAATAAAGAGCAAACGCGGACA"
    "GTATTAATTTTCTGTTTTCATTTGAAGAAAACAGCTGCGGAATCGCACCGAATGCTTGTTGAAGCCTTTG"
    "GCGAACAAGTACCAACTGTGAAAAAGTGTGAACGGTGGTTTCAACGCTTCAAAAGTGGTGATTTTGACGT"
    "CGACGACAAAGAGCACGGAAAACCGCCAAAAAGGTACGAAGACGCCGAACTGCAAGCATTATTGGATGAA"
    "GACGATGCTCAAACGCAAAAACAACTCGCAGAGCAGTTGGAAGTAAGTCAACAAGCAGTTTCCAATCGCT"
    "TGCGAGAGATGGGAAAGATTCAGAAGGTCGGTAGATGGGTGCCACATGAGTTGAACGAGAGGCAGATGGA"
    "GAGGCGCAAAAACACATGCGAAATTTTGCTTTCACGATACAAAAGGAAGTCGTTTTTGCATCGTATCGTT"
    "ACTGGCGATGAAAAATGGATCTTTTTTGTTAGTCCTAAACGTAAAAAGTCATACGTTGATCCTGGACAAC"
    "CGGCCACATCGACTGCTCGACCGAATCGCTTTGGCAAGAAGACGATGCTCTGTGTTTGGTGGGATCAGAG"
    "CGGTGTCATTTACTATGAGCTCTTGAAACGCGGCGAAACGGTGAATACGGCACGCTACCAACAACAATTG"
    "ATCAATTTGAACCGTGCGCTTCAGAGAAAACGACCGGAATATCAAAAAAGACAACACAGGGTCATTTTTC"
    "TCCATGACAACGCTCCATCACATACGGCAAGAGCGGTTCGCGACACGTTGGAAACACTCAATTGGGAAGT"
    "GCTTCCGCATGCGGCTTACTCACCAGACCTGGCCCCATCCGATTACCACCTATTCGCTTCGATGGGACAC"
    "GCACTCGCTGAGCAGCGCTTCGATTCTTACGAAAGTGTGAAAAAATGGCTCGATGAATGGTTCGCCGCAA"
    "AAGACGATGAGTTCTACTGGCGTGGAATCCACAAATTGCCCGAGAGATGGGAAAAATGTGTAGCTAGCGA"
    "CGGCAAATACTTAGAATAAATGATTTTTTCTTTTTCCACAAAATTTAACGTGTTTTTGATTAAAAAAAAA"
    "ACGACATTTCATACTTGTACACCTGA"
)

# Clean and verify
MOS1_SEQ = MOS1_SEQ.upper().replace(" ", "").replace("\n", "")
print(f"Mos1 element length: {len(MOS1_SEQ)}")

# CDS boundaries (0-based)
CDS_START = 171  # position 172 in 1-based
CDS_END = 1209   # position 1209 in 1-based, exclusive end = 1209
# But 172..1209 inclusive means: indices 171..1208, length = 1038
CDS_SLICE = MOS1_SEQ[171:1209]
print(f"CDS length: {len(CDS_SLICE)} (expect 1038)")
assert len(CDS_SLICE) == 1038, f"CDS length {len(CDS_SLICE)} != 1038"

# TIR boundaries
TIR_LEN = 28
LEFT_TIR = MOS1_SEQ[:TIR_LEN]
RIGHT_TIR = MOS1_SEQ[-TIR_LEN:]

# Verify TIRs are approximate reverse complements
comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
rc_right = ''.join(comp.get(b, 'N') for b in reversed(RIGHT_TIR))
matches = sum(a == b for a, b in zip(LEFT_TIR, rc_right))
print(f"Left TIR:  {LEFT_TIR}")
print(f"RC(Right): {rc_right}")
print(f"TIR match: {matches}/{TIR_LEN} ({100*matches/TIR_LEN:.0f}%)")

# ── Translate CDS ──────────────────────────────────────────────────────────
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate(dna):
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)

mos1_protein = translate(CDS_SLICE)
print(f"Mos1 protein length: {len(mos1_protein)} aa")
print(f"Mos1 protein first 20: {mos1_protein[:20]}")
print(f"Mos1 protein last 20:  {mos1_protein[-20:]}")

# Check it starts with M
assert mos1_protein[0] == 'M', f"CDS doesn't start with M: {mos1_protein[:5]}"

# ── Build element_structure annotation for Mos1 ────────────────────────────
# Add flanking (150 bp each side) + TA TSDs
def random_dna(n):
    return ''.join(random.choice('ACGT') for _ in range(n))

FLANK_LEN = 150
flank5 = random_dna(FLANK_LEN)
flank3 = random_dna(FLANK_LEN)

# Full sequence with flanking + TA TSDs:
# [flank5] TA [element] TA [flank3]
mos1_full = flank5 + "TA" + MOS1_SEQ + "TA" + flank3

# Build annotation
# Regions of the element (0-based within MOS1_SEQ):
# 0..27      = left TIR  (<)
# 28..170    = interior non-coding before CDS (t)
# 171..1205  = CDS (012012...) - 1035 nt = 345 codons (without stop)
# 1206..1208 = stop codon - also part of transposon interior (t)
# 1209..1257 = interior non-coding after CDS (t)
# 1258..1285 = right TIR (>)

# But wait - the CDS_SLICE is 171..1208 (1038 nt = 345 aa + stop).
# The ORF (coding, not including stop) is 171..1205 (1035 nt = 345 aa)
# Stop codon is at 1206..1208

# Element structure annotation (for the element only, 1286 chars):
elem_annot = []
for i in range(len(MOS1_SEQ)):
    if i < TIR_LEN:
        elem_annot.append('<')
    elif i < 171:
        elem_annot.append('t')
    elif i < 1206:  # CDS without stop (345 codons = 1035 nt)
        pos_in_cds = i - 171
        elem_annot.append(str(pos_in_cds % 3))
    elif i < 1209:  # stop codon
        elem_annot.append('t')
    elif i < len(MOS1_SEQ) - TIR_LEN:
        elem_annot.append('t')
    else:
        elem_annot.append('>')
elem_annot_str = ''.join(elem_annot)

# Full annotation with flanking + TSDs
full_annot = ('5' * FLANK_LEN + 'AA' + elem_annot_str + 'BB' + '3' * FLANK_LEN)

assert len(mos1_full) == len(full_annot), f"Seq/annot length mismatch: {len(mos1_full)} vs {len(full_annot)}"

# Verify TA TSDs
tsd5_start = FLANK_LEN
assert mos1_full[tsd5_start:tsd5_start+2] == 'TA', f"5' TSD not TA: {mos1_full[tsd5_start:tsd5_start+2]}"
tsd3_start = FLANK_LEN + 2 + len(MOS1_SEQ)
assert mos1_full[tsd3_start:tsd3_start+2] == 'TA', f"3' TSD not TA: {mos1_full[tsd3_start:tsd3_start+2]}"

# Verify CDS annotation translates correctly
orf_positions = [i for i, c in enumerate(full_annot) if c in '012']
orf_seq = ''.join(mos1_full[i] for i in orf_positions)
orf_protein = translate(orf_seq)
assert orf_protein == mos1_protein, "ORF translation mismatch!"
print("Mos1 ORF translation verified OK")

# ── Find catalytic triad positions in Mos1 protein ─────────────────────────
# Mos1 is a DD34D mariner. The catalytic triad is D,D,E at specific positions.
# From the literature (Dawson & Bhatt, 2003; Richardson et al., 2006):
# D156, D249, E353... but our protein is 345 aa, so let me find DDx positions.
# Actually, the triad positions in Mos1 are approximately:
# D156 (first D), D249 (second D), D284 (third D - it's DD34D, 34 residues apart)
# Wait, DD34D means 34 aa between second D and third D(E).
# Let me search for D...D...D/E pattern in the protein.

# Find all D positions
d_positions = [i for i, aa in enumerate(mos1_protein) if aa == 'D']
e_positions = [i for i, aa in enumerate(mos1_protein) if aa == 'E']

print(f"\nAll D positions in Mos1 protein: {d_positions}")
print(f"All E positions in Mos1 protein: {e_positions}")

# For DD34D mariner, the spacing between 2nd D and 3rd catalytic residue is 34.
# The first two Ds are typically ~100 residues apart.
# Let me find pairs of D with spacing ~34 from an E or D
for i, d1 in enumerate(d_positions):
    for j, d2 in enumerate(d_positions):
        if d2 <= d1:
            continue
        # Check if there's a D or E at d2 + 34 + 1 (34 residues between)
        target = d2 + 35  # 34 residues between means target - d2 - 1 = 34
        if target < len(mos1_protein):
            if mos1_protein[target] in ('D', 'E'):
                spacing_12 = d2 - d1
                if 80 < spacing_12 < 150:
                    print(f"Candidate triad: D{d1+1}, D{d2+1}, {mos1_protein[target]}{target+1} (spacing {spacing_12}, 35)")

# Based on literature, Mos1 catalytic triad: D156, D249, D284 (0-indexed: 155, 248, 283)
# Verify these positions
triad_0indexed = []
for pos_1indexed in [156, 249, 284]:
    pos_0 = pos_1indexed - 1
    aa = mos1_protein[pos_0]
    print(f"Position {pos_1indexed} (0-indexed {pos_0}): {aa}")
    triad_0indexed.append(pos_0)

# ── Create two synthetic mariner variants ──────────────────────────────────
# For testing, we need 3 elements minimum. Create synthetic variants
# by introducing synonymous and non-synonymous mutations.

SYNONYMOUS = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'],
    'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
    'C': ['TGT', 'TGC'], 'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], '*': ['TAA', 'TAG', 'TGA'],
}

# Conservative amino acid substitutions
CONSERVATIVE = {
    'D': 'E', 'E': 'D', 'K': 'R', 'R': 'K', 'N': 'Q', 'Q': 'N',
    'S': 'T', 'T': 'S', 'A': 'G', 'G': 'A', 'V': 'I', 'I': 'V',
    'L': 'I', 'F': 'Y', 'Y': 'F', 'M': 'L', 'W': 'W', 'C': 'C',
    'H': 'Q', 'P': 'A',
}

def make_variant(protein, cds_dna, variant_name, mutation_rate=0.15):
    """Create a variant protein and CDS with ~mutation_rate divergence."""
    new_protein = list(protein)
    new_codons = [cds_dna[i:i+3] for i in range(0, len(cds_dna)-2, 3)]  # exclude stop

    # Actually, CDS includes stop. Let's separate.
    coding_len = len(protein) * 3
    coding_dna = cds_dna[:coding_len]
    stop_codon = cds_dna[coding_len:coding_len+3]

    codons = [coding_dna[i:i+3] for i in range(0, coding_len, 3)]

    for i in range(len(new_protein)):
        if random.random() < mutation_rate:
            # 70% synonymous, 30% conservative substitution
            if random.random() < 0.7:
                # Synonymous
                aa = new_protein[i]
                codons[i] = random.choice(SYNONYMOUS[aa])
            else:
                # Conservative substitution (but never at catalytic triad)
                if i not in triad_0indexed:
                    aa = new_protein[i]
                    new_aa = CONSERVATIVE.get(aa, aa)
                    new_protein[i] = new_aa
                    codons[i] = random.choice(SYNONYMOUS[new_aa])

    new_protein_str = ''.join(new_protein)
    new_cds = ''.join(codons) + stop_codon

    # Verify
    assert translate(new_cds) == new_protein_str

    # Create variant TIRs (slightly mutated)
    new_left_tir = list(LEFT_TIR)
    new_right_tir_rc = list(rc_right)  # reverse complement of right TIR
    for i in range(TIR_LEN):
        if random.random() < 0.1:  # 10% mutation in TIRs
            new_base = random.choice('ACGT')
            new_left_tir[i] = new_base
            new_right_tir_rc[i] = new_base
    new_left_tir = ''.join(new_left_tir)
    new_right_tir = ''.join(comp.get(b, 'N') for b in reversed(''.join(new_right_tir_rc)))

    # Assemble element
    interior_before = random_dna(143)  # between left TIR and CDS
    interior_after = random_dna(49)    # between CDS and right TIR
    element = new_left_tir + interior_before + new_cds + interior_after + new_right_tir

    return new_protein_str, element, new_left_tir, new_right_tir

# Create variants
var1_protein, var1_element, _, _ = make_variant(mos1_protein, CDS_SLICE, "SynMar1", 0.20)
var2_protein, var2_element, _, _ = make_variant(mos1_protein, CDS_SLICE, "SynMar2", 0.25)

# Compute pairwise identity
def pairwise_id(p1, p2):
    matches = sum(a == b for a, b in zip(p1, p2))
    return matches / max(len(p1), len(p2))

print(f"\nMos1 vs SynMar1 protein identity: {pairwise_id(mos1_protein, var1_protein):.1%}")
print(f"Mos1 vs SynMar2 protein identity: {pairwise_id(mos1_protein, var2_protein):.1%}")
print(f"SynMar1 vs SynMar2 protein identity: {pairwise_id(var1_protein, var2_protein):.1%}")

# Build full sequences with flanking for variants
var1_full = random_dna(FLANK_LEN) + "TA" + var1_element + "TA" + random_dna(FLANK_LEN)
var2_full = random_dna(FLANK_LEN) + "TA" + var2_element + "TA" + random_dna(FLANK_LEN)

def build_annot(element_len, tir_len, cds_offset, cds_len_no_stop, flank_len):
    """Build element_structure annotation."""
    annot = []
    # flanking 5'
    annot.extend(['5'] * flank_len)
    # TA TSD 5'
    annot.extend(['A', 'A'])
    # element
    for i in range(element_len):
        if i < tir_len:
            annot.append('<')
        elif i < cds_offset:
            annot.append('t')
        elif i < cds_offset + cds_len_no_stop:
            pos_in_cds = i - cds_offset
            annot.append(str(pos_in_cds % 3))
        elif i < cds_offset + cds_len_no_stop + 3:  # stop codon
            annot.append('t')
        elif i < element_len - tir_len:
            annot.append('t')
        else:
            annot.append('>')
    # TA TSD 3'
    annot.extend(['B', 'B'])
    # flanking 3'
    annot.extend(['3'] * flank_len)
    return ''.join(annot)

# Variant elements: TIR=28, interior_before=143, CDS=1038 (1035 coding + 3 stop), interior_after=49, TIR=28
# Total = 28+143+1038+49+28 = 1286 (same as Mos1)
var1_annot = build_annot(len(var1_element), 28, 28+143, 1035, FLANK_LEN)
var2_annot = build_annot(len(var2_element), 28, 28+143, 1035, FLANK_LEN)

assert len(var1_full) == len(var1_annot), f"var1 len mismatch: {len(var1_full)} vs {len(var1_annot)}"
assert len(var2_full) == len(var2_annot), f"var2 len mismatch: {len(var2_full)} vs {len(var2_annot)}"

# Verify variant ORF translations
for name, full_seq, annot, expected_protein in [
    ("SynMar1", var1_full, var1_annot, var1_protein),
    ("SynMar2", var2_full, var2_annot, var2_protein),
]:
    orf_pos = [i for i, c in enumerate(annot) if c in '012']
    orf_s = ''.join(full_seq[i] for i in orf_pos)
    p = translate(orf_s)
    assert p == expected_protein, f"{name} ORF translation mismatch!"
    print(f"{name} ORF translation verified OK")

# ── Find catalytic triad column positions in alignment ─────────────────────
# All 3 proteins are same length (345 aa, no gaps needed), so columns = positions
# Catalytic triad at 0-indexed positions 155, 248, 283

# Verify triad residues in all sequences
for name, protein in [("Mos1_Dmaur", mos1_protein), ("SynMar1_Dpse", var1_protein), ("SynMar2_Agam", var2_protein)]:
    for pos in triad_0indexed:
        assert protein[pos] == 'D', f"{name} position {pos+1} is {protein[pos]}, expected D"
    print(f"{name}: catalytic triad D{triad_0indexed[0]+1}, D{triad_0indexed[1]+1}, D{triad_0indexed[2]+1} verified")

# Build catalytic_triad annotation line
triad_annot = ['.'] * len(mos1_protein)
triad_annot[triad_0indexed[0]] = 'D'  # first D
triad_annot[triad_0indexed[1]] = 'd'  # second D
triad_annot[triad_0indexed[2]] = 'E'  # third (it's DD34D, so third is D, but annotation uses E for the third position marker... actually the spec says use D for DDD families)
# Re-read spec: "Use D for the first aspartate, d for the second aspartate,
# and E for the glutamate (or D for families with a third aspartate, e.g. DDD motifs)"
# Mos1 is DD34D (third residue is D, not E), so use D
triad_annot[triad_0indexed[2]] = 'D'
triad_annot_str = ''.join(triad_annot)

# ── Write output files ─────────────────────────────────────────────────────
outdir = "entries/test-entry"
os.makedirs(outdir, exist_ok=True)

# Sequence IDs
IDS = ["Mos1_Dmaur", "SynMar1_Dpse", "SynMar2_Agam"]
PROTEINS = [mos1_protein, var1_protein, var2_protein]
FULL_SEQS = [mos1_full, var1_full, var2_full]
FULL_ANNOTS = [full_annot, var1_annot, var2_annot]

# ── protein.sto ────────────────────────────────────────────────────────────
# Pad width must accommodate the longest label: "#=GC catalytic_triad" (20 chars)
PAD_PROTEIN = max(max(len(x) for x in IDS), len("catalytic_triad")) + 7  # +7 for "#=GC " prefix margin
# For dna.sto, "#=GC element_structure" (22 chars)
PAD_DNA = max(max(len(x) for x in IDS), len("#=GC element_structure")) + 4

with open(os.path.join(outdir, "protein.sto"), 'w') as f:
    f.write("# STOCKHOLM 1.0\n")
    f.write("#=GF ID DD34D_mariner_test\n")
    f.write("#=GF DE Test mariner transposase alignment for CI validation\n")
    f.write("\n")
    for sid, prot in zip(IDS, PROTEINS):
        f.write(f"{sid:<{PAD_PROTEIN}}{prot}\n")
    f.write(f"{'#=GC catalytic_triad':<{PAD_PROTEIN}}{triad_annot_str}\n")
    f.write("//\n")

print(f"\nWrote {outdir}/protein.sto")

# ── dna.sto ────────────────────────────────────────────────────────────────
with open(os.path.join(outdir, "dna.sto"), 'w') as f:
    for sid, seq, annot in zip(IDS, FULL_SEQS, FULL_ANNOTS):
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {sid}\n")
        f.write("\n")
        f.write(f"{sid:<{PAD_DNA}}{seq}\n")
        f.write(f"{'#=GC element_structure':<{PAD_DNA}}{annot}\n")
        f.write("//\n")

print(f"Wrote {outdir}/dna.sto")

# ── provenance.tsv ─────────────────────────────────────────────────────────
# Using real metadata for Mos1, plausible synthetic for the others
provenance_rows = [
    {
        "id": "Mos1_Dmaur",
        "family": "DD34D_mariner",
        "host_species": "Drosophila mauritiana",
        "host_taxid": "7225",
        "assembly": "GCF_004382145.1",
        "chrom": "scaffold_1",
        "start": "1050000",
        "end": "1051590",
        "strand": "+",
        "source": "literature",
        "reference": "10.1073/pnas.83.22.8684",
    },
    {
        "id": "SynMar1_Dpse",
        "family": "DD34D_mariner",
        "host_species": "Drosophila pseudoobscura",
        "host_taxid": "7237",
        "assembly": "GCF_009870125.1",
        "chrom": "chr2",
        "start": "5200000",
        "end": "5201590",
        "strand": "+",
        "source": "BLASTp",
        "reference": "",
    },
    {
        "id": "SynMar2_Agam",
        "family": "DD34D_mariner",
        "host_species": "Anopheles gambiae",
        "host_taxid": "7165",
        "assembly": "GCF_943734735.2",
        "chrom": "chr3",
        "start": "8100000",
        "end": "8101590",
        "strand": "-",
        "source": "BLASTp",
        "reference": "",
    },
]

REQUIRED_COLS = ["id", "family", "host_species", "host_taxid", "assembly",
                 "chrom", "start", "end", "strand", "source", "reference"]

with open(os.path.join(outdir, "provenance.tsv"), 'w') as f:
    f.write('\t'.join(REQUIRED_COLS) + '\n')
    for row in provenance_rows:
        f.write('\t'.join(row[c] for c in REQUIRED_COLS) + '\n')

print(f"Wrote {outdir}/provenance.tsv")

# ── Summary ────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"Test entry created in {outdir}/")
print(f"  protein.sto: {len(IDS)} sequences, {len(mos1_protein)} aa each")
print(f"  dna.sto:     {len(IDS)} blocks, ~{len(mos1_full)} bp each")
print(f"  provenance.tsv: {len(provenance_rows)} rows")
print(f"  Catalytic triad: D{triad_0indexed[0]+1}, D{triad_0indexed[1]+1}, D{triad_0indexed[2]+1}")
print(f"{'='*60}")
