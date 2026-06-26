'''Decision-grade sequence property calculations: the unit and default traps'''
# Reference: biopython 1.83+ | Verify API if version differs
import inspect
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight, MeltingTemp
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# GC content: fraction vs percent, and the ambiguity-default trap
print('=== GC Content ===')
seq = Seq('ATGCGATCGNNNCGATCGATCG')
print(f'gc_fraction (fraction 0-1): {gc_fraction(seq):.4f}')
print(f'gc_fraction * 100 (percent): {gc_fraction(seq) * 100:.2f}')
print(f"ambiguous='remove' (default): {gc_fraction(seq, ambiguous='remove') * 100:.2f}")
print(f"ambiguous='ignore' (legacy GC drop-in): {gc_fraction(seq, ambiguous='ignore') * 100:.2f}")

# Molecular weight: single-strand default vs full duplex
print('\n=== Molecular Weight ===')
dna = Seq('ATGCGATCGATCGATCGATCG')
mw_ss = molecular_weight(dna)
mw_ds = molecular_weight(dna, double_stranded=True)
print(f'single-stranded (DEFAULT): {mw_ss:.2f} Da')
print(f'double_stranded=True: {mw_ds:.2f} Da')
print(f'ratio ds/ss (not exactly 2 for non-self-complementary): {mw_ds / mw_ss:.4f}')

# Melting temperature: accurate nearest-neighbor for a primer
print('\n=== Melting Temperature (primer) ===')
primer = Seq('ACGGTCAGGTCAGGTACGGT')  # 20-mer
print(f'Tm_Wallace (2+4, wrong for primers): {MeltingTemp.Tm_Wallace(primer):.1f} C')
print(f'Tm_NN default salt: {MeltingTemp.Tm_NN(primer, strict=True):.1f} C')
# Mg/dNTPs only honored at saltcorr 6 or 7
tm_no_mg = MeltingTemp.Tm_NN(primer, Mg=1.5, dNTPs=0.2, saltcorr=5)
tm_with_mg = MeltingTemp.Tm_NN(primer, Mg=1.5, dNTPs=0.2, saltcorr=7)
print(f'Mg=1.5 at saltcorr=5 (Mg silently ignored): {tm_no_mg:.1f} C')
print(f'Mg=1.5 at saltcorr=7 (Owczarzy 2008, Mg honored): {tm_with_mg:.1f} C')

# ProtParam summary after sanitizing non-standard residues
print('\n=== Protein Report ===')
raw = 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQ*X'
clean = raw.upper().replace('*', '').replace('X', '')
protein = ProteinAnalysis(clean)
helix, turn, sheet = protein.secondary_structure_fraction()
print(f'length: {len(clean)} aa')
print(f'molecular_weight: {protein.molecular_weight():.2f} Da')
print(f'isoelectric_point: {protein.isoelectric_point():.2f}')
print(f'charge_at_pH 7.0: {protein.charge_at_pH(7.0):+.2f}')
ii = protein.instability_index()
print(f'instability_index: {ii:.2f} ({"unstable" if ii > 40 else "stable"})')
print(f'gravy: {protein.gravy():.3f}')
print(f'aromaticity: {protein.aromaticity():.4f}')
print(f'helix/turn/sheet: {helix * 100:.1f}% / {turn * 100:.1f}% / {sheet * 100:.1f}%')

# Confirm the misspelled default scale literal by introspection
default_scale = inspect.signature(ProteinAnalysis.gravy).parameters['scale'].default
print(f"\ngravy() default scale literal (note one 't'): {default_scale!r}")
