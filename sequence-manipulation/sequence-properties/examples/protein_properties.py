'''Protein sequence property calculations'''
from Bio.SeqUtils.ProtParam import ProteinAnalysis

protein_seq = 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGMKNGKKIELKDLVSGFLAEKQGSPTFFGYMKFLSNSEIVVLPNNVAPNVRYIIQQYGFYHHVGTWNNNSHAKIGLIILYLNKEKTLFNNNVQNKRTSHLLSQMYDPKK'
protein = ProteinAnalysis(protein_seq)

print(f'Protein sequence: {protein_seq[:50]}...')
print(f'Length: {len(protein_seq)} aa')

# Basic properties
print('\n=== Basic Properties ===')
print(f'Molecular weight: {protein.molecular_weight():.2f} Da')
print(f'Isoelectric point (pI): {protein.isoelectric_point():.2f}')

# Stability
print('\n=== Stability ===')
ii = protein.instability_index()
stable = 'stable' if ii < 40 else 'unstable'
print(f'Instability index: {ii:.2f} ({stable})')

# Hydropathicity
print('\n=== Hydropathicity ===')
gravy = protein.gravy()
char = 'hydrophilic' if gravy < 0 else 'hydrophobic'
print(f'GRAVY: {gravy:.3f} ({char})')

# Aromaticity
print('\n=== Aromaticity ===')
print(f'Aromaticity: {protein.aromaticity():.4f}')

# Secondary structure
print('\n=== Secondary Structure Prediction ===')
helix, turn, sheet = protein.secondary_structure_fraction()
print(f'Helix: {helix * 100:.1f}%')
print(f'Turn: {turn * 100:.1f}%')
print(f'Sheet: {sheet * 100:.1f}%')

# Extinction coefficient
print('\n=== Extinction Coefficient (280 nm) ===')
eps = protein.molar_extinction_coefficient()
print(f'Reduced cysteines: {eps[0]} M^-1 cm^-1')
print(f'Cystine bridges: {eps[1]} M^-1 cm^-1')

# Amino acid composition (top 5)
print('\n=== Amino Acid Composition (Top 5) ===')
aa_pct = protein.amino_acids_percent
sorted_aa = sorted(aa_pct.items(), key=lambda x: x[1], reverse=True)[:5]
for aa, pct in sorted_aa:
    print(f'{aa}: {pct * 100:.1f}%')
