'''DNA sequence property calculations'''
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight, MeltingTemp

seq = Seq('ATGCGATCGATCGATCGATCG')
print(f'Sequence: {seq}')
print(f'Length: {len(seq)} bp')

# GC Content
print('\n=== GC Content ===')
gc = gc_fraction(seq)
print(f'GC fraction: {gc:.4f}')
print(f'GC percent: {gc * 100:.2f}%')

# Base composition
print('\n=== Base Composition ===')
seq_str = str(seq)
total = len(seq_str)
for base in 'ATGC':
    count = seq_str.count(base)
    pct = count / total * 100
    print(f'{base}: {count} ({pct:.1f}%)')

# Molecular weight
print('\n=== Molecular Weight ===')
mw_ss = molecular_weight(seq)
mw_ds = molecular_weight(seq, double_stranded=True)
print(f'Single-stranded: {mw_ss:.2f} Da')
print(f'Double-stranded: {mw_ds:.2f} Da')

# Melting temperature
print('\n=== Melting Temperature ===')
primer = Seq('ATGCGATCGATCGATCGATC')  # 20-mer
print(f'Primer: {primer} ({len(primer)} bp)')
tm_wallace = MeltingTemp.Tm_Wallace(primer)
tm_gc = MeltingTemp.Tm_GC(primer)
tm_nn = MeltingTemp.Tm_NN(primer)
print(f'Tm (Wallace rule): {tm_wallace:.1f} C')
print(f'Tm (GC method): {tm_gc:.1f} C')
print(f'Tm (Nearest neighbor): {tm_nn:.1f} C')
