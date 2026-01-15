'''Using reverse complement for primer design'''
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

def design_primers(template, target_start, target_end, primer_length=20):
    '''Design forward and reverse primers for a target region'''
    forward = template[target_start:target_start + primer_length]
    reverse_region = template[target_end - primer_length:target_end]
    reverse = reverse_region.reverse_complement()
    return forward, reverse

def primer_info(name, primer):
    '''Calculate basic primer properties'''
    gc = gc_fraction(primer) * 100
    return f'{name}: 5\'-{primer}-3\' (Tm estimate: {2 * (str(primer).count("A") + str(primer).count("T")) + 4 * (str(primer).count("G") + str(primer).count("C"))}C, GC: {gc:.1f}%)'

# Example template sequence
template = Seq('ATGCGATCGATCGATCGATCGAAAAAAAAAAAGATCGATCGATCGATCGATCG')
print('=== Template ===')
print(f'5\'-{template}-3\'')
print(f'Length: {len(template)} bp')

# Design primers for region 10-40
print('\n=== Primer Design (region 10-40) ===')
fwd, rev = design_primers(template, 10, 40, primer_length=18)
print(primer_info('Forward', fwd))
print(primer_info('Reverse', rev))

# Show amplicon
print('\n=== Amplicon ===')
amplicon = template[10:40]
print(f'Amplicon: {amplicon}')
print(f'Length: {len(amplicon)} bp')

# Verify primers bind correctly
print('\n=== Verification ===')
print(f'Forward binds at 5\' end: {template[10:28]}')
print(f'Reverse binds at 3\' end (shown as revcomp): {rev}')
print(f'Match: {str(fwd) == str(template[10:28])}')
