'''Select enzymes based on criteria'''

from Bio import SeqIO
from Bio.Restriction import Analysis, CommOnly

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq

print(f'Analyzing: {record.id} ({len(seq)} bp)')
print('=' * 50)

analysis = Analysis(CommOnly, seq)

once = analysis.once_cutters()
twice = analysis.twice_cutters()
non = analysis.only_dont_cut()

print(f'\nSingle-cutters ({len(once)}):')
for enzyme, sites in list(once.items())[:10]:
    print(f'  {enzyme}: position {sites[0]}, site {enzyme.site}')
if len(once) > 10:
    print(f'  ... and {len(once) - 10} more')

print(f'\nDouble-cutters ({len(twice)}):')
for enzyme, sites in list(twice.items())[:10]:
    print(f'  {enzyme}: positions {sites}')

print(f'\nNon-cutters: {len(non)} enzymes')

print('\n\nGrouped by overhang type (single-cutters only):')
blunt = [e for e in once.keys() if e.is_blunt()]
five_prime = [e for e in once.keys() if e.is_5overhang()]
three_prime = [e for e in once.keys() if e.is_3overhang()]

print(f"  Blunt ({len(blunt)}): {[str(e) for e in blunt[:5]]}")
print(f"  5' overhang ({len(five_prime)}): {[str(e) for e in five_prime[:5]]}")
print(f"  3' overhang ({len(three_prime)}): {[str(e) for e in three_prime[:5]]}")

print('\n\nGrouped by recognition site length:')
for length in [4, 6, 8]:
    enzymes = [e for e in once.keys() if len(e.site) == length]
    print(f'  {length}-cutters ({len(enzymes)}): {[str(e) for e in enzymes[:5]]}')
