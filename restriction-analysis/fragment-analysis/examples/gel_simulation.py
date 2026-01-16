'''Simulate gel electrophoresis pattern'''

from Bio import SeqIO
from Bio.Restriction import EcoRI, BamHI, HindIII, RestrictionBatch, Analysis

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq

ladder = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250]

def get_fragment_sizes(seq, enzyme, linear=True):
    fragments = enzyme.catalyze(seq, linear=linear)[0]
    return sorted([len(f) for f in fragments], reverse=True)

digests = {
    'EcoRI': get_fragment_sizes(seq, EcoRI),
    'BamHI': get_fragment_sizes(seq, BamHI),
    'HindIII': get_fragment_sizes(seq, HindIII),
}

all_bands = set(ladder)
for sizes in digests.values():
    all_bands.update(sizes)
all_bands = sorted(all_bands, reverse=True)

header = f'{"Size":>6} | {"Ladder":^6}'
for name in digests.keys():
    header += f' | {name:^8}'
print(header)
print('-' * len(header))

for band in all_bands:
    row = f'{band:>6} |'
    row += f' {"---" if band in ladder else "   ":^6} |'
    for name, sizes in digests.items():
        count = sizes.count(band)
        mark = '====' * count if count else ''
        row += f' {mark:^8} |'
    print(row)

print('\n\nFragment summary:')
for name, sizes in digests.items():
    if sizes:
        print(f'  {name}: {sizes} (n={len(sizes)})')
    else:
        print(f'  {name}: No cuts')
