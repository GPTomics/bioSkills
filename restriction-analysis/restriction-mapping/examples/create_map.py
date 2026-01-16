'''Create a restriction map'''

from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Restriction import EcoRI, BamHI, HindIII, XhoI, NotI, XbaI

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq
seq_len = len(seq)

print(f'Restriction Map: {record.id}')
print(f'Sequence length: {seq_len} bp')
print('=' * 60)

enzymes = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI, NotI, XbaI])
analysis = Analysis(enzymes, seq)

print('\nMap format:')
analysis.print_as('map')

print('\n\nLinear format:')
analysis.print_as('linear')

print('\n\nDetailed site list:')
results = analysis.full()

all_cuts = []
for enzyme, sites in results.items():
    for site in sites:
        all_cuts.append((site, str(enzyme)))

all_cuts.sort()

for pos, enz in all_cuts:
    pct = (pos / seq_len) * 100
    print(f'  {pos:6d} bp ({pct:5.1f}%) - {enz}')

print('\n\nFragment sizes between consecutive sites:')
if len(all_cuts) > 1:
    positions = [c[0] for c in all_cuts]
    for i in range(len(positions) - 1):
        size = positions[i + 1] - positions[i]
        print(f'  {positions[i]} -> {positions[i + 1]}: {size} bp')
