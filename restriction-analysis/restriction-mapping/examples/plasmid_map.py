'''Create restriction map for a circular plasmid with feature overlap analysis'''

from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Restriction import EcoRI, BamHI, HindIII, XhoI, NotI, NcoI, NdeI

record = SeqIO.read('plasmid.gb', 'genbank')
seq = record.seq
seq_len = len(seq)

print(f'Plasmid: {record.id}')
print(f'Description: {record.description}')
print(f'Length: {seq_len} bp (circular)')
print('=' * 70)

enzymes = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI, NotI, NcoI, NdeI])
analysis = Analysis(enzymes, seq, linear=False)

results = analysis.full()

print('\nRestriction sites:')
all_cuts = []
for enzyme, sites in results.items():
    if sites:
        print(f'\n{enzyme} ({len(sites)} site{"s" if len(sites) > 1 else ""}):')
        for site in sites:
            all_cuts.append((site, str(enzyme)))
            overlapping = []
            for feature in record.features:
                start = int(feature.location.start)
                end = int(feature.location.end)
                if start <= site <= end:
                    label = feature.qualifiers.get('label', feature.qualifiers.get('gene', [feature.type]))[0]
                    overlapping.append(f'{feature.type}:{label}')
            overlap_str = ', '.join(overlapping) if overlapping else 'intergenic'
            print(f'  {site:6d} bp - {overlap_str}')

all_cuts.sort()

print('\n\nFragment sizes (circular):')
if len(all_cuts) >= 2:
    positions = [c[0] for c in all_cuts]
    enzymes_at_pos = {c[0]: c[1] for c in all_cuts}

    for i in range(len(positions)):
        pos1 = positions[i]
        pos2 = positions[(i + 1) % len(positions)]
        if pos2 > pos1:
            size = pos2 - pos1
        else:
            size = (seq_len - pos1) + pos2
        print(f'  {enzymes_at_pos[pos1]}({pos1}) -> {enzymes_at_pos[pos2]}({pos2}): {size} bp')
