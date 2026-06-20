'''Predict restriction digest fragments'''
# Reference: biopython 1.83+ (API verified on 1.86) | Verify API if version differs
# catalyze() returns the tuple of fragments DIRECTLY -- do NOT index [0]
# (that returns only the first fragment). A linear molecule with n cuts -> n+1 fragments.

from Bio import SeqIO
from Bio.Restriction import EcoRI, BamHI, HindIII

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq

print(f'Sequence: {record.id} ({len(seq)} bp)')
print('=' * 50)

for enzyme in [EcoRI, BamHI, HindIII]:
    sites = enzyme.search(seq)
    if not sites:
        print(f'\n{enzyme}: No cut sites')
        continue

    fragments = enzyme.catalyze(seq, linear=True)
    sizes = sorted([len(f) for f in fragments], reverse=True)

    print(f'\n{enzyme} ({enzyme.site}):')
    print(f'  Cut sites: {sites}')
    print(f'  Fragments: {len(fragments)}')
    print(f'  Sizes: {sizes}')
    print(f'  Total: {sum(sizes)} bp')          # must equal len(seq)

ecori_sites = EcoRI.search(seq)
bamhi_sites = BamHI.search(seq)
all_sites = sorted(set(ecori_sites + bamhi_sites))

if all_sites:
    print('\n\nDouble digest (EcoRI + BamHI):')
    print(f'  Combined sites: {all_sites}')

    sizes = [all_sites[0]]
    sizes += [all_sites[i + 1] - all_sites[i] for i in range(len(all_sites) - 1)]
    sizes.append(len(seq) - all_sites[-1])

    sizes = sorted(sizes, reverse=True)
    print(f'  Fragments: {len(sizes)}')
    print(f'  Sizes: {sizes}')
    print(f'  Total: {sum(sizes)} bp')
