'''Find restriction sites in a sequence'''

from Bio import SeqIO
from Bio.Restriction import EcoRI, BamHI, HindIII, XhoI, RestrictionBatch, Analysis

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq
print(f'Sequence: {record.id}, Length: {len(seq)} bp')

ecori_sites = EcoRI.search(seq)
print(f'\nEcoRI sites: {ecori_sites}')

common_enzymes = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI])
analysis = Analysis(common_enzymes, seq)

print('\nAll cut sites:')
results = analysis.full()
for enzyme, sites in results.items():
    if sites:
        print(f'  {enzyme}: {sites}')

print('\nEnzymes that cut once:')
once_cutters = analysis.once_cutters()
for enzyme, sites in once_cutters.items():
    print(f'  {enzyme}: position {sites[0]}')

print('\nEnzymes that do not cut:')
non_cutters = analysis.only_dont_cut()
for enzyme in non_cutters:
    print(f'  {enzyme}')
