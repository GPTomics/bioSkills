'''Find restriction sites in a sequence'''
# Reference: biopython 1.83+ (API verified on 1.86) | Verify API if version differs

from Bio import SeqIO
from Bio.Restriction import EcoRI, BamHI, HindIII, XhoI, RestrictionBatch, Analysis

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq
print(f'Sequence: {record.id}, Length: {len(seq)} bp')

ecori_sites = EcoRI.search(seq)
print(f'\nEcoRI sites (1-based cut positions): {ecori_sites}')

common_enzymes = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI])
analysis = Analysis(common_enzymes, seq)

print('\nAll cut sites:')
for enzyme, sites in analysis.with_sites().items():
    print(f'  {enzyme}: {sites}')

print('\nEnzymes that cut exactly once:')
for enzyme, sites in analysis.with_N_sites(1).items():
    print(f'  {enzyme}: position {sites[0]}')

print('\nEnzymes that do not cut:')
for enzyme in analysis.without_site():
    print(f'  {enzyme}')
