'''Find restriction sites in a circular plasmid'''
# Reference: biopython 1.83+ (API verified on 1.86) | Verify API if version differs
# linear=False matters: a circular search finds sites spanning the origin and
# changes fragment counts (a circle with n cuts gives n fragments, not n+1).

from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from Bio.Restriction import EcoRI, BamHI, HindIII, XhoI, NotI, XbaI, NcoI, NdeI

record = SeqIO.read('plasmid.gb', 'genbank')
seq = record.seq
print(f'Plasmid: {record.id}, Length: {len(seq)} bp')

cloning_enzymes = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI, NotI, XbaI, NcoI, NdeI])
analysis = Analysis(cloning_enzymes, seq, linear=False)

print('\nEnzymes that cut once (good for linearization):')
for enzyme, sites in analysis.with_N_sites(1).items():
    print(f'  {enzyme}: position {sites[0]}')

print('\nEnzymes that cut twice (good for excision):')
for enzyme, sites in analysis.with_N_sites(2).items():
    print(f'  {enzyme}: positions {sites}')

print('\nEnzymes that do not cut (good for cloning):')
for enzyme in analysis.without_site():
    print(f'  {enzyme}')

print('\nCommercially available enzymes that cut once:')
comm_analysis = Analysis(CommOnly, seq, linear=False)
comm_once = comm_analysis.with_N_sites(1)
print(f'  Found {len(comm_once)} single-cutters')
