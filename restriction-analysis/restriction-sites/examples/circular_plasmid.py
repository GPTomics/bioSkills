'''Find restriction sites in a circular plasmid'''

from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from Bio.Restriction import EcoRI, BamHI, HindIII, XhoI, NotI, XbaI, NcoI, NdeI

record = SeqIO.read('plasmid.gb', 'genbank')
seq = record.seq
print(f'Plasmid: {record.id}, Length: {len(seq)} bp')

cloning_enzymes = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI, NotI, XbaI, NcoI, NdeI])
analysis = Analysis(cloning_enzymes, seq, linear=False)

print('\nEnzymes that cut once (good for linearization):')
once = analysis.once_cutters()
for enzyme, sites in once.items():
    print(f'  {enzyme}: position {sites[0]}')

print('\nEnzymes that cut twice (good for excision):')
twice = analysis.twice_cutters()
for enzyme, sites in twice.items():
    print(f'  {enzyme}: positions {sites}')

print('\nEnzymes that do not cut (good for cloning):')
non_cutters = analysis.only_dont_cut()
for enzyme in non_cutters:
    print(f'  {enzyme}')

print('\nCommercially available enzymes that cut once:')
comm_analysis = Analysis(CommOnly, seq, linear=False)
comm_once = comm_analysis.once_cutters()
print(f'  Found {len(comm_once)} single-cutters')
