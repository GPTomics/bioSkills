'''Extract protein sequences from structure'''

from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

ppb = PPBuilder()
for i, pp in enumerate(ppb.build_peptides(structure)):
    seq = pp.get_sequence()
    print(f'Polypeptide {i + 1}: length={len(seq)}')
    print(f'  {seq[:60]}...' if len(seq) > 60 else f'  {seq}')
