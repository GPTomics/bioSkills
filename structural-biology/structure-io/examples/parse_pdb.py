'''Parse a PDB file and inspect structure contents'''

from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('1abc', '1abc.pdb')

print(f'Structure: {structure.id}')
print(f'Models: {len(list(structure.get_models()))}')
print(f'Chains: {len(list(structure.get_chains()))}')
print(f'Residues: {len(list(structure.get_residues()))}')
print(f'Atoms: {len(list(structure.get_atoms()))}')

for model in structure:
    for chain in model:
        residues = list(chain.get_residues())
        print(f'  Chain {chain.id}: {len(residues)} residues')
