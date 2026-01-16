'''Iterate over structure hierarchy and count components'''

from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

print(f'Structure: {structure.id}')
print(f'Models: {len(list(structure.get_models()))}')
print(f'Chains: {len(list(structure.get_chains()))}')
print(f'Residues: {len(list(structure.get_residues()))}')
print(f'Atoms: {len(list(structure.get_atoms()))}')

print('\nPer-chain breakdown:')
for model in structure:
    for chain in model:
        amino_acids = [r for r in chain if r.id[0] == ' ']
        hetero = [r for r in chain if r.id[0].startswith('H_')]
        waters = [r for r in chain if r.id[0] == 'W']
        print(f'  Chain {chain.id}: {len(amino_acids)} AA, {len(hetero)} ligands, {len(waters)} waters')
