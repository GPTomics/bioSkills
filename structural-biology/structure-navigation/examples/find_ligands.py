'''Find and list all ligands in a structure'''

from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

print('Ligands found:')
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0].startswith('H_'):
                hetfield, resseq, icode = residue.id
                n_atoms = len(list(residue.get_atoms()))
                print(f'  Chain {chain.id}: {residue.resname} at position {resseq} ({n_atoms} atoms)')
