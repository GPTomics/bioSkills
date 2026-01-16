'''Transform structure coordinates'''

from Bio.PDB import PDBParser, PDBIO
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

coords = np.array([a.coord for a in structure.get_atoms()])
center = coords.mean(axis=0)
print(f'Original center: {center}')

for atom in structure.get_atoms():
    atom.coord = atom.coord - center

coords = np.array([a.coord for a in structure.get_atoms()])
new_center = coords.mean(axis=0)
print(f'New center: {new_center}')

io = PDBIO()
io.set_structure(structure)
io.save('centered.pdb')
print('Saved centered structure to centered.pdb')
