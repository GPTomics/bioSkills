'''Transform structure coordinates with the correct rotation convention'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.vectors import rotaxis, Vector
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

center = np.array([a.coord for a in structure.get_atoms()]).mean(axis=0)
print(f'Original center: {center}')

# rotaxis returns a row-convention matrix for Entity.transform (coords @ rot + tran).
rot = rotaxis(np.radians(90), Vector(0, 0, 1))

# Rotate about the center of mass: pick tran so the center is the fixed point.
tran = center - center @ rot
structure.transform(rot, tran)

new_center = np.array([a.coord for a in structure.get_atoms()]).mean(axis=0)
print(f'New center after rotation about COM: {new_center}')

io = PDBIO()
io.set_structure(structure)
io.save('rotated.pdb')
print('Saved rotated structure to rotated.pdb')
