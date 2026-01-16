'''Measure distances between atoms in a structure'''

from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

chain = structure[0]['A']

atom1 = chain[50]['CA']
atom2 = chain[100]['CA']

distance = atom1 - atom2
print(f'Distance between CA 50 and CA 100: {distance:.2f} Angstroms')

ca_atoms = [r['CA'] for r in chain if r.has_id('CA') and r.id[0] == ' ']
n = len(ca_atoms)
print(f'\nDistance matrix for {n} CA atoms:')

dist_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        dist = ca_atoms[i] - ca_atoms[j]
        dist_matrix[i, j] = dist
        dist_matrix[j, i] = dist

print(f'Min distance: {dist_matrix[dist_matrix > 0].min():.2f} A')
print(f'Max distance: {dist_matrix.max():.2f} A')
