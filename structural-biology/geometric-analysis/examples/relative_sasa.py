'''Compute per-residue relative SASA and classify burial'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# Tien et al 2013 theoretical Gly-X-Gly max-ASA (A^2), the current best normalization scale
MAX_ASA = {'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0}

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# probe_radius 1.4 A = water; the absolute SASA is meaningless without recording it
ShrakeRupley(probe_radius=1.40, n_points=100).compute(structure, level='R')

buried, exposed = 0, 0
for res in structure.get_residues():
    if res.resname in MAX_ASA and hasattr(res, 'sasa'):
        rsa = res.sasa / MAX_ASA[res.resname]
        if rsa < 0.20:                               # RSA < 0.20 is the common buried heuristic (rule of thumb, not a law)
            buried += 1
        else:
            exposed += 1

print(f'Buried (RSA < 0.20): {buried}')
print(f'Exposed: {exposed}')
