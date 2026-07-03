'''Overload the B-factor column for coloring after snapshotting the originals'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# This column is a real temperature factor (or AlphaFold pLDDT) until overwritten - snapshot it.
original_bfactors = {atom.get_full_id(): atom.bfactor for atom in structure.get_atoms()}

conservation_scores = {i: (i % 10) for i in range(1, 200)}

for residue in structure.get_residues():
    if residue.id[0] != ' ':
        continue
    score = conservation_scores.get(residue.id[1])
    if score is None:
        continue
    for atom in residue:
        atom.bfactor = score  # set on every atom so per-atom coloring is not patchy

io = PDBIO()
io.set_structure(structure)
io.save('colored_by_conservation.pdb')
print(f'Snapshotted {len(original_bfactors)} original B-factors before overwriting')
print('B-factors set to conservation scores; do not send this file to refinement')
print('Open in PyMOL and color by B-factor: spectrum b')
