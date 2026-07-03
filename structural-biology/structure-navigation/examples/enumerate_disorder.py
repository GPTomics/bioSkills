'''Enumerate disordered atoms instead of silently using one conformer'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# A DisorderedAtom forwards uncaught calls to its highest-occupancy child, so plain
# iteration measures ONE altloc. Enumerate every conformer before any geometry.
n_disordered = 0
for residue in structure.get_residues():
    for atom in residue:
        if atom.is_disordered():
            n_disordered += 1
            hetflag, resseq, icode = residue.id
            for alt in atom.disordered_get_list():
                print(f'{residue.resname}{resseq}{icode} {atom.name} altloc {alt.altloc}: occ={alt.occupancy}')

print(f'\ntotal disordered atoms: {n_disordered}')
