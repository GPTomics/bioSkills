'''Superimpose two structures and calculate RMSD'''

from Bio.PDB import PDBParser, Superimposer, PDBIO

parser = PDBParser(QUIET=True)
ref_structure = parser.get_structure('ref', 'reference.pdb')
mobile_structure = parser.get_structure('mobile', 'mobile.pdb')

ref_atoms = [r['CA'] for r in ref_structure.get_residues() if r.has_id('CA') and r.id[0] == ' ']
mobile_atoms = [r['CA'] for r in mobile_structure.get_residues() if r.has_id('CA') and r.id[0] == ' ']

n = min(len(ref_atoms), len(mobile_atoms))
ref_atoms = ref_atoms[:n]
mobile_atoms = mobile_atoms[:n]

sup = Superimposer()
sup.set_atoms(ref_atoms, mobile_atoms)
print(f'RMSD: {sup.rms:.2f} Angstroms')

sup.apply(mobile_structure.get_atoms())

io = PDBIO()
io.set_structure(mobile_structure)
io.save('aligned.pdb')
print('Aligned structure saved to aligned.pdb')
