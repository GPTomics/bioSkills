'''Pure-Python pairwise structure superposition with Bio.PDB.Superimposer.

Use when residue correspondence is known (e.g. apo vs holo of the same protein,
mutant vs wild-type). For unknown correspondence, prefer TMalign / USalign.
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio.PDB import PDBParser, Superimposer, PDBIO

if __name__ == '__main__':
    parser = PDBParser(QUIET=True)
    mobile_structure = parser.get_structure('mobile', 'mobile.pdb')
    reference_structure = parser.get_structure('reference', 'reference.pdb')

    mobile_atoms = list(mobile_structure.get_atoms())
    reference_atoms = list(reference_structure.get_atoms())

    ca_mobile = [a for a in mobile_atoms if a.get_id() == 'CA']
    ca_reference = [a for a in reference_atoms if a.get_id() == 'CA']

    if len(ca_mobile) != len(ca_reference):
        print(f'Warning: CA-atom counts differ (mobile={len(ca_mobile)}, reference={len(ca_reference)}).')
        print('Truncating to the shorter list assumes residue order matches; if structures are not')
        print('co-numbered (different inserts/deletions), use TMalign / USalign which handle')
        print('unknown correspondence via dynamic programming.')

    n_pairs = min(len(ca_mobile), len(ca_reference))
    if n_pairs == 0:
        raise RuntimeError('No CA atoms found in one or both structures; check input PDB files')

    sup = Superimposer()
    sup.set_atoms(ca_reference[:n_pairs], ca_mobile[:n_pairs])
    sup.apply(mobile_atoms)

    print(f'Superposed {n_pairs} CA atoms')
    print(f'RMSD: {sup.rms:.3f} A')
    print(f'Rotation matrix:\n{sup.rotran[0]}')
    print(f'Translation vector: {sup.rotran[1]}')

    io = PDBIO()
    io.set_structure(mobile_structure)
    io.save('mobile_superposed.pdb')
    print('\nSuperposed structure written to mobile_superposed.pdb')
