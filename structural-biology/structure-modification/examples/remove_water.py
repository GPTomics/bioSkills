'''Remove water molecules from a structure'''

from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

water_count = 0
for model in structure:
    for chain in model:
        residues_to_remove = [r.id for r in chain if r.id[0] == 'W']
        water_count += len(residues_to_remove)
        for res_id in residues_to_remove:
            chain.detach_child(res_id)

print(f'Removed {water_count} water molecules')

io = PDBIO()
io.set_structure(structure)
io.save('no_water.pdb')
print('Saved to no_water.pdb')
