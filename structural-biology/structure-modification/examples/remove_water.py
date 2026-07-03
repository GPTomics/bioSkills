'''Strip crystallographic water by HETFLAG while keeping metals and cofactors'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# 'W' HETFLAG isolates water; r.id[0] != ' ' would also delete Zn/Mg/heme/FAD and MSE.
water_count = 0
for chain in structure[0]:
    water_ids = [r.id for r in chain if r.id[0] == 'W']
    water_count += len(water_ids)
    for res_id in water_ids:
        chain.detach_child(res_id)

print(f'Removed {water_count} water molecules; hetero ligands and metals preserved')

io = PDBIO()
io.set_structure(structure)
io.save('no_water.pdb')
print('Saved to no_water.pdb')
