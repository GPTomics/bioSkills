'''Find residue contacts using NeighborSearch'''

from Bio.PDB import PDBParser, NeighborSearch

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

all_atoms = list(structure.get_atoms())
ns = NeighborSearch(all_atoms)

contact_distance = 4.0
contacts = ns.search_all(contact_distance, level='R')

print(f'Residue contacts within {contact_distance} A:')
print(f'Total contacts: {len(contacts)}')
# Expected output: ~500-2000 contacts for typical small protein (100-300 residues)
# Rule of thumb: ~3-5 contacts per residue at 4A cutoff, fewer at stricter thresholds

print('\nFirst 10 contacts:')
for res1, res2 in contacts[:10]:
    print(f'  {res1.get_parent().id}:{res1.resname}{res1.id[1]} - {res2.get_parent().id}:{res2.resname}{res2.id[1]}')
