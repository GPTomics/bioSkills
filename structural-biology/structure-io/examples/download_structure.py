'''Download structures from RCSB PDB'''
# Reference: biopython 1.83+, scanpy 1.10+ | Verify API if version differs

from Bio.PDB import PDBList

pdbl = PDBList()

pdb_id = '4HHB'
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='mmCif')
print(f'Downloaded {pdb_id}: {file_path}')

# Download biological assembly
assembly_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='mmCif', assembly_num=1)
print(f'Downloaded assembly: {assembly_path}')
