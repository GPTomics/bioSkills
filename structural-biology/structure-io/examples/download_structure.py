'''Download a deposited structure and its biological assembly from RCSB'''
# Reference: biopython 1.85+ | Verify API if version differs

import gzip
import shutil
import urllib.request
from Bio.PDB import PDBList

pdbl = PDBList()

pdb_id = '4HHB'
# pdir=None writes a two-char divided subdir tree; pass pdir to control location.
deposited = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='mmCif')
print(f'Downloaded deposited (asymmetric unit) coords: {deposited}')

# retrieve_pdb_file has no assembly_num; fetch the pre-built assembly directly.
# The deposited coords are the ASU, often not the functional oligomer.
assembly_gz = f'{pdb_id}-assembly1.cif.gz'
url = f'https://files.rcsb.org/download/{assembly_gz}'
urllib.request.urlretrieve(url, assembly_gz)
with gzip.open(assembly_gz, 'rb') as fin, open(f'{pdb_id}-assembly1.cif', 'wb') as fout:
    shutil.copyfileobj(fin, fout)
print(f'Downloaded biological assembly: {pdb_id}-assembly1.cif')
