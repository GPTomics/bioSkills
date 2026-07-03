'''Extract the OBSERVED sequence and compare it to declared SEQRES length'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser(QUIET=True, get_header=True)
structure = parser.get_structure('protein', 'protein.pdb')

# PPBuilder returns the observed (ATOM) sequence; missing-density loops are concatenated away.
ppb = PPBuilder()
observed = ''.join(str(pp.get_sequence()) for pp in ppb.build_peptides(structure))
print(f'observed (ATOM) length: {len(observed)}')

# SEQRES is the full declared sequence, including residues that were never resolved.
for record in SeqIO.parse('protein.pdb', 'pdb-seqres'):
    print(f'{record.id} declared (SEQRES) length: {len(record.seq)}')

missing = structure.header.get('missing_residues', [])
print(f'unmodeled residues (the difference): {len(missing)}')
