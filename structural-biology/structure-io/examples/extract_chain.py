'''Extract a single chain from a structure'''

from Bio.PDB import PDBParser, PDBIO, Select

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', '1abc.pdb')

io = PDBIO()
io.set_structure(structure)
io.save('chain_A.pdb', ChainSelect('A'))
print('Extracted chain A to chain_A.pdb')
