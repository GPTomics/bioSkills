'''Extract a single protein chain, excluding water and hetero residues'''
# Reference: biopython 1.85+ | Verify API if version differs

from Bio.PDB import MMCIFParser, PDBIO, Select

class ProteinChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # id[0] is the hetflag: ' ' standard, 'W' water, 'H_XXX' hetero.
        return residue.id[0] == ' '

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('4hhb', '4hhb.cif')

io = PDBIO()
io.set_structure(structure)
io.save('chain_A_protein.pdb', ProteinChainSelect('A'))
print('Extracted protein chain A to chain_A_protein.pdb')
