'''Parse a structure and inspect contents, with the auth-vs-label numbering trap'''
# Reference: biopython 1.85+ | Verify API if version differs

from Bio.PDB import MMCIFParser

# auth numbering (default) matches the paper; label numbering is gap-free internal.
auth_parser = MMCIFParser(QUIET=True)
structure = auth_parser.get_structure('4hhb', '4hhb.cif')

print(f'Structure: {structure.id}')
print(f'Models: {len(list(structure.get_models()))}')
print(f'Chains: {len(list(structure.get_chains()))}')
print(f'Residues: {len(list(structure.get_residues()))}')
print(f'Atoms: {len(list(structure.get_atoms()))}')

for model in structure:
    for chain in model:
        residues = list(chain.get_residues())
        first = residues[0].id[1] if residues else None
        print(f'  Chain {chain.id}: {len(residues)} residues, first resseq {first}')

label_parser = MMCIFParser(QUIET=True, auth_residues=False, auth_chains=False)
label_structure = label_parser.get_structure('4hhb', '4hhb.cif')
label_first = next(label_structure.get_residues()).id[1]
print(f'Label-scheme first resseq (differs from auth): {label_first}')
