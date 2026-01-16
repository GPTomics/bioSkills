import requests
import os
from Bio.PDB import PDBParser

def download_alphafold(uniprot_id, output_dir='.'):
    os.makedirs(output_dir, exist_ok=True)
    base_url = 'https://alphafold.ebi.ac.uk/files'

    pdb_url = f'{base_url}/AF-{uniprot_id}-F1-model_v4.pdb'
    response = requests.get(pdb_url)
    if response.status_code != 200:
        print(f'Not found: {uniprot_id}')
        return None

    pdb_path = f'{output_dir}/AF-{uniprot_id}-F1-model_v4.pdb'
    with open(pdb_path, 'w') as f:
        f.write(response.text)
    print(f'Downloaded: {pdb_path}')
    return pdb_path

def analyze_plddt(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    plddt_scores = []
    for residue in structure[0].get_residues():
        if residue.id[0] == ' ' and 'CA' in residue:
            plddt_scores.append(residue['CA'].get_bfactor())

    return {
        'mean': sum(plddt_scores) / len(plddt_scores),
        'min': min(plddt_scores),
        'max': max(plddt_scores),
        'n_residues': len(plddt_scores),
        'high_confidence': sum(1 for s in plddt_scores if s > 70) / len(plddt_scores) * 100
    }

if __name__ == '__main__':
    pdb_file = download_alphafold('P04637', 'alphafold_structures')
    if pdb_file:
        stats = analyze_plddt(pdb_file)
        print(f"\npLDDT Statistics:")
        print(f"  Mean: {stats['mean']:.1f}")
        print(f"  Range: {stats['min']:.1f} - {stats['max']:.1f}")
        print(f"  Residues: {stats['n_residues']}")
        print(f"  High confidence (>70): {stats['high_confidence']:.1f}%")
