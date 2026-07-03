'''Fetch an AFDB model via REST metadata and read pLDDT + PAE confidence correctly.'''
# Reference: biopython 1.83+, numpy 1.26+, requests 2.31+ | Verify API if version differs
import os
import requests
import numpy as np
from Bio.PDB import MMCIFParser

def afdb_metadata(accession):
    r = requests.get(f'https://alphafold.ebi.ac.uk/api/prediction/{accession}')
    r.raise_for_status()
    return r.json()  # list; one object per fragment/isoform, empty if no model exists

def fetch_afdb(accession, out_dir='.'):
    entries = afdb_metadata(accession)
    if not entries:
        return None  # >2700-aa non-human proteins and non-UniProt sequences are often absent
    entry = entries[0]  # discover URLs from metadata - version suffix (v6 as of 2025) drifts
    cif_text = requests.get(entry['cifUrl']).text
    pae = requests.get(entry['paeDocUrl']).json()
    os.makedirs(out_dir, exist_ok=True)  # write target may be a nested dir that does not exist yet
    cif_path = f'{out_dir}/AF-{accession}.cif'
    with open(cif_path, 'w') as f:
        f.write(cif_text)
    return cif_path, pae

def extract_plddt(cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('af', cif_file)
    plddt = {}
    for residue in structure[0].get_residues():
        if residue.id[0] == ' ' and 'CA' in residue:
            # pLDDT rides in the B-factor column but is CONFIDENCE (0-100, higher=better),
            # opposite polarity to a thermal B-factor - never read it as motion
            plddt[residue.id[1]] = residue['CA'].get_bfactor()
    return plddt

def plddt_band(score):
    if score > 90:                       # 90: side-chain-trustworthy core (AFDB very-high cut)
        return 'very_high'
    if score >= 70:                      # 70: backbone-reliable fold (AFDB confident cut)
        return 'confident'
    if score >= 50:                      # 50: coarse topology only below this
        return 'low'
    return 'very_low'                    # <50: usually an intrinsically disordered region

def load_pae(pae_json):
    entry = pae_json[0] if isinstance(pae_json, list) else pae_json
    # Compact format (2023+): 2D num_res x num_res array; legacy 1D fields were removed
    return np.array(entry['predicted_aligned_error'])

if __name__ == '__main__':
    result = fetch_afdb('P04637', 'alphafold_structures')  # human p53
    if result:
        cif_path, pae_json = result
        plddt = extract_plddt(cif_path)
        pae = load_pae(pae_json)
        mean_plddt = sum(plddt.values()) / len(plddt)
        disordered = [res for res, s in plddt.items() if s < 50]
        print(f'Residues: {len(plddt)}')
        print(f'Mean pLDDT: {mean_plddt:.1f}')
        print(f'Very-low (candidate IDR) residues: {len(disordered)}')
        print(f'PAE matrix shape: {pae.shape}')
