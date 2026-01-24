import requests
from pathlib import Path

def predict_esmfold(sequence, output_file=None):
    '''Predict protein structure using ESMFold API'''
    url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    # timeout=300: ESMFold can take several minutes for long sequences (>500 residues)
    response = requests.post(url, data=sequence, timeout=300)
    if response.status_code != 200:
        raise Exception(f'ESMFold failed with status {response.status_code}')

    pdb_text = response.text
    if output_file:
        Path(output_file).write_text(pdb_text)
        print(f'Structure saved to {output_file}')
    return pdb_text

def extract_plddt(pdb_text):
    '''Extract per-residue pLDDT scores from ESMFold PDB output'''
    plddt = {}
    for line in pdb_text.split('\n'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            resnum = int(line[22:26])
            bfactor = float(line[60:66])
            plddt[resnum] = bfactor
    return plddt

def analyze_confidence(plddt):
    '''Summarize pLDDT confidence regions'''
    # pLDDT thresholds: >90 very high, 70-90 confident, 50-70 low, <50 very low
    very_high = [r for r, s in plddt.items() if s > 90]
    confident = [r for r, s in plddt.items() if 70 <= s <= 90]
    low = [r for r, s in plddt.items() if 50 <= s < 70]
    very_low = [r for r, s in plddt.items() if s < 50]

    avg = sum(plddt.values()) / len(plddt)
    print(f'Average pLDDT: {avg:.1f}')
    print(f'Very high confidence (>90): {len(very_high)} residues')
    print(f'Confident (70-90): {len(confident)} residues')
    print(f'Low confidence (50-70): {len(low)} residues')
    print(f'Very low (<50, likely disordered): {len(very_low)} residues')

    return {'avg': avg, 'very_high': very_high, 'confident': confident, 'low': low, 'very_low': very_low}

if __name__ == '__main__':
    # Example: Human hemoglobin alpha chain (first 50 residues)
    # For full proteins, use sequences from UniProt
    sequence = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH'
    pdb_text = predict_esmfold(sequence, 'esmfold_prediction.pdb')
    plddt = extract_plddt(pdb_text)
    analyze_confidence(plddt)
