# Reference: fair-esm 2.0+, biopython 1.83+, requests 2.31+ | Verify API if version differs
import requests
from pathlib import Path

def predict_esmfold_local(sequence, output_file=None, device='cuda'):
    '''Predict a single chain locally with ESMFold (no MSA, needs ~16 GB GPU).'''
    import torch, esm
    model = esm.pretrained.esmfold_v1().eval().to(device)
    with torch.no_grad():
        pdb_text = model.infer_pdb(sequence)  # pLDDT rides in the B-factor column
    if output_file:
        Path(output_file).write_text(pdb_text)
    return pdb_text

def predict_esmfold_api(sequence, output_file=None):
    '''Fallback to the hosted endpoint; it is intermittently down (SSL/500). Prefer local.'''
    url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    resp = requests.post(url, data=sequence, timeout=300)  # 300 s: long sequences take minutes
    resp.raise_for_status()
    pdb_text = resp.text
    if output_file:
        Path(output_file).write_text(pdb_text)
    return pdb_text

def extract_plddt(pdb_text):
    '''Extract per-residue pLDDT (B-factor column) from an ESMFold PDB string.'''
    plddt = {}
    for line in pdb_text.split('\n'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            plddt[int(line[22:26])] = float(line[60:66])
    return plddt

def analyze_confidence(plddt):
    '''Band pLDDT by the AlphaFold/EBI convention: >90 very high, 70-90 confident, 50-70 low, <50 very low.'''
    very_high = [r for r, s in plddt.items() if s > 90]
    confident = [r for r, s in plddt.items() if 70 <= s <= 90]
    low = [r for r, s in plddt.items() if 50 <= s < 70]
    very_low = [r for r, s in plddt.items() if s < 50]  # a contiguous band usually flags an IDR, not an error
    avg = sum(plddt.values()) / len(plddt)
    print(f'mean pLDDT {avg:.1f}: very-high {len(very_high)}, confident {len(confident)}, low {len(low)}, very-low {len(very_low)}')
    return {'avg': avg, 'very_high': very_high, 'confident': confident, 'low': low, 'very_low': very_low}

def interface_is_reliable(iptm, ptm=None, floor=0.5):
    '''Gate a complex interface on ipTM (interface pTM), not per-chain pLDDT.

    floor=0.5: below this the interface is unreliable or the chains likely do not interact;
    a confident interface wants ipTM > ~0.8. AF-Multimer ranks by 0.8*ipTM + 0.2*pTM.
    '''
    if iptm is None or iptm < floor:
        print(f'interface NOT reliable (ipTM={iptm}); inspect the inter-chain PAE block before any claim')
        return False
    print(f'interface plausible (ipTM={iptm:.2f}' + (f', pTM={ptm:.2f}' if ptm is not None else '') + ')')
    return True

if __name__ == '__main__':
    # Human hemoglobin alpha chain, first 50 residues. Use full sequences from UniProt in practice.
    sequence = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH'
    pdb_text = predict_esmfold_local(sequence, 'esmfold_prediction.pdb')
    analyze_confidence(extract_plddt(pdb_text))
