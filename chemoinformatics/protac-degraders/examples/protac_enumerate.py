# Reference: RDKit 2024.09+ | Verify API if version differs
# Enumerate PROTAC linkers between target-ligand and E3-ligand fragments

from rdkit import Chem
from rdkit.Chem import AllChem


# Common PROTAC linker building blocks
LINKERS = {
    # Length 2-4 atoms
    'short_alkyl': '[*]CC[*]',
    'short_pegylated': '[*]COC[*]',
    'piperazine_short': '[*]N1CCN([*])CC1',
    # Length 5-9 atoms
    'medium_alkyl': '[*]CCCCC[*]',
    'medium_pegylated': '[*]CCOCC[*]',
    'piperazine_acid': '[*]N1CCN(CC1)C(=O)[*]',
    'triazole': '[*]Cn1cc(C[*])nn1',
    # Length 10-14 atoms
    'long_alkyl': '[*]CCCCCCCCCC[*]',
    'long_pegylated': '[*]COCCOCCOC[*]',
    'rigid_piperidine': '[*]C(=O)CN1CCC(CCC[*])CC1',
    'triazole_extended': '[*]CCCCn1cc(CCCC[*])nn1',
}


def enumerate_linkers(target_fragment_smi, e3_fragment_smi, linker_smarts_list=None):
    '''
    Combine target-ligand fragment + linker + E3-ligand fragment.
    Fragments must have [*] attachment points.
    '''
    if linker_smarts_list is None:
        linker_smarts_list = list(LINKERS.values())

    protac_smiles = []
    target_frag = Chem.MolFromSmiles(target_fragment_smi)
    e3_frag = Chem.MolFromSmiles(e3_fragment_smi)
    if target_frag is None or e3_frag is None:
        return []

    for linker_smarts in linker_smarts_list:
        # Use BRICSBuild-style combination
        # Simplified: replace [*] in target with linker[*][*]; then replace second [*] with e3 fragment
        # Real implementation: use RDKit's BRICS or custom fragment-combination
        try:
            t_smi = Chem.MolToSmiles(target_frag).replace('[*]', '*[N+H3]', 1)
            l_mol = Chem.MolFromSmiles(linker_smarts)
            if l_mol is None:
                continue
            combined = t_smi + '.' + linker_smarts + '.' + Chem.MolToSmiles(e3_frag)
            protac_smiles.append(combined)
        except Exception:
            continue
    return protac_smiles


def estimate_linker_atoms(distance_A, rigidity='flexible'):
    '''Estimate atoms needed to span distance.'''
    # Flexible: ~1.5 A per C-C bond + 20% slack
    if rigidity == 'flexible':
        return int(distance_A / 1.4 * 1.2)
    # Rigid: 1.4 A per sp2 C-C + curvature factor
    return int(distance_A / 1.4 * 1.5)


def compute_protac_size(protac_smi):
    '''Compute MW and TPSA for permeability assessment.'''
    from rdkit.Chem import Descriptors
    mol = Chem.MolFromSmiles(protac_smi)
    if mol is None:
        return None
    return {
        'MolWt': Descriptors.MolWt(mol),
        'TPSA': Descriptors.TPSA(mol),
        'LogP': Descriptors.MolLogP(mol),
        'RotBonds': Descriptors.NumRotatableBonds(mol),
        'permeability_concern': Descriptors.MolWt(mol) > 1200 or Descriptors.TPSA(mol) > 150,
    }


if __name__ == '__main__':
    # Example: BTK ibrutinib fragment + pomalidomide (CRBN) fragment
    target_frag = '[*]N1CCC(c2cnc3nc(...))CC1'  # placeholder ibrutinib-like
    e3_frag = '[*]CC(=O)NC1CCC(=O)NC1=O'  # placeholder pomalidomide-like

    protacs = enumerate_linkers(target_frag, e3_frag)
    for smi in protacs[:5]:
        size = compute_protac_size(smi)
        print(f'PROTAC: {smi[:80]}...; size: {size}')

    # Estimate linker length for 12 A target-E3 distance
    print(f'Atoms needed for 12 A: {estimate_linker_atoms(12, "flexible")}')
