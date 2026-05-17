# Reference: RDKit 2024.09+ | Verify API if version differs
# Classify and score covalent warheads in a compound library

from rdkit import Chem


WARHEAD_SMARTS = {
    'acrylamide': '[CX3](=O)[NH][CX3]=[CX3]',
    'alpha_substituted_acrylamide': '[CX3](=O)[NH][CX3]([!H])=[CX3]',  # GSH-stable variant
    'methacrylamide': '[CX3](=O)[NH][CX3](C)=[CX3]',
    'chloroacetamide': '[CX3](=O)[CH2][Cl]',
    'bromoacetamide': '[CX3](=O)[CH2][Br]',
    'vinyl_sulfone': '[SX4](=O)(=O)[CX3]=[CX3]',
    'sulfonyl_fluoride': '[SX4](=O)(=O)[F]',
    'fluorosulfate_sufex': '[O][SX4](=O)(=O)[F]',
    'aldehyde': '[CX3H1](=O)',
    'boronate': '[B]([O])[O]',
    'nitrile': '[C]#[N]',
    'epoxide': '[C]1[O][C]1',
    'aziridine': '[C]1[N][C]1',
    'maleimide': 'O=C1N(C(=O)/C=C/1)',
    'isothiocyanate': '[NX2]=C=[SX1]',
    'isocyanate': '[NX2]=C=[OX1]',
}

# Reactivity tiers (descending; clinical drugs use 'moderate' or 'low'; ABPP uses 'high')
REACTIVITY_TIER = {
    'chloroacetamide': 'high',
    'bromoacetamide': 'high',
    'maleimide': 'very_high',
    'isothiocyanate': 'high',
    'isocyanate': 'high',
    'epoxide': 'high',
    'aziridine': 'high',
    'acrylamide': 'moderate',
    'alpha_substituted_acrylamide': 'low_moderate',  # GSH-stable; drug-like
    'methacrylamide': 'low',
    'vinyl_sulfone': 'moderate',
    'sulfonyl_fluoride': 'moderate',
    'fluorosulfate_sufex': 'moderate',
    'aldehyde': 'reversible',
    'boronate': 'reversible',
    'nitrile': 'low',
}

# Target residue selectivity
RESIDUE_SELECTIVITY = {
    'chloroacetamide': ['Cys'],
    'bromoacetamide': ['Cys'],
    'acrylamide': ['Cys'],
    'alpha_substituted_acrylamide': ['Cys'],
    'methacrylamide': ['Cys'],
    'vinyl_sulfone': ['Cys'],
    'sulfonyl_fluoride': ['Lys', 'Tyr', 'Ser'],
    'fluorosulfate_sufex': ['Tyr', 'Lys'],
    'aldehyde': ['Cys', 'Lys', 'Ser'],
    'boronate': ['Ser', 'Thr'],
    'nitrile': ['Cys'],
    'epoxide': ['Cys', 'Lys', 'Asp'],
    'aziridine': ['Cys', 'Lys'],
    'maleimide': ['Cys'],
    'isothiocyanate': ['Cys', 'Lys'],
    'isocyanate': ['Lys', 'Ser'],
}


def classify_warheads(smi):
    '''Identify warheads in a SMILES; return per-warhead matches + reactivity tier.'''
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    matches = {}
    for name, smarts in WARHEAD_SMARTS.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            continue
        m = mol.GetSubstructMatches(pat)
        if m:
            matches[name] = {
                'count': len(m),
                'reactivity_tier': REACTIVITY_TIER.get(name, 'unknown'),
                'targets': RESIDUE_SELECTIVITY.get(name, ['unknown']),
            }
    return matches


def filter_for_drug_like_covalent(smi):
    '''Filter to GSH-stable warheads suitable for drug development.'''
    matches = classify_warheads(smi)
    if matches is None:
        return False
    drug_like_warheads = {'alpha_substituted_acrylamide', 'methacrylamide',
                          'acrylamide', 'vinyl_sulfone', 'fluorosulfate_sufex',
                          'nitrile', 'sulfonyl_fluoride'}
    return any(w in matches for w in drug_like_warheads)


if __name__ == '__main__':
    # Sotorasib / KRAS G12C-like acrylamide
    sotorasib_warhead = 'O=C(C=C)N1CCC(c2c...)cc1'
    print(classify_warheads(sotorasib_warhead))

    # Chloroacetamide ABPP probe
    abpp = 'O=C(CCl)NCc1ccccc1'
    print(classify_warheads(abpp))

    print(filter_for_drug_like_covalent(sotorasib_warhead))
    print(filter_for_drug_like_covalent(abpp))
