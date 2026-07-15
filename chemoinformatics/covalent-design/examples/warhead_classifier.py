# Reference: RDKit 2024.09+ | Verify API if version differs
# Classify and score covalent warheads in a compound library

from rdkit import Chem


WARHEAD_SMARTS = {
    'acrylamide': '[CX3](=[OX1])([NX3])[CX3]=[CX3]',
    'alpha_substituted_acrylamide': '[CX3](=[OX1])([NX3])[CX3]([!#1])=[CX3]',
    'methacrylamide': '[CX3](=[OX1])([NX3])[CX3]([#6])=[CX3]',
    'chloroacetamide': '[CX3](=[OX1])([NX3])[CH2][Cl]',
    'bromoacetamide': '[CX3](=[OX1])([NX3])[CH2][Br]',
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

# Nominal warhead classes are qualitative annotations, not compound-level rate predictions.
REACTIVITY_TIER = {
    'chloroacetamide': 'high',
    'bromoacetamide': 'high',
    'maleimide': 'very_high',
    'isothiocyanate': 'high',
    'isocyanate': 'high',
    'epoxide': 'high',
    'aziridine': 'high',
    'acrylamide': 'moderate',
    'alpha_substituted_acrylamide': 'context_dependent',
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


def has_recognized_warhead(smi):
    '''Return whether any catalogued substructure is present; not a suitability filter.'''
    matches = classify_warheads(smi)
    if matches is None:
        return False
    return bool(matches)


if __name__ == '__main__':
    # Valid generic N-substituted acrylamide example
    acrylamide_example = 'C=CC(=O)N1CCCCC1'
    print(classify_warheads(acrylamide_example))

    # Chloroacetamide ABPP probe
    abpp = 'O=C(CCl)NCc1ccccc1'
    print(classify_warheads(abpp))

    print(has_recognized_warhead(acrylamide_example))
    print(has_recognized_warhead(abpp))
