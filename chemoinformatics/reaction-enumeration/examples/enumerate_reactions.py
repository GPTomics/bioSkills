#!/usr/bin/env python3
'''
Reaction enumeration for virtual library generation.
'''

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from itertools import product


REACTION_SMARTS = {
    'amide_coupling': '[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]',
    'reductive_amination': '[C:1]=O.[N:2]>>[C:1][N:2]',
    'suzuki': '[c:1][Br].[c:2][B](O)O>>[c:1][c:2]',
    'buchwald': '[c:1][Br].[N:2]>>[c:1][N:2]',
    'ester_formation': '[C:1](=[O:2])O.[O:3]>>[C:1](=[O:2])[O:3]',
}


def enumerate_reaction(rxn_smarts, reactant_lists, deduplicate=True):
    '''
    Enumerate products from combinatorial reaction.

    Args:
        rxn_smarts: Reaction SMARTS string
        reactant_lists: List of lists of SMILES
        deduplicate: Remove duplicate products
    '''
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)

    if rxn.Validate()[0] != 0:
        raise ValueError('Invalid reaction SMARTS')

    products = []
    seen = set()

    for reactants in product(*reactant_lists):
        mols = [Chem.MolFromSmiles(s) for s in reactants]
        if None in mols:
            continue

        try:
            prods = rxn.RunReactants(tuple(mols))
            for prod_set in prods:
                for prod in prod_set:
                    try:
                        Chem.SanitizeMol(prod)
                        smiles = Chem.MolToSmiles(prod)

                        if deduplicate:
                            if smiles not in seen:
                                seen.add(smiles)
                                products.append(smiles)
                        else:
                            products.append(smiles)
                    except Exception:
                        continue
        except Exception:
            continue

    return products


def validate_products(smiles_list, mw_max=500, logp_max=5):
    '''Validate and filter enumerated products.'''
    valid = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        mw = Descriptors.MolWt(mol)
        if mw > mw_max:
            continue

        logp = Descriptors.MolLogP(mol)
        if logp > logp_max:
            continue

        try:
            Chem.SanitizeMol(mol)
            valid.append(smiles)
        except Exception:
            continue

    return valid


def multi_step_synthesis(building_blocks, reaction_sequence):
    '''
    Enumerate products from multi-step synthesis.

    Args:
        building_blocks: Dict of {step_index: [smiles_list]}
        reaction_sequence: List of reaction SMARTS
    '''
    current = building_blocks[0]

    for step, rxn_smarts in enumerate(reaction_sequence):
        next_bbs = building_blocks.get(step + 1, [])
        if not next_bbs:
            break

        current = enumerate_reaction(rxn_smarts, [current, next_bbs])
        print(f'Step {step + 1}: {len(current)} intermediates')

    return current


def apply_rgroup_decoration(core_smiles, r_groups):
    '''
    Apply R-group decoration to a core scaffold.

    Args:
        core_smiles: Core with * attachment point
        r_groups: List of R-group SMILES
    '''
    products = []

    for rg in r_groups:
        product_smiles = core_smiles.replace('*', f'({rg})', 1)
        mol = Chem.MolFromSmiles(product_smiles)
        if mol:
            try:
                Chem.SanitizeMol(mol)
                products.append(Chem.MolToSmiles(mol))
            except Exception:
                continue

    return products


if __name__ == '__main__':
    acids = ['CC(=O)O', 'c1ccccc1C(=O)O', 'OC(=O)CCc1ccccc1']
    amines = ['CCN', 'c1ccc(N)cc1', 'NCC(C)C']

    print('Amide library enumeration:')
    products = enumerate_reaction(
        REACTION_SMARTS['amide_coupling'],
        [acids, amines]
    )
    print(f'Generated {len(products)} products')

    print('\nSample products:')
    for p in products[:5]:
        print(f'  {p}')

    valid = validate_products(products, mw_max=400)
    print(f'\n{len(valid)} products pass MW < 400 filter')

    print('\nR-group decoration:')
    core = '*c1ccccc1'
    r_groups = ['C', 'CC', 'C(=O)O', 'N']
    decorated = apply_rgroup_decoration(core, r_groups)
    for d in decorated:
        print(f'  {d}')
