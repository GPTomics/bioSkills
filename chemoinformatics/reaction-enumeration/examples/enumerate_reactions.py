#!/usr/bin/env python3
'''
Reaction enumeration for virtual library generation.
'''
# Reference: rdkit 2024.03+ | Verify API if version differs

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

    num_warnings, num_errors = rxn.Validate()
    if num_errors:
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
    Connect a core and R group through one explicit [*] on each input.

    Args:
        core_smiles: Core with exactly one [*] attachment point
        r_groups: List of R-group SMILES, each with exactly one [*]
    '''
    core = Chem.MolFromSmiles(core_smiles)
    if core is None:
        raise ValueError('Core is not valid SMILES')
    core_dummies = [a.GetIdx() for a in core.GetAtoms() if a.GetAtomicNum() == 0]
    if len(core_dummies) != 1:
        raise ValueError('Core must contain exactly one [*] attachment point')

    products = []

    for rg in r_groups:
        rgroup = Chem.MolFromSmiles(rg)
        if rgroup is None:
            raise ValueError(f'R group is not valid SMILES: {rg}')
        rg_dummies = [
            a.GetIdx() for a in rgroup.GetAtoms() if a.GetAtomicNum() == 0
        ]
        if len(rg_dummies) != 1:
            raise ValueError(f'R group must contain exactly one [*]: {rg}')

        core_dummy = core.GetAtomWithIdx(core_dummies[0])
        rg_dummy = rgroup.GetAtomWithIdx(rg_dummies[0])
        if core_dummy.GetDegree() != 1 or rg_dummy.GetDegree() != 1:
            raise ValueError('Each attachment dummy must have exactly one neighbor')

        core_neighbor = core_dummy.GetNeighbors()[0].GetIdx()
        rg_neighbor = rg_dummy.GetNeighbors()[0].GetIdx()
        core_bond = core.GetBondBetweenAtoms(core_dummies[0], core_neighbor)
        rg_bond = rgroup.GetBondBetweenAtoms(rg_dummies[0], rg_neighbor)
        if (
            core_bond.GetBondType() != Chem.BondType.SINGLE
            or rg_bond.GetBondType() != Chem.BondType.SINGLE
        ):
            raise ValueError('Attachment dummy bonds must be single bonds')
        offset = core.GetNumAtoms()
        combined = Chem.RWMol(Chem.CombineMols(core, rgroup))
        combined.AddBond(core_neighbor, offset + rg_neighbor, Chem.BondType.SINGLE)
        for idx in sorted((core_dummies[0], offset + rg_dummies[0]), reverse=True):
            combined.RemoveAtom(idx)
        product_mol = combined.GetMol()
        Chem.SanitizeMol(product_mol)
        products.append(Chem.MolToSmiles(product_mol, canonical=True))

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
    r_groups = ['[*]C', '[*]CC', '[*]C(=O)O', '[*]N']
    decorated = apply_rgroup_decoration(core, r_groups)
    for d in decorated:
        print(f'  {d}')
