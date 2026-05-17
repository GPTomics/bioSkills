# Reference: RDKit 2024.09+, chembl_structure_pipeline 1.2+, pandas 2.2+ | Verify API if version differs

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_structure_pipeline import standardize_mol, get_parent_mol
import pandas as pd


def chembl_standardize(smi):
    '''Apply ChEMBL structure pipeline: standardize + get parent.'''
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, 'parse_failure'
    try:
        std_mol, _ = standardize_mol(mol)
        parent_mol, _ = get_parent_mol(std_mol)
        return Chem.MolToSmiles(parent_mol), 'ok'
    except Exception as e:
        return None, f'standardize_error: {e}'


def rdkit_standardize(smi, keep_isotopes=False):
    '''Full rdMolStandardize pipeline: sanitize -> largest fragment -> normalize -> uncharge -> canon tautomer.'''
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    Chem.SanitizeMol(mol)

    largest = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    mol = largest.choose(mol)

    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)

    # canonicalOrdering=True ensures stable output regardless of input atom ordering
    uncharger = rdMolStandardize.Uncharger(canonicalOrdering=True)
    mol = uncharger.uncharge(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    if not keep_isotopes:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(mol)


def prepare_qsar_data(df, smiles_col='smiles', activity_col='pIC50', pipeline='chembl'):
    '''Standardize, generate InChIKey, deduplicate, mean-aggregate activity for ML training data.'''
    rows = []
    for _, row in df.iterrows():
        smi = row[smiles_col]
        if pipeline == 'chembl':
            std_smi, status = chembl_standardize(smi)
        else:
            std_smi = rdkit_standardize(smi)
            status = 'ok' if std_smi else 'failed'
        if std_smi is None or status != 'ok':
            continue
        mol = Chem.MolFromSmiles(std_smi)
        # Skip inorganic / no-carbon compounds
        if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
            continue
        rows.append({
            'smiles': std_smi,
            'inchikey': Chem.MolToInchiKey(mol),
            'activity': row[activity_col],
            'original_smiles': smi,
        })
    df_std = pd.DataFrame(rows)

    # Deduplicate by InChIKey; mean-aggregate activity
    df_dedup = df_std.groupby('inchikey').agg(
        smiles=('smiles', 'first'),
        activity=('activity', 'mean'),
        activity_std=('activity', 'std'),
        n_replicates=('activity', 'count'),
        original_smiles=('original_smiles', list),
    ).reset_index()
    return df_dedup


if __name__ == '__main__':
    test_smiles = [
        'CC(=O)Oc1ccccc1C(=O)O',           # aspirin
        '[Na+].CC(=O)Oc1ccccc1C(=O)[O-]',  # aspirin sodium salt
        'CCC(=O)Nc1ccc(C(=O)C)cc1',         # different tautomer-ambiguous
        'invalid_smiles',                    # parse failure
    ]
    df = pd.DataFrame({'smiles': test_smiles, 'pIC50': [6.0, 6.2, 5.5, None]})
    result = prepare_qsar_data(df)
    print(result[['inchikey', 'smiles', 'activity', 'n_replicates']])
