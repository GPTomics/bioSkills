# Reference: RDKit 2024.09+ | Verify API if version differs
# Coarse feature-family prefilter for pharmacophore workflow

import os
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures
from rdkit.RDPaths import RDDataDir


def get_feature_factory(fdef_path=None):
    '''Build RDKit feature factory; default uses BaseFeatures.fdef from RDKit.'''
    if fdef_path is None:
        fdef_path = os.path.join(RDDataDir, 'BaseFeatures.fdef')
    return ChemicalFeatures.BuildFeatureFactory(fdef_path)


def molecule_features(mol, factory):
    '''Extract pharmacophore features from a molecule.'''
    feats = factory.GetFeaturesForMol(mol)
    return [(f.GetFamily(), f.GetType(), list(f.GetAtomIds()), f.GetPos())
            for f in feats]


def shared_feature_types_prefilter(active_mols, factory, n_conf=20):
    '''Find shared feature labels; this does not derive a 3D pharmacophore.'''
    # Generate 3D conformers
    embedded = []
    for mol in active_mols:
        if mol is None:
            continue
        m = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(m, numConfs=n_conf, params=AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMoleculeConfs(m)
        embedded.append(m)

    # Extract features per molecule
    all_features = [molecule_features(m, factory) for m in embedded]
    if not all_features:
        return set()

    # Cross-reference: find features appearing in all active compounds
    common_types = set([(f[0], f[1]) for f in all_features[0]])
    for feats in all_features[1:]:
        types = set([(f[0], f[1]) for f in feats])
        common_types &= types
    return common_types


def has_feature_types(query_features, target_mol, factory):
    '''Check feature-family presence only; this is not a 3D distance match.'''
    target = Chem.AddHs(target_mol)
    AllChem.EmbedMolecule(target, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(target)

    target_features = molecule_features(target, factory)

    # Simple feature-type match (loose); production: distance-constrained
    target_types = set([(f[0], f[1]) for f in target_features])
    query_types = set([(qf[0], qf[1]) for qf in query_features])
    return query_types.issubset(target_types)


def feature_family_prefilter(query_mol_list, library_smiles):
    '''Prefilter by feature families before a distance-constrained 3D search.'''
    factory = get_feature_factory()
    active_mols = [Chem.MolFromSmiles(s) if isinstance(s, str) else s
                   for s in query_mol_list]

    common = shared_feature_types_prefilter(active_mols, factory)
    if not common:
        return []

    hits = []
    for smi in library_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        if has_feature_types(list(common), mol, factory):
            hits.append(smi)
    return hits


if __name__ == '__main__':
    queries = ['CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1',
               'CC(=O)Nc1ccc(C(=O)c2ccc(F)cc2)cc1']
    library = ['CCC(=O)Nc1ccc(C(=O)c2ccc(Cl)cc2)cc1',
               'CCCCCC',
               'c1ccncc1']
    hits = feature_family_prefilter(queries, library)
    print(f'Feature-family prefilter matches: {hits}')
