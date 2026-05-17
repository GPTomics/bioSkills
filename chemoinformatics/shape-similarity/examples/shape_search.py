# Reference: RDKit 2024.09+, usrcat 1.2+ (optional) | Verify API if version differs

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit import DataStructs


def prepare_mol_3d(smiles, n_conf=20, seed=42):
    '''Generate 3D conformer ensemble for shape comparison.'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    AllChem.MMFFOptimizeMoleculeConfs(mol)
    return mol


def open3dalign_score(query_mol, target_mol):
    '''Compute best-conformer Open3DAlign score between query and target.'''
    if query_mol is None or target_mol is None:
        return None
    best_score = -1.0
    best_qcid = -1
    best_tcid = -1
    for q_cid in range(query_mol.GetNumConformers()):
        for t_cid in range(target_mol.GetNumConformers()):
            O3A = rdMolAlign.GetO3A(target_mol, query_mol,
                                     prbCid=t_cid, refCid=q_cid)
            score = O3A.Score()
            if score > best_score:
                best_score = score
                best_qcid = q_cid
                best_tcid = t_cid
    return best_score, best_qcid, best_tcid


def ecfp4_tanimoto(smi1, smi2):
    '''Compute ECFP4 Tanimoto for diversity check.'''
    m1 = Chem.MolFromSmiles(smi1)
    m2 = Chem.MolFromSmiles(smi2)
    if m1 is None or m2 is None:
        return None
    fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, 2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def scaffold_hop_candidates(query_smi, library_smiles, shape_threshold=0.5,
                             ecfp_threshold=0.5, n_conf=20):
    '''Find scaffold hops: high shape similarity, low ECFP4 similarity.'''
    query_mol = prepare_mol_3d(query_smi, n_conf=n_conf)
    candidates = []
    for smi in library_smiles:
        target_mol = prepare_mol_3d(smi, n_conf=n_conf)
        if target_mol is None:
            continue
        shape_result = open3dalign_score(query_mol, target_mol)
        if shape_result is None:
            continue
        shape_score, _, _ = shape_result
        ecfp = ecfp4_tanimoto(query_smi, smi)
        if shape_score >= shape_threshold and ecfp is not None and ecfp < ecfp_threshold:
            candidates.append({'smiles': smi, 'shape_score': shape_score, 'ecfp_tanimoto': ecfp})
    return sorted(candidates, key=lambda x: x['shape_score'], reverse=True)


if __name__ == '__main__':
    query = 'CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1'
    library = [
        'CC(=O)Nc1ccc(C(=O)c2ccc(F)cc2)cc1',  # similar 2D + 3D
        'O=S(=O)(c1ccccc1)Nc2ccc(C(=O)c3ccccc3)cc2',  # scaffold hop candidate
        'CCC',
    ]
    hops = scaffold_hop_candidates(query, library, shape_threshold=0.3, ecfp_threshold=0.5)
    for h in hops:
        print(h)
