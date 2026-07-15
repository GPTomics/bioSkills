# Reference: RDKit 2024.09+ | Verify API if version differs

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdFingerprintGenerator
from rdkit import DataStructs


def prepare_mol_3d(smiles, n_conf=20, seed=42):
    '''Generate 3D conformers; 20 is a starting budget to convergence-check.'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not ids:
        raise ValueError(f'Embedding failed for {smiles!r}')
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError(f'MMFF parameters unavailable for {smiles!r}')
    optimization = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=1000)
    failed = [ids[i] for i, (status, _) in enumerate(optimization) if status != 0]
    if failed:
        raise RuntimeError(f'MMFF optimization did not converge for conformers {failed}')
    return mol


def open3dalign_score(query_mol, target_mol):
    '''Align conformers with O3A and return the best normalized shape score.'''
    if query_mol is None or target_mol is None:
        return None
    best_shape_score = -1.0
    best_o3a_score = None
    best_qcid = -1
    best_tcid = -1
    best_coords = None
    for q_cid in range(query_mol.GetNumConformers()):
        for t_cid in range(target_mol.GetNumConformers()):
            probe = Chem.Mol(target_mol)
            O3A = rdMolAlign.GetO3A(probe, query_mol,
                                     prbCid=t_cid, refCid=q_cid)
            O3A.Align()
            o3a_score = O3A.Score()
            shape_score = 1.0 - rdShapeHelpers.ShapeTanimotoDist(
                probe, query_mol, confId1=t_cid, confId2=q_cid
            )
            if shape_score > best_shape_score:
                best_shape_score = shape_score
                best_o3a_score = o3a_score
                best_qcid = q_cid
                best_tcid = t_cid
                best_coords = list(probe.GetConformer(t_cid).GetPositions())
    if best_coords is not None:
        conf = target_mol.GetConformer(best_tcid)
        for atom_idx, point in enumerate(best_coords):
            conf.SetAtomPosition(atom_idx, point)
    return best_shape_score, best_o3a_score, best_qcid, best_tcid


def ecfp4_tanimoto(smi1, smi2):
    '''Compute ECFP4 Tanimoto for diversity check.'''
    m1 = Chem.MolFromSmiles(smi1)
    m2 = Chem.MolFromSmiles(smi2)
    if m1 is None or m2 is None:
        return None
    # Radius 2 / 2048 bits is a repository ECFP4 starting representation;
    # compare only fingerprints built with identical, project-selected settings.
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp1 = generator.GetFingerprint(m1)
    fp2 = generator.GetFingerprint(m2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def scaffold_hop_candidates(query_smi, library_smiles, shape_threshold=0.5,
                             ecfp_threshold=0.5, n_conf=20):
    '''Find scaffold hops using repository defaults that require calibration.'''
    query_mol = prepare_mol_3d(query_smi, n_conf=n_conf)
    candidates = []
    for smi in library_smiles:
        target_mol = prepare_mol_3d(smi, n_conf=n_conf)
        if target_mol is None:
            continue
        shape_result = open3dalign_score(query_mol, target_mol)
        if shape_result is None:
            continue
        shape_score, o3a_score, _, _ = shape_result
        ecfp = ecfp4_tanimoto(query_smi, smi)
        if shape_score >= shape_threshold and ecfp is not None and ecfp < ecfp_threshold:
            candidates.append({
                'smiles': smi,
                'shape_tanimoto': shape_score,
                'o3a_score': o3a_score,
                'ecfp_tanimoto': ecfp,
            })
    return sorted(candidates, key=lambda x: x['shape_tanimoto'], reverse=True)


if __name__ == '__main__':
    query = 'CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1'
    library = [
        'CC(=O)Nc1ccc(C(=O)c2ccc(F)cc2)cc1',  # similar 2D + 3D
        'O=S(=O)(c1ccccc1)Nc2ccc(C(=O)c3ccccc3)cc2',  # scaffold hop candidate
        'CCC',
    ]
    # Demo-only cutoffs are chosen to emit an illustrative result for this fixed
    # toy set; calibrate both thresholds on a task-relevant benchmark.
    hops = scaffold_hop_candidates(query, library, shape_threshold=0.3, ecfp_threshold=0.55)
    if not hops:
        raise RuntimeError('The fixed demonstration set produced no scaffold-hop result')
    for h in hops:
        print(h)
