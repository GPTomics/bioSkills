# Reference: RDKit 2024.09+, numpy 1.26+ | Verify API if version differs

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def gen_conformer_ensemble(smiles, n_conf=20, seed=42, optimize=True):
    '''Generate 3D conformer ensemble via ETKDGv3 + MMFF94 optimization.'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, []
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.maxIterations = 1000
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)

    if not optimize:
        return mol, [(int(cid), None) for cid in ids]

    conf_data = []
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    use_uff = mmff_props is None
    if use_uff and not AllChem.UFFHasAllMoleculeParams(mol):
        raise ValueError('Neither MMFF94s nor UFF covers this molecule')

    for cid in ids:
        ff = (AllChem.UFFGetMoleculeForceField(mol, confId=cid) if use_uff else
              AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid))
        if ff is None:
            raise ValueError(f'Could not construct force field for conformer {cid}')
        status = ff.Minimize(maxIts=1000)
        if status != 0:
            raise RuntimeError(f'Optimization did not converge for conformer {cid}')
        conf_data.append((int(cid), float(ff.CalcEnergy())))
    return mol, conf_data


def macrocycle_conformers(smiles, n_conf=200, seed=42):
    '''Generate conformers with macrocycle torsion preferences for >=12 atom rings.'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smiles!r}')
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.useMacrocycleTorsions = True
    params.useSmallRingTorsions = True
    params.maxIterations = 5000
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not ids:
        raise RuntimeError(f'No macrocycle conformers embedded for {smiles!r}')
    return mol


def prune_rmsd(mol, conf_data, rmsd_cutoff=0.5):
    '''Remove conformers within RMSD cutoff of an already-kept one.'''
    keep = []
    for cid, e in conf_data:
        is_unique = True
        for kept_cid, _ in keep:
            rmsd = AllChem.GetBestRMS(mol, mol, cid, kept_cid)
            if rmsd < rmsd_cutoff:
                is_unique = False
                break
        if is_unique:
            keep.append((cid, e))
    return keep


def filter_energy_window(conf_data, window_kcal=10.0):
    '''Keep only conformers within energy window of the minimum.'''
    if not conf_data:
        return []
    if any(e is None for _, e in conf_data):
        raise ValueError('Energy filtering requires optimized conformers with energies')
    min_e = min(e for _, e in conf_data)
    return [(cid, e) for cid, e in conf_data if (e - min_e) <= window_kcal]


def boltzmann_average(values, energies, T=300.0):
    '''Boltzmann-weighted average of per-conformer properties.'''
    energies = np.array(energies)
    kt = 0.001987 * T  # kcal/mol per K (Boltzmann's constant in chemistry units)
    rel = energies - energies.min()
    w = np.exp(-rel / kt)
    w = w / w.sum()
    return float(np.sum(np.array(values) * w))


if __name__ == '__main__':
    # Example: caffeine, n_conf=20 with MMFF94 optimization
    mol, conf_data = gen_conformer_ensemble('Cn1cnc2c1c(=O)n(C)c(=O)n2C', n_conf=20)
    pruned = prune_rmsd(mol, conf_data, rmsd_cutoff=0.5)
    in_window = filter_energy_window(pruned, window_kcal=10.0)
    print(f'Generated: {len(conf_data)}, after RMSD pruning: {len(pruned)}, in 10 kcal/mol window: {len(in_window)}')
