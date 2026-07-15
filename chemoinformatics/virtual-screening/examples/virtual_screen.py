#!/usr/bin/env python3
'''
Virtual screening with AutoDock Vina.
'''
# Reference: autodock vina 1.2+, rdkit 2024.09+, pandas 2.2+ | Verify API if version differs

import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy


def prepare_ligand(smiles, output_pdbqt):
    '''Prepare a supplied, project-approved protomer for docking.'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'invalid SMILES: {smiles}')
    mol = Chem.AddHs(mol)
    embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if embed_status != 0:
        raise RuntimeError('ETKDGv3 failed to generate a ligand conformer')
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError('MMFF94 parameters are unavailable for this ligand')
    optimization_status = AllChem.MMFFOptimizeMolecule(mol)
    if optimization_status != 0:
        raise RuntimeError('MMFF94 ligand optimization did not converge')

    setups = MoleculePreparation().prepare(mol)
    pdbqt_text, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        raise RuntimeError(f'Meeko PDBQT export failed: {err}')
    Path(output_pdbqt).write_text(pdbqt_text)

    return output_pdbqt


def prepare_receptor(protonated_pdb, output_pdbqt):
    '''Convert an already reviewed/repaired/protonated receptor with Meeko.'''
    output_basename = str(Path(output_pdbqt).with_suffix(''))
    subprocess.run([
        'mk_prepare_receptor.py', '--read_pdb', protonated_pdb,
        '-o', output_basename, '-p'
    ], check=True)
    return output_pdbqt


def dock_single(receptor_pdbqt, ligand_pdbqt, center, box_size, exhaustiveness=8, n_poses=10):
    '''Dock a single ligand using AutoDock Vina.'''
    from vina import Vina

    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    energies = v.energies()
    return energies


def virtual_screen(receptor_pdbqt, ligand_dict, center, box_size, output_dir, exhaustiveness=8):
    '''
    Screen compound library against receptor.

    Args:
        receptor_pdbqt: Prepared receptor
        ligand_dict: Dict of {name: smiles}
        center: (x, y, z) binding site center
        box_size: (x, y, z) search box dimensions
        output_dir: Directory for output
        exhaustiveness: Search thoroughness
    '''
    from vina import Vina
    import pandas as pd

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)

    results = []

    for name, smiles in ligand_dict.items():
        try:
            ligand_pdbqt = f'{output_dir}/{name}.pdbqt'
            prepare_ligand(smiles, ligand_pdbqt)

            v.set_ligand_from_file(ligand_pdbqt)
            v.dock(exhaustiveness=exhaustiveness, n_poses=5)

            energies = v.energies()
            has_poses = len(energies) > 0
            best_affinity = energies[0][0] if has_poses else None

            if has_poses:
                v.write_poses(f'{output_dir}/{name}_poses.pdbqt', n_poses=5)

            results.append({
                'name': name,
                'smiles': smiles,
                'affinity_kcal_mol': best_affinity
            })

        except Exception as e:
            results.append({
                'name': name,
                'smiles': smiles,
                'affinity_kcal_mol': None,
                'error': str(e)
            })

    df = pd.DataFrame(results)
    df = df.sort_values('affinity_kcal_mol')
    df.to_csv(f'{output_dir}/results.csv', index=False)

    return df


def find_binding_site(receptor_pdb, ligand_pdb, padding=5.0):
    '''Define binding site from co-crystallized ligand.'''
    mol = Chem.MolFromPDBFile(ligand_pdb)
    conf = mol.GetConformer()
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]

    x = [c.x for c in coords]
    y = [c.y for c in coords]
    z = [c.z for c in coords]

    center = (sum(x)/len(x), sum(y)/len(y), sum(z)/len(z))
    box_size = (max(x)-min(x)+2*padding, max(y)-min(y)+2*padding, max(z)-min(z)+2*padding)

    return center, box_size


if __name__ == '__main__':
    print('Virtual Screening Pipeline')
    print('=' * 40)
    print('1. Prepare receptor: prepare_receptor(pdb, pdbqt)')
    print('2. Define binding site: find_binding_site(receptor, ligand)')
    print('3. Screen library: virtual_screen(receptor, compounds, center, box)')
    print()
    print('Requirements: pip install vina rdkit meeko')
