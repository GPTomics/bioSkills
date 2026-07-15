# Reference: PoseBusters 0.6+, RDKit 2024.09+, pandas 2.2+ | Verify API if version differs

from posebusters import PoseBusters
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd


def run_posebusters(pred_sdf, receptor_pdb, ref_sdf=None):
    '''Run PoseBusters and return per-pose validity DataFrame.'''
    config = 'redock' if ref_sdf else 'dock'
    bust = PoseBusters(config=config)
    if ref_sdf:
        results = bust.bust(mol_pred=pred_sdf, mol_true=ref_sdf, mol_cond=receptor_pdb)
    else:
        results = bust.bust(mol_pred=pred_sdf, mol_cond=receptor_pdb)

    # Collect boolean check columns (PoseBusters returns metadata + bool checks)
    bool_cols = [
        col for col in results.select_dtypes(include='bool').columns
        if not col.lower().startswith('rmsd')
    ]
    results['pb_valid'] = results[bool_cols].all(axis=1)
    return results


def ligand_strain_mmff(docked_sdf, n_ref_conf=20):
    '''Compute relative MMFF94 strain versus the lowest sampled conformer.'''
    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False, sanitize=True)
    rows = []
    for i, docked in enumerate(suppl):
        if docked is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'parse_fail'})
            continue

        # Generate reference unconstrained ensemble for the same chemistry
        smi = Chem.MolToSmiles(docked)
        ref_mol = Chem.MolFromSmiles(smi)
        if ref_mol is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'reference_parse_fail'})
            continue
        ref_mol = Chem.AddHs(ref_mol)

        mmff_props_ref = AllChem.MMFFGetMoleculeProperties(ref_mol)
        if mmff_props_ref is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'no_mmff_params'})
            continue
        conf_ids = list(AllChem.EmbedMultipleConfs(
            ref_mol, numConfs=n_ref_conf, params=AllChem.ETKDGv3()
        ))
        if not conf_ids:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'embedding_failed'})
            continue
        AllChem.MMFFOptimizeMoleculeConfs(ref_mol)

        ref_energies = []
        for c in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(ref_mol, mmff_props_ref, confId=c)
            if ff is not None:
                ref_energies.append(ff.CalcEnergy())
        if not ref_energies:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'reference_force_field_failed'})
            continue
        min_ref = min(ref_energies)

        # MMFF energies are comparable only for the same explicit atom system.
        # Add any missing H coordinates and relax H atoms while preserving the
        # docked heavy-atom pose.
        docked_h = Chem.AddHs(Chem.Mol(docked), addCoords=True)
        if docked_h.GetNumAtoms() != ref_mol.GetNumAtoms():
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'atom_system_mismatch'})
            continue
        mmff_props_dock = AllChem.MMFFGetMoleculeProperties(docked_h)
        if mmff_props_dock is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'no_dock_mmff_params'})
            continue
        ff_dock = AllChem.MMFFGetMoleculeForceField(docked_h, mmff_props_dock)
        if ff_dock is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'dock_force_field_failed'})
            continue
        for atom in docked_h.GetAtoms():
            if atom.GetAtomicNum() != 1:
                ff_dock.AddFixedPoint(atom.GetIdx())
        ff_dock.Minimize(maxIts=200)
        docked_e = ff_dock.CalcEnergy()

        rows.append({
            'pose_idx': i,
            'strain_kcal': docked_e - min_ref,
            'note': 'ok',
        })
    return pd.DataFrame(rows)


def pose_qc_pipeline(docked_sdf, receptor_pdb, strain_cutoff=None, ref_sdf=None):
    '''Run PoseBusters and report relative strain with an optional user cutoff.'''
    pb_df = run_posebusters(docked_sdf, receptor_pdb, ref_sdf=ref_sdf)
    pb_df['pose_idx'] = range(len(pb_df))

    strain_df = ligand_strain_mmff(docked_sdf)

    combined = pb_df.merge(strain_df, on='pose_idx', how='left')
    if strain_cutoff is not None:
        combined['passes_user_strain_cutoff'] = (
            combined['strain_kcal'].notna()
            & (combined['strain_kcal'] <= strain_cutoff)
        )
    return combined


if __name__ == '__main__':
    # Example usage (requires real SDF + PDB)
    # results = pose_qc_pipeline('docked.sdf', 'receptor.pdb')
    # print(results[['pose_idx', 'pb_valid', 'strain_kcal']])
    print('Example: provide docked.sdf and receptor.pdb to run')
