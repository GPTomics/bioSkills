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
    bool_cols = results.select_dtypes(include='bool').columns
    results['pb_valid'] = results[bool_cols].all(axis=1)
    return results


def ligand_strain_mmff(docked_sdf, n_ref_conf=20):
    '''Compute MMFF94 strain per docked pose vs lowest-energy unconstrained conformer.'''
    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False, sanitize=True)
    rows = []
    for i, docked in enumerate(suppl):
        if docked is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'parse_fail'})
            continue

        # Generate reference unconstrained ensemble for the same chemistry
        smi = Chem.MolToSmiles(docked)
        ref_mol = Chem.MolFromSmiles(smi)
        ref_mol = Chem.AddHs(ref_mol)
        AllChem.EmbedMultipleConfs(ref_mol, numConfs=n_ref_conf,
                                    params=AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMoleculeConfs(ref_mol)

        mmff_props_ref = AllChem.MMFFGetMoleculeProperties(ref_mol)
        if mmff_props_ref is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'no_mmff_params'})
            continue

        ref_energies = []
        for c in range(ref_mol.GetNumConformers()):
            ff = AllChem.MMFFGetMoleculeForceField(ref_mol, mmff_props_ref, confId=c)
            ref_energies.append(ff.CalcEnergy())
        min_ref = min(ref_energies)

        # MMFF94 energy of docked pose
        mmff_props_dock = AllChem.MMFFGetMoleculeProperties(docked)
        if mmff_props_dock is None:
            rows.append({'pose_idx': i, 'strain_kcal': None, 'note': 'no_dock_mmff_params'})
            continue
        ff_dock = AllChem.MMFFGetMoleculeForceField(docked, mmff_props_dock)
        docked_e = ff_dock.CalcEnergy()

        rows.append({
            'pose_idx': i,
            'strain_kcal': docked_e - min_ref,
            'note': 'ok',
        })
    return pd.DataFrame(rows)


def pose_qc_pipeline(docked_sdf, receptor_pdb, strain_cutoff=5.0, ref_sdf=None):
    '''Full pose QC: PoseBusters + strain; flag PB-valid AND strain < cutoff.'''
    pb_df = run_posebusters(docked_sdf, receptor_pdb, ref_sdf=ref_sdf)
    pb_df['pose_idx'] = range(len(pb_df))

    strain_df = ligand_strain_mmff(docked_sdf)

    combined = pb_df.merge(strain_df, on='pose_idx', how='left')
    # Strain cutoff: 5 kcal/mol is realistic upper bound for bioactive conformers (Boström 1998)
    combined['strain_ok'] = combined['strain_kcal'].fillna(999) <= strain_cutoff
    combined['fep_ready'] = combined['pb_valid'] & combined['strain_ok']
    return combined


if __name__ == '__main__':
    # Example usage (requires real SDF + PDB)
    # results = pose_qc_pipeline('docked.sdf', 'receptor.pdb')
    # print(results[['pose_idx', 'pb_valid', 'strain_kcal', 'fep_ready']])
    print('Example: provide docked.sdf and receptor.pdb to run')
