# Reference: OpenFE 1.7+, OpenMM 8.1+, alchemlyb 2.1+, pymbar 4.0+ | Verify API if version differs

# Note: OpenFE setup is complex; this is a sketch of the API call pattern.
# Full reproduction requires GPU + multi-hour run.

import pandas as pd
from rdkit import Chem


def estimate_rbfe_pair(
    receptor_pdb,
    lig1_sdf,
    lig2_sdf,
    n_lambdas=12,
    sim_time_ns=5,
):
    '''
    Sketch of OpenFE RBFE setup. In production, use openfe.setup + openfe.execute.

    Returns delta-delta-G with statistical uncertainty estimate.
    '''
    from openfe import (
        SmallMoleculeComponent,
        ProteinComponent,
        SolventComponent,
        ChemicalSystem,
    )
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

    protein = ProteinComponent.from_pdb_file(receptor_pdb)
    lig1 = SmallMoleculeComponent.from_sdf_file(lig1_sdf)
    lig2 = SmallMoleculeComponent.from_sdf_file(lig2_sdf)
    solvent = SolventComponent()

    system_bound_1 = ChemicalSystem({'protein': protein, 'ligand': lig1, 'solvent': solvent})
    system_bound_2 = ChemicalSystem({'protein': protein, 'ligand': lig2, 'solvent': solvent})
    system_solv_1 = ChemicalSystem({'ligand': lig1, 'solvent': solvent})
    system_solv_2 = ChemicalSystem({'ligand': lig2, 'solvent': solvent})

    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.alchemical_settings.lambda_steps = n_lambdas
    settings.simulation_settings.production_length = f'{sim_time_ns} * nanosecond'

    protocol = RelativeHybridTopologyProtocol(settings)
    return {
        'note': 'sketch_only; run full OpenFE setup script + cluster execution',
        'protocol': protocol,
        'n_lambdas': n_lambdas,
        'sim_time_ns': sim_time_ns,
    }


def analyze_mbar(u_nk_files, T=300.0):
    '''Analyze GROMACS / OpenMM xvg files with alchemlyb + pymbar.'''
    from alchemlyb.parsing import gmx
    from alchemlyb.estimators import MBAR

    u_nks = []
    for f in u_nk_files:
        u_nks.append(gmx.extract_u_nk(f, T=T))
    u_nk = pd.concat(u_nks)
    mbar = MBAR().fit(u_nk)
    delta_g = mbar.delta_f_.iloc[0, -1]
    d_delta_g = mbar.d_delta_f_.iloc[0, -1]
    return delta_g, d_delta_g


def cycle_closure(edges):
    '''
    Compute cycle closure error for a closed cycle of RBFE edges.

    edges: list of (lig_a, lig_b, delta_g, sd)
    Returns total (should be ~0 if closure is good) + propagated SD.
    '''
    total = sum(e[2] for e in edges)
    total_var = sum(e[3] ** 2 for e in edges)
    return total, total_var ** 0.5


if __name__ == '__main__':
    # Demonstration only; full run requires SDF, PDB, GPU cluster
    print('See full OpenFE setup at https://docs.openfree.energy/')
