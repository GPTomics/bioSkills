# Reference: OpenFE 1.7+, OpenMM 8.1+, alchemlyb 2.1+, pymbar 4.0+ | Verify API if version differs

# Note: OpenFE setup is complex; this is a sketch of the API call pattern.
# Full reproduction requires GPU + multi-hour run.

def estimate_rbfe_pair(
    receptor_pdb,
    lig1_sdf,
    lig2_sdf,
    n_lambdas=12,
    sim_time_ns=5,
):
    '''
    Build version-checked OpenFE RBFE transformations.

    Execute the returned transformations with their protocol DAGs or serialize
    them and use ``openfe quickrun``. This setup function does not run MD.
    '''
    from openfe import (
        SmallMoleculeComponent,
        ProteinComponent,
        SolventComponent,
        ChemicalSystem,
        Transformation,
    )
    from openfe.setup import LomapAtomMapper
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

    protein = ProteinComponent.from_pdb_file(receptor_pdb)
    lig1 = SmallMoleculeComponent.from_sdf_file(lig1_sdf)
    lig2 = SmallMoleculeComponent.from_sdf_file(lig2_sdf)
    solvent = SolventComponent()

    system_bound_1 = ChemicalSystem({'protein': protein, 'ligand': lig1, 'solvent': solvent})
    system_bound_2 = ChemicalSystem({'protein': protein, 'ligand': lig2, 'solvent': solvent})
    system_solv_1 = ChemicalSystem({'ligand': lig1, 'solvent': solvent})
    system_solv_2 = ChemicalSystem({'ligand': lig2, 'solvent': solvent})

    # Twelve replicas and 5 ns are repository setup defaults only; adjust both
    # from state-overlap, exchange, and replicate-convergence diagnostics.
    solvent_settings = RelativeHybridTopologyProtocol.default_settings()
    solvent_settings.lambda_settings.lambda_windows = n_lambdas
    solvent_settings.simulation_settings.n_replicas = n_lambdas
    solvent_settings.simulation_settings.production_length = f'{sim_time_ns} ns'
    solvent_protocol = RelativeHybridTopologyProtocol(solvent_settings)

    complex_settings = RelativeHybridTopologyProtocol.default_settings()
    complex_settings.lambda_settings.lambda_windows = n_lambdas
    complex_settings.simulation_settings.n_replicas = n_lambdas
    complex_settings.simulation_settings.production_length = f'{sim_time_ns} ns'
    # OpenFE 1.7 recommends reducing complex-leg padding from its solvent default.
    complex_settings.solvation_settings.solvent_padding = '1 nm'
    complex_protocol = RelativeHybridTopologyProtocol(complex_settings)
    # LOMAP is an explicitly selected supported mapper. OpenFE 1.7's CLI
    # default is Kartograf. Inspect the proposed mapping before execution.
    mapping = next(LomapAtomMapper().suggest_mappings(lig1, lig2))
    bound_transformation = Transformation(
        stateA=system_bound_1,
        stateB=system_bound_2,
        mapping=mapping,
        protocol=complex_protocol,
        name='lig1_to_lig2_bound',
    )
    solvent_transformation = Transformation(
        stateA=system_solv_1,
        stateB=system_solv_2,
        mapping=mapping,
        protocol=solvent_protocol,
        name='lig1_to_lig2_solvent',
    )
    return {
        'note': 'setup_only; inspect mappings, then execute both transformations',
        'bound_transformation': bound_transformation,
        'solvent_transformation': solvent_transformation,
        'n_lambdas': n_lambdas,
        'sim_time_ns': sim_time_ns,
    }


def analyze_mbar(u_nk_files, T=300.0):
    '''Analyze GROMACS XVG files; pass the simulation temperature (300 K is a starting default).'''
    from alchemlyb.parsing import gmx
    from alchemlyb import concat
    from alchemlyb.estimators import MBAR
    from alchemlyb.postprocessors.units import to_kcalmol

    u_nks = []
    for f in u_nk_files:
        u_nks.append(gmx.extract_u_nk(f, T=T))
    u_nk = concat(u_nks)
    mbar = MBAR().fit(u_nk)
    delta_g = to_kcalmol(mbar.delta_f_).iloc[0, -1]
    d_delta_g = to_kcalmol(mbar.d_delta_f_).iloc[0, -1]
    return delta_g, d_delta_g


def cycle_closure_residual(edges):
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
