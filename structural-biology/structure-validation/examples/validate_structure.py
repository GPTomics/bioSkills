'''Validate a deposited structure: read refinement metadata, flag the R-free overfitting
gap, sanity-check B-factors, and screen backbone geometry for Ramachandran and cis-peptide
outliers. The Python screen decides whether to run MolProbity; MolProbity gives the numbers.'''
# Reference: biopython 1.85+ | Verify API if version differs

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, calc_dihedral
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# Thresholds with rationale (all sourced; verify current MolProbity targets before reporting):
RFREE_GAP_CONCERN = 0.05      # R-free minus R-work above this flags overfitting at ~2 A (Kleywegt & Brunger 1996)
BFACTOR_Z_CUT = 2.0           # robust z on per-residue B; high B marks disorder/model error, not portable dynamics
CIS_OMEGA_DEG = 30.0          # |omega| < 30 deg is a cis peptide bond (trans ~180, cis ~0); cis non-Pro is usually an error
RAMA_FAVORED_GOAL = 98.0      # MolProbity goal: >98% of residues in favored phi/psi basins (Williams 2018)
FSC_HALFMAP = 0.143           # cryo-EM half-map FSC resolution criterion (Rosenthal & Henderson 2003)
FSC_MAPMODEL = 0.5            # cryo-EM map-vs-model FSC criterion (distinct from the resolution claim)

# Coarse general-allowed phi/psi basins (deg): alpha, beta/PPII, left-handed. NOT MolProbity's
# rama8000 contours - this only screens gross outliers to decide whether to run phenix.molprobity.
_BASINS = [(-160, -40, -80, 30), (-180, -40, 90, 180), (30, 90, -30, 90)]


def read_refinement_metadata(cif_path):
    meta = MMCIF2Dict(cif_path)
    def first(key):
        val = meta.get(key, ['NA'])[0]
        try:
            return float(val)
        except ValueError:
            return val
    resolution, r_work, r_free = first('_refine.ls_d_res_high'), first('_refine.ls_R_factor_R_work'), first('_refine.ls_R_factor_R_free')
    gap = r_free - r_work if isinstance(r_free, float) and isinstance(r_work, float) else None
    return {'method': meta.get('_exptl.method', ['NA'])[0], 'resolution': resolution,
            'r_work': r_work, 'r_free': r_free, 'gap': gap,
            'overfit_flag': gap is not None and gap > RFREE_GAP_CONCERN}


def bfactor_outliers(structure, chain_id):
    residue_b = {res.id: np.mean([a.get_bfactor() for a in res]) for res in structure[0][chain_id] if res.id[0] == ' '}
    vals = np.array(list(residue_b.values()))
    median, mad = np.median(vals), np.median(np.abs(vals - np.median(vals))) + 1e-9
    return {rid: b for rid, b in residue_b.items() if (b - median) / (1.4826 * mad) > BFACTOR_Z_CUT}


def is_rama_allowed(phi, psi):
    d = np.degrees([phi, psi])
    return any(lo_p <= d[0] <= hi_p and lo_s <= d[1] <= hi_s for lo_p, hi_p, lo_s, hi_s in _BASINS)


def geometry_flags(structure):
    outliers, cis_nonpro = [], []
    for pp in PPBuilder().build_peptides(structure[0]):
        for residue, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            if phi is not None and psi is not None and not is_rama_allowed(phi, psi):
                outliers.append((residue.get_parent().id, residue.id[1], residue.resname))
        residues = list(pp)
        for prev, curr in zip(residues, residues[1:]):
            omega = calc_dihedral(prev['CA'].get_vector(), prev['C'].get_vector(), curr['N'].get_vector(), curr['CA'].get_vector())
            if abs(np.degrees(omega)) < CIS_OMEGA_DEG and curr.resname != 'PRO':
                cis_nonpro.append((curr.get_parent().id, curr.id[1], curr.resname))
    return {'rama_outlier_candidates': outliers, 'cis_nonproline': cis_nonpro}


def validate(path, chain_id, is_cif=True):
    parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
    structure = parser.get_structure('model', path)
    report = {'geometry': geometry_flags(structure), 'high_b_residues': list(bfactor_outliers(structure, chain_id))}
    if is_cif:
        report['refinement'] = read_refinement_metadata(path)
    return report


if __name__ == '__main__':
    import sys
    result = validate(sys.argv[1], sys.argv[2], is_cif=sys.argv[1].endswith('.cif'))
    result
