---
name: bio-structural-biology-structure-validation
description: Judges whether a macromolecular model (or a region of it) is reliable enough to build on, using resolution, R-free, B-factors, MolProbity geometry, and predicted-model confidence with Bio.PDB. Use when deciding if a structure or a specific region is trustworthy before docking/mechanism/measurement; reading resolution, R-work vs R-free and the R-free-minus-R-work overfitting gap; sanity-checking per-residue and mean B-factors; flagging clashscore, Ramachandran and rotamer outliers and cis non-proline peptides; validating a PREDICTED (AlphaFold/ESMFold) model via pLDDT bands and PAE before docking or molecular replacement; and interpreting cryo-EM global-vs-local resolution (FSC 0.143 half-map vs 0.5 map-model) or an NMR ensemble spread. Keywords validation, resolution, R-free, B-factor, MolProbity, clashscore, Ramachandran, rotamer, pLDDT, PAE, wwPDB, cryo-EM local resolution.
tool_type: python
primary_tool: Bio.PDB
---

## Version Compatibility

Reference examples tested with: biopython 1.85+, numpy 1.26+

MolProbity, phenix (`phenix.molprobity`, `phenix.process_predicted_model`), and DSSP (`mkdssp`) are external CLI tools invoked via `subprocess`, not pip packages; install them separately and confirm they are on PATH before use.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structure Validation

**"Is this structure good enough to build on?"** -> Read the data-quality metadata, run geometry validation, and decide per-region rather than per-file.
- Python: `Bio.PDB.MMCIF2Dict` (resolution/R-free), `Bio.PDB.calc_dihedral` (phi/psi/omega), `subprocess` -> `phenix.molprobity` (clashscore, rotamer, Ramachandran).

## Governing Principle

A coordinate file is an INTERPRETED MODEL fit to experimental data, not the data and not ground truth, and its reliability varies PER-ATOM within one file. The single global resolution is a DATA CEILING on what the map can resolve, never a per-region quality certificate: a 1.5 Angstrom structure can still carry a guesswork surface loop, and a 3.2 Angstrom structure can have a rigid, locally excellent active site. The number that answers "can I trust THIS region" is the LOCAL signal - the per-residue B-factor and the real-space fit (RSRZ) for the residues actually used - not the headline resolution. Validate before any geometric interpretation: clashscore, Ramachandran and rotamer outliers, cis-peptides, and bond/angle deviations tell whether the model is even self-consistent before a distance or an angle computed on it means anything.

R-free, not R-work, is the cross-validation statistic. R-work is computed on the reflections used in refinement and can be driven down by adding parameters (waters, alt-confs, loose restraints) that fit noise; R-free is the same metric on a held-out test set that never entered refinement (Brunger 1992 *Nature* 355:472-475). The R-free-MINUS-R-work GAP is the overfitting flag - a gap much larger than expected for the resolution signals a model fitting its own noise, and a suspiciously small gap signals test-set leakage (Kleywegt & Brunger 1996 *Structure* 4:897-904). B-factors conflate genuine thermal motion, static disorder, and MODEL ERROR into one number, so they are a within-structure relative signal ("which atoms are least certain here"), NOT a portable cross-structure dynamics readout; normalize before comparing across structures and exclude high-B atoms (say >60-80 Angstrom^2 at moderate resolution) from precise geometric claims.

For a PREDICTED model the confidence scores bound SELF-CONSISTENCY, NOT correctness. pLDDT is per-residue and PAE is inter-residue confidence in the model's OWN frame (Jumper et al. 2021 *Nature* 596:583-589); a confident model can be confidently wrong when AlphaFold modeled a monomer/apo state where the truth is a complex/holo/alternative state. Trim low-pLDDT residues and split by PAE-defined domains before docking or molecular replacement (Oeffner et al. 2022 *Acta Crystallogr D* 78:1303-1314). Cryo-EM has GLOBAL vs LOCAL resolution - flexible peripheries in a big map are often docked hypotheses at 6-8 Angstrom local resolution inside a "2.5 Angstrom" map - and TWO different FSC thresholds: 0.143 on the gold-standard half-map FSC for the resolution claim, 0.5 for map-vs-model agreement (Rosenthal & Henderson 2003 *J Mol Biol* 333:721-745). An NMR deposit is an ENSEMBLE of models; validate per-model and report the per-residue spread - never average coordinates, which produces a physically impossible structure with distorted bonds and clashes (Montelione et al. 2013 *Structure* 21:1563-1570).

## Decision: what the resolution buys (X-ray / cryo-EM global)

| Resolution (A) | Reliably resolvable | Do NOT over-read |
|---|---|---|
| < 1.2 (atomic) | Individual atoms, H atoms, anisotropic ADPs, alt-confs | - |
| 1.2-1.8 | Side-chain rotamers, ordered waters, alt-confs | Surface/loop atoms with high B still uncertain |
| 1.8-2.5 (typical) | Backbone plus most side chains; fold solid | Long polar rotamers, water networks, weak ligand density are model-dependent |
| 2.5-3.5 | Backbone trace, secondary structure, domain arrangement | Side-chain conformations, exact ligand pose, H-bond geometry are interpretive |
| > 3.5 | Overall fold, gross assembly | Individual side chains largely placed by geometry; treat atomic detail as hypothesis |

## Decision: R-free and the overfitting gap

| Resolution | Typical R-work | Typical R-free | Concern if gap (R-free - R-work) exceeds ~ |
|---|---|---|---|
| ~1.5 A | 0.13-0.17 | 0.15-0.19 | 0.03-0.04 |
| ~2.0 A | 0.16-0.20 | 0.19-0.24 | 0.04-0.05 |
| ~2.5-3.0 A | 0.18-0.24 | 0.23-0.29 | 0.05-0.06 |

Absolute R rises with resolution (weaker high-angle reflections are noisier); a gap wider than the row's threshold flags overfitting, a near-zero gap flags work/test contamination. These are heuristics - verify against contemporaneous PDB statistics for the resolution.

## Decision: geometry validation targets (MolProbity)

| Metric | Measures | Target (goal) | Concern | Source |
|---|---|---|---|---|
| Clashscore | Serious all-atom overlaps (>0.4 A) per 1000 atoms | low, high same-resolution percentile | high vs same-resolution peers | Chen 2010 |
| Ramachandran favored | % residues in favored phi/psi basins | > 98% | < 95% | Williams 2018 |
| Ramachandran outliers | % residues in disallowed phi/psi | < 0.2% (goal 0) | > 0.5% | Williams 2018 |
| Poor rotamers | % side chains in disallowed rotamers | < 0.3% | > 1.5% | Williams 2018 |
| RSRZ outliers (X-ray) | Residues poorly fitting local density (RSRZ > 2) | few, dispersed | many, or clustered in the region of interest | Gore 2017 |
| MolProbity score | Composite mapped to a resolution-equivalent | <= deposited resolution | >> deposited resolution | Chen 2010 |

The poor-rotamer target tightened from ~1.0% to 0.3% with the updated reference distributions (Williams 2018). Read the wwPDB report's PERCENTILE sliders alongside the raw numbers: they rank an entry against the whole archive AND against same-resolution entries (Gore et al. 2017 *Structure* 25:1916-1927).

## Decision: how to validate this model (experimental vs predicted fork)

| The model is | Trust question | What to read | Authoritative route |
|---|---|---|---|
| X-ray | Is the region well-determined? | resolution + local B + RSRZ; R-free and its gap | wwPDB report + `phenix.molprobity` |
| Cryo-EM | Is the region rigid or a docked guess? | LOCAL resolution map, not global; map-model FSC (0.5) | EMDB local-resolution + map-model FSC |
| NMR | How wide is the conformational spread? | per-residue RMSD across all models | validate each model; report spread |
| Predicted (AlphaFold/ESMFold) | Is it self-consistent AND in the right biological context? | pLDDT bands + PAE blocks; is it apo/monomer where truth is holo/complex? | `phenix.process_predicted_model` (trim + PAE split) |

## Read Validation Metadata From the mmCIF Header

**Goal:** Pull resolution, R-work, R-free, and method from a deposited entry and flag the overfitting gap - the numbers Bio.PDB's thin `structure.header` omits.

**Approach:** Read the raw mmCIF categories with `MMCIF2Dict` (which reaches anything in the file), cast the strings, and compare the R-free-minus-R-work gap against a resolution-scaled expectation.

```python
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def read_refinement_metadata(cif_path):
    meta = MMCIF2Dict(cif_path)
    def first(key):
        val = meta.get(key, ['NA'])[0]
        try:
            return float(val)
        except ValueError:
            return val
    method = meta.get('_exptl.method', ['NA'])[0]
    resolution = first('_refine.ls_d_res_high')
    r_work = first('_refine.ls_R_factor_R_work')
    r_free = first('_refine.ls_R_factor_R_free')
    gap = r_free - r_work if isinstance(r_free, float) and isinstance(r_work, float) else None
    # A gap wider than ~0.05 flags overfitting at typical (~2 A) resolution (Kleywegt & Brunger 1996).
    overfit_flag = gap is not None and gap > 0.05
    return {'method': method, 'resolution': resolution, 'r_work': r_work, 'r_free': r_free, 'gap': gap, 'overfit_flag': overfit_flag}
```

## Sanity-Check B-factors Within One Structure

**Goal:** Locate the least-certain atoms of THIS model so precise-distance claims avoid them - a relative, within-structure read, never a cross-structure comparison.

**Approach:** Collect per-residue mean B for a chain, then flag residues whose B sits far above the structure's own median as low-confidence.

```python
import numpy as np
from Bio.PDB import PDBParser

def bfactor_outliers(structure, chain_id, z_cut=2.0):
    residue_b = {}
    for residue in structure[0][chain_id]:
        if residue.id[0] != ' ':  # skip hetero/water; validate the polymer only
            continue
        residue_b[residue.id] = np.mean([a.get_bfactor() for a in residue])
    vals = np.array(list(residue_b.values()))
    median, mad = np.median(vals), np.median(np.abs(vals - np.median(vals))) + 1e-9
    # Robust z on B: |B - median| / (1.4826*MAD); high B marks disorder/model error, not portable dynamics.
    return {rid: b for rid, b in residue_b.items() if (b - median) / (1.4826 * mad) > z_cut}
```

## Flag Ramachandran Outliers and cis Non-Proline Peptides

**Goal:** Catch backbone geometry that is almost always a modeling error - residues in disallowed phi/psi and cis peptide bonds that are not proline.

**Approach:** Get per-residue phi/psi from `PPBuilder`, coarsely classify against the canonical basins, and compute omega (Ca-C-N-Ca) directly to find cis bonds; a cis assignment at a non-proline is a validation flag until proven by density.

```python
import numpy as np
from Bio.PDB import PDBParser, PPBuilder, calc_dihedral

# Coarse general-allowed basins (deg): alpha, beta/PPII, left-handed. Authoritative
# favored/allowed/outlier percentages need MolProbity's rama8000 contours (Williams 2018);
# this only screens gross outliers to decide whether to run MolProbity.
_BASINS = [(-160, -40, -80, 30), (-180, -40, 90, 180), (30, 90, -30, 90)]

def is_rama_allowed(phi, psi):
    d = np.degrees([phi, psi])
    return any(lo_p <= d[0] <= hi_p and lo_s <= d[1] <= hi_s for lo_p, hi_p, lo_s, hi_s in _BASINS)

def geometry_flags(structure):
    outliers, cis_nonpro = [], []
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure[0]):
        for residue, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            if phi is not None and psi is not None and not is_rama_allowed(phi, psi):
                outliers.append((residue.get_parent().id, residue.id[1], residue.resname))
        residues = list(pp)
        for prev, curr in zip(residues, residues[1:]):
            omega = calc_dihedral(prev['CA'].get_vector(), prev['C'].get_vector(), curr['N'].get_vector(), curr['CA'].get_vector())
            # omega ~180 = trans, ~0 = cis; |omega| < 30 deg is cis. cis at non-Pro is rare and usually an error.
            if abs(np.degrees(omega)) < 30 and curr.resname != 'PRO':
                cis_nonpro.append((curr.get_parent().id, curr.id[1], curr.resname))
    return {'rama_outliers': outliers, 'cis_nonproline': cis_nonpro}
```

## Run MolProbity for Authoritative Geometry (the real answer)

**Goal:** Get archive-calibrated clashscore, rotamer, and Ramachandran outlier percentages rather than the coarse Python screen above.

**Approach:** Shell out to `phenix.molprobity` (or the MolProbity web service / `molprobity.molprobity`), which carries the reference contour and rotamer distributions the Python screen cannot reproduce.

```python
import subprocess

def run_molprobity(model_path):
    # phenix.molprobity writes molprobity.out with clashscore, rotamer_outliers,
    # ramachandran_outliers, ramachandran_favored, molprobity_score. Parse that file.
    result = subprocess.run(['phenix.molprobity', model_path], capture_output=True, text=True)
    return result.stdout
```

The coarse Python screen decides WHETHER to run MolProbity; MolProbity (Chen et al. 2010 *Acta Crystallogr D* 66:12-21; Williams et al. 2018 *Protein Sci* 27:293-315) gives the numbers to report. For a deposited entry, prefer the pre-computed wwPDB validation report (`https://files.rcsb.org/pub/pdb/validation_reports/<xy>/<id>/<id>_validation.xml.gz`) - it already carries the percentile sliders and per-residue RSRZ.

## Validate a Predicted Model Before Docking or MR

**Goal:** Turn an AlphaFold/ESMFold model into a trustworthy input by trimming low-confidence residues and reading inter-domain confidence - because pLDDT/PAE bound self-consistency, not correctness.

**Approach:** Read pLDDT from the B-factor column into bands, read the PAE matrix for domain segmentation, then hand the raw file to `phenix.process_predicted_model` (which converts pLDDT to a pseudo-B, trims, and splits by PAE-defined domains).

```python
import json
import numpy as np
import subprocess
from Bio.PDB import MMCIFParser

_PLDDT_BANDS = [(90, 'very_high'), (70, 'confident'), (50, 'low'), (0, 'very_low')]

def plddt_bands(cif_path):
    # pLDDT rides in the B-factor column but is confidence (high = good), OPPOSITE polarity to a real B-factor.
    structure = MMCIFParser(QUIET=True).get_structure('pred', cif_path)
    counts = {label: 0 for _, label in _PLDDT_BANDS}
    for residue in structure[0].get_residues():
        if 'CA' not in residue:
            continue
        score = residue['CA'].get_bfactor()
        counts[next(label for cut, label in _PLDDT_BANDS if score >= cut)] += 1
    return counts  # a long very_low run usually flags an intrinsically disordered region, not an error

def pae_interdomain_confident(pae_json, block_a, block_b, cutoff=5.0):
    # Off-diagonal PAE (A) between two domain blocks; low = relative orientation trusted, high = independent guess.
    pae = np.array(json.load(open(pae_json))[0]['predicted_aligned_error'])
    return pae[np.ix_(block_a, block_b)].mean() < cutoff

def process_for_mr(model_path):
    # Trims below ~0.7 fractional pLDDT, converts pLDDT->pseudo-B, splits into PAE-defined domains (Oeffner 2022).
    return subprocess.run(['phenix.process_predicted_model', model_path], capture_output=True, text=True).stdout
```

## Report an NMR Ensemble Spread (never average coordinates)

**Goal:** Turn a multi-model NMR deposit into a per-residue uncertainty map instead of over-claiming precision from model 1.

**Approach:** Superpose all models on a reference and report per-residue Ca RMSD across the ensemble; wide spread marks flexible or under-restrained regions.

```python
import numpy as np
from Bio.PDB import PDBParser

def ensemble_ca_spread(structure, chain_id):
    coords = []
    for model in structure:
        coords.append(np.array([res['CA'].get_coord() for res in model[chain_id] if 'CA' in res]))
    stack = np.stack(coords)  # (n_models, n_residues, 3); assumes consistent residue set across models
    # Per-residue spread = mean distance of each model's Ca from the ensemble mean position.
    return np.linalg.norm(stack - stack.mean(axis=0), axis=2).mean(axis=0)
```

Secondary-structure validation (DSSP, `Bio.PDB.DSSP(model, path, dssp='mkdssp')`) needs backbone geometry to place the amide H and computes its own H-bond energy; it processes only the first model and its output drifts across the dssp->mkdssp v2->v4 rewrites, so name the version (Kabsch & Sander 1983 *Biopolymers* 22:2577-2637). See geometric-analysis for dihedral and DSSP mechanics.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| "It is a 1.5 A structure so every atom is accurate" | Global resolution read as per-region quality | Read the local B-factor and RSRZ for the specific residues used; resolution is a data ceiling |
| Low R-work reported as proof of a good model | R-work is fit on refinement reflections and rewards overfitting | Report R-free (held-out) and the R-free-minus-R-work gap; a wide gap flags overfitting |
| B-factors of two structures compared at face value | B conflates thermal motion, disorder, and model error and is not portable | Compare within one structure only; normalize (z/percentile) before any cross-structure claim |
| Coloring an AlphaFold model "by B-factor" to infer flexibility | pLDDT rides in the B-factor column with OPPOSITE polarity (high = confident) | Read the column as pLDDT bands; high value means high confidence, not high motion |
| Deleting all low-pLDDT residues as junk | Low pLDDT often marks a real intrinsically disordered region | Distinguish disorder (keep, annotate) from misfold; trim only for MR/geometry pipelines |
| Docking straight into an AlphaFold pocket | Pocket is apo, side-chain rotamers are the least reliable atoms, may be wrong state | Prefer an experimental holo structure; if using the model, ensemble/flexible-side-chain dock and treat hits as hypotheses |
| Confident predicted model trusted for a complex | pLDDT/PAE bound self-consistency, not biological correctness | Ask what context AF could not see (partner, ligand, PTM); a monomer/apo model can be confidently wrong |
| Cryo-EM peripheral domain trusted at the headline resolution | Global resolution hides a low-local-resolution flexible arm | Consult the EMDB local-resolution map; treat low-local-res regions as docked hypotheses |
| Quoting FSC 0.5 as the cryo-EM resolution | 0.143 is the half-map criterion; 0.5 is the map-vs-model curve | Use gold-standard half-map FSC at 0.143 for resolution; 0.5 for checking the model against the map |
| Averaging NMR model coordinates into one structure | The mean of two valid conformers is physically impossible (distorted bonds, clashes) | Compute on each model and report the spread, or pick a representative/medoid model |
| cis peptide flagged everywhere, or missed entirely | omega near 0 (cis) vs 180 (trans) not checked; cis-Pro is common but cis non-Pro is rare | Compute omega; treat cis non-proline as a validation flag pending density, cis-Pro as plausible |
| Python Ramachandran percentages reported as authoritative | Coarse basin boxes are not MolProbity's rama8000 contours | Use the screen to decide whether to run `phenix.molprobity`; report MolProbity's numbers |
| `resolution`/`R-free` come back `None` from `structure.header` | Bio.PDB's header dict is thin | Read `_refine.ls_d_res_high`, `_refine.ls_R_factor_R_free`, `_exptl.method` via `MMCIF2Dict` |

## Related Skills

- structure-io - Read resolution/R-free via MMCIF2Dict and fetch the biological assembly the validation applies to
- structure-navigation - Resolve altlocs, insertion codes, and multi-model NMR files before validating per-model
- geometric-analysis - Compute the dihedrals and DSSP secondary structure this skill validates; measure only after validation passes
- structure-modification - Trim low-pLDDT residues or edit B-factors once a predicted model is validated
- structure-preparation - Add hydrogens, protonation, and missing atoms after validation and before docking/MD
- alphafold-predictions - Download the AlphaFold model plus its PAE JSON that this skill reads for confidence
- modern-structure-prediction - Reconcile a re-run prediction with pLDDT/PAE/pTM when the AFDB entry is untrustworthy
- interface-analysis - Validate the assembly before interpreting an interface that only exists in it
- database-access/uniprot-access - Map validated residues back to a UniProt reference sequence

## References

- Ramachandran GN, Ramakrishnan C, Sasisekharan V (1963). Stereochemistry of polypeptide chain configurations. *J Mol Biol* 7:95-99. DOI 10.1016/S0022-2836(63)80023-6.
- Brunger AT (1992). Free R value: a novel statistical quantity for assessing the accuracy of crystal structures. *Nature* 355(6359):472-475. DOI 10.1038/355472a0.
- Kabsch W, Sander C (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. *Biopolymers* 22(12):2577-2637. DOI 10.1002/bip.360221211.
- Kleywegt GJ, Brunger AT (1996). Checking your imagination: applications of the free R value. *Structure* 4(8):897-904. DOI 10.1016/S0969-2126(96)00097-4.
- Rosenthal PB, Henderson R (2003). Optimal determination of particle orientation, absolute hand, and contrast loss in single-particle electron cryomicroscopy. *J Mol Biol* 333(4):721-745. DOI 10.1016/j.jmb.2003.07.013.
- Chen VB, Arendall WB III, Headd JJ, Keedy DA, Immormino RM, Kapral GJ, Murray LW, Richardson JS, Richardson DC (2010). MolProbity: all-atom structure validation for macromolecular crystallography. *Acta Crystallogr D* 66(1):12-21. DOI 10.1107/S0907444909042073.
- Read RJ, Adams PD, Arendall WB III, Brunger AT, Emsley P, Joosten RP, Kleywegt GJ, Krissinel EB, Luetteke T, Otwinowski Z, Perrakis A, Richardson JS, Sheffler WH, Smith JL, Tickle IJ, Vriend G, Zwart PH (2011). A new generation of crystallographic validation tools for the Protein Data Bank. *Structure* 19(10):1395-1412. DOI 10.1016/j.str.2011.08.006.
- Gore S, Sanz Garcia E, Hendrickx PMS, Gutmanas A, Westbrook JD, Yang H, Feng Z, Baskaran K, Berrisford JM, et al. (2017). Validation of structures in the Protein Data Bank. *Structure* 25(12):1916-1927. DOI 10.1016/j.str.2017.10.009.
- Williams CJ, Headd JJ, Moriarty NW, Prisant MG, Videau LL, Deis LN, Verma V, Keedy DA, Hintze BJ, Chen VB, Jain S, Lewis SM, Arendall WB III, Snoeyink J, Adams PD, Lovell SC, Richardson JS, Richardson DC (2018). MolProbity: more and better reference data for improved all-atom structure validation. *Protein Sci* 27(1):293-315. DOI 10.1002/pro.3330.
- Jumper J, Evans R, Pritzel A, et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596(7873):583-589. DOI 10.1038/s41586-021-03819-2.
- Oeffner RD, Croll TI, Millan C, Poon BK, Schlicksup CJ, Read RJ, Terwilliger TC (2022). Putting AlphaFold models to work with phenix.process_predicted_model and ISOLDE. *Acta Crystallogr D* 78(11):1303-1314. DOI 10.1107/S2059798322010026.
- Montelione GT, Nilges M, Bax A, et al. (2013). Recommendations of the wwPDB NMR Validation Task Force. *Structure* 21(9):1563-1570. DOI 10.1016/j.str.2013.07.021.
