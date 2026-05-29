---
name: bio-admet-prediction
description: Predicts ADMET properties using ADMETlab 3.0 (119 endpoints with uncertainty), ADMET-AI, DeepChem MolNet, and chemprop D-MPNN with explicit handling of OECD QSAR principles, applicability domain assessment, calibration, hERG/CYP/AMES gold-standard endpoints, and PAINS / Lipinski / Ro5 / Veber / BBB druglikeness filters. Use when filtering compounds for drug-likeness, prioritizing leads by predicted safety, or building an in-house ADMET QSAR model.
tool_type: python
primary_tool: ADMETlab
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, requests 2.31+, DeepChem 2.8+, chemprop 2.0+ (note major API change from 1.x), admet-ai 1.3+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# ADMET Prediction

Predict absorption, distribution, metabolism, excretion, and toxicity properties of drug candidates. ADMET prediction underpins lead selection and de-risking; calibrated, applicability-domain-aware predictions distinguish a working filter from a costly false-confidence rejection. Modern best practice combines online services (ADMETlab 3.0 with uncertainty estimates), open-source models (chemprop D-MPNN), and rule-based filters (Lipinski / Veber / BBB heuristics) -- each with known failure modes.

For PAINS / Brenk / structural alerts, see `chemoinformatics/substructure-search`. For QSAR model building from in-house data, see `chemoinformatics/qsar-modeling`.

## ADMET Model Taxonomy

| Tool | Endpoints | Architecture | Uncertainty | Access | Fails when |
|------|-----------|--------------|-------------|--------|------------|
| ADMETlab 3.0 | 119 (A,D,M,E,T + physchem + medchem) | Multi-task DMPNN + descriptors | Per-prediction | REST API (free, no auth) | Outside training distribution; metals; macrocycles |
| ADMET-AI (NVIDIA) | ~50 (focus on safety) | chemprop D-MPNN | Ensemble variance | Python package | Limited endpoints vs ADMETlab |
| DeepChem MolNet | ~30 (tox21, ToxCast, ClinTox) | Various GCN/GAT | Per-task variance | Python package | Models trained on small datasets |
| pkCSM | ~30 (absorption, distribution, toxicity) | Graph signatures + RF | None | Web service | Smaller training data |
| SwissADME | ~30 (filters + physchem) | Hand-curated rules | None | Web service (NO API) | Cannot batch programmatically |
| ProTox-3.0 | ~46 (toxicity end-points) | DT + descriptors | None | Web service | Toxicity only; LD50 categorical |
| ADMETpredictor (Simulations Plus) | ~140 | Proprietary | Per-prediction | Commercial | License cost |
| FAF-Drugs4 | filters | Rule-based | None | Web | Static rules |
| chemprop (in-house) | User-defined | D-MPNN ± descriptors | Bayesian ensemble | Python package | Requires training data |

**Decision:** For batch screening of <10k compounds with no in-house data, **ADMETlab 3.0** (free API, 119 endpoints, calibrated uncertainty) is the modern standard. For in-house QSAR on a specific endpoint with >500 measurements, train a **chemprop D-MPNN + Morgan + MOE descriptors** model (Liu et al. 2024 achieved AUC 0.956 on hERG with this combo).

## Decision Tree by Scenario

| Scenario | Workflow | Reasoning |
|----------|----------|-----------|
| Library triage, no in-house data | ADMETlab 3.0 API batch | Calibrated 119 endpoints |
| Single target, large in-house data (>500 datapoints) | chemprop D-MPNN ensemble | Beats generic models on target-specific data |
| Need calibrated probabilities | chemprop with ensemble + Platt | Native deep learning rarely calibrated |
| FDA / regulatory submission | OECD-compliant QSAR with AD | See OECD principles below |
| Quick filter for VS | Lipinski + Veber + QED >= 0.5 | Rule-based, no model needed |
| BBB penetration | TPSA <= 90, MW <= 500, HBD <= 3 (Pfizer CNS) | Wager 2010 |
| Cardiotox liability | hERG model (ADMETlab + ProTox + literature lit-check) | Triangulate; hERG critical |
| Drug-drug interaction (CYP) | CYP1A2/2C9/2C19/2D6/3A4 inhibitor + substrate | Standard set of 5 CYPs |

## OECD QSAR Principles (5 Pillars)

For regulatory-grade ADMET QSAR (REACH, ECHA, FDA submissions), models must satisfy:

1. **Defined endpoint** -- specific bioassay, units, conditions
2. **Unambiguous algorithm** -- reproducible model + code
3. **Defined applicability domain (AD)** -- where the model is valid
4. **Appropriate statistical validation** -- external test set, cross-validation
5. **Mechanistic interpretation** -- biological / chemical rationale

For non-regulatory work, AD assessment is still critical. The OECD's *applicability domain* is the workhorse: predictions outside the AD are unreliable, but operational AD measures (leverage, kNN, conformal prediction) often disagree.

## Applicability Domain Methods

| Method | Definition | Flags out-of-AD when |
|--------|-----------|----------------------|
| kNN distance | Mean distance to k nearest neighbors in training set | > training-set distribution P95 |
| Leverage (Williams) | Hat-matrix diagonal | > 3p/n (p = features, n = compounds) |
| Density (KDE on PCA) | Density in feature space | < density of training set P5 |
| Conformal prediction | Per-prediction confidence interval | Interval > tolerance |
| Bayesian variance | Ensemble or MC-dropout variance | > training-set variance P95 |

For deep-learning ADMET, **conformal prediction** (Bostrom et al. 2024) is becoming the standard.

## ADMETlab 3.0 API

The current standard for free ADMET prediction. 119 endpoints across 6 categories; per-prediction uncertainty.

**Goal:** Predict 119 ADMET endpoints with uncertainty for a batch of SMILES using a hosted API.

**Approach:** POST batches of <=500 SMILES to ADMETlab 3.0 REST endpoint and parse the returned JSON into a per-compound endpoint DataFrame.

```python
import requests
import pandas as pd

# ADMETlab 3.0 API base: https://admetlab3.scbdd.com
# Endpoints: /api/admet (full 119-endpoint batch), /api/wash (standardization),
# /api/single/admet (single SMILES), /api/render (visualization)
# See https://admetlab3.scbdd.com/apis/ for current API spec.
def admetlab_predict(smiles_list, endpoint='admet'):
    url = f'https://admetlab3.scbdd.com/api/{endpoint}'
    payload = {'smiles': smiles_list}
    response = requests.post(url, json=payload, timeout=120)
    response.raise_for_status()
    return pd.DataFrame(response.json())

smiles = ['CCO', 'c1ccc(C(=O)O)cc1', 'CC(=O)Oc1ccccc1C(=O)O']
results = admetlab_predict(smiles)  # POST batches of <=500 SMILES at a time
```

ADMETlab 3.0 endpoints (sample):
- Absorption: Caco-2 permeability (logPapp), HIA (%), Pgp inhibitor/substrate, MDCK
- Distribution: BBB+, PPB (%), VDss (L/kg), Fu (fraction unbound)
- Metabolism: CYP1A2/2C9/2C19/2D6/3A4 inhibitor / substrate
- Excretion: CL (mL/min/kg), T1/2 (h)
- Toxicity: hERG, AMES, hepatotoxicity (DILI), carcinogenicity, immunotoxicity, mutagenicity, respiratory, skin, eye, cardiotoxicity, mitochondrial, NR-AR, NR-ER, SR-MMP
- Drug-likeness: Lipinski, Veber, Ghose, Egan, Muegge, QED, SAscore

## chemprop D-MPNN for Custom Endpoints

When in-house data is available, train a target-specific model. chemprop's D-MPNN architecture + atom/bond features + optional Morgan / MOE descriptors is the modern open-source SOTA.

**Goal:** Train a target-specific ADMET classifier or regressor on in-house bioassay data.

**Approach:** Use chemprop 2.x CLI with rdkit_2d_normalized descriptor features, 5-fold scaffold-balanced split, and 5-model ensemble for uncertainty estimation.

```python
# chemprop 2.x CLI (current; 'chemprop train' with space + dashed args)
# chemprop train --data-path data.csv --task-type classification \
#                --save-dir model_dir --molecule-featurizers rdkit_2d_normalized \
#                --num-folds 5 --ensemble-size 5

# Or chemprop 2.x programmatic API (full programmatic API documented at chemprop.readthedocs.io)
# See chemoinformatics/qsar-modeling for the full chemprop 2.x training pipeline.
```

**Key:** 5-fold ensemble for calibrated uncertainty; D-MPNN + rdkit_2d_normalized typically outperforms either alone.

## hERG Cardiotoxicity (Gold Standard Endpoint)

hERG (KCNH2) blockade causes QT prolongation, Torsades de Pointes, and is the #1 reason for late-stage drug attrition. ICH S7B and FDA require non-clinical assessment.

| Model | Architecture | Training data | AUC | Reference |
|-------|--------------|---------------|-----|-----------|
| Cai et al. D-MPNN + MOE | chemprop + 206 MOE descriptors | 7,889 compounds | 0.956 | Liu 2024 |
| CardioTox-net | ECFP + RF stacking | ChEMBL hERG | 0.93 | Aniketh 2021 |
| ADMETlab 3.0 hERG | DMPNN multi-task | Internal | 0.92 (reported) | Fu 2024 |
| ProTox-3.0 | DT | ProTox training | 0.86 | Banerjee 2024 |

**Postdoc-grade interpretation:** A single-model probability > 0.5 is NOT a kill signal. Triangulate ADMETlab + ProTox + literature search for cardiotox liability. Consider the **active concentration** vs predicted hERG IC50 (10 uM threshold for "concern"); a hERG-positive compound at 100 nM active concentration may still be safe with sufficient hERG selectivity.

## CYP Inhibition (DDI Risk)

5 CYP isoforms cover most clinically relevant DDIs:

| CYP | Substrates (drugs) | Inhibitor flag if predicted prob | Action |
|-----|--------------------|-----------------------------------|--------|
| CYP3A4 | ~50% of drugs | > 0.5 | Often acceptable if substrate; flag if inhibitor |
| CYP2D6 | beta-blockers, antidepressants | > 0.5 | Higher concern (polymorphism) |
| CYP2C9 | warfarin, NSAIDs | > 0.5 | Caution if patient on warfarin |
| CYP2C19 | PPIs, clopidogrel | > 0.5 | Variable (polymorphism) |
| CYP1A2 | caffeine, theophylline | > 0.5 | Diet/smoking-influenced |

## PAINS, BRENK, REOS Filters

ADMET prediction is separate from structural alerts; combine. See `chemoinformatics/substructure-search` for PAINS/BRENK/REOS pattern catalogs.

```python
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

def alerts(mol, catalogs=('PAINS_A', 'BRENK', 'ZINC')):
    params = FilterCatalogParams()
    for cat in catalogs:
        params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, cat))
    catalog = FilterCatalog(params)
    hits = catalog.GetMatches(mol)
    return [h.GetDescription() for h in hits]
```

## Lipinski / Veber / Drug-Likeness

See `chemoinformatics/molecular-descriptors` for full physchem table. Quick filter:

```python
from rdkit.Chem import Descriptors, Lipinski, QED

def druglike_score(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotbonds = Lipinski.NumRotatableBonds(mol)
    qed = QED.qed(mol)

    lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_pass = rotbonds <= 10 and tpsa <= 140
    bbb_pfizer = tpsa <= 90 and mw <= 500 and hbd <= 3

    return {'MW': mw, 'LogP': logp, 'HBD': hbd, 'HBA': hba,
            'TPSA': tpsa, 'RotBonds': rotbonds, 'QED': qed,
            'Lipinski_violations': lipinski_violations,
            'Veber_pass': veber_pass, 'BBB_pfizer': bbb_pfizer}
```

## Per-Tool Failure Modes

### ADMETlab 3.0 -- out-of-distribution prediction

**Trigger:** Macrocycle, metal-coordinated, oligomer (peptide >5 AA), or PROTAC submitted.

**Mechanism:** ADMETlab training set is drug-like organic molecules. Predictions on PROTACs, macrocycles, peptides extrapolate.

**Symptom:** Uncertainty estimates near max; ambiguous probabilities (0.4-0.6).

**Fix:** Check uncertainty band; if interval is broad, do not trust point estimate. For PROTACs / macrocycles, prefer literature-derived experimental data.

### hERG D-MPNN -- training data bias

**Trigger:** Compound is novel chemotype not in training set (drug-like but in unexplored region).

**Mechanism:** D-MPNN learns local chemical features; for genuinely new scaffolds, extrapolation is unreliable.

**Symptom:** Model predicts hERG- (false negative) for compound that experimentally inhibits.

**Fix:** Use ensemble + applicability-domain assessment (kNN distance, ensemble variance). If kNN distance to training set > P95, treat prediction as low-confidence.

### CYP3A4 inhibitor + substrate ambiguity

**Trigger:** Model trained on either inhibitor OR substrate; predictions confused.

**Mechanism:** CYP3A4 inhibitors and substrates have similar SAR; many compounds are both.

**Symptom:** Both classes report > 0.5.

**Fix:** Two separate models (inhibitor model, substrate model); compounds that score high in both are flagged for in vitro confirmation.

### SwissADME -- no API

**Trigger:** Wanting to batch programmatically.

**Mechanism:** SwissADME explicitly forbids automated access (TOS).

**Symptom:** No programmatic access; manual web upload only.

**Fix:** Use ADMETlab 3.0 API for programmatic batch (free, no TOS restrictions).

### PAINS as a kill filter

**Trigger:** Treating PAINS_A match as a categorical exclusion.

**Mechanism:** PAINS is calibrated against HTS assay-interference; matches do NOT predict failed drug development.

**Symptom:** Library purged of valid leads (curcumin analogs, polyphenol natural products).

**Fix:** Flag PAINS for orthogonal-assay confirmation; do not exclude pre-emptively. See substructure-search for details.

### Class-imbalanced AMES dataset

**Trigger:** Training/predicting AMES mutagenicity.

**Mechanism:** AMES public datasets have ~3:1 negative:positive imbalance; baseline F1 misleading.

**Symptom:** Model reports high accuracy but predicts negative for all.

**Fix:** Use AUC-ROC or balanced accuracy. Apply SMOTE / class-weighted loss in chemprop.

## Reconciliation Across Models

When ADMETlab, ProTox, and chemprop disagree on hERG:
- All three predict hERG+ -> high confidence, deprioritize
- Two predict hERG+, one hERG- -> moderate confidence, plan in vitro patch-clamp
- Single positive (other two negative) -> likely false positive of disagreeing model; verify chemotype is in distribution
- All three hERG- -> low confidence (especially novel chemotype); still consider in vitro screen for clinical candidates

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ADMETlab API 504 timeout | Batch too large | Submit <500 SMILES per request |
| chemprop training overfits | Random split | Use scaffold split (`--split scaffold_balanced`) |
| hERG prediction 50/50 | Out-of-distribution | Check applicability domain |
| QED returns nan | Stereo undefined or invalid | Standardize first |
| ProTox endpoints missing | Web scrape uses CSS selector | Use formal API |
| BBB+ true but TPSA > 90 | Different BBB model | Use Pfizer CNS heuristic for orthogonal check |
| Predictions inconsistent across runs | Random seed for chemprop ensemble | `--seed 42` and reuse model |
| Calibration mismatch | DL native probabilities not calibrated | Apply Platt scaling on validation set |

## References

- Fu et al., *Nucleic Acids Res.* 52:W422 (2024) -- ADMETlab 3.0.
- Shan W, Jiang Y, Chen J, Qin H, Qin H, Cheng F. 2022. RSC Adv 12(7):3935-3942 -- hERG ML benchmarks (DOI 10.1039/D1RA07956E).
- Heid et al., *J. Chem. Inf. Model.* 64:9 (2024) -- chemprop 2.0 redesign.
- OECD, "OECD Principles for the Validation of QSAR Models" (2007).
- OECD, "(Q)SAR Assessment Framework" (2023).
- Capuzzi et al., *J. Chem. Inf. Model.* 57:417 (2017) -- PAINS reality check.
- Wager et al., *ACS Chem. Neurosci.* 1:435 (2010) -- Pfizer CNS MPO.
- Bickerton et al., *Nat. Chem.* 4:90 (2012) -- QED.

## Related Skills

- chemoinformatics/molecular-descriptors - Compute drug-likeness physchem
- chemoinformatics/substructure-search - PAINS / BRENK / REOS filter
- chemoinformatics/qsar-modeling - Build custom QSAR for in-house data
- chemoinformatics/molecular-standardization - Canonicalize before prediction
- machine-learning/biomarker-discovery - Adjacent ML approaches
- clinical-databases/pharmacogenomics - Patient genotype overlay
