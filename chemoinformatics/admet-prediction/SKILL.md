---
name: bio-admet-prediction
description: Predicts ADMET properties using ADMETlab 3.0 (119 platform features, including 77 prediction models with modeled-endpoint uncertainty), ADMET-AI, DeepChem MolNet, and chemprop D-MPNN with explicit handling of OECD QSAR principles, applicability domain assessment, calibration, hERG/CYP/AMES endpoints, and PAINS / Lipinski / Ro5 / Veber / BBB druglikeness filters. Use when filtering compounds for drug-likeness, prioritizing leads by predicted safety, or building an in-house ADMET QSAR model.
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
| ADMETlab 3.0 | 119 reported features: 77 prediction models, 34 computed properties, 8 rules | Multi-task DMPNN + descriptors for modeled endpoints | Evidential uncertainty for modeled endpoints | Web service; hosted API documented by the authors | Outside training distribution; metals; macrocycles |
| ADMET-AI | TDC-derived ADMET tasks; inspect installed model metadata | Chemprop D-MPNN | Inspect version-specific outputs; do not assume calibrated uncertainty | Python package | v2 package predictions differ from the v1 paper/server |
| DeepChem MolNet | Dataset-dependent tasks including Tox21, ToxCast, and ClinTox | Model-dependent | Model-dependent | Python package | Coverage and uncertainty depend on the selected dataset/model |
| pkCSM | Service-defined ADMET endpoints | Graph signatures + ML | Inspect current service output | Web service | Applicability domain and service contract must be checked |
| SwissADME | Physchem, pharmacokinetics, drug-likeness, and medchem outputs | Published models and rules | None advertised | Web service (no public API) | Automated access is restricted by its terms |
| ProTox-3.0 | 61 toxicity models/endpoints | RF/DNN + fingerprints, similarity, and pharmacophore methods | Confidence score | Web service / sample API | Toxicity only; reports LD50 and toxicity class |
| ADMETpredictor (Simulations Plus) | ~140 | Proprietary | Per-prediction | Commercial | License cost |
| FAF-Drugs4 | filters | Rule-based | None | Web | Static rules |
| chemprop (in-house) | User-defined | D-MPNN ± descriptors | Ensemble and other estimators; optional calibration | Python package | Requires suitable training and calibration data |

**Decision:** For batch screening with no in-house data, **ADMETlab 3.0** provides 119 reported platform features and uncertainty for modeled endpoints through its web service and hosted API; verify the live API documentation before automating access. For a sufficiently large, relevant in-house endpoint dataset, benchmark a **chemprop D-MPNN**, descriptors, and simpler baselines under a deployment-relevant split rather than assuming a universal sample-size threshold. Shan et al. (2022) reported an AUC of 0.956 for a D-MPNN combined with 206 MOE descriptors on their random-split hERG benchmark.

## Decision Tree by Scenario

| Scenario | Workflow | Reasoning |
|----------|----------|-----------|
| Library triage, no in-house data | ADMETlab 3.0 API batch | Broad platform coverage plus modeled-endpoint uncertainty |
| Single endpoint, adequate in-house data | Benchmark chemprop D-MPNN, descriptors, and simpler baselines | Select by prospective or deployment-relevant validation |
| Need calibrated probabilities | chemprop with ensemble + Platt | Native deep learning rarely calibrated |
| FDA / regulatory submission | OECD-compliant QSAR with AD | See OECD principles below |
| Quick annotation for VS | Lipinski, Veber, and QED reported separately | Rank or annotate; do not impose a universal QED gate |
| BBB penetration | Simple screen: TPSA <= 90, MW <= 500, HBD <= 3 | Repository heuristic; not the six-factor CNS MPO |
| Cardiotox liability | ADMETlab hERG + ProTox-3.0 cardiotoxicity + literature check | Both hosted endpoints model hERG blockade; compare applicability domains and assay definitions |
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

For deep-learning ADMET, **conformal prediction** can provide calibrated prediction sets or intervals when its exchangeability and calibration assumptions are appropriate (McShane et al. 2024).

## ADMETlab 3.0 API

ADMETlab 3.0 reports 119 platform features: 77 prediction models, 34 computed physicochemical properties, and 8 medicinal-chemistry rules. The modeled endpoints include prediction uncertainty; do not imply that computed properties and rules have model uncertainty.

**Goal:** Obtain the platform's 119 features for a batch of SMILES, including uncertainty for the 77 modeled endpoints, using the hosted service.

**Approach:** Follow the live ADMETlab 3.0 API tutorial to wash molecules, submit batch predictions, and retrieve the returned results. The 2024 paper documents API/batch support and modeled-endpoint uncertainty. Obtain current rate limits, routes, payloads, task identifiers, and output contracts from the live official documentation rather than attributing them to the paper or hard-coding an unofficial example.

```python
import pandas as pd

# After submitting with the current official API example, load its CSV output.
results = pd.read_csv('admetlab3_results.csv')
# Preserve the uncertainty columns and task identifier in downstream reports.
```

ADMETlab 3.0 endpoints (sample):
- Absorption: Caco-2 permeability (logPapp), HIA (%), Pgp inhibitor/substrate, MDCK
- Distribution: BBB+, PPB (%), VDss (L/kg), Fu (fraction unbound)
- Metabolism: CYP1A2/2C9/2C19/2D6/3A4 inhibitor / substrate
- Excretion: CL (mL/min/kg), T1/2 (h)
- Toxicity: hERG, AMES, hepatotoxicity (DILI), carcinogenicity, immunotoxicity, mutagenicity, respiratory, skin, eye, cardiotoxicity, mitochondrial, NR-AR, NR-ER, SR-MMP
- Drug-likeness: Lipinski, Veber, Ghose, Egan, Muegge, QED, SAscore

## chemprop D-MPNN for Custom Endpoints

When in-house data is available, train a target-specific model. chemprop provides a widely used open-source D-MPNN architecture with atom/bond features and optional molecular descriptors; benchmark it against appropriate baselines on the project's data.

**Goal:** Train a target-specific ADMET classifier or regressor on in-house bioassay data.

**Approach:** Use the installed Chemprop 2.x CLI with a scaffold split, a release-supported descriptor featurizer, replicated models, and an explicit prediction-time uncertainty/calibration workflow.

```python
# Chemprop 2.2 CLI; verify flags against the installed release.
# chemprop train --data-path data.csv --task-type classification \
#                --save-dir model_dir --split-type scaffold_balanced \
#                --molecule-featurizers v1_rdkit_2d_normalized \
#                --num-replicates 5 --ensemble-size 5
# chemprop predict --test-path test.csv --model-paths model_dir \
#                  --uncertainty-method ensemble --preds-path predictions.csv

# Or chemprop 2.x programmatic API (full programmatic API documented at chemprop.readthedocs.io)
# See chemoinformatics/qsar-modeling for the full chemprop 2.x training pipeline.
```

**Key:** Replicates and an ensemble estimator produce an uncertainty estimate, not automatic calibration. Fit and evaluate a documented calibrator on a separate calibration set when calibrated probabilities or intervals are required. Descriptor benefit must be demonstrated on the intended endpoint.

## hERG Cardiotoxicity Endpoint

hERG (KCNH2) blockade can contribute to QT prolongation and Torsades de Pointes and is an important non-clinical cardiac-safety endpoint. Follow the current ICH S7B/E14 and regulator-specific guidance applicable to the program rather than treating one model output as a regulatory conclusion.

| Model | Architecture | Training data | AUC | Reference |
|-------|--------------|---------------|-----|-----------|
| Shan et al. D-MPNN + MOE | D-MPNN + 206 MOE descriptors | 7,889 compounds | 0.956 (random split) | Shan 2022 |
| CardioTox-net | Five DL base representations + neural meta-ensemble | BindingDB, ChEMBL, and literature | 0.930 (10-fold meta-validation) | Karim 2021 |
| ADMETlab 3.0 hERG | DMPNN multi-task | Internal | 0.92 (reported) | Fu 2024 |
| ProTox-3.0 cardiotoxicity | RF-based classifier | 5,252 ChEMBL compounds with hERG IC50/Ki | 0.86 CV; 0.95 external | Banerjee 2024 |

**Interpretation:** A single-model probability > 0.5 is NOT a kill signal. Triangulate multiple hERG-specific models and a literature search. ProTox-3.0 calls the endpoint cardiotoxicity, but its model specifically predicts small-molecule hERG blockers; it should not be treated as an independent non-hERG mechanism. Consider exposure relative to measured hERG potency and confirm important decisions experimentally rather than applying a universal safe/unsafe IC50 cutoff.

## CYP Inhibition (DDI Risk)

5 CYP isoforms cover most clinically relevant DDIs:

| CYP | Substrates (drugs) | Inhibitor flag if predicted prob | Action |
|-----|--------------------|-----------------------------------|--------|
| CYP3A4 | many drug classes | Model-specific threshold | Interpret inhibitor and substrate assays separately |
| CYP2D6 | beta-blockers, antidepressants | Model-specific threshold | Include polymorphism and exposure context |
| CYP2C9 | warfarin, NSAIDs | Model-specific threshold | Evaluate clinical substrate/exposure context |
| CYP2C19 | PPIs, clopidogrel | Model-specific threshold | Include polymorphism and assay context |
| CYP1A2 | caffeine, theophylline | Model-specific threshold | Include induction, diet, and smoking context |

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
    bbb_simple_screen = tpsa <= 90 and mw <= 500 and hbd <= 3

    return {'MW': mw, 'LogP': logp, 'HBD': hbd, 'HBA': hba,
            'TPSA': tpsa, 'RotBonds': rotbonds, 'QED': qed,
            'Lipinski_violations': lipinski_violations,
            'Veber_pass': veber_pass, 'BBB_simple_screen': bbb_simple_screen}
```

## Per-Tool Failure Modes

### ADMETlab 3.0 -- out-of-distribution prediction

**Trigger:** Chemistry materially unlike the service's documented training/applicability domain, such as metal-containing complexes, many peptides, PROTACs, or unusual macrocycles.

**Mechanism:** ADMETlab training set is drug-like organic molecules. Predictions on PROTACs, macrocycles, peptides extrapolate.

**Symptom:** High reported uncertainty, disagreement with neighbors or orthogonal models, or unstable conclusions under reasonable preprocessing.

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

**Mechanism:** SwissADME's terms restrict automated crawler/data-retrieval access, and no public API is documented.

**Symptom:** No programmatic access; manual web upload only.

**Fix:** Use a currently documented programmatic service and follow its access policy; for ADMETlab 3.0, verify the live API tutorial before writing a client.

### PAINS as a kill filter

**Trigger:** Treating PAINS_A match as a categorical exclusion.

**Mechanism:** PAINS is calibrated against HTS assay-interference; matches do NOT predict failed drug development.

**Symptom:** Library purged of valid leads (curcumin analogs, polyphenol natural products).

**Fix:** Flag PAINS for orthogonal-assay confirmation; do not exclude pre-emptively. See substructure-search for details.

### Class-imbalanced AMES dataset

**Trigger:** Training/predicting AMES mutagenicity.

**Mechanism:** Public AMES datasets can be imbalanced and differ in assay definition and curation; aggregate accuracy can therefore be misleading.

**Symptom:** Model reports high accuracy but predicts negative for all.

**Fix:** Report class balance and use suitable metrics such as PR-AUC, ROC-AUC, MCC, or balanced accuracy. Compare class weighting or resampling inside training folds without leaking validation/test data.

## Reconciliation Across Models

When ADMETlab, ProTox-3.0, and an independently trained chemprop model disagree on hERG:
- All predict hERG+ -> higher concern; plan in vitro patch-clamp
- Results disagree -> inspect applicability domains, activity thresholds, and assay definitions before deciding
- All predict hERG- -> lower concern, but still consider in vitro screening for clinical candidates and novel chemotypes
- Do not count correlated models as independent evidence merely because they are hosted by different services

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ADMETlab API timeout | Service load, payload, or current quota | Follow live batch limits; retry with backoff and record failures |
| chemprop training overfits | Random split | Use scaffold split (`--split scaffold_balanced`) |
| hERG prediction 50/50 | Out-of-distribution | Check applicability domain |
| QED calculation fails | Molecule is missing, unsanitized, or unsupported | Reject parse failures; sanitize inputs and handle calculation exceptions |
| ProTox endpoints missing | Web scrape uses CSS selector | Use formal API |
| BBB+ true but TPSA > 90 | Different BBB model | Use the simple physicochemical screen as an orthogonal heuristic |
| Predictions inconsistent across runs | Random seed for chemprop ensemble | `--seed 42` and reuse model |
| Calibration mismatch | DL native probabilities not calibrated | Apply Platt scaling on validation set |

## References

- Fu et al., *Nucleic Acids Res.* 52:W422-W431 (2024) -- ADMETlab 3.0 (DOI 10.1093/nar/gkae236).
- Shan M, Jiang C, Chen J, Qin L-P, Qin J-J, Cheng G. *RSC Adv.* 12:3423-3430 (2022) -- D-MPNN/MOE hERG benchmark (DOI 10.1039/D1RA07956E).
- Karim A, Lee M, Balle T, Sattar A. *J. Cheminformatics* 13:60 (2021) -- CardioTox-net (DOI 10.1186/s13321-021-00541-z).
- Banerjee P, Kemmler E, Dunkel M, Preissner R. *Nucleic Acids Res.* 52:W513-W520 (2024) -- ProTox-3.0 (DOI 10.1093/nar/gkae303).
- McShane SA et al. *J. Cheminformatics* 16:75 (2024) -- conformal prediction for molecular-property models (DOI 10.1186/s13321-024-00870-9).
- Heid E et al., *J. Chem. Inf. Model.* 64:9-17 (2024) -- Chemprop redesign (DOI 10.1021/acs.jcim.3c01250).
- Chemprop documentation, training, descriptors, and uncertainty: https://chemprop.readthedocs.io/
- ADMET-AI official repository and version notes: https://github.com/swansonk14/admet_ai
- SwissADME Terms of Use: https://www.swissadme.ch/termsofuse.php
- OECD, "Principles for the Validation, for Regulatory Purposes, of (Q)SAR Models" (agreed 2004).
- OECD, *Guidance Document on the Validation of (Q)SAR Models*, OECD Series on Testing and Assessment No. 69 (2007).
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
