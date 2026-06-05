---
name: bio-immunoinformatics-mhc-binding-prediction
description: Predict peptide-MHC class I binding and natural presentation with MHCflurry, NetMHCpan-4.1, and MixMHCpred to nominate candidate CD8 T-cell epitopes. Covers the binding-affinity (BA) vs eluted-ligand (EL/presentation) distinction, why %Rank beats raw nM for cross-allele work, the MS abundance bias that misranks low-expression neoantigens, allele-coverage inequity, and length bias. Use when scanning a protein or peptide set for class I epitopes, scoring neoantigen candidates, or choosing a binding predictor. For CD4/HLA class II see mhc-class-ii-prediction.
tool_type: python
primary_tool: mhcflurry
---

## Version Compatibility

Reference examples tested with: MHCflurry 2.1+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: MHCflurry 2.2.0+ switched its backend from TensorFlow to PyTorch (Python 3.10+) — confirm `mhcflurry-downloads fetch` succeeded and the backend imports before scoring. NetMHCpan-4.1 and MixMHCpred are standalone academic binaries, not pip-installable; the IEDB REST API wraps NetMHCpan if a local install is unavailable. Tool versions move fast — re-verify supported-allele lists and default %Rank thresholds against current docs.

# MHC Binding Prediction

**"Predict which peptides bind/are presented by MHC class I"** -> Score peptide-HLA class I binding affinity and natural-presentation likelihood to nominate candidate CD8 epitopes.
- Python: `mhcflurry.Class1PresentationPredictor.load().predict()` (pip-installable, forgiving allele parser)
- CLI: `netMHCpan` (field default; EL score by default, `-BA` adds affinity) or `MixMHCpred` (MS-deconvolution, EL-only)

## The Single Most Important Modern Insight -- a strong predicted binder is a candidate for the next experiment, not an epitope

Binding to MHC is necessary but nowhere near sufficient for immunogenicity. The real path is a funnel: expression -> proteasomal processing -> TAP transport and loading -> stable surface display -> a cognate T cell that survived thymic selection and activates. Binding prediction addresses essentially one stage. Each downstream stage discards a large fraction of binders, so the precision of "predicted binder -> validated epitope" is low even when the binding model itself is excellent. Two operational corollaries follow. First, never report a presentation score to a collaborator as an "immunogenicity" or "epitope" probability — that is a different, far weaker prediction (immunoinformatics/immunogenicity-scoring). Second, the modern EL/MS models that now define the field learned natural presentation from mass-spec immunopeptidomes, which over-represent peptides from highly expressed proteins; the model therefore partly learns "comes from an abundant protein" as a proxy for "is presented." That bias is exactly backwards for neoantigen discovery, where the targets are mutated and often lowly expressed, living in the under-detected tail the model systematically under-ranks.

## Tool Taxonomy (Class I)

| Tool | Citation | Score type | Form | Use when |
|------|----------|-----------|------|----------|
| NetMHCpan-4.1 | Reynisson 2020 | EL (default) + BA (`-BA`) | standalone/web/IEDB | Field default; broadest allele coverage; presentation discovery |
| MHCflurry 2.0 | O'Donnell 2020 | BA + processing + presentation | pip Python | Scripting, messy allele strings, integrated presentation score |
| MixMHCpred 3.0 | Tadros 2025 | EL only (MS motifs) | standalone | MS-grounded presentation; cross-allele/species extrapolation study |
| NetMHC-4.0 | Andreatta 2016 | BA only | standalone/web | Legacy reproducibility; allele-specific, data-rich common alleles |
| MHCnuggets | Shao 2020 | BA (IC50) | pip Python | High-throughput TCGA-scale screens; rare-allele transfer learning |

## BA vs EL -- the conceptual axis that determines which score to read

BA (binding affinity) models train on in-vitro competitive-binding IC50 assays and measure only whether the groove can hold the peptide thermodynamically. EL (eluted-ligand / presentation) models train on mass-spec immunopeptidomics — peptides actually eluted from MHC on real cells — so the label implicitly folds in processing, transport, editing, and surface stability. Read BA for "could this bind the groove if delivered there"; read EL/presentation for "is this likely naturally presented," which is the default and recommended output of NetMHCpan-4.1, MHCflurry's presentation predictor, and MixMHCpred. IEDB codifies the split: `netmhcpan_ba` = recommended-binding, `netmhcpan_el` = recommended-epitope. Pick by intent.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Quick scriptable scan, messy allele names | MHCflurry presentation predictor | pip-installable, normalizes `A*02:01`/`A0201`/`HLA-A0201` |
| Maximal accuracy / broadest alleles | NetMHCpan-4.1 EL | Field default; concurrent MS motif deconvolution |
| Neoantigen candidate scoring | EL/presentation + expression check | EL under-ranks low-expression mutants (abundance bias) |
| "Can it physically bind" (engineered/delivered peptide) | BA mode (`-BA`, MHCflurry affinity) | Question is thermodynamic, not presentation |
| Rare / non-European allele | NetMHCpan-4.1, but verify training support | Pan-models extrapolate; confidence drops off the manifold |
| Cross-allele ranking in a multi-HLA patient | %Rank, never raw nM | nM scales differ per allele; nM cutoffs are allele-biased |

## Predict Presentation with MHCflurry

**Goal:** Score peptides against a patient genotype and report the best-presenting allele per peptide.

**Approach:** Load the presentation predictor; pass `alleles` as a sample->genotype dict so the model reports `best_allele`, `affinity` (nM), `affinity_percentile` (%Rank), and `presentation_score` (0-1, higher = more likely presented). Supply real `n_flanks`/`c_flanks` only if the genomic context is known.

```python
from mhcflurry import Class1PresentationPredictor

predictor = Class1PresentationPredictor.load()
df = predictor.predict(
    peptides=['SIINFEKL', 'GILGFVFTL', 'NLVPMVATV'],
    alleles={'patient1': ['HLA-A*02:01', 'HLA-A*24:02', 'HLA-B*07:02']},
    include_affinity_percentile=True,   # required: %Rank column is off by default
    verbose=0,
)
# columns: peptide, sample_name, affinity, best_allele, processing_score,
#          presentation_score, and affinity_percentile (only with the flag above)
# affinity nM: LOWER is stronger. presentation_score: HIGHER is more likely presented.
```

## Interpret with %Rank, Not Raw nM

**Goal:** Classify binding strength in a way that is comparable across alleles.

**Approach:** Threshold on %Rank (percentile of the score against random peptides for that same allele), not on absolute IC50. The 500 nM convention is allele-biased — it over-calls permissive alleles and under-calls restrictive ones, skewing a multi-HLA patient's epitope list toward a subset of the genotype.

```python
def classify_by_percentile(affinity_percentile):
    '''Class I %Rank cutoffs (NetMHCpan convention). LOWER percentile = stronger.
    Strong binder <= 0.5%; weak binder <= 2.0%. Use %Rank for any cross-allele
    comparison; raw nM is only meaningful within a single allele.'''
    if affinity_percentile <= 0.5:
        return 'strong'
    elif affinity_percentile <= 2.0:
        return 'weak'
    return 'non-binder'
```

## Scan a Protein for Class I Epitopes

**Goal:** Enumerate candidate epitopes across a protein for a patient genotype.

**Approach:** Tile 8-11mers (9mers dominate real ligands), score all windows in one batched call, keep windows under the 2% weak-binder cutoff. See examples/mhc_binding.py for the full tiling-and-rank script.

```python
def scan_protein(protein_seq, genotype, lengths=(8, 9, 10, 11)):
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    peptides = [protein_seq[i:i + k] for k in lengths for i in range(len(protein_seq) - k + 1)]
    df = predictor.predict(peptides=peptides, alleles={'patient': list(genotype)},
                           include_affinity_percentile=True, verbose=0)
    return df[df['affinity_percentile'] <= 2.0].sort_values('affinity_percentile')
```

## Per-Method Failure Modes

### Pan-model extrapolation on rare alleles
**Trigger:** scoring an allele with little/no training support (much of HLA-C, many non-European alleles). **Mechanism:** pan-models emit a confident %Rank for any allele sequence — there is no built-in "I don't know." **Symptom:** a flat/mushy predicted motif; calls that don't validate. **Fix:** check the allele is in the trained/supported list and that close pseudosequence neighbors had real ligands; downgrade confidence when extrapolating.

### EL abundance bias misranks neoantigens
**Trigger:** ranking mutated, low-expression peptides by EL/presentation score alone. **Mechanism:** MS immunopeptidomes over-represent abundant proteins; EL partly learns expression as a presentation proxy. **Symptom:** housekeeping-gene peptides float to the top; real low-expression neoantigens sink. **Fix:** combine EL with measured expression (TPM) and judge within-target, not against the proteome.

### Placeholder flanks corrupt the processing score
**Trigger:** passing dummy `n_flanks`/`c_flanks` to get a presentation/processing number. **Mechanism:** the processing model reads flanking context; wrong flanks inject noise. **Symptom:** processing_score that tracks nothing biological. **Fix:** supply the true genomic flanks, or omit flanks and read affinity/EL only.

### Allele in the list != well-trained on that allele
**Trigger:** trusting a number because the allele appears in `-listMHC`/`supported_alleles`. **Mechanism:** coverage is not training support. **Fix:** treat coverage and data depth as separate questions.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Class I strong binder <= 0.5% Rank | NetMHCpan-4.1 default (`-rth 0.5`) | Percentile normalizes per-allele score scales |
| Class I weak binder <= 2.0% Rank | NetMHCpan-4.1 default (`-rlt 2.0`) | Standard recall/precision balance for candidate lists |
| Peptide length 8-11mers (9 dominant) | Immunopeptidome composition | 9mers dominate training; non-9mers thinner evidence |
| IC50 <= 500 nM "strong" (legacy) | Pre-pan-allele convention | Allele-biased; AVOID for cross-allele work, use %Rank |
| 2-field (4-digit) HLA resolution | IMGT/HLA, groove determinants | Higher fields are synonymous/intronic; serotype is insufficient |
| Evaluate by PPV@top-N, not bare AUC | Zhao & Sher 2018; imbalance | True ligands ~1 in 10,000+; AUC is computed on an unreal balance |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Epitope list skewed to one HLA in a patient | Thresholded on nM not %Rank | Use affinity_percentile / %Rank |
| MHCflurry import/backend error | TF->PyTorch backend change (2.2.0+) | Use Python 3.10+; re-run `mhcflurry-downloads fetch` |
| Confident number on an untrained allele | Pan-model extrapolation | Verify supported allele + training depth |
| Low-expression neoantigen under-ranked | EL/MS abundance bias | Integrate expression; rank within-target |
| Reporting presentation as "immunogenicity" | Conflating funnel stages | Defer to immunogenicity-scoring; caveat the report |
| Class II call trusted like class I | Different maturity regime | Use mhc-class-ii-prediction; treat II as hypothesis |

## References

- Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. 2020. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research* 48(W1):W449-W454.
- O'Donnell TJ, Rubinsteyn A, Laserson U. 2020. MHCflurry 2.0: improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. *Cell Systems* 11(1):42-48.e7.
- Tadros DM, Racle J, Gfeller D, et al. 2025. Predicting MHC-I ligands across alleles and species: how far can we go? *Genome Medicine* 17:25.
- Andreatta M, Nielsen M. 2016. Gapped sequence alignment using artificial neural networks: application to the MHC class I system (NetMHC-4.0). *Bioinformatics* 32(4):511-517.
- Shao XM, Bhattacharya R, Huang J, et al. 2020. High-throughput prediction of MHC class I and II neoantigens with MHCnuggets. *Cancer Immunology Research* 8(3):396-408.
- Zhao W, Sher X. 2018. Systematically benchmarking peptide-MHC binding predictors: from synthetic to naturally processed epitopes. *PLOS Computational Biology* 14(11):e1006457.
- Trolle T, Metushi IG, Greenbaum JA, et al. 2015. Automated benchmarking of peptide-MHC class I binding predictions. *Bioinformatics* 31(13):2174-2181.

## Related Skills

- immunoinformatics/mhc-class-ii-prediction - CD4/HLA class II binding (the harder, less-reliable regime; open groove, register, DQ pairing)
- immunoinformatics/neoantigen-prediction - applies class I binding to tumor mutations; where the EL abundance bias bites
- immunoinformatics/immunogenicity-scoring - the separate, weaker prediction of T-cell response (binding != immunogenicity)
- immunoinformatics/epitope-prediction - T-cell epitope mapping reduces to MHC presentation; B-cell epitopes are a different problem
- clinical-databases/hla-typing - determine the patient genotype that conditions every prediction
