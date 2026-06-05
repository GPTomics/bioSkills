---
name: bio-immunoinformatics-epitope-prediction
description: Predict B-cell and T-cell epitopes for vaccine antigen design and epitope mapping with BepiPred-3.0, DiscoTope-3.0, the IEDB tools, and EL-mode MHC presentation. Encodes the load-bearing asymmetry that T-cell epitope prediction is mature (it reduces to MHC presentation, AUC>0.9) while B-cell prediction is unreliable (linear predictors ~AUC 0.6 because ~90% of real epitopes are conformational) — so structure-based DiscoTope-3.0 on AlphaFold models is the only defensible B-cell path, propensity scales are obsolete, and NetChop is largely redundant on EL-trained models. Use when mapping epitopes or selecting vaccine antigens. MHC binding lives in mhc-binding-prediction.
tool_type: python
primary_tool: BepiPred
---

## Version Compatibility

Reference examples tested with: BepiPred-3.0, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: BepiPred-3.0 ships as the `bepipred3` package and auto-downloads ESM-2 weights on first run; its default threshold is 0.1512 (NOT 0.5). DiscoTope-3.0, ElliPro, SEPPA, NetChop, and NetCTLpan are standalone/web (IEDB or DTU). The IEDB classic and next-generation REST APIs wrap most predictors. Re-verify thresholds and the supported-method list against current docs.

# Epitope Prediction

**"Predict the B-cell and T-cell epitopes in my antigen"** -> Identify antibody-binding (B-cell) and MHC-presented (T-cell) immunogenic regions, with appropriately different confidence for each.
- Python: `bepipred3` for linear B-cell epitopes; IEDB REST API for B-cell/T-cell tools
- CLI/web: DiscoTope-3.0 for conformational B-cell epitopes (structure-based); NetMHCpan/MHCflurry (EL) for T-cell epitopes

## The Single Most Important Modern Insight -- "epitope prediction" is two fields at different maturity, wrongly conflated

T-cell epitope prediction is mature and trustworthy because it reduces to MHC binding/presentation — a sharply constrained problem (a peptide fits the groove or it does not) with an enormous mass-spec eluted-ligand training corpus; NetMHCpan-4.1 and MHCflurry routinely exceed AUC 0.9 for class I. B-cell epitope prediction is unreliable: linear sequence-based predictors land around AUC 0.6, and even the ESM-2-based BepiPred-3.0 falls to AUC 0.663 on the real IEDB external test set. This is structural, not a tuning problem the next network will fix: ~90% of natural B-cell epitopes are conformational/discontinuous — residues clustered in 3D but far apart in sequence — which a sequence-only model is by construction blind to. The single most damaging mistake in this domain is letting the well-deserved confidence in MHC/T-cell prediction leak into unwarranted confidence in B-cell prediction. Write down which problem is being solved before running anything.

## Tool Taxonomy

| Tool | Citation | Target | Input | When |
|------|----------|--------|-------|------|
| NetMHCpan-4.1 EL / MHCflurry | Reynisson 2020; O'Donnell 2020 | T-cell (MHC-I presentation) | sequence + HLA | Default T-cell path; EL encodes processing |
| NetMHCIIpan / NetCTLpan | Nilsson 2023; Stranzl 2010 | T-cell (CD4 / integrated CTL) | sequence + HLA | CD4 epitopes; integrated cleavage+TAP+MHC |
| DiscoTope-3.0 | Høie 2024 | B-cell (conformational) | 3D structure (AlphaFold OK) | The only defensible B-cell method when a structure exists |
| BepiPred-3.0 | Clifford 2022 | B-cell (linear) | sequence | Linear/denatured-target reagents; misses ~90% native |
| ElliPro / SEPPA 3.0 | Ponomarenko 2008; Zhou 2019 | B-cell (conformational) | 3D structure | Fast geometric baseline; SEPPA for glycoproteins |
| Propensity scales | Kolaskar 1990 etc. | B-cell (linear) | sequence | Obsolete; decoration, not data |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| T-cell (CD8) epitopes | NetMHCpan-4.1 EL / MHCflurry | Mature; defer to mhc-binding-prediction |
| T-cell (CD4) epitopes | NetMHCIIpan-4.3 | Defer to mhc-class-ii-prediction; less reliable |
| B-cell, structure available or foldable | DiscoTope-3.0 on AlphaFold model | Conformational; ~no penalty for predicted structures |
| B-cell glycoprotein (Env/S/HA) | SEPPA 3.0 | Models glycan shielding |
| B-cell, sequence only, peptide/denatured target | BepiPred-3.0 (linear/top-X%) | Legitimate narrow use; state the conformational caveat |
| B-cell, sequence only, native antibody response | Fold a structure first, then DiscoTope-3.0 | Linear prediction structurally cannot see native epitopes |
| Broadly-protective vaccine | + conservation + HLA population coverage | A high-scoring epitope in a hypervariable loop is worthless |

## Predict Linear B-Cell Epitopes (BepiPred-3.0)

**Goal:** Score per-residue linear B-cell epitope probability from sequence, for a linear/denatured-target use case.

**Approach:** Run the `bepipred3` CLI (or package) on a FASTA; it emits per-residue probabilities, a binary FASTA (upper = epitope), and top-X% selections. Use the default threshold 0.1512 or the top-X% mode; treat output as a hypothesis that misses most native conformational epitopes.

```bash
# bepipred3 auto-downloads ESM-2 weights on first run; default threshold 0.1512 (NOT 0.5)
python bepipred3_CLI.py -i antigen.fasta -o bp3_out/ -pred vt_pred -t 0.1512
# or select the top 20% scoring residues per sequence instead of a fixed cutoff:
python bepipred3_CLI.py -i antigen.fasta -o bp3_out/ -pred vt_pred -top 20
```

## Predict Conformational B-Cell Epitopes (DiscoTope-3.0)

**Goal:** Identify antibody-accessible surface patches from a 3D structure (the defensible B-cell path).

**Approach:** Provide a single antigen chain (experimental or AlphaFold). DiscoTope-3.0 scores per-residue conformational propensity and was trained on predicted structures, so AF2 models incur essentially no penalty (AUC 0.799 vs 0.807). Gate trust by pLDDT — accuracy drops ~5 percentile points per 10-point pLDDT decrease — and remember AUC-PR is only ~0.22 (low precision, many false positives).

```python
def gate_discotope_by_plddt(df, plddt_col='pLDDT', score_col='DiscoTope-3.0 score', min_plddt=70):
    '''Keep DiscoTope-3.0 calls only in confidently-folded regions; low-pLDDT loops
    (where antibodies often bind) are exactly where structure-based calls are least
    reliable. df: per-residue DiscoTope-3.0 output joined with model pLDDT.'''
    return df[df[plddt_col] >= min_plddt].sort_values(score_col, ascending=False)
```

## T-Cell Epitopes Reduce to MHC Presentation

**Goal:** Nominate CD8/CD4 epitopes from an antigen.

**Approach:** Tile the antigen and score with EL-mode MHC presentation (class I: mhc-binding-prediction; class II: mhc-class-ii-prediction). Do NOT add NetChop by default — EL models are trained on eluted ligands that already survived proteasomal cleavage and TAP, so the processing signal is implicit; explicit cleavage prediction is largely redundant and can double-penalize. Reserve NetChop/NetCTLpan for long source proteins as a cleavage sanity check or alleles lacking EL coverage.

## Per-Method Failure Modes

### Linear predictor used for native antibody response
**Trigger:** running BepiPred on a folded viral spike to predict neutralizing epitopes. **Mechanism:** native epitopes are conformational; sequence models cannot see them. **Symptom:** "predicted epitopes" that no native antibody targets. **Fix:** fold a structure and use DiscoTope-3.0; reserve linear predictors for peptide/denatured targets.

### Predicting epitopes of a wrong model
**Trigger:** DiscoTope on a low-confidence AlphaFold surface loop or a monomer of an oligomeric antigen. **Mechanism:** a subtly wrong surface moves the predicted epitope; an oligomer interface looks exposed in the monomer. **Symptom:** false-positive epitopes at buried/flexible sites. **Fix:** gate by pLDDT; model the biological assembly when the antigen oligomerizes.

### Propensity-scale cargo cult
**Trigger:** reporting Kolaskar-Tongaonkar/Parker/Emini "antigenic regions" as data. **Mechanism:** these are coarse 1980s physicochemical descriptors at/near random. **Symptom:** confident-looking but uninformative B-cell calls. **Fix:** treat as obsolete decoration; everything they encode is subsumed by BepiPred/structure methods.

### Confusing presentation with immunodominance
**Trigger:** ranking vaccine epitopes purely by binding/presentation score. **Mechanism:** immunodominance depends on repertoire, competition, processing kinetics, immune history — none modeled. **Symptom:** a strong predicted binder that is subdominant or ignored in vivo. **Fix:** treat presentation as necessary-not-sufficient; validate by ELISpot/tetramer.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| BepiPred-3.0 default 0.1512 | Clifford 2022 | Balances sens/spec on their benchmark; NOT 0.5 |
| Linear B-cell AUC ~0.6 | Field benchmarks | Barely above random; report as hypothesis |
| DiscoTope-3.0 AUC-ROC ~0.80, AUC-PR ~0.22 | Høie 2024 | Moderate ranking, low precision (minority class) |
| pLDDT >= 70 to trust DiscoTope calls | Høie 2024 | ~5 percentile-point drop per 10-point pLDDT loss |
| ~90% of B-cell epitopes conformational | B-cell literature | Why sequence-only prediction has a low ceiling |
| Skip NetChop on EL-mode predictions | Reynisson 2020 | EL training already encodes cleavage/TAP |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Over-trusting B-cell predictions | Conflated with mature T-cell prediction | State the maturity asymmetry; treat B-cell as hypothesis |
| Few/no BepiPred epitopes | Applied 0.5 threshold | Use default 0.1512 or top-X% mode |
| False epitopes in flexible loops | Low-pLDDT AlphaFold model | Gate by pLDDT; assess model quality |
| Epitope worthless across strains | No conservation analysis | Add IEDB Epitope Conservancy + MSA |
| Redundant/over-penalized T-cell calls | NetChop stacked on EL model | Use EL presentation as the primary filter |
| Vaccine "designed" in silico | Over-trusting reverse-vaccinology scores | Treat VaxiJen/Vaxign as candidate funnels; validate experimentally |

## References

- Clifford JN, Høie MH, Deleuran S, Peters B, Nielsen M, Marcatili P. 2022. BepiPred-3.0: improved B-cell epitope prediction using protein language models. *Protein Science* 31(12):e4497.
- Høie MH, Gade FS, Johansen JM, et al. 2024. DiscoTope-3.0: improved B-cell epitope prediction using inverse folding latent representations. *Frontiers in Immunology* 15:1322712.
- Jespersen MC, Peters B, Nielsen M, Marcatili P. 2017. BepiPred-2.0: improving sequence-based B-cell epitope prediction using conformational epitopes. *Nucleic Acids Research* 45(W1):W24-W29.
- Kringelum JV, Lundegaard C, Lund O, Nielsen M. 2012. Reliable B cell epitope predictions: impacts of method development and improved benchmarking (DiscoTope-2.0). *PLoS Computational Biology* 8(12):e1002829.
- Ponomarenko J, Bui HH, Li W, et al. 2008. ElliPro: a new structure-based tool for the prediction of antibody epitopes. *BMC Bioinformatics* 9:514.
- Stranzl T, Larsen MV, Lundegaard C, Nielsen M. 2010. NetCTLpan: pan-specific MHC class I pathway epitope predictions. *Immunogenetics* 62(6):357-368.
- Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. 2020. NetMHCpan-4.1 and NetMHCIIpan-4.0. *Nucleic Acids Research* 48(W1):W449-W454.
- Calis JJA, Maybeno M, Greenbaum JA, et al. 2013. Properties of MHC class I presented peptides that enhance immunogenicity. *PLoS Computational Biology* 9(10):e1003266.
- Bui HH, Sidney J, Li W, Fusseder N, Sette A. 2007. Development of an epitope conservancy analysis tool. *BMC Bioinformatics* 8:361.

## Related Skills

- immunoinformatics/mhc-binding-prediction - T-cell (CD8) epitope prediction reduces to class I presentation
- immunoinformatics/mhc-class-ii-prediction - T-cell (CD4) epitopes; the class II presentation regime
- immunoinformatics/immunogenicity-scoring - ranking epitope candidates by likely T-cell response
- structural-biology/alphafold-predictions - fold an antigen with AlphaFold to enable DiscoTope-3.0
- database-access/entrez-fetch - retrieve antigen sequences/structures for epitope mapping
