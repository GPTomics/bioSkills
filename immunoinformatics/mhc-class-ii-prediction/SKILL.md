---
name: bio-immunoinformatics-mhc-class-ii-prediction
description: Predict peptide-MHC class II (HLA-DR/DQ/DP) binding and presentation for CD4 T-cell epitopes with NetMHCIIpan-4.3 and MixMHC2pred-2.0. Covers why class II is far less reliable than class I (open binding groove, 9-mer register ambiguity, sparse noisy training data, DR>DP>DQ accuracy asymmetry), the DQ/DP heterodimer alpha/beta pairing trap, and the looser 1%/5% %Rank thresholds. Use when predicting CD4 epitopes for vaccine help, mapping class II neoantigens, or scoring long peptides against DR/DQ/DP. For CD8/class I see mhc-binding-prediction.
tool_type: cli
primary_tool: NetMHCIIpan
---

## Version Compatibility

Reference examples tested with: NetMHCIIpan 4.3+, MixMHC2pred 2.0+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: NetMHCIIpan and MixMHC2pred are standalone academic binaries (not pip-installable); the IEDB MHC-II REST API wraps NetMHCIIpan if a local install is unavailable. Allele nomenclature differs sharply between tools and between isotypes (DR single-chain vs DQ/DP heterodimer) — confirm the exact string format against the installed tool before scoring. Class II %Rank thresholds (1%/5%) are LOOSER than class I (0.5%/2.0%); do not copy class I cutoffs.

# MHC Class II Prediction

**"Predict which long peptides bind/are presented by HLA class II"** -> Score peptide-HLA class II (DR/DQ/DP) presentation to nominate candidate CD4 T-cell epitopes, inferring the 9-mer binding core within each long peptide.
- CLI: `NetMHCIIpan` (field default; EL score by default, `-BA` adds affinity; pan-DR/DQ/DP)
- CLI: `MixMHC2pred` (MS-deconvolution motifs; models the reverse DP binding mode)

## The Single Most Important Modern Insight -- class II is basically broken, and that must be stated plainly

For class I, modern pan-allele EL predictors recover most true ligands at high precision and the field has hit diminishing returns. The same architectures, on the same conceptual pipeline, produce dramatically weaker class II predictions. The honest one-line summary to give a collaborator is: "trust a class I strong-binder call; treat a class II call as a ranked hypothesis, not a fact." Four compounding reasons, not one, cause this. The groove is open at both ends, so a 12-25mer can sit in multiple registers and the model must infer which latent 9-residue core is the true binding frame — an error-prone latent-variable problem class I (closed groove, defined termini) never faces. The training data are smaller and noisier: in-vitro class II binding assays are notoriously irreproducible, and class II immunopeptidomics yields fewer, longer, more heterogeneous peptides. The three isotypes are unequally tractable — historically DR >> DP > DQ, because DR was studied first and most while DQ was data-starved (NetMHCIIpan-4.3's headline 2023 contribution was finally closing this gap with tailored data acquisition, a sign of how recent and data-driven the fix is). And DP/DQ are obligate alpha/beta heterodimers whose chains are independently polymorphic, so the effective number of distinct molecules is the combinatorial product of alpha and beta alleles.

## The DQ/DP heterodimer pairing trap

A donor's DQA1 and DQB1 alleles pair both in cis (same haplotype) and in trans (across haplotypes), so a heterozygous individual can express up to four DQ heterodimers — and some trans-pairs are non-functional or rare. Mechanically feeding all DQA1 x DQB1 combinations to NetMHCIIpan generates molecules that do not biologically exist; taking only cis pairs may miss real trans-dimers. There is no fully automated, universally agreed resolution. The expert move is to be explicit about the pairing assumption, prefer documented haplotype pairings, and flag DQ (and to a lesser extent DP) epitope calls as lower-confidence than DR. DR is single-chain (the alpha is effectively invariant), so it carries none of this combinatorial burden and is the most trustworthy isotype.

## Tool Taxonomy (Class II)

| Tool | Citation | Score type | Loci | Use when |
|------|----------|-----------|------|----------|
| NetMHCIIpan-4.3 | Nilsson 2023 | EL (default) + BA (`-BA`) | DR, DQ, DP (+ mouse H-2, BoLA) | Field default; broadest coverage; closes DQ gap; reverse-mode binders |
| MixMHC2pred-2.0 | Racle 2023 | EL only (MS motifs) | DR, DQ, DP | MS-grounded motifs; models reverse (C->N) DP binding mode |
| MHCnuggets | Shao 2020 | BA (IC50) | class I + II | High-throughput screens; rare-allele transfer learning |
| NetMHCIIpan-4.0 | Reynisson 2020 | EL + BA | DR, DQ, DP | Reproducing 2020-era results; superseded by 4.3 |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Default CD4 epitope screen | NetMHCIIpan-4.3 EL | Broadest, pan-allele, current accuracy |
| DR-restricted, want highest confidence | NetMHCIIpan-4.3 (DRB1_*) | DR is the most reliable isotype (single-chain) |
| DQ or DP restriction | NetMHCIIpan-4.3 + explicit pairing | Heterodimer combinatorics; flag as lower-confidence |
| MS-grounded motif / DP reverse binders | MixMHC2pred-2.0 | Built from deconvolved immunopeptidomes; models reverse mode |
| Class II neoantigens (CD4 help) | NetMHCIIpan-4.3 EL + expression | CD4 help boosts vaccine efficacy; EL still abundance-biased |
| No local install | IEDB MHC-II REST API | Wraps NetMHCIIpan, always-current versions |

## Running the Predictions (the fragile commands - run as written)

NetMHCIIpan reads a peptide list or FASTA and scores against one or more alleles. EL %Rank is the default output; `-BA` adds an affinity prediction. The model reports the inferred 9-mer core and its offset.

```bash
# DR (single-chain: beta allele names the molecule)
netMHCIIpan -f peptides.txt -inptype 1 -a DRB1_0101 -BA -xls -xlsfile out.tsv

# DQ heterodimer (BOTH chains, hyphen-joined) and DP
netMHCIIpan -f antigen.fasta -a HLA-DQA10501-DQB10201,HLA-DPA10103-DPB10401 -length 15
```
Key flags: `-a` allele(s, comma-separated), `-f` input, `-inptype` (0=FASTA, 1=peptide list), `-length` peptide length(s) to consider, `-BA` add affinity, `-xls`/`-xlsfile` tab output, `-list` dump supported alleles.

MixMHC2pred uses chain-underscore allele names with a DOUBLE underscore between heterodimer chains, and alleles are space-separated:
```bash
MixMHC2pred -i peptides.txt -o out.txt -a DRB1_15_01 DRB5_01_01 DPA1_02_01__DPB1_01_01
```

## Per-Method Failure Modes

### Register ambiguity in the open groove
**Trigger:** any class II prediction on a long peptide. **Mechanism:** the 9-mer binding core can sit in several frames; the model infers the latent core. **Symptom:** the reported core shifts with small input changes; unstable rankings. **Fix:** treat the call as a hypothesis; corroborate with MixMHC2pred and check core consistency; never over-interpret a single offset.

### DQ/DP heterodimer mis-pairing
**Trigger:** expanding a genotype to all DQA1 x DQB1 (or DPA1 x DPB1) combinations. **Mechanism:** not all alpha/beta pairs form stable functional dimers; trans-pairs may be rare. **Symptom:** epitope calls against molecules that do not exist in the donor. **Fix:** restrict to documented/cis pairings, state the assumption, flag DQ/DP as lower-confidence than DR.

### Copying class I thresholds
**Trigger:** applying 0.5%/2.0% %Rank to class II. **Mechanism:** class II distributions and recommended cutoffs differ. **Symptom:** over-stringent filtering, missed real binders. **Fix:** use class II cutoffs (strong <= 1%, weak <= 5%).

### EL abundance bias (shared with class I)
**Trigger:** ranking class II neoantigens by EL alone. **Mechanism:** MS immunopeptidomes over-represent abundant proteins. **Symptom:** low-expression CD4 neoantigens under-ranked. **Fix:** integrate expression; judge within-target.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Class II strong binder <= 1% Rank | NetMHCIIpan-4.x default | Looser than class I; reflects class II score distributions |
| Class II weak binder <= 5% Rank | NetMHCIIpan-4.x default | Standard recall/precision balance for class II |
| Peptide length 12-25mers (core = 9) | Open-groove biology | Class II ligands are long with ragged termini; core always 9 |
| MixMHC2pred input 12-21mers | Racle 2023 | Outside this range or non-standard residues -> NA |
| 2-field (4-digit) typing for both chains | IMGT/HLA | Both alpha and beta needed for DQ/DP heterodimers |
| Isotype confidence DR > DP > DQ | Nilsson 2023 | Reflects historical training-data depth per isotype |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Allele not recognized | Wrong nomenclature for the tool | NetMHCIIpan `DRB1_0101`/`HLA-DQA10501-DQB10201`; MixMHC2pred `DRB1_15_01`/`DPA1_02_01__DPB1_01_01` |
| Calls against non-existent molecules | Naive DQ/DP combinatorial expansion | Use cis/documented pairings; flag DQ/DP |
| Over-stringent, few binders | Class I cutoffs applied | Use 1%/5% class II thresholds |
| Unstable core/offset | Register ambiguity | Corroborate across tools; treat as hypothesis |
| `NA` scores from MixMHC2pred | Peptide outside 12-21mer / non-standard residue | Filter input length and alphabet first |
| Class II trusted like class I | Different maturity regime | Report as ranked hypotheses, not facts |

## References

- Nilsson JB, Kaabinejadian S, Yari H, et al. 2023. Accurate prediction of HLA class II antigen presentation across all loci using tailored data acquisition and refined machine learning (NetMHCIIpan-4.3). *Science Advances* 9(47):eadj6367.
- Racle J, Guillaume P, Schmidt J, et al. 2023. Machine learning predictions of MHC-II specificities reveal alternative binding mode of class II epitopes (MixMHC2pred-2.0). *Immunity* 56(6):1359-1375.e13.
- Racle J, Michaux J, Rockinger GA, et al. 2019. Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes (MixMHC2pred-1.0). *Nature Biotechnology* 37:1283-1286.
- Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. 2020. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research* 48(W1):W449-W454.
- Shao XM, Bhattacharya R, Huang J, et al. 2020. High-throughput prediction of MHC class I and II neoantigens with MHCnuggets. *Cancer Immunology Research* 8(3):396-408.

## Related Skills

- immunoinformatics/mhc-binding-prediction - CD8/HLA class I binding (the solved regime; closed groove, 0.5%/2.0% cutoffs)
- immunoinformatics/neoantigen-prediction - class II neoantigens for CD4 help; pVACseq runs both classes
- immunoinformatics/immunogenicity-scoring - CD4 immunogenicity is even less solved than CD8
- immunoinformatics/epitope-prediction - T-cell epitope prediction reduces to MHC presentation
- clinical-databases/hla-typing - resolve DR/DQ/DP alleles for both chains
