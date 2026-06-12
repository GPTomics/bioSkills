---
name: bio-phasing-imputation-foundations
description: Frames the phasing/imputation pipeline before any tool runs: phasing and imputation are one Li-Stephens copying HMM (recombination is the transition, mutation the emission, the genetic map and Ne set the rates), imputation's honest output is a dosage with a self-estimated quality (INFO/R2/DR2) not a hard genotype, and the stages are ordered and each fails silently (QC, align build and strand to the panel, phase, impute per chromosome, filter by INFO/R2 plus a MAF floor, carry dosages to GWAS). Covers the strategy fork (array vs low-coverage WGS plus genotype-likelihood imputation), why the panel ancestry is the prior, and why a flipped strand or build mismatch destroys accuracy without an error. Use when deciding a genotyping strategy, sequencing the pipeline, choosing array vs low-coverage WGS, or diagnosing silently-wrong imputation. Mechanics route to reference-panels, haplotype-phasing, genotype-imputation, imputation-qc; read-backed phasing is long-read-sequencing/haplotype-phasing.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

This skill is the conceptual spine for the category; it sequences the pipeline and routes the per-stage work to the sibling skills rather than running the engines, so the binding versions are the bcftools/Beagle/SHAPEIT5/Minimac4/GLIMPSE stack the siblings share. A reference panel is a dated release on a specific genome build (GRCh37 vs GRCh38); record the panel name, build, and release date with every result, because the same code against a different panel release imputes different sites.

# Phasing and Imputation Foundations -- Sequencing the Pipeline and Carrying the Uncertainty

**"Phase and impute my genotypes"** -> Decide the genotyping strategy, order the stages, and carry a dosage (not a hard call) downstream - because phasing and imputation are one statistical machine whose output is a posterior, and every stage of the pipeline fails silently when misordered or misaligned.
- CLI: confirm the genome build and strand, normalize to biallelic, phase, then impute per chromosome against an ancestry-matched panel, then filter on the imputer's own quality field plus a MAF floor

Scope: the meta-decisions - which genotyping strategy (array vs low-coverage WGS), the stage ordering, why dosages carry the uncertainty, why the panel is the prior, and where each stage silently fails. Panel acquisition/selection -> reference-panels. Statistical phasing engines -> haplotype-phasing. Imputation engines and dosage fields -> genotype-imputation. Quality metrics and filtering -> imputation-qc. Read-backed / molecular single-sample phasing (long reads, Hi-C, 10x) -> long-read-sequencing/haplotype-phasing. The input VCF and its QC -> variant-calling/vcf-basics. The GWAS test on the dosages -> population-genetics/association-testing.

## The Single Most Important Modern Insight -- Phasing and Imputation Are One Li-Stephens Copying HMM, and Its Honest Output Is a Dosage, Not a Genotype

There is no separate "phasing algorithm" and "imputation algorithm" at the core: both model a target haplotype as an imperfect mosaic of reference haplotypes, where switching which reference is copied is a recombination (an HMM transition whose rate is set by the genetic map and Ne) and a site where the copy disagrees is a mutation (an emission whose rate shrinks as the panel grows) (Li & Stephens 2003 *Genetics* 165:2213). Phasing asks which haplotype pair explains the diploid genotypes; imputation aligns those haplotypes onto the panel and reads the copied allele at untyped sites. Three facts drive every decision:

1. **The honest output of imputation is a DOSAGE, not a hard call.** The dosage is the posterior expected alternate-allele count `E[#alt]` in [0,2] (`= P(het) + 2*P(hom-alt)`), and the per-variant INFO/R2/DR2 is the model's self-estimate of how much real information that dosage carries (1 = perfect, 0 = nothing beyond allele frequency). Hard-thresholding the dosage into 0/1/2 and analyzing it as observed throws away the uncertainty imputation exists to quantify - it loses power at exactly the rare, low-R2 variants and, when imputation quality differs across batches, manufactures false associations that fail to replicate. Carry dosage (or the full genotype-probability triple) into a dosage-aware regression.
2. **The pipeline is ORDERED and every stage fails silently.** QC the array -> align genome build and strand to the panel -> phase -> impute per chromosome -> filter on INFO/R2 plus a MAF floor -> carry dosages to the GWAS. Skipping or misordering (imputing unphased input, a GRCh37/GRCh38 build mismatch, an unflipped A/T or C/G strand, no INFO filter) produces output that runs clean and is wrong. There is no exception or warning; the wrongness is only visible downstream as deflated accuracy or non-replicating hits.
3. **The panel is the prior, and its self-graded quality cannot see ancestry mismatch.** Imputation can only copy alleles the panel contains, so panel ancestry-match plus size sets the ceiling (Marchini & Howie 2010 *Nat Rev Genet* 11:499). When the target's ancestry is absent from the panel, the HMM still finds some template and reports a confident-looking R2 about a copy that is systematically wrong - high R2 in an under-represented ancestry is not reassurance. Pre-phasing (Howie 2012 *Nat Genet* 44:955) commits the study haplotypes once and propagates their switch errors into imputation, so phase well; the upside is that re-imputation against a newer panel needs no re-phasing.

## The Copying Model in One Paragraph

A new haplotype is modeled as copying a panel of K reference haplotypes, switching templates along the chromosome. The hidden state at each marker is which template is being copied; the transition probability of switching rises with the genetic distance between markers (in centimorgans, NOT base pairs) and with Ne through the `rho ~ 4*Ne` recombination scaling; the emission allows a mismatch with a probability `~ theta/(theta+K)` that shrinks as K grows (a bigger panel almost always has a near-perfect template, which is the quantitative root of "bigger panels impute better"). The forward-backward algorithm gives, at every site, the posterior over which template is copied; imputation sums that posterior to a dosage and reports its variance ratio as INFO/R2. IMPUTE/IMPUTE2, MaCH/minimac, Beagle, SHAPEIT, Eagle, and GLIMPSE are all this model with different approximations (which haplotypes to condition on, how PBWT makes the K-quadratic match-finding near-linear, EM vs MCMC, genotype-likelihood vs hard-call emission) - the differences are engineering, not different theories.

## The Strategy Fork -- Array + Imputation vs Low-Coverage WGS + Imputation

The upstream design decision, made before any engine is chosen, is how to generate the genotypes that will be imputed. It is not primarily an accuracy question; it is an ascertainment question.

| | SNP array + pre-phase + impute | Low-coverage WGS (~0.5-4x) + impute from genotype likelihoods |
|---|---|---|
| Input to the HMM | hard genotype calls (array error is tiny) | genotype LIKELIHOODS (PL/GL); a hard call at 1x is mostly noise |
| Ascertainment | FIXED - only the SNPs the array was designed for, biased to the design population | UNBIASED - whatever is in the genome is observed |
| Common variants | excellent and cheap | excellent |
| Rare variants | limited by the array scaffold and panel | matches or beats dense arrays (Rubinacci 2021 *Nat Genet* 53:120) |
| Under-represented ancestry | poor (no good array, panel-mismatched) | the main route around array/panel bias |
| Tools | Beagle / Minimac4 / IMPUTE5 -> genotype-imputation | GLIMPSE2 (panel) / STITCH (no panel) -> genotype-imputation |

The judgment: common-variant GWAS in a well-paneled ancestry -> array plus imputation is cheap and adequate; rare variants, under-represented ancestry, or a need for unbiased genome-wide ascertainment -> low-coverage WGS plus genotype-likelihood imputation, the direction the field is moving as sequencing costs fall. Low-coverage WGS is not magic: it is only as good as its panel and its likelihoods (bad mapping, contamination, or damage produce garbage GLs that impute garbage).

## Tool Taxonomy

| Tool / panel | Citation | Role | Owns / route to |
|--------------|----------|------|------------------|
| Beagle | Browning 2018 *Am J Hum Genet* 103:338 | one tool that BOTH phases and imputes; bref3 panel format | haplotype-phasing, genotype-imputation |
| SHAPEIT5 | Hofmeister 2023 *Nat Genet* 55:1243 | scalable phasing incl. rare/singleton; split into phase_common/phase_rare/ligate | haplotype-phasing |
| Eagle2 | Loh 2016 *Nat Genet* 48:1443 | PBWT reference-based phasing; the Michigan-server phaser | haplotype-phasing |
| Minimac4 | Das 2016 *Nat Genet* 48:1284 | imputation engine; msav/m3vcf; the Michigan/TOPMed server engine | genotype-imputation |
| IMPUTE5 | Rubinacci 2020 *PLoS Genet* 16:e1009049 | PBWT-based imputation for large panels | genotype-imputation |
| GLIMPSE2 | Rubinacci 2023 *Nat Genet* 55:1088 | low-coverage WGS imputation from genotype likelihoods | genotype-imputation |
| STITCH | Davies 2016 *Nat Genet* 48:965 | low-coverage imputation with NO reference panel (latent panel via EM) | genotype-imputation |
| PBWT | Durbin 2014 *Bioinformatics* 30:1266 | the data structure making biobank-scale haplotype matching near-linear | (underlies the phasers/imputers) |
| 1000G / HRC / TOPMed | 1000G 2015 *Nature* 526:68; McCarthy 2016 *Nat Genet* 48:1279; Taliun 2021 *Nature* 590:290 | the reference panels (the prior) | reference-panels |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| SNP array, common-variant GWAS, well-paneled ancestry | array -> pre-phase -> impute against a matched panel | mature, cheap, accurate for common variants |
| Rare variants, under-represented ancestry, or unbiased ascertainment wanted | low-coverage WGS -> GLIMPSE2 (or STITCH if no panel) -> genotype-imputation | removes array ascertainment bias; wins at rare variants |
| Moderate cohort, want one tool for phase and impute | Beagle -> haplotype-phasing, genotype-imputation | does both with the bref3 panel format |
| Biobank-scale cohort (100k+), phasing | SHAPEIT5 -> haplotype-phasing | scales; phases rare variants on a common scaffold |
| Choosing/preparing the panel for a target ancestry | -> reference-panels | panel ancestry-match is the first-order accuracy lever |
| Filtering imputed output before analysis | -> imputation-qc | INFO/R2/DR2 + MAF floor on the imputer's own field |
| Single individual, long reads / Hi-C / linked reads | -> long-read-sequencing/haplotype-phasing | read-backed phasing, a physically different signal |
| The association test on the dosages | -> population-genetics/association-testing | the GWAS test is downstream |
| Establishing target ancestry for panel choice | -> population-genetics/population-structure | PCA, not phasing |

Default when uncertain: for a common-variant GWAS in a well-served ancestry, array plus imputation against a large matched panel; otherwise consider low-coverage WGS plus GLIMPSE2.

## The Three Physically Distinct Sources of Phase

Phase information comes from three non-interchangeable physical sources; they resolve different things and the modern answer is to combine, not rank, them. This skill and its phasing sibling own the first; the third is a mode of the same tools; the second is a different category routed out.

- **Statistical (population LD).** The Li-Stephens copying model. Genome-wide reach from one sample, but it cannot phase a truly private variant (no population shares it, so there is nothing to copy) -> haplotype-phasing.
- **Read-backed (molecular).** Two heterozygous sites observed on the same physical DNA molecule (long reads, Hi-C, linked reads). The only way to phase a private variant, but reach is the molecule length -> long-read-sequencing/haplotype-phasing.
- **Pedigree / trio (Mendelian).** Deterministic and genome-wide where parents are informative; the gold standard for validating the others, useless where parents are uninformative -> trio mode of SHAPEIT/Beagle.

## Per-Method Failure Modes

### Genome build or strand not aligned to the panel
**Trigger:** imputing GRCh37 array data against a GRCh38 panel, or leaving A/T and C/G (palindromic) SNPs unflipped. **Mechanism:** positions or alleles silently disagree with the panel, so the HMM copies the wrong templates. **Symptom:** imputation accuracy near zero over whole regions, with no error and a plausible-looking VCF. **Fix:** confirm the build, lift over if needed, and align strand/alleles to the panel (`bcftools +fixref` or the imputation-server check tool) before phasing -> reference-panels.

### Hard-calling the dosage
**Trigger:** thresholding DS into 0/1/2 and analyzing those as observed genotypes. **Mechanism:** thresholding discards the posterior uncertainty, worst at low-R2 rare variants. **Symptom:** lost power read as a true null, or false hits from batch-differential imputation quality. **Fix:** regress on dosage or the genotype-probability triple -> imputation-qc, population-genetics/association-testing.

### Trusting R2 in an under-represented ancestry
**Trigger:** reporting a confident INFO/R2 for a sample whose ancestry the panel barely contains. **Mechanism:** R2 is the model grading its own posterior; it cannot see that no panel haplotype resembled the target. **Symptom:** high self-graded quality on systematically wrong copies. **Fix:** match the panel ancestry, prefer low-coverage WGS for under-served ancestries, and validate with masked/held-out truth stratified by frequency -> reference-panels, imputation-qc.

### Switch errors invisible to genotype QC
**Trigger:** reporting phasing by tool name with no switch-error rate. **Mechanism:** a switch error changes which haplotype an allele sits on without changing any genotype, so per-site QC is structurally blind to it. **Symptom:** genotype-concordant data that is phased garbage; a single switch inverts a compound-het (cis vs trans) conclusion. **Fix:** report a switch-error rate against a trio or read-backed truth set for any phase-dependent claim -> haplotype-phasing.

### Expecting imputation to recover a private variant
**Trigger:** assuming a large panel "recovers" ultra-rare or sample-private variants. **Mechanism:** a private allele is carried by no panel haplotype, so there is no template to copy - it is un-imputable by construction. **Symptom:** the variant is absent or imputed at near-zero R2. **Fix:** sequence (do not impute) private variants; use larger panels only for the merely-rare.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Dosage range E[#alt] in [0,2], not {0,1,2} | Li & Stephens 2003 *Genetics* 165:2213; Marchini & Howie 2010 *Nat Rev Genet* 11:499 | the dosage is the posterior mean alt-count; analyze on dosage to propagate uncertainty |
| INFO/R2 filter at 0.3 (common variants), 0.7-0.8 stricter | Marchini & Howie 2010 *Nat Rev Genet* 11:499 (framing) | R2<0.3 means the dosage carries little real information; the exact cutoff is a study choice, not a law -> imputation-qc |
| Genetic distance in centimorgans (1 cM ~ 1% recombination/meiosis; 1 Morgan = 100 cM) | Li & Stephens 2003 *Genetics* 165:2213 | LS transitions use genetic, not physical, distance; a wrong-build or flat map silently degrades phasing |
| Low-coverage WGS sweet spot ~0.5-4x (often ~1x) with a large panel | Rubinacci 2021 *Nat Genet* 53:120 | GLIMPSE operating range; below ~0.5x accuracy degrades, above ~4x hard-calling becomes viable |
| HRC accurate to MAF ~0.1% | McCarthy 2016 *Nat Genet* 48:1279 | 64,976 haplotypes; larger panels (TOPMed) push the frequency floor lower |
| Pre-phase the cohort once, re-impute against any panel | Howie 2012 *Nat Genet* 44:955 | imputation collapses from the diploid (K^2) to the haploid (K) problem; phasing is committed, so phase well |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Imputation accuracy near zero across a region | genome build or strand mismatch to the panel | confirm build, lift over, align strand before phasing -> reference-panels |
| Hits do not replicate across batches | hard-called dosages plus batch-differential imputation quality | analyze dosages; impute batches together -> imputation-qc |
| All R2 low for a cohort | panel ancestry mismatch or wrong-build genetic map | use a matched panel; check the map build -> reference-panels |
| Rare-variant burden differs across cohorts | a single global R2 cutoff removes more rare variants than common, differently per panel | use frequency-stratified accuracy; match panels -> imputation-qc |
| Engine errors or mis-imputes at a site | multiallelic / non-biallelic input | `bcftools norm -m -any` to split and left-align first -> variant-calling/variant-normalization |
| Expected variant missing after imputation | private / panel-absent allele | un-imputable by construction; sequence it |

## References

- Li N, Stephens M. 2003. Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data. *Genetics* 165:2213-2233.
- Marchini J, Howie B. 2010. Genotype imputation for genome-wide association studies. *Nat Rev Genet* 11:499-511.
- Howie B, Fuchsberger C, Stephens M, Marchini J, Abecasis GR. 2012. Fast and accurate genotype imputation in genome-wide association studies through pre-phasing. *Nat Genet* 44:955-959.
- Browning BL, Zhou Y, Browning SR. 2018. A one-penny imputed genome from next-generation reference panels. *Am J Hum Genet* 103:338-348.
- Hofmeister RJ, Ribeiro DM, Rubinacci S, Delaneau O. 2023. Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. *Nat Genet* 55:1243-1249.
- Loh PR, Danecek P, Palamara PF, et al. 2016. Reference-based phasing using the Haplotype Reference Consortium panel. *Nat Genet* 48:1443-1448.
- Das S, Forer L, Schonherr S, et al. 2016. Next-generation genotype imputation service and methods. *Nat Genet* 48:1284-1287.
- Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O. 2021. Efficient phasing and imputation of low-coverage sequencing data using large reference panels. *Nat Genet* 53:120-126.
- Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. 2023. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. *Nat Genet* 55:1088-1090.
- Davies RW, Flint J, Myers S, Mott R. 2016. Rapid genotype imputation from sequence without reference panels. *Nat Genet* 48:965-969.
- Durbin R. 2014. Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT). *Bioinformatics* 30:1266-1272.
- McCarthy S, Das S, Kretzschmar W, et al. 2016. A reference panel of 64,976 haplotypes for genotype imputation. *Nat Genet* 48:1279-1283.
- Taliun D, Harris DN, Kessler MD, et al. 2021. Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program. *Nature* 590:290-299.

## Related Skills

- reference-panels - Select and prepare the panel (the prior); ancestry-match, build, strand alignment, server access
- haplotype-phasing - Statistical phasing engines and the switch-error metric
- genotype-imputation - Imputation engines, dosage fields, servers, and low-coverage WGS
- imputation-qc - INFO/R2/DR2 quality metrics, MAF-stratified accuracy, dosage filtering
- long-read-sequencing/haplotype-phasing - Read-backed / molecular single-sample phasing (a physically different signal)
- variant-calling/vcf-basics - Source and basic handling of the input VCF
- variant-calling/variant-normalization - Split multiallelics and left-align before phasing/imputation
- population-genetics/association-testing - GWAS test on the imputed dosages
- population-genetics/population-structure - PCA to establish target ancestry for panel selection
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate -> visualize
