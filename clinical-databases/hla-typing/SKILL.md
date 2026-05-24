---
name: bio-clinical-databases-hla-typing
description: Calls HLA class I and class II alleles at 2/4/6/8-field resolution from WGS/WES/RNA-seq/long-read data using OptiType, HLA-LA, T1K, Polysolver, HLA-HD, arcasHLA, StarPhase, or HIBAG imputation. Use when typing for HSCT, solid-organ transplant, neoantigen prediction, PGx screening (B*57:01, B*15:02, etc.), or disease-association studies, with reconciliation across tools and IPD-IMGT/HLA version mismatch handling.
tool_type: cli
primary_tool: T1K
---

## Version Compatibility

Reference examples tested with: OptiType 1.3.5, HLA-LA 1.0.4, T1K 1.0.6 (Song 2023), Polysolver 4.0, HLA-HD 1.7.1, arcasHLA 0.6.0, StarPhase 1.0+ (PacBio), HIBAG 1.40+, samtools 1.19+, bwa-mem 0.7.17+. IPD-IMGT/HLA database release frequency is quarterly; tools must be re-bundled with the current release to capture new alleles (~38,000 alleles at Jan 2024; ~43,000+ by Jul 2025).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. Tool reference-bundle vintage matters more than algorithm choice for non-European cohorts; a 2022-bundled HLA-LA will silently miss thousands of post-2022 alleles dominant in African and South Asian ancestry.

# HLA Typing for Clinical Applications

**'Determine HLA genotype for HSCT / neoantigen prediction / PGx screening'** -> Call HLA class I (A, B, C) and class II (DRB1, DRB3/4/5, DQA1, DQB1, DPA1, DPB1) alleles at the resolution required by the downstream application.

- CLI (general-purpose all-rounder): `t1k --preset hla -1 R1.fq -2 R2.fq -f hla_reference.fa`
- CLI (class I gold standard from WES/WGS): `OptiTypePipeline.py -i R1.fq R2.fq -d`
- CLI (class I + II with PRG): `HLA-LA.pl --BAM input.bam --graph PRG_MHC_GRCh38_withIMGT`
- CLI (RNA-seq): `arcasHLA extract sample.bam -o out && arcasHLA genotype out/sample.extracted.fq.gz`
- CLI (long-read transplant-grade): PacBio HiFi StarPhase
- R (imputation from SNP arrays): `HIBAG::predict()` with ancestry-stratified reference panel

## Resolution Levels and What Each Application Requires

HLA nomenclature: **HLA-A\*02:01:01:01** = family : protein-changing : synonymous : intronic/UTR. **Expression suffixes:** `N` (null; DNA present, no protein expressed); `L` (low expression); `S` (secreted); `Q` (questionable); `A` (aberrant). A serologically apparent DR4-positive donor carrying `DRB4*01:03:01:02N` is functionally DR53-negative; a classic HSCT donor-selection failure.

| Application | Min resolution | Why |
|-------------|---------------|-----|
| **HSCT (unrelated donor)** | 6-field (12/12 match) | Null alleles + permissive DPB1 + Bw4/Bw6 + TCE3 core/non-core |
| **Solid organ transplant** | 4-field (2-digit:2-digit) | Eplet-level epitope match (HLAMatchmaker, PIRCHE-II) |
| **ICI neoantigen prediction** | 4-field class I + II | NetMHCpan-4.1 minimum |
| **HLA-disease association** | 4-field | Standard for GWAS HLA fine-mapping |
| **HLA-B\*57:01 abacavir screen** | 4-field, specific | Other \*57 alleles (\*57:03) do NOT cause HSS |
| **HLA-B\*15:02 carbamazepine** | 4-field, specific | \*15:02 only; \*15:01 (NFE-common) is not the risk allele |

## G-Groups vs P-Groups: Routinely Confused

- **G-groups** collapse alleles with identical DNA sequence across the antigen-recognition exons (class I exons 2-3; class II exon 2). Use for sequence-level lab QC.
- **P-groups** collapse alleles encoding identical mature protein across class I positions 1-90 (or class II beta1 domain positions 1-94). Use for epitope-based matching and neoantigen prediction.

## DRB1 + DRB3/4/5 Linkage: The Mandatory Sanity Check

DR haplotype linkage is fixed and is the canonical sanity check on any DR typing:

| DRB1 allele family | Linked DRB3/4/5 |
|---------------------|-----------------|
| DR1 (\*01), DR8 (\*08), DR10 (\*10) | None |
| DR3 (\*03), DR11 (\*11), DR12 (\*12), DR13 (\*13), DR14 (\*14) | DRB3 |
| DR4 (\*04), DR7 (\*07), DR9 (\*09) | DRB4 |
| DR15 (\*15), DR16 (\*16) | DRB5 |

Any caller reporting DRB4 with `DRB1*15:01` is broken or has a chimera. Use this as a routine QC check on automated pipelines.

## Algorithmic Taxonomy: Short-Read Tools

| Tool | Class I | Class II | KIR | Resolution | Approach | Fails when |
|------|---------|----------|-----|-----------|----------|-----------|
| **OptiType** (Szolek 2014 *Bioinformatics* 30:3310) | Yes (~98% 4-digit) | No | No | 4-field | ILP on exons 2-3 | Class II needed; very deep contamination |
| **Polysolver** (Shukla 2015 *Nat Biotechnol* 33:1152) | Yes (~95% 4-digit) | No | No | 4-field | Allele-specific ref alignment | Class II; non-European ancestry under-typing |
| **HLA-LA** (Dilthey 2019 *Bioinformatics* 35:4054) | Yes (~94% class I) | Yes (best class II of WES tools) | No | 4-field | Graph-based PRG | High RAM/disk (~30-100 GB scratch) |
| **T1K** (Song 2023 *Genome Res*) | Yes (~99% 4-digit) | Yes (~99%) | Yes (KIR + KIR3DL2 ligand) | 4-field | EM on consensus reference | Newer; less benchmarking on edge cases |
| **HLA-HD** (Kawaguchi 2017 *Hum Mutat* 38:788) | Yes (~98%) | Yes (~95%) | No | 4-field | Bowtie2 against IPD-IMGT | License required for commercial use |
| **arcasHLA** (Orenbuch 2020 *Bioinformatics* 36:33) | Yes (~100% 2-field) | Yes (>99% 2-field) | No | 4-field from RNA-seq | EM on STAR alignment | DNA-seq; population prior bias in non-EUR |
| **PHLAT, HLAforest, HLAminer, seq2HLA, HLAreporter** | Yes | Some | No | Mostly 2-4 field | Various | Older; superseded |

**Operational benchmark consensus (Claeys 2023 *BMC Genomics*; Matey-Hernandez 2018):** T1K is currently the best general-purpose all-rounder; HLA-LA is the class-II reference; OptiType is the class-I anchor for WES. For full coverage of class I + II + KIR on WGS/WES, T1K is the 2024-2026 recommendation.

## Long-Read and Ultra-High-Resolution

| Tool | Platform | Resolution | Use case |
|------|----------|-----------|----------|
| **StarPhase** (PacBio official 2024+) | PacBio HiFi | 8-field (full-field) | Transplant-grade typing |
| **HLA*ASM** | PacBio HiFi | 8-field | Assembly-based |
| **FuFiHLA** (2025 bioRxiv) | PacBio HiFi + ONT R10 | 8-field | Platform-agnostic |
| **HLAminer streaming** (Warren 2025) | ONT long-read | 4-field | Streaming nanopore |
| **pbaa + StarPhase** | PacBio amplicon | 8-field | Cost-effective targeted typing |
| **IGenotyper** (Roe 2021) | PacBio long-read | 8-field | Immunogenetics-focused |

ONT R9 was historically unreliable for null-allele discrimination due to homopolymer errors; R10.4 with duplex closes the gap for class I and is competitive with PacBio HiFi for class II. PacBio HiFi remains the gold standard for DPB1 4-field typing.

## SNP-Based HLA Imputation: The Ancestry Footgun

When only SNP-array genotypes are available (GWAS cohorts), use imputation:

| Tool | Approach | Reference panel | Best for |
|------|----------|----------------|----------|
| **HIBAG** (Zheng 2014 *Pharmacogenomics J* 14:192) | Random forest from SNP-array | Pre-fit per-ancestry classifiers (EUR, AS, AFR, HIS) | Population-stratified GWAS |
| **HLA-TAPAS** (Luo 2021 *Nat Genet* 53:1504) | Multi-ancestry imputation | 21,546 multi-ancestry reference | Cross-ancestry GWAS |
| **HLA*IMP:02** (Dilthey 2013) | Hidden Markov | EUR-only | Legacy; EUR-only |
| **SNP2HLA** (Jia 2013) | Beagle-based | Type 1 Diabetes / EUR | Older; EUR-only |
| **CookHLA** (Cook 2021) | Hybrid SNP2HLA + supplementary | Multi-ancestry refs | Modern alternative to SNP2HLA |
| **Multi-Ethnic Reference Panel** (Degenhardt 2019) | Multi-ancestry imputation | Cross-population samples | Cross-ancestry GWAS |

**Critical caveat:** imputation panel quality is the limiting factor, NOT the imputation algorithm. EUR-trained HIBAG on East-Asian SNP-array data produces confidently wrong calls. African-ancestry imputation accuracy drops 10-20 percentage points without an ancestry-matched panel (Douillard 2024 *HLA*). For populations underrepresented in IPD-IMGT/HLA itself, imputation is fundamentally limited regardless of method.

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| WGS/WES, class I only, max speed | OptiType | Best class-I accuracy, ILP-based, fast |
| WGS/WES, class I + II, general-purpose | T1K | Best all-rounder; class I + II + KIR co-typing |
| WGS/WES, class II reference grade | HLA-LA | Highest class-II accuracy in benchmarks |
| RNA-seq tumor/normal for ICI | arcasHLA | RNA-seq native; expressed-allele-aware |
| Transplant 6+ field resolution | StarPhase (PacBio HiFi) | 8-field native; reference standard |
| Cost-effective targeted typing | pbaa + StarPhase amplicons | Lower cost than WGS |
| TCGA-style cancer cohort | Polysolver | TCGA convention; reproduces published values |
| SNP array (e.g., UKB) | HIBAG with population-matched panel | No sequencing data |
| Multi-ancestry GWAS | HLA-TAPAS | Cross-ancestry reference |
| Class II DPB1 4-field certainty | StarPhase or HiFi | Pre-2021 WES kits under-cover DPB1 |
| ONT-only data | T1K or HLAminer streaming for class I; ONT R10.4+ duplex for class II | R9 unreliable for nulls |

## HLA and Pharmacogenomics

| HLA allele | Drug | Reaction | Population enrichment | OR |
|------------|------|----------|----------------------|-----|
| **B\*57:01** | Abacavir | Hypersensitivity syndrome | All ancestries (5-8% NFE) | ~100 |
| **B\*15:02** | Carbamazepine, oxcarbazepine | SJS/TEN | Han Chinese, Thai, Malay (>=5%) | ~2500 |
| **B\*58:01** | Allopurinol | SJS/TEN | Han Chinese, Korean, Thai | ~580 |
| **A\*31:01** | Carbamazepine | DRESS, MPE | Europeans, Japanese | ~12 |
| **B\*13:01** | Dapsone | DDS | Han Chinese, SE Asian | -- |
| **B\*35:02** (NOT \*35:01) | Minocycline | DILI | All ancestries | -- |
| **B\*35:01** | TMP-SMX | DILI | Mixed | -- |
| **B\*14:01** | TMP-SMX | DILI | African | -- |
| **A\*33:01/03** | Terbinafine | DILI | Multi-ancestry | -- |
| **DRB1\*15:01 + DQB1\*06:02 haplotype** | Amoxicillin-clavulanate | DILI | Europeans | -- |
| **B\*15:13** | Phenytoin | SJS | Malaysian | -- |

**Operational rule:** Pharmacogenomic HLA screening requires 4-field resolution; 2-field (e.g., "B*15") misses the specific allele.

## Standard Workflow: T1K on WGS/WES

**Goal:** Type HLA class I, class II, KIR from short-read sequencing with KIR3DL1 Bw4/Bw6 ligand prediction.

**Approach:** Extract MHC-region reads, run T1K with IPD-IMGT/HLA reference; T1K outputs allele-pair calls + class II haplotype + KIR.

```bash
# Extract chr6:28-34 Mb plus alt contigs (alt-aware alignment is critical)
samtools view -b -h input.bam chr6:28000000-34000000 chr6_GL000250v2_alt chr6_GL000251v2_alt \
                              chr6_GL000252v2_alt chr6_GL000253v2_alt chr6_GL000254v2_alt \
                              chr6_GL000255v2_alt chr6_GL000256v2_alt > hla_region.bam

samtools sort -n hla_region.bam -o hla_sorted.bam
samtools fastq -1 hla_R1.fq -2 hla_R2.fq -s singletons.fq -0 /dev/null hla_sorted.bam

# Run T1K (preset hla; includes class I + II).
# Some releases ship the entry point as `run-t1k` (a wrapper script) rather than `t1k`;
# verify with `which run-t1k` / `which t1k` before scripting.
t1k --preset hla \
    -1 hla_R1.fq -2 hla_R2.fq \
    -f hla_idx/hlaidx_rna_seq.fa \
    -o sample_hla \
    --threads 8

# Output: sample_hla_genotype.tsv with HLA-A, B, C, DRB1, DRB3/4/5, DQA1, DQB1, DPA1, DPB1
```

## OptiType for Class I (TCGA-Compatible)

**Goal:** Type HLA-A, B, C at 4-field from WES with high accuracy.

**Approach:** Razers3-based alignment to IMGT class-I reference; ILP optimization to assign reads to allele pairs.

```bash
samtools view -h input.bam chr6:28000000-34000000 | samtools fastq -1 R1.fq -2 R2.fq -
OptiTypePipeline.py -i R1.fq R2.fq -d -o optitype_out -c config.ini
```

```ini
# config.ini
[mapping]
razers3=/usr/bin/razers3
threads=8
[ilp]
solver=glpk
threads=8
[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
```

## HLA-LA for Class II (PRG-Based)

**Goal:** Type both class I and class II at 4-field with the highest class-II accuracy of any WES tool.

**Approach:** Population reference graph (PRG) covering the MHC; HLA-LA maps reads to the PRG and infers the most likely paths.

```bash
HLA-LA.pl \
    --BAM input.bam \
    --graph PRG_MHC_GRCh38_withIMGT \
    --workingDir hla_la_out \
    --sampleID sample_name \
    --maxThreads 8

# Output: hla_la_out/sample_name/hla/R1_bestguess_G.txt
# Format: Locus, Allele1, Allele2, AverageCoverage
```

## arcasHLA for RNA-seq

**Goal:** Type HLA class I + II directly from RNA-seq for ICI neoantigen prediction.

**Approach:** Extract HLA-mapped reads from STAR BAM, EM-based genotype call against IMGT.

```bash
# Update reference to current IPD-IMGT/HLA release
arcasHLA reference --update

# Extract and genotype
arcasHLA extract sample.bam -o arcas_out --threads 8
arcasHLA genotype arcas_out/sample.extracted.fq.gz -o arcas_out --threads 8 --population prior

# Output: arcas_out/sample.genotype.json
```

## SNP-Array Imputation (HIBAG): For GWAS Cohorts

**Goal:** Impute HLA from SNP array genotypes when sequencing is unavailable.

**Approach:** HIBAG random-forest classifier with population-matched reference panel.

```r
library(HIBAG)

# Population-matched panel is critical; mismatch causes systematic errors
# Available panels: EUR, ASN, AFR, HIS (download from HIBAG release page)
load('European-HLA4-hg19.RData')

# Load PLINK genotype (.bed/.bim/.fam)
gen <- hlaBED2Geno(bed.fn='cohort.bed', fam.fn='cohort.fam', bim.fn='cohort.bim')

# Predict each locus
hla_A <- predict(model.list[['A']], gen, type='response+prob')
hla_B <- predict(model.list[['B']], gen, type='response+prob')
hla_DRB1 <- predict(model.list[['DRB1']], gen, type='response+prob')

# Filter on probability >= 0.5 for downstream use; lower for exploratory
```

## Per-Operation Failure Modes

**1. Alt-aware alignment missing**
- Trigger: BAM was aligned with bwa-mem against GRCh38 *without* `--alt-aware`; HLA reads are coerced to chr6 primary contigs.
- Mechanism: GRCh38 has ~8 alternate HLA contigs (chr6_GL000250v2_alt, etc.); without alt-aware alignment, reads from these regions get assigned to suboptimal positions on the primary chr6.
- Symptom: HLA typing accuracy drops 5-10 percentage points; high read-coverage variants get miscalled.
- Fix: Re-align the HLA region with bwa-mem-alt or use the original cDNA reference for HLA typing (extract reads to FASTQ first).

**2. Stale IPD-IMGT/HLA bundle**
- Trigger: Tool was installed in 2022 with the corresponding IPD-IMGT/HLA release; never updated.
- Mechanism: ~5000+ new alleles added between 2022 and 2025; new alleles dominant in under-represented ancestries.
- Symptom: Non-European samples get common alleles reported as ambiguous or as the closest legacy match.
- Fix: Update the tool's reference bundle (HLA-LA: rebuild PRG; T1K: re-run `t1k-build`; OptiType: update `data/hla_reference_dna.fasta`).

**3. EUR-trained imputation on non-EUR samples**
- Trigger: Use HIBAG European panel on East-Asian or African ancestry samples.
- Mechanism: Random forest was trained on EUR allele frequencies; non-EUR alleles missing from training set.
- Symptom: Confidently wrong calls; high probability assigned to incorrect alleles.
- Fix: Use ancestry-matched HIBAG panel; or switch to HLA-TAPAS multi-ancestry; or fall back to sequencing.

**4. Cross-mapping DRB-related loci**
- Trigger: Naive bwa-mem alignment without read-grouping at DRB1/DRB3/DRB4/DRB5.
- Mechanism: DRB1, DRB3, DRB4, DRB5 share extensive sequence identity; reads map ambiguously.
- Symptom: DR3/DR4/DR5 paralog reads contaminate DRB1 calls; haplotype linkage rule (e.g., DRB1*15:01 + DRB5) violated.
- Fix: Use HLA-LA or T1K which model paralogous loci jointly; verify DRB1+DRB3/4/5 haplotype rule.

**5. DPB1 under-coverage in pre-2021 WES kits**
- Trigger: Used SureSelect v5 or Nextera Rapid Capture WES; DPB1 reports homozygous typing.
- Mechanism: Pre-2021 capture kits under-covered DPB1 exon 2.
- Symptom: Heterozygous DPB1 reported as homozygous; affects HSCT matching.
- Fix: Confirm capture coverage at DPB1; if insufficient, supplement with targeted amplicon or use WGS/long-read.

**6. Class II expression-allele confusion**
- Trigger: Report `DRB4*01:03:01:02N` as functional DR53.
- Mechanism: N-suffix = null allele (DNA present but no protein expressed).
- Symptom: Functionally DR53-negative donor reported as DR53-positive; transplant matching failure.
- Fix: Parse 4-field suffix (`N`, `L`, `S`, `Q`, `A`); treat N as null in functional analysis; preserve full nomenclature for typing report.

**7. Specific allele vs allele family confusion**
- Trigger: PGx screen reports "B*57" carrier as abacavir-risk-positive.
- Mechanism: HLA-B*57 family includes \*57:01 (abacavir HSS risk), \*57:02, \*57:03 (no HSS risk).
- Symptom: False-positive abacavir contraindication; patient denied effective therapy.
- Fix: Report at 4-field minimum; B\*57:01 specifically, not B\*57.

**8. KIR co-typing mistaken for HLA**
- Trigger: Report KIR allele as HLA.
- Mechanism: KIR (chromosome 19) and HLA (chromosome 6) are functionally paired (KIR3DL1 binds HLA-Bw4) but are distinct loci.
- Symptom: Wrong locus annotation; downstream tools fail.
- Fix: Use T1K which co-types HLA + KIR + KIR3DL2 ligand and labels output correctly.

## Reconciliation: When Tools Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| OptiType vs HLA-LA class I disagree | Stale reference bundle in one; non-EUR ancestry | Update both; rerun; prefer the one with current reference |
| HLA-LA vs T1K class II disagree | DRB1+DRB3/4/5 haplotype rule violated in one | Check haplotype linkage; the consistent caller is correct |
| HIBAG vs sequencing disagree | EUR-trained model on non-EUR sample | Trust sequencing; use ancestry-matched HIBAG panel |
| Tumor vs normal HLA differ | Tumor LOH at HLA locus (frequent in NSCLC, HNSCC) | Run LOHHLA / DASH to confirm somatic loss; report germline + somatic |
| DPB1 homozygous on WES, het on WGS | WES kit under-covers DPB1 exon 2 | Trust WGS; flag WES result as low confidence |
| Class I 4-field stable across tools, class II differs | Class II is fundamentally harder | Prefer HLA-LA or StarPhase for class II |
| arcasHLA vs OptiType for tumor RNA | arcasHLA returns expressed-allele only (may miss silenced allele due to LOH) | Confirm with DNA-based typing for transplant context |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| IPD-IMGT/HLA quarterly release | Updates Jan/Apr/Jul/Oct | IPD-IMGT/HLA database |
| Current allele count | ~43,000+ at Jul 2025 | IPD-IMGT/HLA database release notes (Robinson J et al, *NAR* DB issue) |
| HLA region coordinates | chr6:28000000-34000000 (GRCh38) | Standard |
| HLA-LA RAM requirement | ~30-100 GB scratch | HLA-LA documentation |
| OptiType class I 4-digit accuracy | ~98% (1000G benchmark) | Claeys 2023 |
| Polysolver class I 4-digit accuracy | ~95% | Matey-Hernandez 2018 |
| HLA-LA class II accuracy | Best of WES tools | Claeys 2023 |
| T1K class I + II accuracy | ~99% / ~99% | Song 2023 |
| HIBAG probability cutoff | >=0.5 for clinical-grade; >=0.3 for exploratory | HIBAG documentation |
| 1000G allele coverage | ~60-70% of African-ancestry alleles still under-represented in IPD-IMGT/HLA | Robinson 2024 |
| HSCT matching standard | 10/10 or 12/12 at 6-field | NMDP/WMDA guidelines |
| TCE3 core alleles | DPB1\*02:01, \*04:01, \*04:02, \*23:01 | Meurer 2024 *Blood* 144:1659 |

## CIWD v3.0.0 Ambiguity Catalogue

Hurley 2020 *HLA* 95:516; compiled from >8M unrelated HSCT donors across 7 geographic/ancestral groups. Categories: Common (18%, n=545), Intermediate (17%, n=513), Well-Documented (65%, n=1,997) at 2-field. Replaces legacy CWD 2.0 (Mack 2013); many older pipelines still hardcode CWD 2.0; a quiet quality failure.

## TCE3 Core vs Non-Core (Meurer 2024 *Blood*)

DPB1 mismatch GvHD/relapse risk depends on TCE3 group:
- **Core (DPB1\*02:01, \*04:01, \*04:02, \*23:01):** GvHD reduction with permissive mismatch in the GvH direction.
- **Non-core:** Relapse-protection effects predominate.

Now operational in NMDP donor selection algorithms; legacy TCE3 frameworks (Crocchiolo 2009) lack this stratification.

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| HLA-DRA in output (DRB1 expected) | Tool confused paralogs | Use HLA-LA or T1K which model paralog loci correctly |
| Class II reports "no call" | Pre-2021 WES kit under-covers class II | Switch to WGS or amplicon |
| Tumor and normal HLA differ | LOH at HLA locus | Confirm with LOHHLA; report germline call as ground truth |
| Imputation reports rare allele with high probability | Reference panel mismatch with cohort ancestry | Switch to ancestry-matched panel |
| 4-field call but only 2-field appears in report | Tool default truncation | Use `--full-field` or equivalent flag |
| Same sample gives different 4-field calls across runs | Stochastic tie-breaking | Pin random seed; report all equally-supported calls |
| DRB4 with DRB1*15 | Linkage rule violated; bug or chimera | Re-run; check for sample swap |
| Null allele not reported in summary | Tool drops N-suffix; output is misleading | Use raw 4-field output; never strip suffixes for clinical reports |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why T1K when HLA-LA is the published reference?" | T1K matches HLA-LA accuracy on class II while also typing class I + KIR in one pass with lower RAM; we cite both. |
| "Your African-ancestry samples have low confidence" | IPD-IMGT/HLA still under-represents African ancestry (~30-40% allele gap); we ran with current 2025 release; for transplant we recommend long-read confirmation. |
| "DRB1 vs DRB3/4/5 reported inconsistently" | We verified DRB1+DRB3/4/5 linkage rule on each sample as routine QC; flagged violations for re-typing. |
| "Why is HLA-B\*15:01 not flagged for carbamazepine?" | \*15:01 (NFE common) is not the SJS risk allele; \*15:02 (Han Chinese) is. PGx requires 4-field specificity. |
| "Imputation results differ from sequencing" | Imputation panel quality is the limiting factor; EUR-trained HIBAG on non-EUR is unreliable; we used ancestry-matched panel. |
| "TCGA pipeline used Polysolver, why T1K?" | TCGA convention is Polysolver; for *current* analysis we use T1K which has better class-II and KIR coverage. We can reproduce Polysolver if back-comparison needed. |

## References

- Robinson J et al. 2024. IPD-IMGT/HLA database: 25 years of evolution. *HLA*.
- Robinson J et al. 2026. IPD-IMGT/HLA recent developments. *Nucleic Acids Res* 54:D1152.
- Szolek A et al. 2014. OptiType: precision HLA typing from NGS data. *Bioinformatics* 30:3310.
- Dilthey AT et al. 2019. HLA*LA; HLA typing from linearly projected graph alignments. *Bioinformatics* 35:4054.
- Song B et al. 2023. T1K: efficient and accurate KIR and HLA genotyping. *Genome Res*.
- Shukla SA et al. 2015. Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. *Nat Biotechnol* 33:1152.
- Kawaguchi S et al. 2017. HLA-HD: An accurate HLA typing algorithm for next-generation sequencing data. *Hum Mutat* 38:788.
- Orenbuch R et al. 2020. arcasHLA: high-resolution HLA typing from RNAseq. *Bioinformatics* 36:33.
- Claeys A et al. 2023. Benchmark of tools for in silico prediction of MHC class I and class II genotypes from NGS data. *BMC Genomics* 24:9351.
- Matey-Hernandez ML et al. 2018. Benchmarking the HLA typing performance of Polysolver and Optitype in 50 Danish parental trios. *BMC Bioinformatics* 19:239.
- Zheng X et al. 2014. HIBAG; HLA genotype imputation with attribute bagging. *Pharmacogenomics J* 14:192.
- Luo Y et al. 2021. A high-resolution HLA reference panel capturing global population diversity. *Nat Genet* 53:1504.
- Hurley CK et al. 2020. Common, intermediate and well-documented HLA alleles in world populations: CIWD version 3.0.0. *HLA* 95:516.
- Meurer T et al. 2024. New HLA-DPB1 T-cell epitope model for mismatched URD transplant. *Blood* 144:1659.
- Douillard V et al. 2024. Population-specific HLA imputation reference panels. *HLA*.

## Related Skills

- clinical-databases/pharmacogenomics - HLA-drug interactions, abacavir/carbamazepine screening
- immunoinformatics/mhc-binding-prediction - Downstream HLA-peptide binding for neoantigen
- workflows/neoantigen-pipeline - HLA typing as upstream step
- clinical-databases/clinvar-lookup - HLA disease associations
- population-genetics/population-structure - Ancestry-aware imputation context
