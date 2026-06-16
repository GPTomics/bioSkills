---
name: bio-phasing-imputation-genotype-imputation
description: Imputes untyped genotypes against a phased reference panel with Beagle, Minimac4, or IMPUTE5 (array data) or from genotype likelihoods with GLIMPSE2, QUILT2, or STITCH (low-coverage WGS), producing per-variant dosages (DS) with a self-estimated quality (Beagle DR2, Minimac R2, IMPUTE INFO). Covers why the honest output is a dosage posterior not a hard call, why GWAS regresses on DS, why the quality metric is an ESTIMATE of r2 from posterior spread (not validation against truth), the DS/GP/HDS fields, the phasing prerequisite, chunking, chrX ploidy, the Michigan/TOPMed servers (the only access to HRC/TOPMed), and low-coverage WGS as the modern array replacement. Use when increasing variant density for GWAS, harmonizing arrays, inferring untyped variants, or imputing low-coverage sequence. Phase first with haplotype-phasing; prepare the panel with reference-panels; filter with imputation-qc; the GWAS test is population-genetics/association-testing; end-to-end orchestration is workflows/gwas-pipeline.
tool_type: cli
primary_tool: Beagle
---

## Version Compatibility

Reference examples tested with: Beagle 5.4 (22Jul22), Minimac4 4.1+, IMPUTE5 1.2, GLIMPSE2, bcftools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Minimac4 (4.x) uses POSITIONAL arguments (`minimac4 panel.msav target.vcf.gz`); the old `--refHaps`/`--haps`/`--prefix`/`--cpus` style is Minimac3 and obsolete. Beagle 5.x emits `DR2`, `AF`, and `IMP` only (`AR2` is a legacy 4.x field) and its default `ne=100000` (not 1,000,000). The panel build (GRCh37 vs GRCh38) must match the data; record the panel name, version, and build with every result.

# Genotype Imputation -- Inferring Untyped Genotypes as Dosages

**"Fill in the variants I did not directly measure"** -> Align the (phased or low-coverage) sample to a reference panel of phased haplotypes and infer the untyped alleles via the Li-Stephens HMM - because the output is a posterior over genotypes summarized as a dosage with a self-estimated quality, not a measured call, so the uncertainty must be carried downstream.
- CLI: `java -jar beagle.jar gt=phased.vcf.gz ref=panel.bref3 map=plink.chr20.map out=imputed` (or `minimac4 panel.msav phased.vcf.gz`, or GLIMPSE2 for low-coverage WGS)

Scope: imputing untyped genotypes from a panel (array data) or from genotype likelihoods (low-coverage WGS), the dosage/quality output, chunking, chrX, and the servers. Phasing the input -> haplotype-phasing. Panel selection/preparation/strand -> reference-panels. Quality metrics and filtering thresholds -> imputation-qc. The GWAS test on the dosages -> population-genetics/association-testing. The genotype likelihoods that low-coverage imputation consumes -> variant-calling/vcf-basics. End-to-end orchestration -> workflows/gwas-pipeline.

## The Single Most Important Modern Insight -- An Imputed Genotype Is a Posterior, and the Deliverable Is a Dosage Plus a Self-Estimated Quality, Not a Hard Call

Imputation aligns a sparsely-genotyped (or low-coverage-sequenced) sample to a densely-typed reference panel of phased haplotypes and infers, via a Li-Stephens HMM, the alleles at positions the sample never observed (Browning 2018 *Am J Hum Genet* 103:338). The output at each untyped variant is a distribution, summarized as an expected allelic dosage in [0,2]. Three facts define the field:

1. **Downstream analysis uses dosages, not hard genotypes.** The dosage DS is the conditional expectation E[genotype | data, panel], the minimum-variance summary; hard-calling forces an uncertain 0.5 dosage to 0 or 1, injecting genotype error that attenuates effects and inflates standard errors. GWAS regresses the trait on DS -> population-genetics/association-testing.
2. **The quality metric (Beagle DR2, Minimac R2, IMPUTE INFO) is an ESTIMATE of r2 from the posterior spread, computed without ever seeing the truth.** Poorly-imputed dosages shrink toward the allele-frequency mean 2p, so low posterior variance relative to the binomial expectation 2p(1-p) flags a low-confidence site. This is NOT a validation against held-out genotypes (that is empirical r2 / EmpRsq, a masked-site quantity). Say "DR2/R2/INFO is an estimate of imputation quality," never "the imputation accuracy was 0.9" as if measured. The metric also cannot detect panel-ancestry mismatch -> imputation-qc.
3. **Low-coverage WGS (0.5-4x) plus GLIMPSE2 has become a credible array replacement.** Because it samples the whole genome rather than a fixed ascertained SNP set, it imputes rare variants and under-represented ancestries better than a dense array at comparable cost (Rubinacci 2023 *Nat Genet* 55:1088). The input is genotype likelihoods, not calls; the array-vs-low-coverage-WGS choice is an ascertainment decision (see Array vs Low-Coverage WGS below).

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| Minimac4 | Das 2016 *Nat Genet* 48:1284 | array imputation; msav/m3vcf panel; the server engine; positional-arg CLI | server-style imputation; meta-imputation |
| Beagle 5.x | Browning 2018 *Am J Hum Genet* 103:338 | Java; phases unphased input AND imputes; bref3 panel | one tool for phase + impute, no compile |
| IMPUTE5 | Rubinacci 2020 *PLoS Genet* 16:e1009049 | PBWT pre-selection then LS HMM; sub-linear in panel size | very large reference panels; local speed |
| GLIMPSE2 | Rubinacci 2023 *Nat Genet* 55:1088 | low-coverage WGS imputation from genotype likelihoods; chunk/split/phase/ligate | 0.5-4x WGS with a panel |
| QUILT2 | Davies 2021 *Nat Genet* 53:1104 | low-coverage, panel-based, read-aware | long-read / haplotagged / ancient DNA / cfDNA |
| STITCH | Davies 2016 *Nat Genet* 48:965 | low-coverage, REFERENCE-FREE; learns ancestral haplotypes by EM | no panel exists (non-model organisms) |
| Michigan / TOPMed servers | Das 2016 *Nat Genet* 48:1284 | Eagle2 phasing + Minimac4; the only access to HRC/TOPMed | turnkey, access-controlled panels |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Array data, want HRC/TOPMed and a turnkey pipeline | TOPMed or Michigan Imputation Server | the only sanctioned access to those panels; runs Eagle2 + Minimac4 |
| Array data, local run, very large panel, want speed | IMPUTE5 (PBWT) or Minimac4 | sub-linear scaling in panel size |
| Array data, local, one tool for phase + impute | Beagle 5.x | phases unphased gt= input itself; bref3 panel |
| Low-coverage WGS (0.5-4x), have a panel | GLIMPSE2 (chunk -> split-reference -> phase -> ligate) | the standard; imputes from genotype likelihoods |
| Low-coverage, read-aware / long-read / ancient DNA / cfDNA | QUILT2 | per-read, base-quality-aware |
| Low-coverage, NO reference panel (non-model organism) | STITCH | learns ancestral haplotypes reference-free |
| Need the panel selected/prepared first | -> reference-panels | the panel is the prior |
| Need the input phased first (Minimac4, IMPUTE5) | -> haplotype-phasing | those engines require a phased target |
| Filter the imputed output before analysis | -> imputation-qc | DR2/R2/INFO + MAF floor |
| The GWAS test on the dosages | -> population-genetics/association-testing | downstream |

## Array vs Low-Coverage WGS: the Imputation-Input Fork

The upstream decision is how to generate the genotypes that will be imputed, and it is an ascertainment question, not just an accuracy one. An array assays a fixed, designed SNP set (biased to its design population); low-coverage WGS samples whatever is in the genome.

| | SNP array + pre-phase + impute | Low-coverage WGS (~0.5-4x) + impute from genotype likelihoods |
|---|---|---|
| Input to the HMM | hard genotype calls (array error is tiny) | genotype LIKELIHOODS (PL/GL); a hard call at 1x is mostly noise |
| Ascertainment | FIXED - only the designed SNPs, biased to the design population | UNBIASED - whatever is in the genome is observed |
| Rare variants | limited by the array scaffold and panel | matches or beats dense arrays (Rubinacci 2021 *Nat Genet* 53:120) |
| Under-represented ancestry | poor (no good array, panel-mismatched) | the main route around array/panel bias |
| Tools | Beagle / Minimac4 / IMPUTE5 | GLIMPSE2 (panel) / STITCH (no panel) |

The judgment: common-variant GWAS in a well-paneled ancestry -> array plus imputation is cheap and adequate; rare variants, under-represented ancestry, or a need for unbiased genome-wide ascertainment -> low-coverage WGS plus genotype-likelihood imputation, the direction the field is moving as sequencing costs fall. Low-coverage WGS is only as good as its panel and its likelihoods (bad mapping, contamination, or damage produce garbage GLs that impute garbage).

## Output Formats and Why Dosages

The central object is the posterior genotype distribution; everything else summarizes it. Request the fields up front (Minimac4 `-f GT,DS,HDS,GP`; Beagle `gp=true ap=true`).

| FORMAT | Meaning | Shape |
|--------|---------|-------|
| GP | genotype probabilities P(0/0),P(0/1),P(1/1); the full posterior | 3 values summing to 1 |
| DS | allelic dosage = P(0/1) + 2*P(1/1) = E[genotype]; the GWAS field | 1 value in [0,2] |
| HDS | haploid (phased per-haplotype) dosage; DS = HDS1 + HDS2 (Minimac4/GLIMPSE) | 2 values, each [0,1] |
| AP1/AP2 | Beagle allele probabilities (P(ALT) per haplotype); DS = AP1 + AP2 (with ap=true) | 1 value each [0,1] |
| GT | hard best-guess genotype (argmax); lossy, discards uncertainty | 0/0, 0/1, 1/1 |

GP is the distribution; DS is its mean - two variants with different GP spreads can share a DS. Use DS for association (it propagates the uncertainty); use HDS/AP for phased/allele-specific analyses. Beagle computes GP from allele probabilities assuming Hardy-Weinberg and sets GT from the per-haplotype argmax, so its GT can occasionally disagree with the argmax of its own GP.

## The Phasing Prerequisite

The reference panel is phased haplotypes; the target must align to that haplotype structure two ways:
- **Pre-phase then impute** (Minimac4, IMPUTE5): phase the target FIRST (Eagle2 or SHAPEIT) into haplotypes, then impute. The server default (Eagle2 -> Minimac4) and the fast local pattern -> haplotype-phasing.
- **Phase-and-impute together** (Beagle, GLIMPSE2): the tool phases internally. Low-coverage tools MUST do this, because there is no confident genotype to phase up front; GLIMPSE2 alternates haploid imputation and phasing, and gains accuracy by imputing all target samples jointly.

## Low-Coverage WGS Workflow (GLIMPSE2)

The input is genotype likelihoods (PL/GL), not calls, because at 0.5-4x no genotype is certain. GLIMPSE2 can read BAM/CRAM directly (computing GLs internally) or a GL BCF made with `bcftools mpileup ... -T panel_sites.vcf.gz | bcftools call -Aim -C alleles -T panel_sites.tsv.gz` (the `-C alleles` constraint needs the panel sites supplied to `call` via `-T`; note the two `-T` files differ in format - a VCF for mpileup, a tab-delimited sites file for call). The pipeline:

1. `GLIMPSE2_chunk` defines windows with buffers.
2. `GLIMPSE2_split_reference` precomputes a binary panel per chunk (the speed innovation that made UK Biobank-scale imputation feasible).
3. `GLIMPSE2_phase` imputes and phases per chunk (`--bam-list` or `--input-gl`; `--ne` default 100000).
4. `GLIMPSE2_ligate` stitches chunks using the overlap buffers to keep phase. Output FORMAT: GT, DS, GP, HS plus a per-variant INFO score.

For chrX with GLIMPSE2, declare each sample's ploidy with `--samples-file` (sample and copy number) and run the PAR/nonPAR split as for the array tools (male nonPAR is haploid) -> reference-panels.

## Imputation Servers

The Michigan (now MIS2) and TOPMed servers run Eagle2 phasing + Minimac4 imputation server-side and are the ONLY sanctioned access to HRC and TOPMed (those panels are controlled-access, not downloadable). Upload a per-chromosome VCF, select the panel, build, and population; the server runs allele-frequency QC and strand-flip detection, phases, imputes in chunks, and returns per-chromosome VCFs in GT,DS,GP plus a Minimac info file with R2 and a QC report. Results are encrypted with a one-time password and auto-deleted after a few days. The reproducibility cost: the panel version (HRC r1.1 vs TOPMed r2 vs r3), tool version, and build can change between runs, so record exactly which server/panel/version produced a result.

## Per-Method Failure Modes

### Obsolete Minimac4 syntax
**Trigger:** `minimac4 --refHaps panel.m3vcf --haps study.vcf --prefix out`. **Mechanism:** that is Minimac3; Minimac4 4.x takes positional args. **Symptom:** the command errors or is not recognized. **Fix:** `minimac4 panel.msav target.phased.vcf.gz -o imputed.vcf.gz -f GT,DS,HDS,GP -t 8`; build the panel with `minimac4 --compress-reference`.

### Imputing unphased input to a pre-phase engine
**Trigger:** feeding unphased genotypes to Minimac4 or IMPUTE5. **Mechanism:** those engines assume a phased target aligned to the panel haplotypes. **Symptom:** garbage or refused input. **Fix:** phase first (Eagle2/SHAPEIT) -> haplotype-phasing, or use Beagle/GLIMPSE2 which phase internally.

### Imputing cases and controls separately
**Trigger:** running imputation per batch (cases, then controls, or per cohort). **Mechanism:** batch-differential imputation quality at a variant creates artifactual genotype structure correlated with phenotype. **Symptom:** genome-wide-significant hits that fail to replicate; every single-batch QC metric passes. **Fix:** impute all samples together (or harmonize panels/versions and check that quality does not differ by batch) -> imputation-qc.

### Hard-calling the dosage
**Trigger:** thresholding DS to 0/1/2 for association. **Mechanism:** discards the posterior uncertainty, worst at low-R2 rare variants. **Symptom:** lost power read as a true null. **Fix:** regress on DS (PLINK2 `dosage=DS`, SNPTEST, REGENIE, BOLT-LMM all accept dosages).

### Missing DS field
**Trigger:** a downstream tool cannot find dosages. **Mechanism:** the FORMAT fields were not requested. **Symptom:** only GT or GP present. **Fix:** request `-f GT,DS,HDS,GP` (Minimac4) or `gp=true ap=true` (Beagle) at run time.

### Genome build or strand not aligned to the panel
**Trigger:** GRCh37 data against a GRCh38 panel, or unflipped palindromic SNPs. **Mechanism:** positions/alleles disagree with the panel; the HMM copies wrong templates. **Symptom:** near-zero accuracy across regions, no error. **Fix:** align build and strand before imputing -> reference-panels.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Regress on DS (dosage), not hard GT | Browning 2018 *Am J Hum Genet* 103:338 | DS = E[genotype | data, panel] is the minimum-variance estimator; hard-calling injects error |
| Beagle `ne=100000` (default) | Beagle 5.x default | effective population size for the HMM; not 1,000,000 |
| Beagle `window=40.0` / `overlap=2.0` cM | Beagle 5.x defaults | window must be >= 1.1x overlap; rarely tuned |
| Impute all samples together | Browning 2018 *Am J Hum Genet* 103:338 (framing) | separate case/control imputation manufactures false associations -> imputation-qc |
| Low-coverage sweet spot ~0.5-4x | Rubinacci 2023 *Nat Genet* 55:1088 | GLIMPSE2 accuracy range; ~1x is array-competitive |
| Request DS explicitly (Minimac4 default is GT,DS) | Minimac4 docs | HDS/GP for phased/probabilistic uses must be named |
| Post-imputation R2/DR2/INFO filter (a QC decision, not a default) | -> imputation-qc | the imputer's number is the INPUT to filtering, not a tool default |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `minimac4 --refHaps` not recognized | Minimac3 syntax | use positional args: `minimac4 panel.msav target.vcf.gz -o out` |
| Beagle OutOfMemoryError | JVM heap too small / whole genome one job | raise `-Xmx`; impute per chromosome |
| No DS in output | fields not requested | `-f GT,DS,HDS,GP` (Minimac4) / `gp=true ap=true` (Beagle) |
| Imputation accuracy near zero across a region | build/strand mismatch to the panel | align build and strand first -> reference-panels |
| Hits do not replicate | cases/controls imputed separately, or hard-called | impute together; regress on dosages -> imputation-qc |
| Engine errors on multiallelic sites | non-biallelic input | `bcftools norm -m -any` first -> variant-calling/variant-normalization |
| Cannot download HRC/TOPMed | controlled-access panels | use the imputation server |

## References

- Das S, Forer L, Schonherr S, et al. 2016. Next-generation genotype imputation service and methods. *Nat Genet* 48:1284-1287.
- Browning BL, Zhou Y, Browning SR. 2018. A one-penny imputed genome from next-generation reference panels. *Am J Hum Genet* 103:338-348.
- Browning BL, Tian X, Zhou Y, Browning SR. 2021. Fast two-stage phasing of large-scale sequence data. *Am J Hum Genet* 108:1880-1890.
- Rubinacci S, Delaneau O, Marchini J. 2020. Genotype imputation using the Positional Burrows-Wheeler Transform. *PLoS Genet* 16:e1009049.
- Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O. 2021. Efficient phasing and imputation of low-coverage sequencing data using large reference panels. *Nat Genet* 53:120-126.
- Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. 2023. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. *Nat Genet* 55:1088-1090.
- Davies RW, Kucka M, Su D, et al. 2021. Rapid genotype imputation from sequence with reference panels. *Nat Genet* 53:1104-1111.
- Davies RW, Flint J, Myers S, Mott R. 2016. Rapid genotype imputation from sequence without reference panels. *Nat Genet* 48:965-969.
- McCarthy S, Das S, Kretzschmar W, et al. 2016. A reference panel of 64,976 haplotypes for genotype imputation. *Nat Genet* 48:1279-1283.
- Taliun D, Harris DN, Kessler MD, et al. 2021. Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program. *Nature* 590:290-299.

## Related Skills

- haplotype-phasing - Pre-phasing the target (required by Minimac4 and IMPUTE5)
- reference-panels - Select and prepare the panel (the prior) and align build/strand
- imputation-qc - Filter by DR2/R2/INFO and MAF; the metric is an estimate, not truth
- variant-calling/vcf-basics - Genotype likelihoods (PL/GL) for low-coverage imputation
- variant-calling/variant-normalization - Split multiallelics before imputation
- population-genetics/association-testing - GWAS test on the imputed dosages
- clinical-databases/polygenic-risk - Polygenic scores from imputed dosages
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate
