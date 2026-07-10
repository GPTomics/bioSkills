---
name: bio-workflows-somatic-variant-pipeline
description: Chains a somatic (tumor-normal) SNV/indel and structural-variant pipeline end to end with GATK Mutect2 (or Strelka2), wiring the somatic-specific machinery - panel-of-normals and gnomAD germline-resource priors, GetPileupSummaries/CalculateContamination, and LearnReadOrientationModel FFPE/oxoG orientation-bias filtering fed into FilterMutectCalls. Use when calling somatic mutations from a tumor-normal pair (or tumor-only with PoN caveats), deciding which artifact filter removes which class of false positive, reasoning about VAF/purity/ploidy and clonal-vs-subclonal detection, adding somatic SV/CNV or TMB/MSI/signatures, or routing variants to AMP/ASCO/CAP tier and oncogenicity interpretation (never germline ACMG).
tool_type: cli
primary_tool: GATK Mutect2
workflow: true
depends_on:
  - read-alignment/bwa-alignment
  - variant-calling/gatk-variant-calling
  - variant-calling/filtering-best-practices
  - variant-calling/structural-variant-calling
  - variant-calling/variant-annotation
  - variant-calling/clinical-interpretation
  - copy-number/cnvkit-analysis
qc_checkpoints:
  - after_alignment: "Tumor + normal mapping rate >95%, tumor coverage adequate for the target VAF"
  - after_contamination: "CalculateContamination estimate low (<~0.02); high contamination inflates false positives"
  - after_filtering: "FilterMutectCalls PASS fraction sane; FFPE/oxoG orientation-bias artifacts removed via --ob-priors"
  - after_interpretation: "Variants tiered by AMP/ASCO/CAP + oncogenicity (never germline ACMG); drivers vs passengers separated"
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, Strelka2 2.9+, Manta 1.6+, Ensembl VEP 111+, bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: this is a WORKFLOW skill - it wires the somatic-specific chain and the DECISIONS between steps. Component mechanism (caller internals, filter thresholds, annotation, CNV) lives in the cross-referenced component skills; interpretation lives in variant-calling/clinical-interpretation and uses the AMP/ASCO/CAP tier system, NOT germline ACMG.

# Somatic Variant Pipeline

**"Call somatic mutations from my tumor-normal pair"** -> Orchestrate somatic SNV/indel calling (Mutect2 or Strelka2) with panel-of-normals + germline-resource priors, contamination and orientation-bias filtering, then somatic SV/CNV, then tier/oncogenicity interpretation.
- CLI: `gatk Mutect2` (+ FilterMutectCalls chain), `configureStrelkaSomaticWorkflow.py` (Strelka2), `configManta.py` (Manta SV), VEP/Funcotator for annotation.

## The governing principle

A somatic callset is not a fact about the tumor; it is a joint property of (tumor, matched normal, caller, filters, reference). Three consequences drive every decision in this pipeline:

1. **There is no universal somatic truth set.** The ICGC-TCGA DREAM Mutation Calling challenge (Alioto 2015 *Nat Commun* 6:10001) gave expert teams the SAME tumor-normal WGS and got widely varying call rates and low SNV concordance - indels and SVs far worse. Somatic reproducibility is intrinsically worse than germline (low VAF, subclonality, purity/ploidy, normal contamination) and there is still no GIAB-grade generalizable somatic reference. Practical corollary: pin the pipeline, benchmark against orthogonal validation for the assay, and prefer multi-caller consensus where feasible - do not treat one caller's VCF as ground truth.
2. **Somatic variants are sub-1.0 VAF and depth/purity-limited.** Unlike a germline 0/0.5/1.0 genotype, a somatic variant sits at a continuous allele fraction set by tumor purity, local copy number, and clonal fraction. Detecting a low-VAF subclonal or low-purity variant is a depth problem (the somatic regime is the one place more coverage genuinely helps). This is why Mutect2 uses a somatic likelihood model, not a diploid genotyper.
3. **Somatic interpretation uses a DIFFERENT framework.** Never apply germline ACMG (PVS1/PM2/PP3) to a tumor variant - it is a category error. Somatic variants are classified by clinical actionability (AMP/ASCO/CAP tiers, Li 2017) and oncogenicity (ClinGen/CGC/VICC, Horak 2022), and the tier is tumor-type-specific. See the interpretation section and variant-calling/clinical-interpretation.

## Pipeline map

```
Tumor BAM + matched Normal BAM   (aligned, dedup'd, BQSR - see read-alignment/*)
    |
    ├── SNV/indel calling
    │     Mutect2 (GATK)    - somatic likelihood, tumor+normal in one command
    │     Strelka2          - faster; pair with Manta candidateSmallIndels
    │
    ├── Somatic-specific filtering  (the four artifact/germline removers, below)
    │     PoN + germline-resource -> Mutect2 call
    │     LearnReadOrientationModel  (FFPE/oxoG)
    │     GetPileupSummaries -> CalculateContamination
    │     FilterMutectCalls (consumes all three) -> PASS somatic VCF
    │
    ├── Structural variants     -> variant-calling/structural-variant-calling (Manta somatic mode)
    ├── Copy number + purity/ploidy -> copy-number/cnvkit-analysis (purity feeds VAF reasoning)
    ├── Annotation (normalize FIRST) -> variant-calling/variant-annotation (VEP/Funcotator)
    │
    └── Interpretation -> variant-calling/clinical-interpretation
          AMP/ASCO/CAP Tier I-IV + oncogenicity (Horak) + TMB/MSI/signatures
```

## The four somatic-specific filters - what each removes

The machinery that separates somatic calling from germline is four independent artifact/germline removers. Knowing WHICH false-positive class each addresses is the core decision - do not treat them as interchangeable boilerplate.

| Filter | Removes | Built from | When it is critical |
|--------|---------|------------|---------------------|
| Panel of Normals (PoN) | Recurrent technical/site artifacts + common germline that reproduce across normals | 40+ unrelated normals from the SAME assay/platform, `CreateSomaticPanelOfNormals` | Always; the ONLY artifact defense in tumor-only mode |
| Germline resource (gnomAD AF-only) | Germline variants, via a population-AF prior in the somatic model | `af-only-gnomad.vcf.gz` (AF field only) | Always; carries the germline burden alone in tumor-only mode |
| Contamination estimate | Cross-individual sample contamination masquerading as low-VAF somatic | `GetPileupSummaries` on common biallelic SNPs -> `CalculateContamination` | Any sample with suspected swap/contamination; low-VAF calls |
| Orientation-bias model | FFPE deamination (C>T/G>A) and oxoG (C>A/G>T) library artifacts, detected as F1R2/F2R1 strand imbalance | `--f1r2-tar-gz` from Mutect2 -> `LearnReadOrientationModel` | FFPE, archival, or oxidatively damaged input; ALWAYS for FFPE |

PoN and germline-resource act at CALL time (priors passed to Mutect2); contamination and orientation-bias are learned separately and injected at FILTER time (`FilterMutectCalls`). The tumor-only trap: without a matched normal, the PoN and germline resource are the ONLY things removing germline and artifacts, so both must be assay-matched and current, and the false-positive rate is materially higher.

## Mutect2 tumor-normal workflow

The end-to-end chained script is in `examples/run_mutect2.sh`; the steps and their decisions:

### Step 1: Build a Panel of Normals (do once per assay)

```bash
# Each normal called in tumor-only mode; --max-mnp-distance 0 is REQUIRED for GenomicsDBImport
for normal in normal1.bam normal2.bam normal3.bam; do
    s=$(basename "$normal" .bam)
    gatk Mutect2 -R reference.fa -I "$normal" --max-mnp-distance 0 -O "${s}.vcf.gz"
done

gatk GenomicsDBImport -R reference.fa --genomicsdb-workspace-path pon_db \
    -V normal1.vcf.gz -V normal2.vcf.gz -V normal3.vcf.gz -L intervals.bed

gatk CreateSomaticPanelOfNormals -R reference.fa -V gendb://pon_db -O pon.vcf.gz
```

A PoN needs 40+ normals from the same platform/chemistry to capture recurrent artifacts; a PoN from a different assay imports the wrong artifact profile and misses real ones. Do NOT build a PoN from tumor-adjacent normals if they may carry tumor-in-normal contamination.

### Step 2: Call somatic variants (tumor + normal in one command)

```bash
gatk Mutect2 -R reference.fa \
    -I tumor.bam -I normal.bam -normal normal_sample_name \
    --germline-resource af-only-gnomad.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    --f1r2-tar-gz f1r2.tar.gz \
    -O unfiltered.vcf.gz
```

`-normal` takes the normal read-group SM name (not the filename). `--f1r2-tar-gz` collects the read-orientation counts needed in Step 3 - omit it and orientation-bias filtering is impossible. Mutect2 also writes `unfiltered.vcf.gz.stats`, which `FilterMutectCalls` reads automatically.

### Step 3: Learn the orientation-bias model

```bash
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
```

Models the strand-orientation artifacts (oxoG C>A/G>T from oxidative shearing damage; FFPE cytosine-deamination C>T/G>A). These masquerade as low-VAF somatic SNVs; the model lets FilterMutectCalls down-weight them by their F1R2/F2R1 imbalance.

### Step 4: Estimate contamination

```bash
gatk GetPileupSummaries -I tumor.bam \
    -V small_exac_common_3.vcf.gz -L small_exac_common_3.vcf.gz -O tumor_pileups.table
gatk GetPileupSummaries -I normal.bam \
    -V small_exac_common_3.vcf.gz -L small_exac_common_3.vcf.gz -O normal_pileups.table

gatk CalculateContamination -I tumor_pileups.table -matched normal_pileups.table \
    -O contamination.table --tumor-segmentation segments.table
```

`GetPileupSummaries` uses a COMMON biallelic-SNP sites resource (e.g. `small_exac_common_3.vcf.gz`), not the af-only-gnomAD used for the germline prior - it needs sites with a known population AF where reference/alt read counts reveal foreign DNA. `--tumor-segmentation` also captures allelic-copy segments that feed the filter.

### Step 5: Filter, then extract PASS

```bash
gatk FilterMutectCalls -R reference.fa -V unfiltered.vcf.gz \
    --contamination-table contamination.table \
    --tumor-segmentation segments.table \
    --ob-priors read-orientation-model.tar.gz \
    -O filtered.vcf.gz

bcftools view -f PASS filtered.vcf.gz -Oz -o somatic_final.vcf.gz
bcftools index -t somatic_final.vcf.gz
```

`FilterMutectCalls` applies a single joint model (contamination + orientation + segmentation + the built-in weak-evidence/germline/strand filters) and sets a per-variant FILTER. Never hand-tune individual thresholds first - the filter is calibrated to balance them together.

## Tumor-only mode and its caveats

When no matched normal exists (archival FFPE, cell lines, legacy cohorts):

```bash
gatk Mutect2 -R reference.fa -I tumor.bam \
    --germline-resource af-only-gnomad.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    -O tumor_only.vcf.gz
```

Decision framing: with no normal, every germline variant is a candidate somatic call, and only the PoN (artifacts) and gnomAD prior (germline) remove them - so both must be assay-matched and the false-positive rate rises sharply. High-VAF (~50% or ~100%) calls are especially suspect for germline. Tumor-only cannot cleanly separate somatic from germline, which creates a disclosure problem (a germline pathogenic finding surfaced as "somatic") - see the interpretation section. Paired tumor-normal is the defensible design; tumor-only requires explicit germline-subtraction logic and reporting caveats.

## VAF, purity, ploidy, and clonality

Somatic VAF is not a genotype - it is `(clonal_fraction x mutation_copies) / local_total_copies`, scaled by tumor purity. Reasoning consequences:

- **Purity sets sensitivity.** At 30% purity a truly clonal heterozygous mutation in diploid regions sits near ~15% VAF; at 10% purity near ~5%. Low-purity samples need higher depth to call the same variant - the one regime where "sequence deeper" is the correct fix (contrast germline, which saturates ~30-35x).
- **Clonal vs subclonal.** High-VAF (adjusted for purity/copy number) variants are clonal (present in most cells, likely early drivers); low-VAF are subclonal (later, spatially/temporally heterogeneous). Clonality informs driver-vs-passenger and treatment-resistance reasoning.
- **Copy number confounds VAF.** A mutation on an amplified allele reads high; on a deleted/LOH background it reads near 1.0. Interpret VAF only alongside local copy number and purity - run copy-number/cnvkit-analysis (or PURPLE) to get purity/ploidy and correct VAF to cancer-cell fraction before calling something subclonal.

## Strelka2 (faster alternative / consensus arm)

```bash
# Run Manta first; its small-indel candidates sharpen Strelka2 indel calls
configManta.py --normalBam normal.bam --tumorBam tumor.bam \
    --referenceFasta reference.fa --runDir manta_run
manta_run/runWorkflow.py -m local -j 16

configureStrelkaSomaticWorkflow.py \
    --normalBam normal.bam --tumorBam tumor.bam --referenceFasta reference.fa \
    --indelCandidates manta_run/results/variants/candidateSmallIndels.vcf.gz \
    --runDir strelka_run
strelka_run/runWorkflow.py -m local -j 16

bcftools concat \
    strelka_run/results/variants/somatic.snvs.vcf.gz \
    strelka_run/results/variants/somatic.indels.vcf.gz \
    -a -Oz -o strelka_somatic.vcf.gz
```

Strelka2 is faster and strong on indels (Kim 2018 *Nat Methods* 15:591); passing Manta's `candidateSmallIndels.vcf.gz` via `--indelCandidates` is the documented coupling that improves Strelka2 indel recall.

## Multi-caller consensus and reproducibility

Given the Alioto/DREAM discordance, running independent callers and requiring agreement raises precision - the concordant core is the reproducible callset:

```bash
# Intersect PASS calls from callers with uncorrelated error modes (2/3 agreement)
bcftools isec -n+2 -p consensus_dir \
    mutect2_pass.vcf.gz strelka2_pass.vcf.gz muse_pass.vcf.gz
```

Strict all-agree intersection sacrifices too much recall; union admits too many false positives; majority voting (e.g. 2 of 3) balances the two. Normalize every caller's VCF to the same representation first (`bcftools norm -f ref.fa -m-`) or the intersection undercounts because an indel left-aligned differently in two callers will not match. Consensus is a precision tool, not a substitute for orthogonal validation on the assay.

## Somatic SV, CNV, and genomic biomarkers

A complete somatic profile is more than SNVs/indels - hand each off to its component skill:

- **Structural variants:** run Manta in tumor-normal (somatic) mode; it emits `somaticSV.vcf.gz`. For complex rearrangements/fusions use GRIDSS -> GRIPSS -> LINX. See variant-calling/structural-variant-calling.
- **Copy number + purity/ploidy:** CNVkit (targeted/exome), or PURPLE for allele-specific CN with purity/ploidy fit. Purity/ploidy is a required input for correct VAF-to-cancer-cell-fraction. See copy-number/cnvkit-analysis.
- **Genomic biomarkers (brief):** Tumor Mutational Burden (TMB) = eligible somatic mutations per Mb of covered target (filter out germline and low-VAF artifacts first - inflated TMB is usually residual germline in tumor-only). Microsatellite instability (MSI) is called from indel patterns at microsatellite loci (e.g. MSIsensor). Mutational signatures (SBS/ID/CNV, COSMIC catalogue) attribute the mutation spectrum to processes (APOBEC, UV, MMR-deficiency, platinum) and cross-check artifacts - a dominant C>A/G>T signature can be residual oxoG, not biology. Verify current tool choices and thresholds against the latest docs.

## Annotation (normalize FIRST)

```bash
bcftools norm -f reference.fa -m- somatic_final.vcf.gz -Oz -o somatic_norm.vcf.gz  # left-align + split multiallelics

gatk Funcotator -R reference.fa -V somatic_norm.vcf.gz -O annotated.vcf.gz \
    --output-file-format VCF --data-sources-path funcotator_dataSources.v1.7 --ref-version hg38

vep -i somatic_norm.vcf.gz -o annotated_vep.vcf --vcf --cache --offline \
    --assembly GRCh38 --everything \
    --custom cosmic.vcf.gz,COSMIC,vcf,exact,0,CNT --fork 4
```

Normalize BEFORE annotating (annotate-then-normalize is an order error): un-normalized records attach annotations to a non-canonical representation and fail to match COSMIC/gnomAD. Full annotation mechanics (transcript choice, MANE, HGVS 3'-shift) live in variant-calling/variant-annotation.

## Interpretation: somatic tiers, NOT germline ACMG

This is the critical decision of the whole pipeline and the most common category error. Route the annotated somatic VCF to variant-calling/clinical-interpretation, which applies TWO orthogonal cancer frameworks:

**AMP/ASCO/CAP four-tier clinical actionability (Li 2017 *J Mol Diagn* 19:4-23):**

| Tier | Meaning | Example |
|------|---------|---------|
| I | Strong clinical significance - FDA-approved therapy or in professional guidelines for THIS tumor type | BRAF V600E in melanoma |
| II | Potential significance - therapy in a different tumor type, or clinical-trial/multi-study evidence | same variant in a non-approved tumor type |
| III | Unknown clinical significance (the somatic "VUS") | rare novel missense, no actionability |
| IV | Benign / likely benign - high population frequency, no oncogenic role | common polymorphism |

Tier is tumor-type-specific - the SAME variant can be Tier I in one cancer and Tier II/III in another (no germline-ACMG analog).

**Oncogenicity (Horak 2022 *Genet Med* 24:986):** a SEPARATE points-based ClinGen/CGC/VICC axis (Oncogenic / Likely Oncogenic / VUS / Likely Benign / Benign) from hotspot recurrence, functional data, and tumor-type frequency. Oncogenicity (is it a driver) and actionability (is there a drug) are complementary: an oncogenic driver may still be Tier III if no therapy exists.

- **Driver vs passenger:** a tumor carries thousands of somatic mutations; only a few drive it. Hotspot recurrence (COSMIC, Tate 2019), presence in known oncogenes/tumor suppressors, functional evidence, and clonality (clonal drivers vs subclonal passengers) distinguish them.
- **Knowledgebase evidence levels:** OncoKB therapeutic levels 1-4 + R1/R2 (Chakravarty 2017); CIViC evidence items graded A (validated) to E (preclinical), read the item not just the letter (Griffith 2017).
- **Tumor-only germline leak:** in tumor-only assays a high-VAF variant may be germline - reporting it as somatic (or applying a somatic tier to a germline pathogenic finding requiring genetic counseling) is a dual danger. Disclose and filter explicitly.

## Quality metrics

```bash
bcftools query -f '%FILTER\n' filtered.vcf.gz | sort | uniq -c   # counts by filter status

# Substitution spectrum: excess C>A/G>T flags residual oxoG; excess C>T/G>A flags FFPE deamination
bcftools query -f '%REF>%ALT\n' somatic_final.vcf.gz | sort | uniq -c

# VAF distribution (a spike near 0.5/1.0 in tumor-only suggests residual germline).
# -s takes the read-group SM name, NOT the filename -- read it from the BAM header, as Step 2 does.
TUMOR_SM=$(samtools view -H tumor.bam | awk '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/){sub(/^SM:/,"",$i); print $i; exit}}')
# -s restricts to the tumor sample: a bare [%AF] iterates BOTH samples and concatenates
# them with no separator (0.25 + 0.01 -> "0.250.01"), which awk then silently truncates to 0.25.
bcftools query -s "$TUMOR_SM" -f '[%AF]\n' somatic_final.vcf.gz | awk '{print int($1*100)/100}' | sort -n | uniq -c
```

Do NOT judge somatic SNVs by a fixed germline-like Ti/Tv (~2-3); the somatic spectrum is signature-dependent, and a collapse toward transversion excess signals artifact contamination, not a target value.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Flood of germline variants in output | Missing/wrong germline-resource, or tumor-only without a good PoN | Pass `--germline-resource af-only-gnomad.vcf.gz`; use an assay-matched 40+ normal PoN |
| Excess C>A/G>T (or C>T/G>A) low-VAF calls | oxoG (or FFPE) artifacts, orientation-bias filter not applied | Emit `--f1r2-tar-gz`, run `LearnReadOrientationModel`, pass `--ob-priors` to FilterMutectCalls |
| `FilterMutectCalls` errors on missing stats | `unfiltered.vcf.gz.stats` not alongside the VCF | Keep the `.stats` Mutect2 wrote next to the VCF, or pass `--stats` |
| GetPileupSummaries gives nonsense contamination | Used af-only-gnomAD instead of a common biallelic-SNP resource | Use `small_exac_common_3.vcf.gz` (common SNP sites with AF) |
| Consensus intersection drops real shared calls | VCFs not normalized before `bcftools isec` | `bcftools norm -f ref.fa -m-` every caller's VCF first |
| Low-purity tumor: expected drivers missing | VAF below detection at that purity/depth | Get purity from copy-number; increase depth; correct VAF to cancer-cell fraction |
| Applied ACMG PVS1/PM2 to a tumor variant | Germline framework used for somatic | Use AMP/ASCO/CAP tiers + oncogenicity via variant-calling/clinical-interpretation |

## Related Skills

- variant-calling/gatk-variant-calling - Mutect2 mechanism and germline HaplotypeCaller context
- variant-calling/filtering-best-practices - FilterMutectCalls internals, normalization, hard filters
- variant-calling/variant-annotation - VEP/Funcotator/SnpEff, transcript choice, HGVS, COSMIC
- variant-calling/clinical-interpretation - AMP/ASCO/CAP tiers and oncogenicity (somatic interpretation)
- variant-calling/structural-variant-calling - Somatic SV detection (Manta somatic mode, GRIDSS/LINX)
- copy-number/cnvkit-analysis - Somatic CNV, purity/ploidy for VAF-to-cancer-cell-fraction
- read-alignment/bwa-alignment - Upstream alignment/dedup/BQSR of tumor and normal BAMs

## References

- Li MM, Datto M, Duncavage EJ, et al. Standards and guidelines for the interpretation and reporting of sequence variants in cancer: a joint consensus recommendation of AMP, ASCO, and CAP. *Journal of Molecular Diagnostics*. 2017;19(1):4-23. doi:10.1016/j.jmoldx.2016.10.002 (four-tier actionability).
- Horak P, Griffith M, Danos AM, et al. Standards for the classification of pathogenicity of somatic variants in cancer (oncogenicity): joint recommendations of ClinGen, CGC, and VICC. *Genetics in Medicine*. 2022;24(5):986-998. doi:10.1016/j.gim.2022.01.001.
- Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C, Lichtenstein L. Calling Somatic SNVs and Indels with Mutect2. *bioRxiv* 861054 (2019). doi:10.1101/861054 (PREPRINT - the Mutect2 somatic workflow reference; never formally journal-published).
- Alioto TS, Buchhalter I, Derdak S, et al. A comprehensive assessment of somatic mutation detection in cancer using whole-genome sequencing. *Nature Communications*. 2015;6:10001. doi:10.1038/ncomms10001 (ICGC-TCGA DREAM; somatic reproducibility).
- Kim S, Scheffler K, Halpern AL, et al. Strelka2: fast and accurate calling of germline and somatic variants. *Nature Methods*. 2018;15(8):591-594. doi:10.1038/s41592-018-0051-x.
- Chen X, Schulz-Trieglaff O, Shaw R, et al. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics*. 2016;32(8):1220-1222. doi:10.1093/bioinformatics/btv710.
- Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*. 2020;581(7809):434-443. doi:10.1038/s41586-020-2308-7 (gnomAD; germline-resource / AF prior).
- Chakravarty D, Gao J, Phillips SM, et al. OncoKB: a precision oncology knowledge base. *JCO Precision Oncology*. 2017;2017:PO.17.00011. doi:10.1200/PO.17.00011 (therapeutic levels 1-4, R1/R2).
- Griffith M, Spies NC, Krysiak K, et al. CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer. *Nature Genetics*. 2017;49:170-174. doi:10.1038/ng.3774 (evidence levels A-E).
- Tate JG, Bamford S, Jubb HC, et al. COSMIC: the Catalogue Of Somatic Mutations In Cancer. *Nucleic Acids Research*. 2019;47(D1):D941-D947. doi:10.1093/nar/gky1015 (hotspot recurrence, signatures).
