---
name: bio-variant-calling-filtering-best-practices
description: Filters germline and somatic variant callsets at the site and genotype level with GATK VQSR (VQSLOD, truth-sensitivity tranches), VETS/ScoreVariantAnnotations, NVScoreVariants, hard filters with per-annotation thresholds, and bcftools/cyvcf2 expressions, plus Ti/Tv-based QC. Use when deciding between VQSR, hard filtering, and ML recalibration by cohort size and platform, setting SNP vs indel thresholds, replicating the missing-annotation-passes rule so hom-alt sites survive, applying genotype-level GQ/DP filters, or validating filter impact. Not for VCF normalization (see variant-calling/variant-normalization) or summary statistics (see variant-calling/vcf-statistics).
tool_type: mixed
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: GATK 4.6+, bcftools 1.19+, cyvcf2 0.30+

Note: CNNScoreVariants is deprecated as of GATK 4.6.1.0 (replaced by NVScoreVariants, a PyTorch drop-in); VETS (ExtractVariantAnnotations/TrainVariantAnnotationsModel/ScoreVariantAnnotations) is BETA. Confirm tool availability with `gatk --list` on the installed build before scripting a pipeline.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Variant Filtering Best Practices

**"Filter my variant calls"** -> Flag or remove low-quality variant sites, and separately null out untrustworthy per-sample genotypes, using a model matched to cohort size, platform, and organism.
- CLI: GATK VariantRecalibrator/ApplyVQSR (large cohorts), VariantFiltration (hard filters), ScoreVariantAnnotations (VETS), NVScoreVariants (single-sample DL); bcftools filter/view
- Python: cyvcf2 for custom per-variant logic

## The governing principle

Filtering decides WHICH errors a callset keeps, not whether they exist. The two site-level paradigms fail in OPPOSITE regimes: VQSR (a ratio of two learned Gaussian-mixture densities over annotation space) collapses on small or exome cohorts; static hard thresholds discard real variants at scale. Two rules follow. First, site-level filtering ("is this SITE real?") and genotype-level filtering ("is this SAMPLE's genotype trustworthy?") are orthogonal -- both are needed, site first. Second, SNPs and indels have different error processes (base-calling/strand vs alignment-ambiguity-in-repeats) and different truth resources, so they are ALWAYS filtered separately then merged. None of these mistakes throws an error: the VCF stays structurally valid while the numbers are silently wrong.

## Site-Level Filter Method Selection

Somatic data is a separate track: use GATK FilterMutectCalls, never VQSR or germline hard filters (the annotations and error model differ).

| Method | Best when | Fails when |
|--------|-----------|------------|
| Hard filters (VariantFiltration) | Single sample, exome, targeted panel, non-model organism, or any callset lacking truth resources | Precision-critical work at scale -- static cutoffs leave real variants on the table |
| VQSR (VariantRecalibrator/ApplyVQSR) | Human, a single deep WGS OR ~30+ jointly-genotyped exomes, HapMap/Omni/Mills truth sets available | A single exome/panel (too few variants): the GMM is non-identifiable and VQSLOD is noise |
| Allele-specific VQSR (`-AS`, `AS_*` annotations) | Very large cohorts (biobank/gnomAD scale) where one bad allele at a multiallelic must not sink the site | Small cohorts; adds nothing over site-level VQSR |
| VETS (ScoreVariantAnnotations, BETA) | Modern GATK replacement for VQSR; scikit-learn isolation-forest on site annotations, more robust than GMM, works down to smaller cohorts | Still BETA -- validate against a truth set before production use |
| NVScoreVariants (deep learning) | A single sample, especially a single exome/panel where VQSR has too few variants to train; PyTorch CNN scores reads+reference, then FilterVariantTranches applies tranches | Needs a GPU-friendly env for the 2D model; replaced deprecated CNNScoreVariants |
| DL-native caller output (DeepVariant, DRAGEN ML) | The caller already emits calibrated QUAL / vendor FILTER flags | Do NOT re-apply GATK hard filters on top -- annotation distributions differ; filter on the caller's own fields |

Methodology is evolving (VETS is displacing VQSR). Verify the current recommended path against the installed GATK version's "How to Filter variants" article before committing a pipeline.

## VQSR -- Mechanism and Why It Breaks

**Goal:** Recalibrate a large jointly-genotyped human callset with a data-driven quality score.

**Approach:** Fit a Gaussian mixture model (GMM) to the annotation profile of known-true sites (positive model), bootstrap a second GMM on the low-probability-tail artifact sites (negative model), and score each variant by VQSLOD = log( P(annotations | positive) / P(annotations | negative) ). Then choose a truth-sensitivity TRANCHE rather than thresholding VQSLOD directly: a "99.7 tranche" is the VQSLOD cutoff that RETAINS 99.7% of the truth-set sites. Tranches are truth-set SENSITIVITIES, not FDRs.

Three load-bearing consequences the agent must respect:
- VQSR estimates full covariance matrices in ~6-8 annotation dimensions, so it needs tens of thousands of variants -- as practical GATK convention, a single deep WGS (which alone supplies millions of sites) OR ~30+ jointly-genotyped exomes. On a single exome or panel the model is non-identifiable or wildly overfit -- it may report "converged" while VQSLOD is garbage. This is the single most common real-world VQSR misuse.
- On exomes, DP must NOT be supplied as a VQSR annotation: capture depth tracks bait design, not truth, and injects a spurious signal.
- SNPs and indels are recalibrated in separate runs (`-mode SNP`, `-mode INDEL`) because indels are ~10x rarer and their GMM fails first (see governing principle).

```bash
# SNP recalibration: fit the GMM in annotation space against truth resources
gatk VariantRecalibrator \
    -R reference.fa -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O snp.recal --tranches-file snp.tranches
# For exomes: OMIT -an DP (capture depth is uninformative of truth), add -an QD -an FS etc. only

# Apply the chosen truth-sensitivity tranche (keeps 99.7% of truth-set SNPs)
gatk ApplyVQSR \
    -R reference.fa -V cohort.vcf.gz \
    -mode SNP --recal-file snp.recal --tranches-file snp.tranches \
    --truth-sensitivity-filter-level 99.7 \
    -O snp.recalibrated.vcf.gz
```

Run the identical pair with `-mode INDEL` and the Mills/1000G gold-indel resource, then merge the recalibrated SNP and indel callsets. For a single exome or panel (too few variants for VQSR), replace this whole block with hard filters or NVScoreVariants.

## GATK Hard Filters (SNPs and indels separately)

**Goal:** Flag artifacts with static, per-annotation thresholds when VQSR is inapplicable.

**Approach:** Split the callset by type (`SelectVariants`), apply type-appropriate OR-combined fail conditions with `VariantFiltration`, then merge. Each annotation targets an independent error mode; a variant fails if it violates ANY one.

**"Filter my variants using GATK best practices"** -> Apply GATK's recommended annotation cutoffs, separately for SNPs and indels.

```bash
# SNPs
gatk VariantFiltration -R reference.fa -V raw_snps.vcf -O filtered_snps.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3"

# Indels: FS loosened to 200, ReadPosRankSum tightened to -20, MQ/MQRankSum DROPPED
gatk VariantFiltration -R reference.fa -V raw_indels.vcf -O filtered_indels.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10"
```

The SNP/indel threshold difference is the point, not an inconsistency: real indels have messier local alignments in repeats, so strand bias (FS) is naturally higher and the gate is loosened to 200; spurious indels cluster at read ends, so ReadPosRankSum is tightened to -20. Mapping-quality metrics (MQ, MQRankSum) are dropped for indels because they are less diagnostic there and the truth model is weaker. Values are GATK-recommended lenient starting points -- verify against the installed version's docs and tune to the annotation histograms of the dataset.

### The hom-alt "missing => PASS" trap

MQRankSum, ReadPosRankSum, and BaseQRankSum are rank-sum tests comparing ref- vs alt-supporting reads, so they are only DEFINED at heterozygous sites. At hom-alt sites there are no ref reads and the annotation is missing (`.`). GATK's `VariantFiltration` fires a filter only when the value is PRESENT and violates the cutoff -- a missing value PASSES. Anyone hand-writing the equivalent in bcftools MUST replicate this: guard every RankSum term with an explicit `|| INFO/X = "."`, or every hom-alt variant silently fails and vanishes.

```bash
# GATK SNP hard filter (plus a QUAL>=30 floor, which is not part of GATK's canonical set)
# -- the "|| = \".\"" guard on each RankSum term lets hom-alt sites (undefined RankSums) pass
bcftools filter -i '
    QUAL >= 30 && (INFO/QD >= 2.0 || INFO/QD = ".") &&
    (INFO/FS <= 60.0 || INFO/FS = ".") && (INFO/MQ >= 40.0 || INFO/MQ = ".") &&
    (INFO/MQRankSum >= -12.5 || INFO/MQRankSum = ".") &&
    (INFO/ReadPosRankSum >= -8.0 || INFO/ReadPosRankSum = ".") &&
    (INFO/SOR <= 3.0 || INFO/SOR = ".")' raw_snps.vcf.gz -Oz -o snps_filtered.vcf.gz
```

## Quality Metric Rationale

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| QD (QualByDepth) | <2.0 | QUAL normalized by alt-supporting depth. Raw QUAL grows with coverage, so a 500x artifact can post a huge QUAL; QD removes that inflation. Bimodal in practice -- real variants ~12-35, artifacts near 0. The workhorse, not QUAL. |
| FS (FisherStrand) | >60 (SNP), >200 (indel) | Phred-scaled Fisher's-exact p-value for strand bias. Real variants are strand-symmetric; many artifacts are strand-specific. Breaks down at exon/read ends where SOR takes over. |
| SOR (StrandOddsRatio) | >3.0 (SNP); >10.0 (indel) is a commonly-added community/WDL convention, not part of GATK's canonical indel set (QD/QUAL/FS/ReadPosRankSum) | Symmetric-odds strand-bias metric that tolerates the legitimate strand imbalance at exon/read ends where FS false-positives. Complements FS, does not replace it. |
| MQ (RMSMappingQuality) | <40.0 | RMS mapping quality of reads at the site. Low MQ => reads map ambiguously (repeats, paralogs, segdups) => likely mapping artifact. |
| MQRankSum | <-12.5 | Rank-sum of mapping quality, alt- vs ref-supporting reads. Strongly negative => alt reads map worse => probable mismapping. Missing at hom-alt sites. |
| ReadPosRankSum | <-8.0 (SNP), <-20.0 (indel) | Rank-sum of within-read position, alt vs ref bases. Strongly negative => alt clusters at read ends (highest error, least reliable alignment). Missing at hom-alt sites. |
| DP (depth) | context-specific | Extreme depth (>2x or <0.3x mean) suggests collapsed repeats or poor capture. Filtering on DP alone removes real variants in duplicated regions -- always combine with MQ/MQRankSum. Never a VQSR annotation on exomes. |
| GQ (genotype quality) | <20 | Genotype-level, not site-level. Phred confidence in the called genotype; GQ 20 = 99%. |

## Site-Level vs Genotype-Level Filtering

The two are orthogonal and both are required, in order. Site filters (above) decide whether a SITE is real. Genotype filters set an individual sample's genotype to no-call (`./.`) when it is untrustworthy at an otherwise-passing site:
- GQ < 20 => set `./.` (genotype confidence below 99%).
- DP < 8-10 => set `./.` (too few reads for a confident diploid call, for WGS).
- Allele balance far from 0.5 at hets (e.g. alt fraction <0.2 or >0.8) => suspect mapping artifact, CNV, or contamination; derive from AD (GATK does not emit AB directly).

Ordering matters: apply genotype-level no-calls BEFORE computing cohort metrics (missingness, HWE, allele frequency). Computing HWE on a matrix full of low-GQ garbage genotypes manufactures spurious deviation. Pipeline: site filter -> genotype filter -> recompute cohort QC.

```bash
# Genotype-level: null out low-confidence genotypes, keeping the site
bcftools filter -S . -e 'FMT/GQ<20 | FMT/DP<8' passing_sites.vcf.gz -Oz -o gt_filtered.vcf.gz
```

## bcftools filter -- Soft vs Hard

**Goal:** Flag (soft) or remove (hard) variants by expression on QUAL, INFO, and FORMAT fields.

**Approach:** `-e` excludes, `-i` includes; `-s NAME` writes a named FILTER label instead of dropping; `bcftools view -f PASS` extracts survivors at the end.

```bash
bcftools filter -e 'QUAL<30' input.vcf.gz -o filtered.vcf          # hard: drop failing
bcftools filter -s 'LowQual' -e 'QUAL<30' input.vcf.gz -o marked.vcf  # soft: label failing
bcftools view -f PASS marked.vcf -o passed.vcf                      # extract PASS survivors
```

Operators: `< <= > >=  = == !=  && ||  !`. Aggregate over samples with `MIN() MAX() AVG() SUM()`. Guard against missing values explicitly (`INFO/DP!="."`), for the same hom-alt reason as above.

## Somatic Variant Filtering

**Goal:** Filter tumor-normal somatic calls with the caller's own model, not germline thresholds.

**Approach:** Run GATK FilterMutectCalls with contamination and segmentation tables, then layer additional thresholds on TLOD and VAF.

```bash
gatk FilterMutectCalls -R reference.fa -V mutect2_raw.vcf \
    --contamination-table contamination.table \
    --tumor-segmentation segments.table \
    -O mutect2_filtered.vcf
bcftools filter -i 'INFO/TLOD>6.3 && FMT/AF[0]>0.05 && FMT/DP[0]>20' \
    mutect2_filtered.vcf -o somatic_final.vcf
```

## Python Filtering (cyvcf2)

**Goal:** Apply custom multi-metric per-variant logic in Python.

**Approach:** Iterate with cyvcf2, read QUAL/INFO fields, write survivors with Writer. `INFO.get` returns None for missing tags -- treat None as pass to avoid the hom-alt trap.

```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('filtered.vcf', vcf)
for variant in vcf:
    qual = variant.QUAL or 0
    dp = variant.INFO.get('DP') or 1e9      # missing depth => do not fail on depth
    fs = variant.INFO.get('FS') or 0.0      # missing strand bias => pass (None -> 0)
    mq = variant.INFO.get('MQ') or 1e9      # missing MQ => pass
    if qual >= 30 and dp >= 10 and fs <= 60.0 and mq >= 40.0:
        writer.write_record(variant)
writer.close(); vcf.close()
```

## Validate Filtering

**Goal:** Confirm filtering removed artifacts without stripping true variants.

**Approach:** Compare before/after `bcftools stats`; check Ti/Tv and Het/Hom against expected ranges and known-variant recovery. A filter that improves one metric while degrading another is miscalibrated.

| Metric | WGS | WES | Interpretation |
|--------|-----|-----|----------------|
| Ti/Tv | 2.0-2.1 | 3.0-3.3 | Below range => excess false positives (random errors have Ti/Tv ~0.5, diluting the signal); a WES set at ~2.1 signals too-loose filtering. WES is higher from CpG-transition-rich coding enrichment. |
| Het/Hom | 1.5-2.0 | 1.5-2.0 | Strongly ancestry-dependent. Elevated => contamination; depressed => inbreeding/ROH. Stratify by ancestry before flagging outliers. |
| Known (dbSNP) % | >99% | >99% | Low known-variant recovery indicates over-filtering. |

```bash
bcftools stats input.vcf > before.txt
bcftools stats filtered.vcf | grep '^TSTV'                    # Ti/Tv after filtering
bcftools query -f '%FILTER\n' filtered.vcf | sort | uniq -c   # counts per FILTER label
```

If Ti/Tv drops after filtering, the filters are preferentially removing true transitions -- relax them. See variant-calling/vcf-statistics for the full QC panel (het/hom by ancestry, contamination, relatedness).

## Region-Based Filtering

Stratify by genomic context; artifact-prone regions dominate false positives. Exclude with `bcftools view -T ^regions.bed`:
- ENCODE exclusion list (github.com/Boyle-Lab/Blacklist) -- anomalous-signal regions (centromeres, satellites).
- GIAB stratification BEDs -- low-complexity, segdups, tandem repeats; essential for honest benchmarking (`bcftools isec` against a GIAB truth set).
- LCR-hs38 (Heng Li) -- homopolymers and simple repeats where indel calling is unreliable.

## Common Filtering Pitfalls

- Applying SNP thresholds to indels: distributions differ (FS, SOR, ReadPosRankSum). Always split by type first.
- Treating missing RankSum as failing: silently deletes every hom-alt site (see the missing => PASS trap).
- Running VQSR on a single exome or panel: too few variants, the GMM is non-identifiable (a single deep WGS is fine); use hard filters, VETS, or NVScoreVariants.
- Re-applying GATK hard filters on DeepVariant/DRAGEN output: their calibrated fields already encode quality; the GATK annotations may be absent or differently distributed.
- Filtering on depth alone: removes real variants in collapsed segdups; combine DP with MQ/MQRankSum.
- Choosing thresholds without looking: plot each annotation stratified by known TP (HapMap) vs likely FP and cut at the valley; GATK defaults are population-level starting points.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `no such INFO tag` | Tag absent from VCF | Check header: `bcftools view -h in.vcf` |
| `syntax error` in expression | Invalid operator | Use `\|\|` not `or`; quote missing as `= "."` |
| Every hom-alt site removed | RankSum missing not guarded | Add `\|\| INFO/X = "."` to each RankSum term |
| VQSR "converged" but nonsense | Too few samples/variants | Switch to hard filters, VETS, or NVScoreVariants |
| empty output | Filter too strict | Relax thresholds; inspect annotation histograms |

## Related Skills

- variant-calling/variant-calling - Variant calling with bcftools to generate VCF files
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller and joint genotyping upstream of VQSR
- variant-calling/deepvariant - Deep-learning caller whose output needs no separate site filter
- variant-calling/variant-annotation - Functional annotation after filtering
- variant-calling/variant-normalization - Left-align and decompose before filtering for consistent comparisons
- variant-calling/vcf-statistics - Ti/Tv, het/hom, contamination, and relatedness QC of filter effects
- variant-calling/vcf-basics - VCF field interpretation, PL/GQ/QUAL, and Number=A/R/G subsetting

## References

- DePristo MA, Banks E, Poplin R, et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics.* 2011;43(5):491-498. -- VQSR foundational description.
- Van der Auwera GA, Carneiro MO, Hartl C, et al. From FastQ Data to High-Confidence Variant Calls: the GATK Best Practices pipeline. *Current Protocols in Bioinformatics.* 2013;43:11.10.1-11.10.33. -- hard-filtering and VQSR usage.
- Poplin R, Chang P-C, Alexander D, et al. A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology.* 2018;36(10):983-987. -- DeepVariant (DL-native calls need no separate filter).
- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. -- bcftools filter/view/stats.
