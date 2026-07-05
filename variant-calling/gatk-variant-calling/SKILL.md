---
name: bio-gatk-variant-calling
description: Call germline SNPs and indels with GATK HaplotypeCaller and the GVCF joint-genotyping workflow. Covers the local-reassembly + PairHMM mechanism (why HC beats pileup callers on indels), the -ERC GVCF reference-confidence model and <NON_REF> allele, BQSR-vs-DRAGSTR and --dragen-mode error modeling, allele-specific (AS_) annotations, and edge cases (ploidy, Mutect2 mitochondria mode, sex chromosomes/PAR, contamination gating). Use when deciding whether to use HaplotypeCaller vs a pileup or DRAGEN caller, whether BQSR still earns its place, whether to call per-sample GVCFs for a cohort, or how to handle non-diploid, mitochondrial, sex-chromosome, or contaminated samples. Not for post-calling filtering depth (see variant-calling/filtering-best-practices) or cohort joint-genotyping scaling (see variant-calling/joint-calling).
tool_type: cli
primary_tool: gatk
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: flag DEFAULTS (e.g. `--max-alternate-alleles`, `--heterozygosity`, `--standard-min-confidence-threshold-for-calling`) drift across GATK 4.x releases and DRAGEN-mode ports models into the caller. Confirm any default with `gatk <Tool> --help` rather than trusting a memorized value.

# GATK Variant Calling

**"Call germline variants from my BAM with GATK"** -> Detect SNPs/indels by locally reassembling candidate haplotypes in active regions and genotyping reads against them.
- CLI: `gatk HaplotypeCaller` (germline diploid), `gatk Mutect2` (somatic / mitochondria / mosaic)

## Why HaplotypeCaller, not a pileup caller

A pileup/position-based caller (bcftools mpileup, the old UnifiedGenotyper) genotypes each column of the alignment independently, trusting the mapper's per-read placement. Near indels and clustered variation that placement is a per-read greedy optimum, not a locus-consistent one, so the same indel gets (mis)placed differently in different reads and is systematically misrepresented. HaplotypeCaller (HC) discards the local alignment and re-derives it: in any region showing signal it performs **local de-novo reassembly** of candidate haplotypes, then realigns every read against those haplotypes. The indel is then represented ONCE, in an assembled haplotype, not independently in each read. This is why HC (and DRAGEN, and DeepVariant) beat pileup callers on indels and complex loci, and why HC is the reference implementation everything else is benchmarked against.

### How HaplotypeCaller works (decision-relevant mechanism)

1. **Active-region determination.** HC computes a per-locus activity score (a fast reference-vs-non-reference genotype-likelihood contrast on the pileup), smooths it with a Gaussian kernel, and thresholds it (`--active-probability-threshold`, default 0.002) to seed active regions, padded so the assembler sees flanking reference. Non-active bases still get a reference-confidence emission in GVCF mode. Pathologically high-depth or highly repetitive regions can fail to assemble -- exactly where DeepVariant/DRAGEN pull ahead.
2. **Local reassembly.** For each active region HC builds a de-Bruijn-like graph from k-mers of the reference plus overlapping reads (default k = 10 and 25, merged), prunes low-weight (error) edges, and enumerates best-supported haplotypes. Each haplotype is Smith-Waterman-realigned to the reference to translate assembled sequence into concrete SNP/indel events.
3. **PairHMM read-vs-haplotype likelihoods.** For every (read, haplotype) pair HC computes P(read | haplotype) with a Pair Hidden Markov Model that integrates over ALL alignments (forward algorithm), not just the best one -- the statistically correct way to weight support when the alignment itself is uncertain. This is the dominant compute cost; production runs use vectorized (AVX/AVX-512) kernels via `--pair-hmm-implementation` and `--native-pair-hmm-threads`. DRAGEN moves this same kernel onto an FPGA.
4. **Genotype likelihoods.** Per-haplotype likelihoods are marginalized to per-allele likelihoods, then genotype likelihoods (PLs) are computed under the assumed ploidy via the Bayesian DePristo/GATK model, emitting GT/AD/DP/GQ/PL plus site annotations (QD, FS, MQ, MQRankSum, ReadPosRankSum, SOR).
5. **The diploid assumption.** Genotyping defaults to diploid (`--sample-ploidy 2`), which hard-codes true alleles at fraction 0, 0.5, or 1.0. Anything violating that (pooled, polyploid, CNV, mosaic, hemizygous chrX/Y) needs an explicit `--sample-ploidy` or a somatic caller -- see Edge Cases.

## Pipeline decision tree

```
What is the analysis context?
├── Single sample, want DRAGEN-like accuracy, open-source -> HaplotypeCaller --dragen-mode (hard-filter on QUAL)
├── Cohort < ~2000, human -> per-sample -ERC GVCF -> joint genotype -> VQSR/VETS (see joint-calling)
├── Cohort > ~2000, human -> ReblockGVCF + GnarlyGenotyper "Biggest Practices", or DeepVariant + GLnexus (see joint-calling)
├── Non-human / non-model organism -> hard filtering (no VQSR training resources)
├── Targeted panel / small exome -> hard filtering (too few variants for VQSR)
├── Non-diploid / pooled / sex chromosomes / mitochondria -> set --sample-ploidy or use Mutect2 (see Edge Cases)
└── Somatic / mosaic variants -> Mutect2 (not HaplotypeCaller)
```

## Single-sample calling

**Goal:** Call germline SNPs and indels from one sample.

**Approach:** Run HaplotypeCaller directly to VCF; add annotations, intervals, or a calling-confidence floor as needed.

```bash
# Basic call
gatk HaplotypeCaller -R reference.fa -I sample.bam -O sample.vcf.gz

# Exome / panel: restrict to capture targets (much faster, fewer off-target artifacts)
gatk HaplotypeCaller -R reference.fa -I sample.bam -L targets.interval_list -O sample.vcf.gz

# Add standard annotations explicitly (usually emitted by default; force when a downstream filter needs them)
gatk HaplotypeCaller -R reference.fa -I sample.bam -O sample.vcf.gz \
    -A Coverage -A QualByDepth -A FisherStrand -A StrandOddsRatio \
    -A MappingQualityRankSumTest -A ReadPosRankSumTest
```

MarkDuplicates is required before calling for all modes. Whether BQSR precedes it is a real decision -- see below.

## The GVCF reference-confidence model (`-ERC GVCF`)

**Goal:** Make per-sample calls that can later be combined into a cohort without re-visiting BAMs.

**Approach:** Emit a GVCF that records confidence at EVERY position, then joint-genotype separately.

`-ERC GVCF` records, at every position (variant and non-variant), how confidently the site is homozygous reference. Two properties make it the backbone of scalable cohort calling:
- **The `<NON_REF>` symbolic allele.** Every record carries a symbolic ALT `<NON_REF>` ("any allele not yet observed") with AD/PL computed against it. At joint-genotyping time a variant discovered in ANOTHER sample can therefore be evaluated in THIS sample even though this sample looked reference -- the evidence against `<NON_REF>` supplies it. This is what makes per-sample GVCFs forward-compatible with alleles found later in the cohort, and lets joint genotyping distinguish confident hom-ref from no-data (the "squaring-off" of a ragged genotype matrix).
- **GQ banding.** Contiguous non-variant sites with similar GQ collapse into homRef blocks so a GVCF is not one line per base. `-ERC BP_RESOLUTION` disables banding (one line per base, larger files).

```bash
# Per-sample GVCF (do this once per sample; reusable when the cohort grows -- the N+1 win)
gatk HaplotypeCaller -R reference.fa -I sample.bam -O sample.g.vcf.gz -ERC GVCF

# Single-sample genotyping straight from one GVCF
gatk GenotypeGVCFs -R reference.fa -V sample.g.vcf.gz -O sample.vcf.gz
```

The **N+1 problem**: naive joint calling re-processes all N samples whenever the cohort changes; the GVCF captures the expensive assembly/likelihood work once per sample, so adding sample N+1 only re-runs the cheap consolidation + GenotypeGVCFs. Cohort consolidation (`GenomicsDBImport` vs `CombineGVCFs`), scaling to tens of thousands (`ReblockGVCF`, `GnarlyGenotyper`), and joint-genotyping mechanics live in variant-calling/joint-calling -- not duplicated here.

## BQSR, DRAGSTR, and DRAGEN-GATK mode

**Does BQSR still earn its place? (honestly unsettled).** BaseRecalibrator/ApplyBQSR builds an empirical error model over base-call covariates (reported quality, read group, cycle, sequence context) to correct systematic, instrument-specific miscalibration. On older continuous-quality 4-color instruments (HiSeq, MiSeq) it mattered. On modern 2-color chemistry (NovaSeq/NextSeq) qualities are emitted in only ~4 coarse bins, leaving little smooth structure to recalibrate; empirically the callset is largely unchanged with vs without BQSR, with the delta concentrated in borderline-GQ sites. Broad keeps BQSR in Best Practices for pipeline consistency; several large pipelines drop it on binned data. Treat it as caller/instrument-dependent, not mandatory -- there is no consensus universal recommendation.

```bash
gatk BaseRecalibrator -R reference.fa -I sample.bam \
    --known-sites dbsnp.vcf.gz --known-sites Mills_and_1000G_gold_standard.indels.vcf.gz \
    -O recal.table
gatk ApplyBQSR -R reference.fa -I sample.bam --bqsr-recal-file recal.table -O sample.recal.bam
```

**Where the indel-accuracy lever actually moved: DRAGSTR.** Indel errors scale with tandem-repeat context, so the bigger gain is STR-aware indel modeling, not base-quality recalibration. DRAGEN-GATK adds **DRAGSTR**: a per-sample auto-calibration (`CalibrateDragstrModel` -> `--dragstr-params-path`) that models a-priori indel error/variant probability as a function of STR period (repeat-unit length) and length (copy number), adjusting the PairHMM indel gap priors before genotyping.

**`--dragen-mode`.** DRAGEN is Illumina's FPGA-accelerated map-align-call engine (a genome in ~20-25 min; wins the difficult-to-map precisionFDA V2 regions using alt-aware mapping). Illumina and Broad co-developed **DRAGEN-GATK**, porting DRAGEN's error models into open-source GATK to a *functionally equivalent* pipeline (Regier et al. 2018: pipelines are functionally equivalent when their call differences are smaller than sequencing-replicate differences). `HaplotypeCaller --dragen-mode` enables DRAGSTR plus **BQD** (Base Quality Dropout, systematic local quality collapse) and **FRD** (Foreign Read Detection, contaminating/mismapped reads), and **replaces classic BQSR** -- error modeling moves inside the caller. QUAL is well-calibrated, so hard-filtering on QUAL is sufficient without VQSR.

```bash
# Single-sample DRAGEN mode (no separate BQSR); GVCF variant for cohorts
gatk HaplotypeCaller -R reference.fa -I sample.markdup.bam -O sample.g.vcf.gz -ERC GVCF --dragen-mode

# DRAGEN-mode hard filter: GATK-recommended QUAL floor for DRAGEN-GATK output
gatk VariantFiltration -R reference.fa -V sample.vcf.gz -O sample.filtered.vcf.gz \
    --filter-expression "QUAL < 10.4139" --filter-name "DRAGENHardQUAL"  # Broad's documented DRAGEN-mode QUAL cutoff
```

DRAGEN-ML's advertised FP/FN reductions and speed figures are vendor-reported (Illumina), trained on GIAB truth; treat GIAB benchmark dominance with the overfitting caveat that DL/ML callers may not transfer identically to non-GIAB ancestries.

## Allele-specific annotations (AS_)

Standard annotations (QD, FS, MQ, ...) lump all reads at a site together, so at a multiallelic site a real allele co-located with an error-driven allele shares one site-level pass/fail. **AS_ annotations** (request with `-G AS_StandardAnnotation` during GVCF calling / genotyping) compute each metric per allele (AS_QD, AS_FS, AS_SOR, AS_MQ, AS_MQRankSum, AS_ReadPosRankSum), letting AS_VQSR (`-AS`) filter each allele independently. The benefit grows with cohort size, because multiallelic sites -- where a true allele and an artifact collide at one position -- proliferate as sample count rises. The GVCF workflow propagates the raw per-allele data through joint genotyping to enable this.

## Filtering: choose the method, then see filtering-best-practices

Filtering separates real variants from artifacts AFTER calling; SNPs and indels always filter separately (different annotation distributions). Pick the method here; the full VariantRecalibrator/ApplyVQSR, VETS, and hard-filter recipes with threshold rationale live in variant-calling/filtering-best-practices.

| Context | Method | Why |
|---|---|---|
| Human WGS cohort, many variants | VQSR (or its successor VETS: ExtractVariantAnnotations -> TrainVariantAnnotationsModel -> ScoreVariantAnnotations) | Enough variants + truth resources to fit the model |
| Large cohort with many multiallelics | AS_VQSR (`-AS`) | Per-allele filtering at colliding sites |
| Single exome / gene panel | Hard filtering | Too few variants for a stable GMM (a single WGS has enough; the floor is exome/panel-specific) |
| Non-model organism | Hard filtering | No HapMap/1000G/Mills truth resources |
| DRAGEN-mode output | Hard filter on QUAL | QUAL is well-calibrated |
| Somatic (Mutect2) | FilterMutectCalls | Dedicated somatic filtering, not VQSR |

VQSR needs many variants overlapping the truth resources to fit a stable multivariate density; it is unreliable on single exomes or panels -- those must hard-filter. GATK is deprecating VQSR toward VETS (isolation-forest backend); verify the current recommended path in the GATK release notes before committing a pipeline.

## Edge cases

**Ploidy (`--sample-ploidy`).** The number of genotypes is the multiset coefficient C(ploidy + alleles - 1, ploidy), so PL vectors blow up in both ploidy and allele count -- the reason high-ploidy and pooled calling are memory-heavy.

| Case | Setting | Why |
|---|---|---|
| Pooled samples (n individuals) | `--sample-ploidy 2n` | Estimate an allele count, not an individual genotype; diploid collapses intermediate frequencies |
| Polyploid organism | `--sample-ploidy 4` (etc.) | Dosage genotypes (AAAB=0.25, AABB=0.5) cannot be represented as het |
| Non-PAR chrX/Y in a 46,XY sample | `--sample-ploidy 1` (or `--ploidy-regions` BED) | Hemizygous; diploid calling emits impossible "hets" from error/paralog mismap |
| PAR1/PAR2 | Diploid (mask PAR on Y, call X-PAR as diploid) | PARs recombine and are diploid in both sexes |

**Mitochondria -> Mutect2, not HaplotypeCaller.** mtDNA heteroplasmy is a continuous VAF (mathematically identical to subclonal somatic variation) that a diploid genotype model cannot express, so use a somatic caller. Run `gatk Mutect2 --mitochondria-mode` (raises low-AF sensitivity). Because rCRS (NC_012920.1, 16,569 bp circular) is linearized in the control region, align twice -- to the standard reference and to one shifted ~8,000 bp (`ShiftFasta`) that moves the artificial breakpoint out of the D-loop -- call control-region variants on the shifted reference, `LiftoverVcf` back, and merge. NUMT-derived reads inflate false low-heteroplasmy calls, so distrust calls below ~5% AF (this underpins the gnomAD mtDNA callset, Laricchia et al. 2022).

**Contamination is a gate before trusting any call.** Even 1-3% cross-sample contamination injects minority alleles that push allele balance far enough from 0/0.5/1 to be scored as low-fraction hets. Estimate it first: **VerifyBamID2** (Zhang et al. 2020, ancestry-agnostic, genotype-free -- models sample ancestry via a PCA/SVD panel, avoiding v1's population-mismatch bias), or **CHARR** (Lu et al. 2023, variant-level only -- needs just a gVCF/VCF, ~$0.0003/sample, from reference-allele leakage at hom-alt sites). Convention: FREEMIX >= 0.03 (3%) flags a probable contaminated/swapped sample (a guide, not a hard constant). For somatic work, feed `CalculateContamination`'s table into `FilterMutectCalls` (`--contamination-table`) so genuine low-fraction variants are separated from contamination artifacts.

## Parallelization

```bash
# Scatter HaplotypeCaller by contig, then gather
for interval in chr{1..22} chrX chrY; do
    gatk HaplotypeCaller -R reference.fa -I sample.bam -L $interval \
        -O sample.${interval}.g.vcf.gz -ERC GVCF &
done
wait
gatk GatherVcfs $(for c in chr{1..22} chrX chrY; do echo "-I sample.${c}.g.vcf.gz"; done) -O sample.g.vcf.gz

# PairHMM is the compute bottleneck: give it native SIMD threads
gatk HaplotypeCaller -R reference.fa -I sample.bam -O sample.vcf.gz --native-pair-hmm-threads 4
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Spurious heterozygous calls across non-PAR chrX/Y in a male | Called as diploid | `--sample-ploidy 1` for non-PAR X/Y (split intervals or `--ploidy-regions`) |
| mtDNA low-heteroplasmy variants missed or all filtered | HaplotypeCaller's diploid model cannot express continuous VAF | Use `Mutect2 --mitochondria-mode` + shifted reference |
| Excess false hets genome-wide, allele balance off 0.5 | Sample contamination / swap | Estimate with VerifyBamID2 or CHARR before trusting calls |
| Real variant not called in a repeat/high-depth region | Assembly failed (cyclic graph at all k, or region too complex/deep) | Compare against DeepVariant/DRAGEN; check `--max-assembly-region-size`, downsampling |
| VariantRecalibrator fails to converge / errors | Too few variants or too little truth-resource overlap | Hard-filter instead (single sample, exome, panel, non-model organism) |
| Indel mis-genotyped near a homopolymer/STR | Base-quality recalibration does not model repeat-context indel error | Use `--dragen-mode` (DRAGSTR STR-aware indel model) |
| A variant appears as `*/A` or `*/*` after joint genotyping | Spanning-deletion `*` allele: this position is inside an upstream deletion in some samples | Expected; `bcftools norm`-decompose and let annotators special-case `*` |

## Related Skills

- variant-calling/joint-calling - Cohort consolidation (GenomicsDBImport/CombineGVCFs), GenotypeGVCFs, and scaling to tens of thousands
- variant-calling/filtering-best-practices - Full VQSR/VETS and hard-filter recipes with threshold rationale
- variant-calling/deepvariant - CNN-based alternative; wins indels/difficult regions, ships platform-specific models
- variant-calling/variant-calling - bcftools pileup alternative (faster, less accurate on indels)
- variant-calling/variant-normalization - Left-align/decompose before annotation or benchmarking
- variant-calling/variant-annotation - Annotate final calls with VEP/SnpEff
- variant-calling/vcf-basics - View and query the resulting VCF
- read-alignment/bwa-alignment - Produce the aligned BAM HaplotypeCaller consumes

## References

- McKenna A, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research* 20:1297-1303 (2010). DOI 10.1101/gr.107524.110.
- DePristo MA, et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics* 43:491-498 (2011). DOI 10.1038/ng.806.
- Poplin R, et al. Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv* 201178 (2018). DOI 10.1101/201178. PREPRINT (GATK's recommended cite for the GVCF reference-confidence + joint-genotyping methodology).
- Van der Auwera GA, et al. From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline. *Current Protocols in Bioinformatics* 43:11.10.1-11.10.33 (2013). DOI 10.1002/0471250953.bi1110s43.
- Regier AA, et al. Functional equivalence of genome sequencing analysis pipelines enables harmonized variant calling across human genetics projects. *Nature Communications* 9:4038 (2018). DOI 10.1038/s41467-018-06159-4.
- Olson ND, et al. PrecisionFDA Truth Challenge V2: Calling variants from short and long reads in difficult-to-map regions. *Cell Genomics* 2:100129 (2022). DOI 10.1016/j.xgen.2022.100129.
- Behera S, et al. Comprehensive genome analysis and variant detection at scale using DRAGEN. *Nature Biotechnology* 43:1177-1191 (2025). DOI 10.1038/s41587-024-02382-1.
- Zhang F, et al. Ancestry-agnostic estimation of DNA sample contamination from sequence reads. *Genome Research* 30:185-194 (2020). DOI 10.1101/gr.246934.118. (VerifyBamID2.)
- Lu W, et al. CHARR efficiently estimates contamination from DNA sequencing data. *American Journal of Human Genetics* 110:2068-2076 (2023). DOI 10.1016/j.ajhg.2023.10.011.
- Laricchia KM, et al. Mitochondrial DNA variation across 56,434 individuals in gnomAD. *Genome Research* 32:569-582 (2022). DOI 10.1101/gr.276013.121.
