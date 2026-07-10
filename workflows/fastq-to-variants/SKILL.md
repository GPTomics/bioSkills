---
name: bio-workflows-fastq-to-variants
description: Orchestrates the end-to-end germline short-variant pipeline from FASTQ to a filtered, normalized, benchmarked VCF, chaining QC/trim, BWA-MEM2 alignment, duplicate marking, optional BQSR, calling (bcftools/GATK HaplotypeCaller/DeepVariant/DRAGEN), normalization, site+genotype filtering, annotation, and hap.py/vcfeval benchmarking. Use when deciding the pipeline-wide reference-genome commitment (GRCh38 analysis set vs T2T, ALT/decoy handling), sequencing the steps in the defensible order (normalize BEFORE annotate, filter site- then genotype-level), choosing the calling engine and single-sample vs cohort joint-calling, picking a filtering strategy by cohort size, or benchmarking stratified within GIAB confident regions. Hands off mechanism to the variant-calling and read-alignment component skills; not a re-teach of any single step.
tool_type: cli
primary_tool: bcftools
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - read-alignment/bwa-alignment
  - alignment-files/alignment-sorting
  - alignment-files/duplicate-handling
  - variant-calling/variant-calling
  - variant-calling/joint-calling
  - variant-calling/variant-normalization
  - variant-calling/filtering-best-practices
  - variant-calling/variant-annotation
  - variant-calling/vcf-statistics
qc_checkpoints:
  - after_qc: "Q30 >85%, adapter content <1%"
  - after_alignment: "Mapping rate >95%, properly paired >90%"
  - after_dedup: "Duplication rate <30% for WGS, <50% for exome"
  - after_calling: "Ti/Tv ratio ~2.0-2.1 for WGS, ~3.0-3.3 for exome; dbSNP overlap >95% only after annotating the ID column"
---

## Version Compatibility

Reference examples tested with: BWA-MEM2 2.2.1+, GATK 4.5+, bcftools 1.19+, samtools 1.19+, fastp 0.23+, DeepVariant 1.6+, hap.py 0.3.15+, Ensembl VEP 111+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: GenotypeGVCFs defaults (`--max-alternate-alleles`, `--heterozygosity`, `--stand-call-conf`), the GATK hard-filter thresholds, and DRAGEN speed/accuracy figures drift by version/vendor; confirm in-tool and against current GIAB benchmarks before quoting.

# FASTQ to Variants Workflow

**"Call variants from my whole-genome or exome FASTQ files"** -> Chain QC/trim, alignment, duplicate marking, an engine-appropriate caller, normalization, filtering, annotation, and benchmarking into one filtered germline VCF.
- CLI: fastp -> bwa-mem2 -> samtools markdup -> (bcftools | gatk HaplotypeCaller | DeepVariant) -> bcftools norm -> filter -> VEP -> hap.py

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A germline pipeline is a chain of commitments, and the two that decide whether the callset is trustworthy are made at the seams between steps, not inside them.

1. **The reference is a pipeline-wide commitment made once and inherited by everything downstream.** The build and analysis set chosen at alignment (step 2) fix the coordinates of every later comparison: the dbSNP/ClinVar/gnomAD records annotation matches against, the truth BED benchmarking scores within, and the cohort other samples are joint-called with must all be the SAME build. Changing it later means redoing alignment, calling, normalization, and annotation. Decide before aligning a single read.
2. **The step order is not arbitrary; two orderings sink reviews.** (a) **Normalize BEFORE annotate.** An indel that is not left-aligned to the database's canonical position is one base off, so the annotator silently misses the ClinVar/gnomAD record and reports a pathogenic variant as novel-absent -- a patient-safety failure that throws no error. (b) **Never `bcftools merge` single-sample VCFs into a cohort.** Absence of a record is then read as homozygous reference, fabricating genotypes; joint-genotype per-sample gVCFs instead so "confident hom-ref" is distinguished from "no data."
3. **Filter site-level first, then genotype-level, then recompute cohort QC.** Setting genotypes to no-call (`./.`) changes missingness, HWE, and allele frequencies, so those metrics must be computed on the genotype-filtered matrix, not before.
4. **A single genome-wide F1 is nearly meaningless.** Callers agree on easy SNPs (F1 > 0.999); they diverge in indels-in-repeats, segmental duplications, and MHC. Benchmark stratified, within the GIAB confident region, or do not claim accuracy (Krusche 2019 *Nat Biotechnol* 37:555-560).

## Workflow overview

```
FASTQ
  | [1] QC & trim -------------------> fastp            (read-qc/fastp-workflow)
  v
  | [2] Align ----------------------> bwa-mem2          (read-alignment/bwa-alignment)
  v     ^-- reference commitment: GRCh38 analysis set / T2T, ALT/decoy handling
  | [3] Mark duplicates ------------> samtools markdup  (alignment-files/duplicate-handling)
  v
  | [4] (BQSR? optional on modern binned-quality instruments)
  v
  | [5] Call -----------------------> bcftools | GATK HaplotypeCaller | DeepVariant | DRAGEN
  v     ^-- single-sample OR per-sample gVCF -> joint-genotype (variant-calling/joint-calling)
  | [6] Normalize (BEFORE annotate) -> bcftools norm -m-any -f ref  (variant-calling/variant-normalization)
  v
  | [7] Filter: site-level THEN genotype-level         (variant-calling/filtering-best-practices)
  v
  | [8] Recompute cohort QC on the genotype-filtered matrix  (missingness, Ti/Tv, het/hom, excess-het HWE)
  v
  | [9] Annotate -------------------> VEP / SnpEff      (variant-calling/variant-annotation)
  v
  | [10] Benchmark & QC ------------> hap.py/vcfeval, bcftools stats  (variant-calling/vcf-statistics)
  v
Filtered, normalized, benchmarked VCF
```

## Reference genome: the pipeline-wide commitment

**Decision made once, before alignment; everything downstream inherits it.** Mechanism of building/indexing the reference lives in read-alignment/bwa-alignment; the deeper reasoning below is what a reviewer expects justified.

| Choice | Commit to it when | Consequence inherited downstream |
|--------|-------------------|----------------------------------|
| GRCh38 analysis set + decoys (hs38DH), ALT-aware (bwa-postalt) | Human germline, research or most clinical | Decoys soak up off-target reads; ALT-aware mapping recovers reads in MHC/segdup loci that ALT-unaware mapping force-fits to the primary, inflating false positives |
| GRCh38 analysis set, ALT-unaware (primary only) | Frozen clinical pipeline needing deterministic simplicity | Simpler and validated, but loses signal in ~5 Mb of ALT-bearing loci |
| Masked GRCh38 (false-duplication fix) | Calling in CBS, U2AF1, KCNE1B, KCNJ18 and other affected genes | Recovers reads whose mapQ collapsed across the phantom duplicate copy |
| T2T-CHM13 | Research needing segdups/centromeres/dark genes; maximum accuracy | Reveals variants in newly resolved regions and removes GRCh38 false-duplication artifacts, but no lossless liftover to GRCh37/38, so the entire annotation/interpretation stack must be revalidated (Nurk 2022 *Science* 376:44-53; Aganezov 2022 *Science* 376:eabl3533) |

The reason to fix this first: annotation databases, panel BEDs, benchmark truth sets, and any cohort a sample is joint-called with are all coordinate-specific. Mixing builds (e.g. normalizing to GRCh38 then annotating against a GRCh37 dbSNP) is a guaranteed silent miss. "We used GRCh38" is under-specified -- plain vs masked vs analysis-set-with-decoys materially changes results in named clinical genes.

## The canonical order and why

Each step assumes the previous; the order is defensible under review (canonical preprocessing order, filtering/representation practice).

1. **QC/trim** -- remove adapters and low-quality tails before they corrupt alignment and duplicate detection.
2. **Align** -- to the committed reference, with read groups (SM/ID/PL/LB); read groups are a hard GATK requirement.
3. **Mark duplicates** -- PCR/optical duplicates are not independent evidence; marking (not removing) lets the caller down-weight them. Skip for amplicon/UMI data.
4. **BQSR -- honestly optional on modern instruments.** BQSR corrected context/cycle-dependent miscalibration on 2010-era continuous-quality Illumina. NovaSeq/NovaSeq X emit ~4 quality bins, leaving little to recalibrate; callsets are largely unchanged with vs without it. DeepVariant explicitly recommends NOT running BQSR (its CNN learned the raw-quality error model); DRAGEN handles quality internally (and uses DRAGSTR for STR/indel error modeling). Keep it for GATK-HaplotypeCaller consistency if a frozen pipeline demands it; otherwise the modern indel-accuracy lever is STR-aware error modeling (`--dragen-mode`), not BQSR.
5. **Call** -- per-sample VCF, or per-sample gVCF (`-ERC GVCF`) if a cohort will be joint-genotyped.
6. **Normalize BEFORE annotate/compare** -- `bcftools norm -m-any -f ref.fa` (split multiallelics, then left-align + parsimony), against the SAME reference used for annotation. Add `-a` (atomize) only when the downstream database is decomposed. This is the most common real ordering bug: annotate-then-normalize attaches consequences to a non-canonical representation that fails to match the database. The binding constraint is normalize-before-ANNOTATE/COMPARE, not normalize-before-filter: GATK convention runs `VariantFiltration` on the raw multiallelic records (its annotations are computed on that representation), then normalizes -- `bwa_gatk_workflow.sh` does exactly that, while `bwa_bcftools_workflow.sh` normalizes first. Both are correct; neither annotates before normalizing.
7. **Filter site-level, then genotype-level** -- site filters (VQSR/hard/ML) decide whether a *site* is real; genotype filters (`GQ`/`DP`/allele-balance) decide whether an *individual genotype* is trustworthy. SNPs and indels are filtered separately (different error processes and truth resources).
8. **Recompute cohort QC on the genotype-filtered matrix** -- missingness, Ti/Tv, het/hom, excess-het HWE. Doing HWE before genotype filtering lets low-GQ garbage drive spurious deviation.
9. **Annotate** on the normalized (and, where consequence matters, haplotype-resolved) representation.
10. **Benchmark/validate** after transforms (`hap.py`, `bcftools stats`).

## Choosing the calling engine

Pipeline-level selection only; the mechanism and full decision table live in variant-calling/variant-calling. Hand off there to pick, then return here for chaining.

| Situation | Lean toward | Hand off to |
|-----------|-------------|-------------|
| Auditable open-source, large cohort, joint calling | GATK HaplotypeCaller GVCF -> GenomicsDBImport -> GenotypeGVCFs | variant-calling/gatk-variant-calling, variant-calling/joint-calling |
| Best indel/difficult-region accuracy, single or cohort | DeepVariant (+ GLnexus for cohorts) | variant-calling/deepvariant |
| Quick/exploratory, non-model organism, limited compute | bcftools mpileup + call | variant-calling/variant-calling |
| Maximum throughput on Illumina, hardware available | DRAGEN (or GATK `--dragen-mode` for the open equivalent) | variant-calling/variant-calling |

**Single-sample vs cohort is a chaining decision, not a caller feature.** For a cohort, emit per-sample gVCFs and joint-genotype them so a variant seen in one sample is evaluated in all (cohort rescue of low-coverage hets, squared-off genotype matrix). This is what makes the pipeline forward-compatible with new samples (the N+1 problem). Full mechanism: variant-calling/joint-calling.

## Primary path: BWA-MEM2 + bcftools

Fast, dependency-light, good for exploratory work and non-model organisms; weaker on indels in homopolymers than reassembly callers.

### Step 1: QC/trim with fastp

```bash
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o trimmed/sample_R1.fq.gz -O trimmed/sample_R2.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html qc/sample_fastp.html
```

### Step 2: Align with BWA-MEM2

Read groups are mandatory; add `-Y` (soft-clip supplementary) if structural-variant calling is downstream, and `-K 100000000` for thread-count-invariant output. Reference/analysis-set choice: read-alignment/bwa-alignment.

```bash
bwa-mem2 index reference.fa   # once
bwa-mem2 mem -t 8 -K 100000000 \
    -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA\tLB:lib1" \
    reference.fa trimmed/sample_R1.fq.gz trimmed/sample_R2.fq.gz \
  | samtools view -bS - > aligned/sample.bam
```

### Step 3: Mark duplicates

Strict order (samtools convention): collate (name) -> fixmate `-m` -> sort (coordinate) -> markdup. Detail: alignment-files/duplicate-handling.

```bash
# collate groups mates by name; fixmate -m adds the ms/MC tags markdup needs; markdup needs coordinate order.
# Do NOT coordinate-sort before fixmate, and do NOT markdup amplicon/PCR data (use UMIs there).
samtools collate -@ 8 -O -u aligned/sample.bam \
  | samtools fixmate -m -@ 8 -u - - \
  | samtools sort -@ 8 -u - \
  | samtools markdup -@ 8 - aligned/sample.markdup.bam
samtools index aligned/sample.markdup.bam
```

### Step 4: Call and normalize

```bash
# Single sample (mpileup passes MQ/BQ filters into the pileup)
# -a FORMAT/DP,FORMAT/AD + call -f GQ emit the per-sample DP/GQ the Step 5 genotype filter needs.
bcftools mpileup -Ou -f reference.fa -a FORMAT/DP,FORMAT/AD --max-depth 250 --min-MQ 20 --min-BQ 20 \
    aligned/sample.markdup.bam \
  | bcftools call -mv -f GQ -Oz -o variants/sample.vcf.gz

# Normalize BEFORE any annotation or cross-callset comparison, against the SAME reference
bcftools norm -m-any -f reference.fa -Oz -o variants/sample.norm.vcf.gz variants/sample.vcf.gz
bcftools index variants/sample.norm.vcf.gz
```

For multi-sample cohorts, bcftools can call several BAMs jointly, but the GATK/DeepVariant gVCF path is preferred at scale (variant-calling/joint-calling).

### Step 5: Filter (site then genotype)

```bash
# Site-level (bcftools flags rather than removes, so failures stay auditable)
bcftools filter -Oz -s LowQual \
    -e 'QUAL<20 || INFO/DP<10 || MQ<30' \
    -o variants/sample.siteflt.vcf.gz variants/sample.norm.vcf.gz

# Genotype-level: set low-confidence genotypes to no-call (NOT 0/0)
bcftools filter -Oz -S . \
    -e 'FMT/GQ<20 | FMT/DP<8' \
    -o variants/sample.filtered.vcf.gz variants/sample.siteflt.vcf.gz
bcftools index variants/sample.filtered.vcf.gz
```

## Alternative path: BWA-MEM2 + GATK HaplotypeCaller

Local reassembly + PairHMM; the auditable reference implementation, strong on indels. Full mechanism: variant-calling/gatk-variant-calling.

```bash
gatk CreateSequenceDictionary -R reference.fa
samtools faidx reference.fa

# DRAGEN mode: no BQSR, STR-aware indel model (DRAGSTR), improved QUAL calibration
gatk HaplotypeCaller -R reference.fa -I aligned/sample.markdup.bam \
    -O gvcf/sample.g.vcf.gz -ERC GVCF --dragen-mode

# Cohort: consolidate gVCFs then joint-genotype (NEVER bcftools merge single-sample VCFs)
gatk GenomicsDBImport --sample-name-map gvcf/map.txt \
    --genomicsdb-workspace-path genomicsdb -L intervals.bed
gatk GenotypeGVCFs -R reference.fa -V gendb://genomicsdb -O variants/cohort.vcf.gz
```

## Filtering strategy depends on cohort size

The site-level filter is chosen by cohort size, platform, and organism; genotype-level filtering is always applied on top. Full mechanism and thresholds: variant-calling/filtering-best-practices.

| Cohort / data | Site-level filter | Why |
|---------------|-------------------|-----|
| Large WGS cohort (~30+ jointly genotyped) | VQSR (or AS_VQSR for huge cohorts) | The Gaussian-mixture model needs tens of thousands of variants and truth-resource overlap to fit; unreliable below that |
| Single sample / small cohort | GATK hard filters or VETS/NVScoreVariants | VQSR is non-identifiable on few variants; a "converged" model on one exome is filtering on noise |
| Exome specifically | Hard filters (do NOT use DP as a VQSR annotation) | Capture-boundary coverage cliffs break the annotation manifold |
| Non-model organism | Hard filters or a bootstrapped truth set | No HapMap/Omni/Mills truth resources exist |
| DeepVariant / DRAGEN output | Use the caller's own calibration; do NOT re-apply GATK hard filters | Their error modes differ; classic annotations do not describe them |

SNPs and indels are filtered separately (different error processes, truth resources, abundance). Hard-filter starting points (SNPs `QD<2, FS>60, MQ<40, MQRankSum<-12.5, ReadPosRankSum<-8, SOR>3`; indels loosen `FS>200`, tighten `ReadPosRankSum<-20`) are lenient heuristics to tune, not universal truth. RankSum annotations are only defined at het sites -- a hand-written filter must treat a missing annotation as PASS, or every hom-alt site vanishes.

## Benchmarking the pipeline

**Goal:** a defensible accuracy statement, not a single number.

The only rigorous way to compare a callset to truth is haplotype-aware, stratified, and confined to the truth set's confident region. Two VCFs can encode the identical haplotype with different records, so a naive `bcftools isec`/line-diff overcounts errors; use `hap.py` wrapping the `vcfeval` engine, which replays variants onto the reference and matches at the haplotype level (Krusche 2019 *Nat Biotechnol* 37:555-560; GIAB truth, Zook 2019 *Nat Biotechnol* 37:561-566).

```bash
# Only meaningful when the sample IS a GIAB genome (HG001-HG007) with a truth VCF + confident BED.
# -f = confident/callable region BED (TP/FP/FN counted ONLY inside it; calls outside are UNK, not FP)
hap.py truth.vcf.gz query.norm.vcf.gz \
    -f HG002_confident.bed \
    -r reference.fa \
    -o bench/hg002 \
    --engine=vcfeval \
    --stratification stratification.tsv   # GIAB region BEDs: low-complexity, segdup, MHC, GC-extreme
```

Discipline that separates a senior benchmark from a naive one:
- **Refuse the global F1.** Report SNP and INDEL separately, and show the low-complexity/segmental-duplication/MHC rows explicitly -- hiding them behind an all-regions average is the most common soft cheat.
- **Do not benchmark an ML caller only on its training genome.** DeepVariant and DRAGEN-ML train on GIAB coordinates, so scoring on HG002 alone partly measures memorization; score on a held-out or semi-blinded sample (HG003/HG004).
- **"We found variants GIAB missed" is almost always a category error** -- calls outside the confident region, which GIAB declined to adjudicate, not accuracy. Real hard-region claims cite CMRG or an assembly-based benchmark.
- **When the sample is not a GIAB genome** (the usual case), there is no truth VCF; fall back to proxy QC -- Ti/Tv (~2.0-2.1 WGS, ~3.0-3.3 exome), dbSNP overlap, het/hom by ancestry, and trio Mendelian concordance if a family is available. These are proxies for the absence of a benchmark, not a substitute for one (variant-calling/vcf-statistics).

## QC checkpoints between steps

| After | Gate | Interpretation |
|-------|------|----------------|
| QC/trim | Q30 >85%, adapter <1% | DNA is typically higher quality than RNA |
| Alignment | Mapped >95%, properly paired >90% (`samtools flagstat`) | Low mapping rate: wrong reference or contamination |
| Dedup | Duplicates <30% WGS, <50% exome | High duplication: PCR over-amplification, low input |
| Calling | Ti/Tv ~2.0-2.1 WGS, ~3.0-3.3 exome; dbSNP overlap >95% | Ti/Tv sliding toward 0.5 (random) signals false-positive inflation -- filters too loose. `bcftools stats` counts known sites from the ID column, and callers leave it as `.`: run `bcftools annotate -c ID -a dbsnp.vcf.gz` first or the dbSNP-overlap line reads 0% regardless of quality |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Annotation reports a known pathogenic variant as novel/absent | Annotated before normalizing; indel one base off the database coordinate | `bcftools norm -m-any -f ref` against the SAME build as the annotation DB, BEFORE annotation |
| Cohort has impossible all-hom-ref genotypes at variant sites | Built the cohort by `bcftools merge` of single-sample VCFs | Emit per-sample gVCFs and joint-genotype (variant-calling/joint-calling) |
| Every hom-alt site filtered out | Hand-written filter treats missing RankSum as failing | Treat missing annotation as PASS; RankSum is defined only at het sites |
| VQSR "converged" on one exome but the callset is garbage | VQSR needs tens of thousands of variants across ~30+ samples | Use hard filters or VETS/NVScoreVariants for single samples/exomes |
| Spurious variants in CBS/U2AF1/KCNE1B | GRCh38 false duplications collapse mapQ | Use a masked GRCh38 or T2T-CHM13; commit the reference before calling |
| GATK error "sample ... has no read group" | Read groups omitted at alignment | Re-run `bwa-mem2 mem -R "@RG\t..."` (SM/ID/PL/LB) |
| Different variant counts from vt vs bcftools on the same data | vt decomposes MNPs by default, bcftools does not | Standardize ONE normalization tool + flags across every cohort compared (variant-calling/variant-normalization) |

## Pipeline map (hand-offs)

- read-qc/fastp-workflow -- QC/trim options and report interpretation
- read-alignment/bwa-alignment -- BWA-MEM2 parameters, read groups, ALT/decoy analysis set, reference indexing
- alignment-files/duplicate-handling -- the collate/fixmate/sort/markdup order and UMI cases
- variant-calling/variant-calling -- engine selection (bcftools vs GATK vs DeepVariant vs DRAGEN) and ploidy
- variant-calling/gatk-variant-calling -- HaplotypeCaller, GVCF, BQSR/DRAGSTR, edge cases
- variant-calling/deepvariant -- CNN calling, platform models, DeepTrio, GLnexus cohorts
- variant-calling/joint-calling -- per-sample gVCF -> cohort joint genotyping, the N+1 problem, scaling
- variant-calling/variant-normalization -- left-align/parsimony, multiallelic split, MNP decomposition
- variant-calling/filtering-best-practices -- VQSR vs hard vs ML by cohort size; site vs genotype filters
- variant-calling/variant-annotation -- VEP/SnpEff on the normalized representation
- variant-calling/vcf-statistics -- Ti/Tv, het/hom, contamination/relatedness identity QC

The complete runnable scripts for both paths are in this skill's examples/ (`bwa_bcftools_workflow.sh`, `bwa_gatk_workflow.sh`).

## Related Skills

- database-access/sra-data - Pull public FASTQ for reanalysis (ENA mirror or STRIDES cloud)
- database-access/ncbi-datasets-cli - Pull reference genome assembly via Datasets v2 CLI
- read-qc/fastp-workflow - Detailed QC options
- sequence-io/fastq-quality - Confirm the FASTQ quality encoding (Phred+33 vs legacy Phred+64/Solexa) before trimming public or pre-2011 data
- sequence-io/paired-end-fastq - Keep R1/R2 mates synchronized; independent per-mate filtering silently desyncs pairs
- read-alignment/bwa-alignment - BWA-MEM2 parameters, read groups, ALT/decoy analysis set, the dedup ordering
- alignment-files/duplicate-handling - Duplicate marking details
- variant-calling/variant-calling - Engine selection and bcftools calling options
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller and DRAGEN mode
- variant-calling/deepvariant - Deep-learning calling and GLnexus cohorts
- variant-calling/joint-calling - Cohort joint genotyping and scaling
- variant-calling/variant-normalization - Normalize before annotate/compare
- variant-calling/filtering-best-practices - VQSR, hard filters, VETS
- variant-calling/variant-annotation - Annotate variants with VEP
- variant-calling/vcf-statistics - Ti/Tv, het/hom, and identity QC

## References

- Krusche P, Trigg L, Boutros PC, et al. (GA4GH Benchmarking Team). Best practices for benchmarking germline small-variant calls in human genomes. *Nature Biotechnology* 37:555-560 (2019). DOI 10.1038/s41587-019-0054-x. Stratified haplotype-aware benchmarking (hap.py/vcfeval).
- Zook JM, McDaniel J, Olson ND, et al. An open resource for accurately benchmarking small variant and reference calls. *Nature Biotechnology* 37:561-566 (2019). DOI 10.1038/s41587-019-0074-6. GIAB truth set + confident regions.
- DePristo MA, Banks E, Poplin R, et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics* 43:491-498 (2011). DOI 10.1038/ng.806. GATK framework.
- Van der Auwera GA, Carneiro MO, Hartl C, et al. From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline. *Current Protocols in Bioinformatics* 43:11.10.1-11.10.33 (2013). DOI 10.1002/0471250953.bi1110s43.
- Poplin R, Ruano-Rubio V, DePristo MA, et al. Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv* 201178 (2018). DOI 10.1101/201178. Preprint only (never journal-published); the GVCF/joint-genotyping reference.
- Poplin R, Chang P-C, Alexander D, et al. A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology* 36:983-987 (2018). DOI 10.1038/nbt.4235. DeepVariant.
- Yun T, Li H, Chang P-C, et al. Accurate, scalable cohort variant calls using DeepVariant and GLnexus. *Bioinformatics* 36:5582-5589 (2020). DOI 10.1093/bioinformatics/btaa1081.
- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008 (2021). DOI 10.1093/gigascience/giab008.
- Nurk S, Koren S, Rhie A, et al. The complete sequence of a human genome. *Science* 376:44-53 (2022). DOI 10.1126/science.abj6987. T2T-CHM13.
- Aganezov S, Yan SM, et al. A complete reference genome improves analysis of human genetic variation. *Science* 376:eabl3533 (2022). DOI 10.1126/science.abl3533. Reference-choice variant-calling payoff.
