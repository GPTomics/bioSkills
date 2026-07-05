---
name: bio-variant-calling
description: Call germline SNPs and indels from a BAM/CRAM with bcftools mpileup and call, and select the right calling engine for the job. Use when generating a VCF from aligned reads, choosing between bcftools, GATK HaplotypeCaller, DeepVariant, and DRAGEN, setting ploidy for haploid/organelle/polyploid/sex-chromosome calling, or deciding whether pileup-based calling is good enough versus a local-reassembly caller for indels and difficult regions. Not for cohort joint genotyping (see variant-calling/joint-calling), GATK-specific workflows (see variant-calling/gatk-variant-calling), deep-learning calling (see variant-calling/deepvariant), or somatic/low-VAF detection.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: bcftools mpileup applies BAQ (per-Base Alignment Quality) by default; this is a real behavior that changes calls, not a nuisance flag (see The Governing Principle).

# Variant Calling from a BAM

**"Call SNPs and indels from my aligned reads"** -> Compute per-position genotype likelihoods from a BAM/CRAM against the reference, then call variant sites under a Bayesian model at the assumed ploidy.
- CLI (fast, position-based): `bcftools mpileup -f ref.fa in.bam | bcftools call -mv`
- CLI (reassembly, higher indel accuracy): GATK HaplotypeCaller (variant-calling/gatk-variant-calling)
- CLI (deep learning): DeepVariant (variant-calling/deepvariant)

This skill does the bcftools calling and is the engine-selection hub: it tells the agent when pileup calling is the right tool and when to hand off to a reassembly or deep-learning caller.

## The Governing Principle

There are two families of short-variant caller, and the choice between them is the single most consequential decision here.

- **Position-based genotype-likelihood callers** (bcftools mpileup|call, the old samtools/UnifiedGenotyper lineage) trust the aligner's per-read placement. At each reference position they tally the pileup, compute P(reads | genotype) per site, and call under a Bayesian model. Fast, transparent, no training data. But the mapper places each read greedily and independently, so an indel near a read end or inside a repeat is placed inconsistently across reads, and the per-position model cannot repair that. bcftools mitigates it with **BAQ** (per-Base Alignment Quality: base qualities near a likely misalignment are downweighted so a shaky column does not produce a confident false SNP), but BAQ suppresses false positives rather than reconstructing the true indel.
- **Local-reassembly / haplotype callers** (GATK HaplotypeCaller, DeepVariant, DRAGEN) discard the local alignment in an active region and re-derive it: assemble candidate haplotypes, realign every read to them (PairHMM or a learned model), then genotype. The indel is represented once, on the assembled haplotype, instead of (mis)placed per read. This is precisely why they beat pileup callers on indels, clustered variants, and difficult regions.

Consequence: bcftools is fine-to-excellent for **simple germline SNPs** and quick genome-wide scans, materially **weaker on indels and in low-complexity / segmental-duplication / MHC regions**, and not built for somatic low-VAF detection or scalable cohort joint calling. Pick the engine from the analysis, not from habit.

## Engine Selection (the decision that comes before any command)

Guidance, not dogma; on a production human pipeline, validate against current GIAB/GA4GH benchmarks (hap.py + vcfeval) before committing.

| Engine | Best when | Fails / weak when | Hand off to |
|--------|-----------|-------------------|-------------|
| **bcftools mpileup\|call** | Simple germline SNPs; non-model/organelle/microbial genomes (no training data, any ploidy); quick exploratory scans; low compute; small multi-sample sets | Indels in homopolymers/STRs; segdups, MHC, low-mappability; low-VAF somatic/mosaic; cohorts beyond ~100 samples | this skill |
| **GATK HaplotypeCaller** | Auditable open-source human WGS/WES; every parameter inspectable; the joint-calling/best-practices orthodoxy (GVCF -> GenomicsDBImport -> GenotypeGVCFs) | Lower indel/difficult-region accuracy than DeepVariant/DRAGEN; local assembly can abort in pathological high-depth/repeat regions | variant-calling/gatk-variant-calling; variant-calling/joint-calling |
| **DeepVariant** | Best open-source accuracy on **indels and difficult regions**; PacBio HiFi / ONT (platform-specific trained models); generalizes off one training sample | Needs the correct platform model (wrong model degrades accuracy); GPU helps; cohort merge needs GLnexus, not GenotypeGVCFs | variant-calling/deepvariant |
| **DRAGEN** | Maximum throughput on Illumina (FPGA, ~20-25 min/genome); leads difficult-to-map benchmarks (alt-aware mapping) | Proprietary/hardware- or license-gated; ML recalibrator trained on GIAB truth (benchmark-overfitting caveat) | vendor pipeline; `HaplotypeCaller --dragen-mode` for an open-source approximation |

Honest state of the field: **DeepVariant and DRAGEN lead on indels and difficult regions**; **GATK is the joint-calling and best-practices reference** everyone else is measured against; **bcftools wins on speed, simplicity, non-model organisms, and organelle/haploid calling**. On easy SNPs every modern caller exceeds F1 0.999, so a caller's headline SNP number is rarely the deciding factor - indels and hard regions are.

## bcftools mpileup + call

**Goal:** Detect germline SNPs and indels from aligned reads with the pileup-and-call pipeline.

**Approach:** Generate per-position genotype likelihoods with mpileup (BAQ on by default), pipe as uncompressed BCF into the multiallelic caller.

### Basic calling
```bash
bcftools mpileup -f reference.fa input.bam | bcftools call -mv -Oz -o variants.vcf.gz
bcftools index variants.vcf.gz
```

### Recommended single-sample pipeline
```bash
# -Ou between steps avoids VCF (de)serialization; -q/-Q drop poorly-supported reads/bases;
# -a requests the FORMAT tags downstream filtering needs (DP, allelic depths, strand-bias p)
bcftools mpileup -Ou -f reference.fa \
    -q 20 -Q 20 \
    -a FORMAT/DP,FORMAT/AD,FORMAT/SP \
    input.bam | \
bcftools call -mv -Oz -o variants.vcf.gz
bcftools index variants.vcf.gz
```

### Region-restricted and multi-sample calling
```bash
# Single region / BED targets
bcftools mpileup -f reference.fa -r chr1:1000000-2000000 input.bam | bcftools call -mv -Oz -o region.vcf.gz
bcftools mpileup -f reference.fa -R targets.bed input.bam | bcftools call -mv -Oz -o targets.vcf.gz

# Multiple BAMs (small cohorts only; see The Governing Principle for the scaling limit)
bcftools mpileup -f reference.fa sample1.bam sample2.bam sample3.bam | bcftools call -mv -Oz -o cohort.vcf.gz

# BAM list file: one path per line
bcftools mpileup -f reference.fa -b bams.txt | bcftools call -mv -Oz -o cohort.vcf.gz
```

## The mpileup / call flags that change results

| Stage | Flag | Effect |
|-------|------|--------|
| mpileup | `-f ref.fa` | Reference FASTA (required); must be the exact one used for alignment |
| mpileup | `-q INT` | Min mapping quality; `-q 20` drops ambiguously placed reads (paralog mismapping) |
| mpileup | `-Q INT` | Min base quality; `-Q 20` drops low-confidence base calls |
| mpileup | `-a LIST` | Extra FORMAT/INFO tags: `FORMAT/AD` (allelic depths), `FORMAT/DP`, `FORMAT/SP` (Phred strand-bias p), `FORMAT/ADF`/`ADR` (per-strand), `INFO/AD` |
| mpileup | `-d INT` | Max per-file depth (default 250); set to 3-4x expected mean coverage to avoid truncating high-coverage sites |
| mpileup | `-B` / `-E` | `-B` disables BAQ (more raw indel signal, more false SNPs near indels); `-E` recomputes BAQ on the fly (more sensitive, slower) |
| call | `-m` | Multiallelic caller - default, recommended for all new work |
| call | `-c` | Consensus caller - legacy; only for reproducing old pipelines |
| call | `-v` | Emit variant sites only (omit to emit all sites, e.g. for hom-ref confidence) |
| call | `-O z\|b\|u\|v` | Output: `z` bgzipped VCF, `b` BCF, `u` uncompressed BCF (piping), `v` VCF |
| call | `--ploidy` / `--ploidy-file` | Sample/region ploidy (below) |
| call | `-P FLOAT` | Mutation-rate prior (default 1.1e-3, human); lower for inbred lines, raise for diverse/outbred populations |

The multiallelic caller (`-m`) handles sites with several ALT alleles natively and is statistically superior; the consensus caller (`-c`) exists only for backward reproducibility.

## Ploidy: sample, organelle, and sex-chromosome calling

**Goal:** Match the caller's ploidy to the biology so genotypes are representable.

**Approach:** Set a scalar ploidy for uniform samples, or a ploidy file (or built-in preset) to vary ploidy by region and sex.

Wrong ploidy silently corrupts calls: calling a diploid as haploid halves heterozygous sensitivity; calling a haploid/hemizygous region as diploid manufactures false heterozygous calls from every error and paralog mismap.

```bash
# Haploid: bacteria, mitochondria (nuclear germline heteroplasmy caveat below), non-PAR chrX/chrY in a male
bcftools mpileup -f reference.fa input.bam | bcftools call -m --ploidy 1 -Oz -o haploid.vcf.gz

# Built-in human preset applies karyotype-aware sex-chromosome ploidy
bcftools call -m --ploidy GRCh38 ...

# Ploidy file: CHROM  FROM  TO  SEX  PLOIDY  (chrY absent in females -> 0)
#   chrX  1  -1  M  1
#   chrX  1  -1  F  2
#   chrY  1  -1  M  1
#   chrY  1  -1  F  0
#   *     1  -1  *  2
bcftools mpileup -f reference.fa input.bam | bcftools call -m --ploidy-file ploidy.txt -Oz -o sexaware.vcf.gz
```

Scope notes: true **mitochondrial heteroplasmy** is continuous-VAF (not 0/0.5/1) and is a somatic-shaped signal - a diploid or haploid genotype model cannot express it; use a somatic caller (GATK Mutect2 `--mitochondria-mode`) for real heteroplasmy work. **Polyploid/pooled** samples need `--ploidy N` set to the true copy number so dosage/allele-count is preserved rather than collapsed to het.

## After calling: the pipeline map

A raw caller VCF is not a finished callset. The standard downstream order:

1. **Normalize** - left-align and split multiallelics so identical variants have identical records: `bcftools norm -f reference.fa -m -any variants.vcf.gz -Oz -o norm.vcf.gz`. Do this before ANY comparison, annotation, or merge. See variant-calling/variant-normalization.
2. **Filter** - bcftools produces no VQSR/DL score, so apply quality/depth/strand hard filters (e.g. `QUAL`, `FORMAT/DP`, `SP`) suited to the depth and platform. See variant-calling/filtering-best-practices.
3. **Inspect / query** - counts, Ti/Tv, per-sample stats. See variant-calling/vcf-basics and variant-calling/vcf-statistics.

## Comparing callers honestly

If the point of choosing bcftools vs a reassembly caller is accuracy, compare them correctly - this is where naive analyses go wrong:

- **Normalize both callsets first** (`bcftools norm -f ref.fa -m -any`). Two VCFs can encode the identical haplotype with different records (indel placement in repeats, MNP vs split SNVs); un-normalized records mismatch spuriously.
- **Use haplotype-aware benchmarking, not `bcftools isec`.** A line-diff / `isec` on raw records overcounts both false positives and false negatives from representation alone. Score against a GIAB truth set with **hap.py + vcfeval** inside the confident-region BED (Krusche 2019), reporting SNVs and indels separately.
- **Stratify.** A genome-wide F1 hides the differences that matter - they live in indels-in-repeats, segdups, and MHC. Report per-region, not one headline number.

## Performance

**Goal:** Speed up calling on large inputs.

**Approach:** Pipe uncompressed BCF between stages, thread both tools, and shard by chromosome.

```bash
# Threaded, uncompressed-BCF pipe
bcftools mpileup -Ou -f reference.fa --threads 4 input.bam | \
    bcftools call -mv --threads 4 -Oz -o variants.vcf.gz

# Parallel by chromosome, then concatenate
for chr in chr1 chr2 chr3; do
    bcftools mpileup -Ou -f reference.fa -r "$chr" input.bam | \
        bcftools call -mv -Oz -o "${chr}.vcf.gz" &
done
wait
bcftools concat -Oz -o all.vcf.gz chr*.vcf.gz
bcftools index all.vcf.gz
```

## Difficult regions (know where pileup calling breaks)

- **Homopolymers / STRs** - the dominant indel false-positive source; slippage + mapping ambiguity + representation ambiguity all concentrate here. Validate indels in homopolymers >6 bp with a reassembly caller or manual review, or switch engines.
- **Segmental duplications / low mappability** - paralog reads pile up and manufacture false SNPs; sites with mean `MQ` <40 signal ambiguous mapping. `-q 20` helps; a reassembly/alt-aware caller helps more.
- **MHC and other hyper-polymorphic loci** - extreme divergence from the reference; expect reduced recall from any linear-reference caller.
- **High-depth regions** - set `-d` to 3-4x expected mean coverage; the default 250 truncates deep targeted panels and can bias likelihoods.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `no FASTA reference` | `-f` omitted | Add `-f reference.fa` |
| `[E::faidx] ... different number of sequences` / reference mismatch | mpileup reference != alignment reference | Use the exact FASTA the BAM was aligned to; compare `@SQ` in `samtools view -H` against `grep '^>' ref.fa` |
| No variants called | Coverage too low, `-q`/`-Q` too strict, empty/wrong BAM | Check `samtools depth`; relax `-q`/`-Q`; confirm reference build |
| False heterozygous calls everywhere on chrX/chrY (male) | Non-PAR sex chromosome called as diploid | Set `--ploidy 1` for non-PAR, or use a `--ploidy-file` / `--ploidy GRCh38` |
| Excess indel false positives in repeats | Position-based limitation, not a bug | Normalize + hard-filter; validate or recall indels with a reassembly caller |
| Downstream tools disagree on the same variant | Records not normalized | `bcftools norm -f ref.fa -m -any` before comparing/merging/annotating |

## Related Skills

- variant-calling/vcf-basics - View and query the resulting VCF
- variant-calling/variant-normalization - Left-align and split multiallelics before comparison
- variant-calling/filtering-best-practices - Hard-filter a bcftools callset (no VQSR/DL score)
- variant-calling/vcf-statistics - Ti/Tv, counts, and callset QC
- variant-calling/gatk-variant-calling - Local-reassembly calling with HaplotypeCaller and DRAGEN-GATK mode
- variant-calling/deepvariant - Deep-learning caller; best indel/difficult-region accuracy, long-read models
- variant-calling/joint-calling - Scalable cohort genotyping (GVCF workflow, GLnexus)
- alignment-files/pileup-generation - Alternative pileup generation
- read-alignment/bwa-alignment - Upstream mapping that determines calling quality

## References

- Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics* 27(21):2987-2993 (2011). DOI 10.1093/bioinformatics/btr509. (The mpileup genotype-likelihood model.)
- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. *GigaScience* 10(2):giab008 (2021). DOI 10.1093/gigascience/giab008. (bcftools mpileup/call/norm implementation.)
- DePristo MA, Banks E, Poplin R, Garimella KV, Maguire JR, Hartl C, et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics* 43(5):491-498 (2011). DOI 10.1038/ng.806. (Local-reassembly genotyping framework - the reassembly contrast.)
- Poplin R, Chang P-C, Alexander D, Schwartz S, Colthurst T, Ku A, et al. A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology* 36(10):983-987 (2018). DOI 10.1038/nbt.4235. (DeepVariant.)
- Krusche P, Trigg L, Boutros PC, Mason CE, De La Vega FM, Moore BL, et al. Best practices for benchmarking germline small-variant calls in human genomes. *Nature Biotechnology* 37:555-560 (2019). DOI 10.1038/s41587-019-0054-x. (hap.py/vcfeval, confident-region model, normalize-before-compare.)
