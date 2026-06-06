# Clair3 Variant Calling - Usage Guide

## Overview
Clair3 calls germline small variants (SNPs and small indels) from Oxford Nanopore and PacBio HiFi long reads using a two-stage deep-learning design: a fast pileup model handles most sites and a slow full-alignment model re-evaluates the hard ones on phased reads. The single most important operational fact is that the model is hand-picked and must match how the reads were basecalled - there is no auto-detection, and a mismatch silently degrades accuracy. ONT indels in homopolymers and short tandem repeats are the residual error mode that a single genome-wide F1 hides, so benchmarking must be stratified. Clair3 is germline-only; somatic, mosaic, trio, and RNA calling belong to the ClairS / Clair3-Trio / Clair3-RNA family.

## Prerequisites
```bash
conda install -c bioconda clair3 whatshap bcftools samtools
# hap.py (with RTG vcfeval) for benchmarking; GIAB truth + stratification BEDs
# The full ONT model set (every version, hac/fast, _with_mv) is in the rerio repo:
#   https://github.com/nanoporetech/rerio  ->  clair3_models/
```

## Quick Start
Tell your AI agent what you want to do:
- "Call germline variants from my ONT R10 BAM with the matching Clair3 model"
- "Which Clair3 model matches reads basecalled with Dorado sup v5.0.0?"
- "Phase my Clair3 calls and produce a haplotagged BAM"
- "Benchmark my Clair3 calls against GIAB HG002 with stratification"

## Example Prompts

### Model matching
> "My reads were basecalled with Dorado sup v5.0.0 on R10.4.1. Pick the matching Clair3 model, call germline variants from the BAM, and confirm the model is not a version mismatch."

### Non-human reference
> "Call variants from my bacterial Nanopore BAM against a draft assembly. Make sure Clair3 does not silently return nothing because the contigs are not named like human chromosomes."

### Phasing
> "Run Clair3 with phasing enabled using LongPhase, and give me a haplotagged BAM for allele-specific methylation."

### Stratified benchmarking
> "Benchmark my ONT Clair3 VCF against GIAB HG002 with hap.py and GIAB stratifications, and tell me the indel F1 inside homopolymer and low-complexity regions specifically."

### Boundary
> "I want to find low-VAF somatic variants in a tumor. Is Clair3 the right tool?"

## What the Agent Will Do
1. Read the basecaller model from the run metadata and pick the matching Clair3 model (chemistry, tier, version, `_with_mv` if applicable).
2. Run `run_clair3.sh` with the correct `--platform` and `--model_path`, adding `--include_all_ctgs` for non-human references.
3. Optionally enable phasing (WhatsHap default, LongPhase for speed/SV-awareness) and produce a haplotagged BAM.
4. Index and summarize `merge_output.vcf.gz` (or `phased_merge_output.vcf.gz` when phasing is enabled; `merge_output.vcf.gz` itself stays unphased).
5. Benchmark against GIAB with hap.py/vcfeval and GIAB stratifications, reporting homopolymer/STR strata separately.
6. Redirect somatic/mosaic/trio/RNA requests to the appropriate Clair3-family tool.

## Tips
- The model string is the experiment: match chemistry (r941 vs r1041), tier (fast/hac/sup), and basecaller version; pick the model version closest to but not above the basecaller version.
- There is no generic `ont`/`hifi` model folder - point `--model_path` at a specific model subfolder.
- `--include_all_ctgs` is mandatory for non-human/draft references, or Clair3 returns near-empty output.
- `--min_coverage` defaults to 2 - that is a floor, not a recommendation; aim for ~20-60x.
- Never report only genome-wide F1; stratify with GIAB BEDs and use CMRG for clinical genes.
- Clair3 is germline-only and not VAF-aware; for tumors use ClairS (paired) or ClairS-TO (tumor-only).
- Clair3 v2 needs PyTorch `.pt` models; v1 TensorFlow models will not load.
- `bcftools merge` on gVCFs is not joint genotyping; use GLnexus for cohort calling.

## Related Skills

- basecalling - Basecaller model+version the Clair3 model must match
- long-read-alignment - Produces the input BAM
- haplotype-phasing - Phasing/haplotagging Clair3 uses internally and can emit
- medaka-polishing - medaka diploid calling is deprecated in favor of Clair3
- structural-variants - SVs (out of Clair3 scope)
- variant-calling/deepvariant - DeepVariant native ONT/HiFi models
- variant-calling/vcf-statistics - Summarize/filter the VCF

## Resources
- [Clair3 GitHub](https://github.com/HKU-BAL/Clair3)
- [Clair3 models (rerio)](https://github.com/nanoporetech/rerio)
- [GIAB benchmarks and stratifications](https://www.nist.gov/programs-projects/genome-bottle)
