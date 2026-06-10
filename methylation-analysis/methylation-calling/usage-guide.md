# Methylation Calling - Usage Guide

## Overview
Methylation calling extracts a per-cytosine methylation level, the beta value M/(M+U), from an aligned bisulfite or EM-seq BAM. The extractor counts reads that show a retained C (methylated) versus a converted T (unmethylated) at each reference cytosine and divides. Two extractor families exist and the choice is dictated by the aligner: bismark_methylation_extractor reads a Bismark BAM (it consumes Bismark's own XM call string), while MethylDackel and BISCUIT recompute calls from BAM plus reference and work on bwa-meth BAMs, which carry no XM tag. The output is a coverage file, a bedGraph track, or a genome-wide cytosine report spanning CpG (and, for plants, CHG/CHH) context. Because a beta is a count ratio rather than a measurement, the trustworthiness of every number depends on conversion efficiency, variant-awareness, overlap handling, and the 5mC vs 5hmC conflation.

## Prerequisites
```bash
conda install -c bioconda bismark methyldackel samtools
# Optional, aligner-dependent:
conda install -c bioconda biscuit   # SNP-aware / allele-specific methylation

bismark_methylation_extractor --version
MethylDackel --version
```

Conceptual prerequisites:
- A deduplicated, aligned bisulfite/EM-seq BAM (WGBS or PBAT deduplicated; RRBS NOT deduplicated). See bismark-alignment.
- The genome FASTA used for alignment (mandatory for the Bismark cytosine report and for MethylDackel/BISCUIT, which recompute calls against it).
- An unmethylated lambda spike-in is strongly recommended at library prep and is mandatory in plants, ESCs, and neurons, where non-CpG methylation is real biology and the CHH-as-conversion proxy fails.
- The extractor must match the aligner: a bwa-meth BAM has no XM tag and cannot be read by bismark_methylation_extractor.

## Quick Start
Tell your AI agent what you want to do:
- "Extract CpG methylation from my Bismark BAM and make a methylKit coverage file"
- "My WGBS was aligned with bwa-meth, extract methylation with MethylDackel"
- "Check the bisulfite conversion rate before I trust these betas"
- "Call methylation SNP-aware so C/T variants are not read as unmethylation"
- "Extract all contexts for my plant methylome"

## Example Prompts

### Conversion-rate QC first
> "Before extracting CpG betas, report my bisulfite conversion efficiency from the CHH rate, and tell me whether I need a lambda spike-in for this tissue."

### Bismark extraction
> "Run methylation extraction on my deduplicated paired-end Bismark BAM, merge the symmetric CpG dyads, produce a genome-wide cytosine report, and trim the R2 5' end-repair artifact based on the M-bias plot."

### bwa-meth / MethylDackel
> "This BAM was aligned with bwa-meth. Extract per-CpG methylation with MethylDackel, merge contexts, and exclude likely C/T SNP positions."

### Variant-aware and allele-specific
> "I am analyzing a human cohort. Call methylation jointly with SNPs so polymorphic CpGs are not miscalled, and phase the reads to test allele-specific methylation at known imprinted loci."

### Plant / non-CpG
> "Extract CpG, CHG, and CHH methylation for my Arabidopsis sample and report the conversion rate from the lambda spike-in, not from CHH."

## What the Agent Will Do
1. Confirm the aligner (Bismark vs bwa-meth) and pick the matching extractor.
2. Measure conversion efficiency from the CHH rate (mammals) or a lambda spike-in (plants/ESC/neurons) and gate on >=99% before trusting any beta.
3. Inspect the M-bias plot and set end-trim bounds from the flat region, not a fixed level.
4. Extract per-cytosine calls in the requested context, merging the symmetric CpG dyad where appropriate.
5. Apply variant-awareness for human/population data so a C/T SNP is not scored as unmethylation.
6. Produce coverage/bedGraph/cytosine-report output, noting the coordinate base (cov and cytosine report 1-based; bedGraphs 0-based), and hand count data (not pre-thresholded betas) to downstream testing.

## Output Files

| File | Content | Coordinate base | Downstream use |
|------|---------|-----------------|----------------|
| *.bismark.cov(.gz) | per-CpG M/U summary | 1-based | methylKit input |
| *.bedGraph(.gz) / *_CpG.bedGraph | methylation track | 0-based | IGV / UCSC |
| *.CpG_report.txt(.gz) | every genome cytosine | 1-based | bsseq input |
| *_splitting_report.txt | per-context % methylated | n/a | conversion-rate QC |
| BISCUIT pileup VCF / epiBED | joint methylation+SNP / per-read epialleles | per spec | SNP-aware betas, ASM |

## Tips
- The single most important calling-stage number is conversion efficiency. Require >=99%; in mammals read it from the CHH rate, but use a lambda spike-in in plants/ESC/neurons where CHH is real methylation.
- A standard WGBS/EM-seq beta is 5mC+5hmC, not 5mC. Label it as such; resolving 5hmC requires a second chemistry (oxBS/TAB/ACE), not an extractor flag.
- For a human cohort, call SNP-aware (BISCUIT/Bis-SNP or MethylDackel --maxVariantFrac); a C/T SNP otherwise masquerades as unmethylation and can fabricate DMRs.
- Symmetric CpG collapse is coverage2cytosine --merge_CpG (Bismark) or --mergeContext (MethylDackel). Note that --merge_non_CpG merges the CHG+CHH files instead and does NOT collapse the CpG dyad.
- --no_overlap is on by default for paired-end Bismark; the real risk is running single-end mode on paired data. MethylDackel handles overlap automatically.
- Read the M-bias plot for FLATNESS, not for a target level; trim the unstable read ends with --ignore/--ignore_r2 (Bismark) or --OT/--OB bounds (MethylDackel). EM-seq/PBAT usually needs --ignore_r2 2.
- Do not pre-threshold to a single beta and t-test it; hand the count data forward so the DMR callers can model the beta-binomial.
- Collapsing to a coverage matrix discards read-level epiallele, strand (hemimethylation), and allele (ASM) information. Keep reads (epiBED/raw BAM) when heterogeneity, clonality, cfDNA deconvolution, or allele-specific methylation is the question.

## Related Skills

- bismark-alignment - Produces the aligned BAM consumed here
- differential-cpg-testing - Per-CpG statistical testing on the counts
- dmr-detection - Region-level methods downstream
- methylkit-analysis - methylKit import of the coverage/cytosine report
- long-read-sequencing/nanopore-methylation - Native long-read MM/ML modification calling
- causal-genomics/mendelian-randomization - mQTL / causal follow-up of allele-specific methylation
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
