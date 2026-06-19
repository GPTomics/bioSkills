# RNA-seq QC - Usage Guide

## Overview
RNA-seq QC runs on the ALIGNED BAM plus a gene model and assesses metrics with no DNA analogue: library strandedness, gene-body 5'-3' coverage, exonic/intronic/intergenic distribution, rRNA/globin/mitochondrial rate, transcript integrity (TIN), and saturation. Two facts dominate: getting strandedness wrong silently halves or zeros counts, and non-UMI bulk RNA-seq must NOT be deduplicated (high-expression genes produce genuine identical-coordinate fragments).

## Prerequisites
```bash
conda install -c bioconda rseqc qualimap rna-seqc picard salmon sortmerna multiqc samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Determine the strandedness of my RNA-seq library before quantifying"
- "Assess gene-body coverage and RNA integrity (TIN) across my samples"
- "Measure the rRNA rate in my libraries"
- "Run full post-alignment RNA-seq QC and aggregate with MultiQC"

## Example Prompts

### Strandedness (do this first)
> "Infer the strandedness of my library with RSeQC and tell me the featureCounts and salmon settings"

> "My counts came out near zero, check whether I used the wrong strand setting"

### Integrity
> "Compute gene-body coverage and TIN to check for 3' degradation"

> "My samples are FFPE; which integrity metric should I use?"

### Distribution and contamination
> "Break down my reads into exonic/intronic/intergenic and tell me if there is gDNA contamination"

> "Measure the rRNA rate to check my depletion worked"

## What the Agent Will Do
1. Infer strandedness first and report the matching downstream tool flags
2. Run the RSeQC suite (read_distribution, geneBody_coverage, tin) on the BAM with a BED12 model
3. Run Picard CollectRnaSeqMetrics with STRAND_SPECIFICITY set to the inferred protocol
4. Interpret 3' bias, intronic/intergenic rate, and rRNA rate against the protocol (not absolute cutoffs)
5. Aggregate with MultiQC and flag cohort outliers; report duplication as a diagnostic, never dedup non-UMI data

## Strandedness decode
| Protocol | infer_experiment dominant | salmon -l | featureCounts -s | htseq |
|----------|---------------------------|-----------|------------------|-------|
| Unstranded | both ~0.5 | IU / U | 0 | no |
| Forward (fr-secondstrand) | 1++,1--,2+-,2-+ | ISF / SF | 1 | yes |
| Reverse (fr-firststrand, dUTP) | 1+-,1-+,2++,2-- | ISR / SR | 2 | reverse |

Single-end infer_experiment drops the read-number prefix: forward = "++,--", reverse = "+-,-+".

## Tips
- Infer strandedness empirically (infer_experiment.py or salmon -l A) before quantifying; never assume from the kit name. dUTP (the common modern protocol) is reverse / -s 2 / ISR.
- Do NOT mark or remove duplicates in non-UMI bulk RNA-seq; high-expression genes create genuine duplicate-coordinate fragments. Report duplication as a diagnostic only.
- RIN is a pre-prep estimate; gene-body coverage and TIN are the post-hoc truth. Use DV200 (not RIN) for FFPE.
- In a cohort with variable RNA quality, regress medTIN out as a covariate in the DE design rather than discarding samples.
- High intronic rate is contamination in bulk poly-A RNA-seq but SIGNAL in single-nucleus RNA-seq; do not apply a bulk gate to snRNA.
- RNA-SeQC requires a collapsed GTF; Picard PCT_* metrics are fractions (0-1), not percentages.

## Resources
- [RSeQC Documentation](http://rseqc.sourceforge.net/)
- [Qualimap](http://qualimap.conesalab.org/)
- [RNA-SeQC 2](https://github.com/getzlab/rnaseqc)
- [Picard RNA Metrics](https://broadinstitute.github.io/picard/)

## Related Skills
read-qc/quality-reports - Raw-FASTQ QC before alignment
read-qc/umi-processing - Molecule-accurate dedup for UMI RNA-seq
read-qc/contamination-screening - rRNA and cross-species contamination
read-alignment/star-alignment - Aligner that emits ReadsPerGene strandedness columns
rna-quantification/featurecounts-counting - Strand-aware quantification after QC
differential-expression/deseq2-basics - Use medTIN as a covariate in the design
