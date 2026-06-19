# Quality Filtering - Usage Guide

## Overview
Quality filtering removes low-quality bases (trimming) and low-quality reads (filtering). The modern lesson is that local aligners soft-clip, so quality trimming is usually unnecessary before alignment and aggressive trimming actively distorts expression estimates and assemblies; trim lightly or not at all, and always pair any trimming with a minimum-length filter. Sliding-window trimming (Trimmomatic SLIDINGWINDOW, fastp --cut_right) cuts decayed 3' tails; per-read filtering (fastp -q/-u/-n, AVGQUAL) drops a junk read subpopulation.

## Prerequisites
```bash
# Trimmomatic
conda install -c bioconda trimmomatic
# fastp
conda install -c bioconda fastp
# Cutadapt
conda install -c bioconda cutadapt
```

## Quick Start
Tell your AI agent what you want to do:
- "Apply a light sliding-window quality trim and a minimum-length filter"
- "Drop reads whose average quality is too low"
- "Remove poly-G tails from my NovaSeq data"
- "Do I need to quality-trim before GATK variant calling?"

## Example Prompts

### Standard filtering
> "Apply sliding-window quality trimming at Q20 to my paired-end reads and drop reads below 36 bp"

> "Filter out the low-quality read subpopulation in my FASTQ files"

### Workflow-specific
> "I am running BWA-MEM then GATK with BQSR, should I quality-trim?"

> "Clean my reads for de novo assembly without over-trimming"

### Platform-specific
> "Remove poly-G artifacts from my NovaSeq reads with fastp"

## What the Agent Will Do
1. Read the FastQC quality plots to decide whether the problem is bad ends (trim) or bad reads (filter)
2. Choose trimming intensity by workflow (light/none for soft-clipping aligners and BQSR; moderate for assembly)
3. Apply window trimming and/or per-read filtering with the right tool
4. Always add a minimum-length filter so over-trimmed fragments do not mis-map
5. Confirm improvement with post-filter FastQC

## Choosing Parameters

### Trimming intensity by workflow
| Workflow | Quality trimming | Notes |
|----------|------------------|-------|
| Alignment-based DNA/RNA | Light or none | Aligner soft-clips; aggressive trim distorts results |
| GATK variant calling (BQSR) | None | BQSR recalibrates; trimming interferes |
| De novo assembly | Moderate (~Q20) | Low-Q errors corrupt the graph; avoid Q>30 |
| Variant calling without BQSR | Moderate | No recalibration safety net |

### Minimum length (always applied with trimming)
| Read length | Min length |
|-------------|------------|
| 150 bp | 36-50 |
| 100 bp | 30-36 |
| Long inserts | 50+ |

## Tips
- Aligners soft-clip low-quality tails, so trim lightly or not at all for alignment-based DNA/RNA; never quality-trim before GATK BQSR.
- Always pair quality trimming with a minimum-length filter; a short over-trimmed read mis-maps (Williams 2016).
- Quality filtering (drop whole reads) and quality trimming (cut bases) are different operations; pick by whether the issue is some bad reads or bad ends.
- On NovaSeq/NextSeq, binned quality scores make window thresholds coarse; do not port HiSeq-era thresholds blindly.
- Poly-G tails are high quality on 2-color chemistry, so a quality filter will not remove them; use poly-G trimming.
- Do not HEADCROP the first ~12 bp of RNA-seq to "fix" the wavy base-content plot; that is random-hexamer bias, not adapter.
- Use fastp `--cut_right` for Trimmomatic-style sliding-window behavior in a single pass.

## Resources
- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic)
- [fastp Documentation](https://github.com/OpenGene/fastp)
- [Cutadapt Documentation](https://cutadapt.readthedocs.io/)

## Related Skills
read-qc/adapter-trimming - Remove adapter before quality filtering
read-qc/quality-reports - Read the quality plots that motivate filtering
read-qc/fastp-workflow - All-in-one preprocessing in a single pass
read-alignment/bwa-alignment - Soft-clipping aligner that absorbs low-quality tails
read-alignment/star-alignment - Soft-clipping RNA aligner (light trimming preferred)
