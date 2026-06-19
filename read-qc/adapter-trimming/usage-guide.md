# Adapter Trimming - Usage Guide

## Overview
Sequencing adapters appear in reads only when the insert is shorter than the read length (read-through), so adapter content is a direct readout of the insert-size distribution. Adapter is foreign sequence with genuine base quality, so aligners may mis-anchor on it -- which is why adapter trimming is the one near-universal preprocessing step (quality trimming usually is not, because aligners soft-clip). Cutadapt gives precise, error-tolerant control; Trimmomatic adds palindrome mode for paired read-through; fastp trims paired adapters from the read overlap with no adapter sequence needed.

## Prerequisites
```bash
# Cutadapt
conda install -c bioconda cutadapt
# Trimmomatic
conda install -c bioconda trimmomatic
# fastp (all-in-one alternative)
conda install -c bioconda fastp
```

## Quick Start
Tell your AI agent what you want to do:
- "Remove TruSeq adapters from my paired-end reads"
- "Trim small-RNA adapters and keep only 18-30 nt inserts"
- "My FastQC adapter-content curve climbs toward the 3' end, fix it"
- "Remove amplicon primers from my reads"

## Example Prompts

### Standard trimming
> "Remove TruSeq adapters from my paired-end FASTQ files and discard reads shorter than 20 bp"

> "Trim Nextera adapters from my reads with Cutadapt"

### Short-insert / small-RNA
> "Trim the small-RNA 3' adapter, gate inserts to 18-30 nt, and drop reads with no adapter"

> "My inserts are short (cfDNA) so most reads read through; trim the adapter from both mates"

### Read-through and troubleshooting
> "Detect paired read-through without knowing the adapter sequence"

> "FastQC still shows adapter after trimming, raise the error tolerance and use the shared stem"

## What the Agent Will Do
1. Identify the adapter from the kit, FastQC overrepresented sequences, or read-through detection
2. Choose the tool (cutadapt for precision/small-RNA/amplicon, fastp for bulk PE overlap, Trimmomatic for palindrome)
3. Trim 3' adapter from both mates together, keeping pairs synchronized
4. Apply a minimum-length filter so adapter dimers (which collapse to ~0 bp) are dropped
5. Verify removal with post-trim FastQC

## Common Adapter Sequences

| Library | Adapter |
|---------|---------|
| Illumina TruSeq (shared stem) | `AGATCGGAAGAGC` |
| Nextera / Tn5 | `CTGTCTCTTATACACATCT` |
| Small RNA | `TGGAATTCTCGGGTGCCAAGG` |

## Tips
- The adapter-content curve in FastQC is an insert-size readout; a climb toward the 3' end means short inserts reading through.
- Trim adapter even when skipping quality trimming -- soft-clipping aligners handle low-quality tails but mis-anchor on adapter.
- Small-RNA inverts the rule: use `--discard-untrimmed` and a tight length gate, because a no-adapter read is junk.
- Always pair adapter trimming with a minimum-length filter so adapter dimers are removed.
- Trimmomatic steps run in command-line order; put ILLUMINACLIP first and MINLEN last.
- Trimmomatic palindrome mode drops the redundant R2 by default (keepBothReads False); add keepBothReads if a downstream tool needs both mates.
- On 2-color instruments, a high-quality poly-G tail is not adapter; use a poly-G-aware trim (cutadapt --nextseq-trim or fastp).

## Resources
- [Cutadapt Documentation](https://cutadapt.readthedocs.io/)
- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic)
- [Illumina Adapter Sequences](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html)

## Related Skills
read-qc/quality-reports - Read the adapter-content panel that triggers trimming
read-qc/quality-filtering - Quality and length filtering after adapter removal
read-qc/fastp-workflow - All-in-one adapter + quality trim with auto poly-G
read-qc/contamination-screening - k-mer removal of PhiX/vector/contaminant sequence
small-rna-seq/smrna-preprocessing - Full small-RNA adapter + length workflow
