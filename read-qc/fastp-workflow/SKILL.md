---
name: bio-read-qc-fastp-workflow
description: Runs all-in-one FASTQ preprocessing with fastp in a single pass - adapter trimming via paired-end overlap analysis, quality/length filtering, 2-color poly-G removal, base correction, optional dedup/UMI/merge, and HTML/JSON reports. Use when preprocessing bulk Illumina data and wanting one fast tool instead of separate Cutadapt, Trimmomatic, and FastQC steps. For precise small-RNA/amplicon adapters use adapter-trimming; for molecule-accurate UMI dedup use umi-processing.
tool_type: cli
primary_tool: fastp
---

## Version Compatibility

Reference examples tested with: fastp 0.23+, FastQC 0.12+, MultiQC 1.21+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# fastp Workflow -- one C++ pass for adapter, quality, poly-G, and QC

Run adapter trimming, quality/length filtering, poly-G removal, and reporting in a single fast pass.

**"Preprocess my reads with fastp"** -> Trim adapters from the read overlap, filter low-quality reads, remove 2-color poly-G, and emit an HTML/JSON report.
- CLI: `fastp -i R1.fq.gz -I R2.fq.gz -o c_R1.fq.gz -O c_R2.fq.gz -h report.html -j report.json`

Scope: this skill OWNS general-purpose single-pass Illumina preprocessing. Precise small-RNA/amplicon/anchored adapters -> read-qc/adapter-trimming. Molecule-accurate UMI dedup/consensus -> read-qc/umi-processing. DNA coordinate dedup -> alignment-files/duplicate-handling. OUT OF SCOPE: transcriptome QC (read-qc/rnaseq-qc).

## The Single Most Important Modern Insight

1. **fastp trims paired-end adapters by OVERLAP ANALYSIS, needing no adapter sequence at all.** It aligns R1 against the reverse complement of R2, finds the insert-derived overlap, and trims whatever extends past it (the read-through region) -- so it can trim adapter down to a SINGLE trailing base, where sequence-matching tools need at least 3. The same overlap drives `--correction` (`-c`): where the mates disagree and one base is high-quality and the other very low, fastp overwrites the low-quality base with the high-quality call. This overlap machinery is why fastp is the default bulk PE preprocessor and why it needs no `--adapter_sequence` for standard libraries.

2. **`--dedup` is SEQUENCE-identity deduplication at the FASTQ level -- no coordinates, no UMI -- so it removes BIOLOGICAL duplicates too.** It cannot tell a PCR duplicate from a highly expressed transcript's fragment or a targeted amplicon. NEVER use `--dedup` for RNA-seq quantification, amplicon, or any assay where identical reads are genuine signal. For molecule-accurate removal use UMIs (read-qc/umi-processing); for DNA variant calling use coordinate-based dedup AFTER alignment (alignment-files/duplicate-handling). fastp `--dedup` is for the narrow case of removing exact-duplicate reads from a non-UMI library where that is known to be safe.

3. **Poly-G trimming auto-enables for 2-color instruments (NextSeq/NovaSeq) from the machine ID, because G is the no-signal call.** Leave it on; a high-quality poly-G tail is invisible to the quality filter. One fast pass does adapter + quality + poly-G + filtering + QC report, and the JSON feeds MultiQC -- but fastp does NOT replace cutadapt's precision for small-RNA 3' adapters, amplicon primers, or anchored/linked adapters.

## Tool Positioning

| Need | Use fastp? | Alternative |
|------|-----------|-------------|
| Bulk PE WGS/WES/RNA/cfDNA preprocessing | Yes (default) | -- |
| One pass: trim + filter + poly-G + QC report | Yes | -- |
| Small-RNA 3' adapter + tight length gate | No | cutadapt (read-qc/adapter-trimming) |
| Amplicon / anchored / linked primers | No | cutadapt |
| Molecule counting / ctDNA consensus | Extract only | umi_tools / fgbio (read-qc/umi-processing) |
| RNA-seq molecule dedup | No (`--dedup` is wrong) | UMIs, or do not dedup |

## Core Operations

```bash
# Single-end and paired-end basics
fastp -i in.fq.gz -o out.fq.gz
fastp -i R1.fq.gz -I R2.fq.gz -o c_R1.fq.gz -O c_R2.fq.gz

# Adapter: PE overlap is automatic; --detect_adapter_for_pe ADDS sequence-based detection on top
fastp -i R1.fq.gz -I R2.fq.gz -o c_R1.fq.gz -O c_R2.fq.gz --detect_adapter_for_pe
# Manual adapter sequences (SE auto-detects from data by default)
fastp -i in.fq.gz -o out.fq.gz --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# Quality FILTER (per-read): base <Q20 unqualified; drop if >40% unqualified or >5 Ns
fastp -i in.fq.gz -o out.fq.gz -q 20 -u 40 -n 5
# Quality TRIM (sliding window from 3', SLIDINGWINDOW analogue) + length gate
fastp -i in.fq.gz -o out.fq.gz --cut_right --cut_window_size 4 --cut_mean_quality 20 -l 36

# 2-color poly-G (auto for NextSeq/NovaSeq); poly-X for 3' poly-A etc.
fastp -i in.fq.gz -o out.fq.gz --trim_poly_g          # --poly_g_min_len 10 default
fastp -i in.fq.gz -o out.fq.gz --trim_poly_x

# Overlap base correction (PE only; high-Q mate fixes low-Q base)
fastp -i R1.fq.gz -I R2.fq.gz -o c_R1.fq.gz -O c_R2.fq.gz --correction

# Merge overlapping pairs (short inserts: cfDNA, small-RNA, aDNA). Produces THREE streams:
# the merged file is single-end (full insert) and un_R1/un_R2 stay paired -- align them separately
# (merged as SE, un_R1/un_R2 as PE) and combine the BAMs.
fastp -i R1.fq.gz -I R2.fq.gz --merge --merged_out merged.fq.gz -o un_R1.fq.gz -O un_R2.fq.gz

# UMI extraction: fastp moves the inline UMI out of the read before trimming; molecule-accurate
# dedup/consensus still happens AFTER alignment (umi_tools/fgbio), not in fastp
fastp -i R1.fq.gz -I R2.fq.gz -o c_R1.fq.gz -O c_R2.fq.gz --umi --umi_loc read1 --umi_len 8
```

Key flags: `-q` qualified quality (default 15), `-u` unqualified percent limit (40), `-n` N limit (5), `-e` average-quality filter (0=off), `-l` length required (15), `--length_limit` (0=off), `--cut_right/--cut_front/--cut_tail` window cut modes (off by default), `--cut_window_size` (4), `--cut_mean_quality` (Q20), `--thread/-w` (default 3), `-h/-j` HTML/JSON report.

## Complete Workflows

```bash
# Standard Illumina PE (4-color: HiSeq/MiSeq)
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz \
      --detect_adapter_for_pe --cut_right --cut_window_size 4 --cut_mean_quality 20 \
      -q 20 -l 36 -w 8 -h sample.html -j sample.json

# NovaSeq / NextSeq (2-color): add poly-G (auto, but explicit for clarity)
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz \
      --detect_adapter_for_pe --trim_poly_g \
      --cut_right --cut_window_size 4 --cut_mean_quality 20 -q 20 -l 36 -w 8 \
      -h sample.html -j sample.json

# RNA-seq: light trim only (aligner soft-clips; do NOT --dedup), longer min length
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o clean_R1.fq.gz -O clean_R2.fq.gz \
      --detect_adapter_for_pe -q 20 -l 50 -w 8 -h sample.html -j sample.json
```

## Parsing the JSON report

```python
import json

with open('sample.json') as f:
    report = json.load(f)

after = report['summary']['after_filtering']
print(f"reads kept: {after['total_reads']}, Q30: {after['q30_rate']:.2%}")
print(f"duplication: {report['duplication']['rate']:.2%}")    # diagnostic only -- do not auto-dedup
```

MultiQC parses fastp JSON directly: `multiqc .` over a directory of `*.json` builds the cohort report (read-qc/quality-reports).

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| RNA-seq counts deflated after fastp | Used `--dedup` (sequence dedup removes biological dups) | Drop `--dedup` for RNA-seq; never sequence-dedup expression data |
| Adapter not trimmed (SE) | SE has no overlap; relies on data auto-detect | Pass `--adapter_sequence` explicitly for SE |
| Poly-G remains | 4-color run, or auto-detect missed the instrument | Add `--trim_poly_g` explicitly |
| Small-RNA results poor | fastp overlap is not precise enough for ~22 nt inserts | Use cutadapt with `--discard-untrimmed` (read-qc/adapter-trimming) |
| Over-trimmed RNA-seq | Aggressive `--cut_right` quality | Light trim only; aligner soft-clips (read-qc/quality-filtering) |
| UMI dedup expected but none happened | `--umi` only EXTRACTS; dedup is post-alignment | Extract here, dedup with umi_tools/fgbio after mapping (read-qc/umi-processing) |

## References

Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17):i884-i890.
Chen S. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2(2):e107.
Ewels P, Magnusson M, Lundin S, Kaller M. 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32(19):3047-3048.

## Related Skills

read-qc/adapter-trimming - Precise adapter/primer control for small-RNA and amplicon
read-qc/quality-filtering - Detailed quality/length filtering options and the trim-light evidence base
read-qc/quality-reports - Aggregate fastp JSON across samples with MultiQC
read-qc/umi-processing - Molecule-accurate UMI dedup and consensus after alignment
alignment-files/duplicate-handling - Coordinate-based duplicate marking for DNA variant calling
