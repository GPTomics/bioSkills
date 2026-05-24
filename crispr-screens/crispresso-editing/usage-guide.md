# CRISPResso2 Editing Analysis - Usage Guide

## Overview

Decision-grade quantification of CRISPR editing outcomes with CRISPResso2 across Cas9 nuclease (indels + HDR), cytosine and adenine base editors (target conversion + bystander), and prime editor (templated edits) modes. Covers single-amplicon (`CRISPResso`), batch (`CRISPRessoBatch`), pooled (`CRISPRessoPooled`), WGS off-target (`CRISPRessoWGS`), and comparison (`CRISPRessoCompare`) workflows. Defines the quantification window math, the substitution-vs-indel ratio for distinguishing BE from Cas9 contamination, and the MMEJ deletion signature in allele tables.

## Prerequisites

```bash
pip install CRISPResso2          # primary tool
# or via conda
conda install -c bioconda crispresso2
# Visualization helpers
pip install pandas matplotlib seaborn
```

Required inputs:
- FASTQ files (paired-end recommended, single-end accepted)
- Amplicon reference sequence (primer-trimmed; do NOT include primer regions)
- Guide protospacer sequence (20 nt; PAM not included)
- For HDR: expected edited amplicon sequence
- For BE: target nucleotide conversion direction (C->T for CBE; A->G for ABE)
- For PE: pegRNA spacer + extension (PBS+RTT) + scaffold sequences

## Quick Start

Tell the AI agent what to analyze:
- "Quantify indel rate at my Cas9 cut site from amplicon sequencing"
- "Compute target editing % and bystander rate from my CBE experiment using quantification window size 10"
- "Run CRISPRessoBatch on 24 samples (timecourse: 0/6/12/24/48 hours) targeting BRCA1 exon 11"
- "Diagnose why my CRISPResso run has 35% alignment rate"
- "Distinguish BE-mediated edits from Cas9-contamination indels using substitution-vs-indel ratio"
- "Run CRISPRessoPooled on 96 arrayed-validation amplicons in a single MiSeq library"

## Example Prompts

### Cas9 Editing

> "Run CRISPResso on a single Cas9 sample at MLH1 exon 12. Amplicon: ACGT... (primer-trimmed). Guide: GUIDE_SEQ. Report % unmodified vs NHEJ; flag if indel rate <70% (incomplete KO)."

> "I have 50 samples from a Cas9 timecourse experiment in K562. Use CRISPRessoBatch with batch_settings.txt; parallelize across 16 processes. Output the aggregated CRISPRessoBatch_quantification_of_editing_frequency.txt for downstream plotting."

### Base Editor Analysis

> "Analyze my CBE sample at BRCA1 c.5135C>T target. CRISPResso with `--base_editor_output`, `--conversion_nuc_from C --conversion_nuc_to T`, quantification window size 10, center -10. Report: target conversion %, bystander rate at adjacent Cs, indel byproduct %, substitution-vs-indel ratio."

> "Compare ABE7.10 vs ABE8.20 editing efficiency at the same target site across 5 cell lines. Use CRISPRessoBatch; output a heatmap of % A->G at the target position."

> "Distinguish my BE3 sample's edits from background Cas9 indels: confirm substitution-vs-indel ratio >10 (clean BE) or <3 (Cas9-like cut activity); recommend whether to repeat with a tighter nCas9-BE3 vector."

### Prime Editor Analysis

> "Quantify PE-2 editing for installation of MLH1 c.677A>G. Specify spacer + extension (PBS+RTT) + scaffold sequences. Report intended-edit %, scaffold incorporation %, indel %; flag if scaffold incorporation >5% (RTT design issue)."

> "Compare 8 candidate pegRNAs for the same intended edit. Use CRISPRessoBatch to run on saturating-dose samples; report which pegRNA achieves highest intended/(scaffold+indel) ratio."

### Pooled-Amplicon Mode

> "Run CRISPRessoPooled on a 96-amplicon arrayed-validation pool. Verify amplicon-misassignment is <2% by checking the primer-overlap diagnostic."

### Comparisons

> "Compare CRISPResso outputs for control vs DNA-damage-treated samples using CRISPRessoCompare. Quantify shift in indel size distribution and MMEJ-like deletion enrichment."

### Diagnostics

> "My CRISPResso alignment rate is 38%. Diagnose: wrong amplicon sequence, primer-dimer contamination, low-quality FASTQ, or strand orientation issue?"

> "Allele table shows a -7 bp deletion in 35% of reads at the cut site. Distinguish MMEJ-mediated repair from random NHEJ; recommend follow-up validation."

## What the Agent Will Do

1. Identify experimental design: nuclease vs base editor vs prime editor; single vs multi-sample vs pooled
2. Pick the correct CRISPResso mode from the decision tree
3. Verify amplicon sequence is primer-trimmed (primers excluded)
4. Verify guide sequence matches amplicon orientation
5. Set quantification window: size 1 for precision Cas9 (default); size 10 for base editor coverage; widened for prime editor template
6. For BE: set `--conversion_nuc_from`/`--conversion_nuc_to` per chemistry
7. For PE: include pegRNA spacer + extension + scaffold parameters
8. Run with `--min_average_read_quality 30` for stringent quality filtering
9. Verify mapping rate >85%; if not, diagnose amplicon sequence / contamination
10. Parse output: editing quantification, allele-frequency table, per-position nucleotide table
11. For BE: report target conversion + bystander rate separately; check substitution-vs-indel ratio
12. For PE: report intended-edit / scaffold-incorporation / indel %
13. Detect MMEJ patterns in allele table (recurring same-size deletions)
14. Output editing-efficiency summary, allele-table flagged hits, recommended follow-up validation

## Tips

- The single most common silent failure is including primers in `--amplicon_seq`. CRISPResso aligns reads to the amplicon; if primers are included, the alignment region differs from the sequenced region. Always trim primers from the amplicon parameter.
- Quantification window size 1 (Cas9) vs 10 (BE) vs wider (PE) is not optional. Mis-sized windows under- or over-count edits; consult the CRISPResso2 paper (Clement 2019) for canonical settings.
- BE samples showing high indel rates (>5%) usually indicate Cas9 / nCas9 expression mismatch; the substitution-vs-indel ratio is the key diagnostic. A ratio <3 = use Cas9-like analysis; >10 = clean BE.
- For BE variant-function screens, always report target-conversion AND bystander rates per position. Bystander edits are not noise -- they are real edits at adjacent bases that may have biological consequences.
- Pooled-amplicon mode (`CRISPRessoPooled`) silently misassigns reads when amplicons share primer regions. Always design pooled amplicons with ≥3 bp distinguishing flanks.
- For prime editor screens, scaffold incorporation indicates RTT design failure. Re-derive pegRNA with PRIDICT2 (see [[prime-editing-screens]]).
- Allele tables reveal MMEJ patterns: look for recurring deletions of specific sizes (e.g., 7 bp, 14 bp) at the cut site that share microhomology at junctions. These deletions are biologically distinct from random NHEJ and may have different functional consequences.
- For WGS off-target analysis, supply both the matched-tumor BAM and GUIDE-seq / CIRCLE-seq predicted sites; randomly chosen regions yield no useful comparison.

## Mode Decision Cheat Sheet

| Design | Mode | Key flags |
|--------|------|-----------|
| Single amplicon | `CRISPResso` | `--amplicon_seq --guide_seq` |
| Same amplicon, many samples | `CRISPRessoBatch` | `--batch_settings` |
| Many amplicons in pool | `CRISPRessoPooled` | `--amplicon_file` |
| Off-target from BAM | `CRISPRessoWGS` | `--bam --regions_file` |
| Compare two runs | `CRISPRessoCompare` | `--crispresso_output_folder_1 _2` |
| HDR | Add `--expected_hdr_amplicon_seq` to CRISPResso | |
| CBE | Add `--base_editor_output --conversion_nuc_from C --conversion_nuc_to T` | |
| ABE | Add `--base_editor_output --conversion_nuc_from A --conversion_nuc_to G` | |
| Prime editor | Add `--prime_editing_pegRNA_*` | |

## Thresholds

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| Functional Cas9 KO | >70% indels | Below this, KO incomplete |
| Clean BE indel byproduct | <5% | Above this = Cas9-like activity |
| Substitution-vs-indel ratio (BE) | >10 = clean BE; <3 = Cas9 contamination | Diagnostic for BE purity |
| Target conversion (BE) | >30% for screen power | Variable by target site and chemistry |
| Intended-edit (PE) | >5% per-edit; >20% at favorable sites | Anzalone 2019 baseline |
| Scaffold incorporation (PE) | <2% | Above = RTT design issue |
| Mapping rate | >85% | Below indicates amplicon design issue |
| Read quality | Phred 30 minimum | Joung 2017 |
| Read depth per sample | 1,000+ for reliable allele table | Higher for low-frequency variants |

## Related Skills

- crispr-screens/base-editing-analysis - Variant-function analysis using CRISPResso2 BE outputs
- crispr-screens/prime-editing-screens - PRIDICT2 pegRNA design + PE-tiling
- crispr-screens/library-design - sgRNA / pegRNA design rules
- crispr-screens/screen-qc - Editing-efficiency QC threshold for variant interpretation
- variant-calling/variant-annotation - Annotate edited variants downstream
- read-alignment/bwa-alignment - For WGS off-target alignment input
