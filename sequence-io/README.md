# sequence-io

## Overview

Sequence file input/output operations using Biopython's Bio.SeqIO module. Handles reading, writing, and converting biological sequence files in 40+ formats including FASTA, FASTQ, GenBank, and specialized formats like ABI traces.

**Tool type:** python | **Primary tools:** Bio.SeqIO

## Skills

| Skill | Description |
|-------|-------------|
| read-sequences | Parse FASTA, FASTQ, GenBank, and 40+ formats; choose streaming vs in-memory vs on-disk index (parse/read/to_dict/index/index_db) |
| write-sequences | Write SeqRecord objects to sequence files, with the per-format field requirements (molecule_type, phred_quality) |
| format-conversion | Convert between formats with SeqIO.convert and re-encode FASTQ offsets, flagging lossy conversions (GenBank to FASTA drops features) |
| compressed-files | Read, write, and index gzip, bzip2, xz, and BGZF files; the BGZF-vs-gzip seekability asymmetry and virtual offsets |
| fastq-quality | Access Phred scores, filter/trim by quality, and decide the Sanger/Solexa/Illumina encoding before parsing |
| filter-sequences | Stream-filter sequences by length, ID, GC content, N content, or regex without desyncing paired reads |
| batch-processing | Process many files (count, merge, split, convert) memory-safely; choose SeqIO vs index_db vs pysam/pyfastx |
| sequence-statistics | Calculate N50/L50, auN, NG50/NGA50, length/GC distributions, and summary statistics |
| paired-end-fastq | Handle R1/R2 pairs, interleave/deinterleave, and filter both mates together with orphan routing |

## Example Prompts

- "Parse my FASTA file and show each sequence ID and length"
- "Read this GenBank file and extract the sequence"
- "Save these modified sequences to a new FASTA file"
- "Convert sequences.gb to FASTA format"
- "Read my gzipped FASTQ and count the reads"
- "Convert my FASTA to BGZF so I can index it"
- "Filter FASTQ reads with mean quality below 25"
- "What's the quality score distribution in my FASTQ?"
- "Keep only sequences longer than 500 bp"
- "Extract sequences matching the IDs in my list"
- "Count sequences in each FASTA file in the data folder"
- "Combine all FASTA files in this directory into one"
- "Calculate N50, L50, and auN for my assembly"
- "Compute NG50 against a 3.1 Gb genome size to compare two assemblies"
- "Which FASTQ encoding does this file use before I parse the qualities?"
- "Re-encode this old Illumina 1.3 FASTQ to standard Phred+33"
- "Show me the GC content distribution, excluding ambiguous bases"
- "Filter my paired FASTQ files, keeping pairs where both pass Q30 and routing orphans"
- "Interleave my R1 and R2 files into a single file"
- "Convert my plain .gz FASTA to BGZF so samtools faidx will index it"
- "Create a persistent index across all FASTA files in this directory"
- "Read my ABI trace file and extract the trimmed sequence"

## Requirements

```bash
pip install biopython
# Optional, for fast/indexed access to huge or BGZF files:
pip install pysam pyfastx   # plus samtools/bgzip on PATH for faidx/BGZF
```

## Related Skills

- **sequence-manipulation** - Work with sequences after reading (transcription, translation, GC content)
- **database-access** - Fetch sequences from NCBI before local processing
- **read-qc** - Quality control and preprocessing before alignment
- **alignment-files** - Process aligned reads after running an aligner
