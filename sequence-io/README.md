# sequence-io

## Overview

Sequence file input/output operations using Biopython's Bio.SeqIO module. Handles reading, writing, and converting biological sequence files in 40+ formats including FASTA, FASTQ, GenBank, and specialized formats like ABI traces.

**Tool type:** python | **Primary tools:** Bio.SeqIO

## Skills

| Skill | Description |
|-------|-------------|
| read-sequences | Parse FASTA, FASTQ, GenBank, ABI, and 40+ formats with SeqIO.parse/read/index |
| write-sequences | Write SeqRecord objects to sequence files with SeqIO.write |
| format-conversion | Convert between formats with SeqIO.convert |
| compressed-files | Handle gzip, bzip2, and BGZF compressed files |
| fastq-quality | Analyze quality scores, filter by quality, convert encodings |
| filter-sequences | Filter sequences by length, ID, GC content, or regex patterns |
| batch-processing | Process multiple files, merge, split, batch convert |
| sequence-statistics | Calculate N50, length/GC distributions, summary statistics |
| paired-end-fastq | Handle R1/R2 pairs, interleave/deinterleave, synchronized filtering |

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
- "Calculate N50 and other statistics for my assembly"
- "Show me the GC content distribution for these sequences"
- "Filter my paired FASTQ files, keeping pairs where both pass Q30"
- "Interleave my R1 and R2 files into a single file"
- "Create a persistent index for my 50GB FASTA file"
- "Read my ABI trace file and extract the trimmed sequence"

## Requirements

```bash
pip install biopython
```

## Related Skills

- **sequence-manipulation** - Work with sequences after reading (transcription, translation, GC content)
- **database-access** - Fetch sequences from NCBI before local processing
- **read-qc** - Quality control and preprocessing before alignment
- **alignment-files** - Process aligned reads after running an aligner
