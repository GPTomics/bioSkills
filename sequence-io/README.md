# sequence-io

Sequence file input/output operations using Biopython's Bio.SeqIO module.

## Skills

| Skill | Description | Key Functions |
|-------|-------------|---------------|
| [read-sequences](read-sequences/) | Parse sequence files | `SeqIO.parse()`, `read()`, `index()`, `to_dict()` |
| [write-sequences](write-sequences/) | Write sequences to files | `SeqIO.write()`, `record.format()` |
| [format-conversion](format-conversion/) | Convert between formats | `SeqIO.convert()` |
| [compressed-files](compressed-files/) | Handle .gz/.bz2 files | `gzip.open()`, `bz2.open()` + SeqIO |
| [fastq-quality](fastq-quality/) | Quality score operations | `letter_annotations['phred_quality']` |
| [filter-sequences](filter-sequences/) | Filter by criteria | Generator expressions + SeqIO |
| [batch-processing](batch-processing/) | Multi-file operations | `pathlib.Path.glob()` + SeqIO |
| [sequence-statistics](sequence-statistics/) | N50, length/GC stats | `statistics` module + SeqIO |
| [paired-end-fastq](paired-end-fastq/) | R1/R2 pair handling | `zip()` + SeqIO for paired iteration |

## Supported Formats

Bio.SeqIO supports 40+ formats. Most common:

| Format | Extension | Read | Write |
|--------|-----------|------|-------|
| FASTA | .fasta, .fa | Yes | Yes |
| FASTQ | .fastq, .fq | Yes | Yes |
| GenBank | .gb, .gbk | Yes | Yes |
| EMBL | .embl | Yes | Yes |
| Swiss-Prot | .dat | Yes | Yes |
| PHYLIP | .phy | Yes | Yes |
| Clustal | .aln | Yes | Yes |
| Stockholm | .sto | Yes | Yes |

## Example Prompts

Ask your AI agent naturally:

| Task | Example Prompt |
|------|----------------|
| Read sequences | "Parse my FASTA file and show each sequence ID and length" |
| Write sequences | "Save these modified sequences to a new FASTA file" |
| Convert format | "Convert sequences.gb to FASTA format" |
| Filter by length | "Keep only sequences longer than 500 bp" |
| Filter by quality | "Filter FASTQ reads with mean quality below 25" |
| Compressed files | "Read my gzipped FASTQ and count the reads" |
| Batch processing | "Count sequences in each FASTA file in the data folder" |
| Assembly stats | "Calculate N50 and other statistics for my assembly" |
| Paired-end data | "Filter my paired FASTQ files, keeping pairs where both pass Q30" |
| GC content | "Show me the GC content distribution for these sequences" |
| Extract by ID | "Extract sequences matching the IDs in my list" |
| Merge files | "Combine all FASTA files in this directory into one" |

## Requirements

```bash
pip install biopython
```
