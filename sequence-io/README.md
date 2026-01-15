# sequence-io

Sequence file input/output operations using Biopython's Bio.SeqIO module.

**Tool type:** python
**Primary tool:** Biopython Bio.SeqIO

## Workflow Context

```
Raw Reads (FASTQ)
    |
    v
[sequence-io] <--- YOU ARE HERE
    |              - Read/write sequence files
    |              - Filter by quality/length
    |              - Format conversion
    |              - Paired-end handling
    v
[Aligner: bwa/bowtie2/STAR]
    |
    v
[alignment-files] --> SAM/BAM processing
    |
    v
[variant-calling] --> VCF generation
```

This category handles the first stage of most bioinformatics pipelines: reading raw sequence data, quality filtering, and format conversion before alignment.

## Skills

| Skill | Description | Key Functions |
|-------|-------------|---------------|
| [read-sequences](read-sequences/) | Parse sequence files | `SeqIO.parse()`, `read()`, `index()`, `index_db()`, `to_dict()` |
| [write-sequences](write-sequences/) | Write sequences to files | `SeqIO.write()`, `record.format()` |
| [format-conversion](format-conversion/) | Convert between formats | `SeqIO.convert()` |
| [compressed-files](compressed-files/) | Handle .gz/.bz2/BGZF files | `gzip.open()`, `bz2.open()`, `bgzf.open()` |
| [fastq-quality](fastq-quality/) | Quality score operations | `letter_annotations['phred_quality']` |
| [filter-sequences](filter-sequences/) | Filter by criteria | Generator expressions + SeqIO |
| [batch-processing](batch-processing/) | Multi-file operations | `pathlib.Path.glob()` + SeqIO |
| [sequence-statistics](sequence-statistics/) | N50, length/GC stats | `statistics` module + SeqIO |
| [paired-end-fastq](paired-end-fastq/) | R1/R2 pair handling | `zip()` + SeqIO for paired iteration |

## Supported Formats

Bio.SeqIO supports 40+ formats. Most common:

| Format | Extension | Read | Write | Index |
|--------|-----------|------|-------|-------|
| FASTA | .fasta, .fa | Yes | Yes | Yes |
| FASTQ | .fastq, .fq | Yes | Yes | Yes |
| GenBank | .gb, .gbk | Yes | Yes | Yes |
| EMBL | .embl | Yes | Yes | Yes |
| Swiss-Prot | .dat | Yes | No | Yes |

### Specialized Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| ABI | .ab1 | Sanger sequencing traces |
| SFF | .sff | 454/Ion Torrent data |
| PDB | .pdb | Protein sequences from structures |
| SnapGene | .dna | Lab plasmid design files |

### Alignment Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| PHYLIP | .phy | Phylogenetic analysis |
| Clustal | .aln | ClustalW alignments |
| Stockholm | .sto | Rfam/Pfam alignments |
| NEXUS | .nex | PAUP/MrBayes |

## High-Performance Options

For large files requiring maximum throughput:

| Function | Use Case | Performance |
|----------|----------|-------------|
| `SimpleFastaParser` | Fast FASTA iteration | 3-6x faster |
| `FastqGeneralIterator` | Fast FASTQ iteration | 3-6x faster |
| `SeqIO.index()` | Random access (single file) | Low memory |
| `SeqIO.index_db()` | Random access (multi-file, persistent) | SQLite-backed |
| BGZF compression | Indexable compressed files | gzip + random access |

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
| Indexable compression | "Convert my FASTA to BGZF so I can index it" |
| Batch processing | "Count sequences in each FASTA file in the data folder" |
| Assembly stats | "Calculate N50 and other statistics for my assembly" |
| Paired-end data | "Filter my paired FASTQ files, keeping pairs where both pass Q30" |
| GC content | "Show me the GC content distribution for these sequences" |
| Extract by ID | "Extract sequences matching the IDs in my list" |
| Merge files | "Combine all FASTA files in this directory into one" |
| Large file index | "Create a persistent index for my 50GB FASTA file" |
| Sanger traces | "Read my ABI trace file and extract the trimmed sequence" |
| Quality encoding | "Convert this old Illumina FASTQ to standard Sanger encoding" |

## Requirements

```bash
pip install biopython
```

## Related Skills

- **sequence-manipulation** - Work with sequences after reading (transcription, translation, GC content)
- **database-access** - Fetch sequences from NCBI before local processing
- **alignment-files** - Process aligned reads after running an aligner
- **variant-calling** - Call and analyze variants from aligned reads
