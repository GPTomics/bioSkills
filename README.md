# bioSkills

A collection of skills that guide AI coding agents (Claude Code, Codex, Gemini) through common bioinformatics tasks.

## Purpose

These skills help AI agents write correct, idiomatic code for bioinformatics workflows. Whether you're an undergrad learning computational biology, a biohacker running analyses, or a PhD processing large-scale data, these skills enable your AI assistant to help effectively.

## Tool Coverage

| Tool | Type | Use Case |
|------|------|----------|
| **Biopython** | Python | Sequence handling, file I/O, database access, phylogenetics |
| **Bioconductor** | R | RNA-seq, single-cell, pathway analysis, statistical methods |
| **samtools** | CLI | SAM/BAM/CRAM alignment file manipulation |
| **bcftools** | CLI | VCF/BCF variant calling and manipulation |
| **NCBI tools** | Mixed | Database queries, SRA downloads, BLAST searches |

## Installation

### Dependencies

**Python tools:**
```bash
pip install biopython pysam cyvcf2
```

**R/Bioconductor:** *(only needed for R-based skills)*
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

**CLI tools:**
```bash
# macOS
brew install samtools bcftools blast

# Ubuntu/Debian
sudo apt install samtools bcftools ncbi-blast+

# conda
conda install -c bioconda samtools bcftools blast
```

### Claude Code

```bash
git clone https://github.com/your-username/bioSkills.git
cd bioSkills

# Install globally (available in all projects)
./install-claude.sh

# Or install to a specific project
./install-claude.sh --project /path/to/your/project

# List available skills
./install-claude.sh --list
```

### Codex CLI

```bash
./install-codex.sh              # Install globally
./install-codex.sh --project /path/to/project
```

### Gemini CLI

```bash
./install-gemini.sh             # Install globally
./install-gemini.sh --project /path/to/project
```

Codex and Gemini installers convert to the Agent Skills standard (`examples/` -> `scripts/`, `usage-guide.md` -> `references/`).

## Available Skills

### sequence-io/ (complete)
Reading, writing, and converting biological sequence files using Biopython Bio.SeqIO.

| Skill | Description |
|-------|-------------|
| read-sequences | Parse FASTA, FASTQ, GenBank, ABI, SFF, PDB and 40+ formats |
| write-sequences | Write SeqRecord objects to sequence files |
| format-conversion | Convert between file formats |
| compressed-files | Handle gzip/bzip2/BGZF compressed files |
| fastq-quality | Analyze quality scores, filter by quality, convert encodings |
| filter-sequences | Filter sequences by length, ID, GC content, patterns |
| batch-processing | Process multiple files, merge, split, batch convert |
| sequence-statistics | Calculate N50, length/GC distributions, summary reports |
| paired-end-fastq | Handle R1/R2 pairs, interleave/deinterleave, sync filtering |

### sequence-manipulation/ (complete)
Working with sequence data programmatically using Biopython Bio.Seq and Bio.SeqUtils.

| Skill | Description |
|-------|-------------|
| seq-objects | Create and modify Seq, MutableSeq, SeqRecord objects |
| transcription-translation | DNA to RNA to protein, codon tables, ORF finding |
| reverse-complement | Reverse complement, palindrome detection, primer design |
| sequence-slicing | Slice, extract, concatenate sequences |
| motif-search | Find patterns with regex, PWM/PSSM, parse motif files |
| sequence-properties | GC content, GC skew, molecular weight, Tm, protein analysis |
| codon-usage | CAI, RSCU, codon optimization for expression |

### database-access/ (complete)
Fetching data from NCBI and other biological databases using Bio.Entrez, SRA toolkit, and BLAST+.

| Skill | Description | Tools |
|-------|-------------|-------|
| entrez-search | Search NCBI databases (ESearch, EInfo, EGQuery) | Bio.Entrez |
| entrez-fetch | Retrieve records by ID (EFetch, ESummary) | Bio.Entrez |
| entrez-link | Cross-database relationships (ELink) | Bio.Entrez |
| batch-downloads | History server, rate limiting, large datasets | Bio.Entrez |
| sra-data | Download sequencing reads (fasterq-dump, prefetch) | SRA toolkit |
| geo-data | Gene expression datasets, link to SRA | Bio.Entrez, GEOparse |
| blast-searches | Remote BLAST via NCBI | Bio.Blast.NCBIWWW |
| local-blast | Local BLAST databases and searches | BLAST+ CLI |

### alignment-files/ (complete)
Working with aligned sequence data (SAM/BAM/CRAM) using samtools and pysam.

| Skill | Description |
|-------|-------------|
| sam-bam-basics | View, convert between SAM/BAM/CRAM |
| alignment-indexing | Create BAI/CSI indices, indexed region access |
| alignment-sorting | Sort by coordinate/name, merge files, collate pairs |
| duplicate-handling | Mark and remove PCR/optical duplicates |
| alignment-statistics | Flagstat, depth, coverage, per-base stats |
| alignment-filtering | Filter by flags, quality, regions |
| reference-operations | Generate consensus, create dict files |
| pileup-generation | Generate pileup for variant calling |

### variant-calling/ (complete)
Variant calling and VCF/BCF manipulation using bcftools and cyvcf2.

| Skill | Description |
|-------|-------------|
| vcf-basics | View, query, understand VCF/BCF format |
| variant-calling | Call SNPs/indels from BAM files |
| vcf-filtering | Apply quality filters, expressions |
| vcf-manipulation | Merge, concat, sort, intersect VCFs |
| variant-normalization | Left-align indels, split multiallelics |
| variant-annotation | Add annotations, predict consequences |
| vcf-statistics | Generate quality metrics, concordance |
| consensus-sequences | Apply variants to reference FASTA |

### alignment/ (complete)
Sequence alignment using Biopython Bio.Align and Bio.AlignIO.

| Skill | Description |
|-------|-------------|
| pairwise-alignment | Global/local alignment (Needleman-Wunsch, Smith-Waterman) |
| alignment-io | Read, write, convert MSA files (Clustal, PHYLIP, Stockholm) |
| msa-parsing | Parse and analyze MSA content: gaps, conservation, consensus |
| alignment-statistics | Calculate identity, conservation, entropy, substitution patterns |

### phylogenetics/ (complete)
Phylogenetic tree analysis using Biopython Bio.Phylo.

| Skill | Description |
|-------|-------------|
| tree-io | Read, write, convert tree files (Newick, Nexus, PhyloXML, NeXML) |
| tree-visualization | Draw trees with matplotlib, customize labels and colors |
| tree-manipulation | Root, prune, ladderize, collapse, modify tree structure |
| distance-calculations | Compute distance matrices, build NJ/UPGMA trees, bootstrap consensus |

### differential-expression/ (complete)
RNA-seq differential expression analysis using R/Bioconductor (DESeq2, edgeR).

| Skill | Description |
|-------|-------------|
| deseq2-basics | DESeq2 workflow: DESeqDataSet, normalization, testing, lfcShrink |
| edger-basics | edgeR workflow: DGEList, TMM normalization, glmQLFit, testing |
| de-visualization | MA plots, volcano plots, PCA, heatmaps with ggplot2/pheatmap |
| de-results | Filter significant genes, add annotations, export results |

### Future Categories

| Category | Description | Tools |
|----------|-------------|-------|
| structural-biology | Protein structure analysis | Bio.PDB |
| single-cell | scRNA-seq workflows | Seurat |
| pathway-analysis | GO/KEGG enrichment | clusterProfiler |
| restriction-analysis | Restriction enzyme operations | Bio.Restriction |

## Typical NGS Workflow

```
Raw Reads (FASTQ)
    |
    v
[sequence-io] -----> QC, filtering, format conversion
    |
    v
[Aligner: bwa/bowtie2/STAR] (external)
    |
    v
[alignment-files] --> Sort, index, mark duplicates, stats
    |
    v
[variant-calling] --> Call variants, filter, annotate
    |
    v
[database-access] --> Compare to known variants, annotate
    |
    v
[differential-expression / pathway-analysis] --> Downstream analysis
```

## Usage

Once skills are deployed, ask your agent naturally:

```
"Parse my FASTA file and show sequence lengths"
"Filter FASTQ reads with mean quality below 25"
"Translate this coding sequence to protein"
"Sort my BAM file by coordinate and index it"
"Call variants from this BAM against the reference"
"Search NCBI for all human BRCA1 sequences"
"Download the SRA run SRR1234567"
"Align these two protein sequences using BLOSUM62"
"Calculate the pairwise identity matrix for this alignment"
"Convert my Clustal alignment to PHYLIP format"
"Build a neighbor joining tree from this alignment"
"Root the tree using Mouse as outgroup"
"Draw this phylogenetic tree and save as PDF"
"Run DESeq2 on my count matrix comparing treated vs control"
"Create a volcano plot of my differential expression results"
"Extract genes with padj < 0.05 and |log2FC| > 1"
```

The agent will use the skill patterns to generate correct code.

## Project Structure

```
bioSkills/
├── README.md              # This file
├── CLAUDE.md              # Project guidance for Claude Code
├── planning.md            # Full roadmap and skill details
├── install-claude.sh      # Installation script for Claude Code
├── install-codex.sh       # Installation script for Codex CLI
├── install-gemini.sh      # Installation script for Gemini CLI
├── sequence-io/           # Sequence file operations (9 skills)
├── sequence-manipulation/ # Sequence manipulation (7 skills)
├── database-access/       # NCBI & databases (8 skills)
├── alignment-files/       # SAM/BAM/CRAM (8 skills)
├── variant-calling/       # VCF/BCF (8 skills)
├── alignment/             # Sequence alignment (4 skills)
├── phylogenetics/         # Phylogenetic trees (4 skills)
├── differential-expression/ # RNA-seq DE analysis (4 skills)
└── ...                    # Additional categories
```

Each skill directory contains:
```
skill-name/
├── SKILL.md           # Agent instructions (YAML frontmatter + markdown)
├── usage-guide.md     # Human documentation
└── examples/          # Sample scripts (.py, .R, .sh)
```

SKILL.md files include:
- `tool_type` - python, r, cli, or mixed
- `primary_tool` - Main package (Bio.SeqIO, samtools, etc.)
- `## Related Skills` - Cross-references to related skills

## Contributing

Skills should be:
- Discrete but not overly granular
- Biopython-first for Python tasks, Bioconductor when R's statistical capabilities are needed
- CLI tools documented with both command-line and Python wrapper alternatives
- Compatible with multiple AI agents
- Include examples and usage guides
