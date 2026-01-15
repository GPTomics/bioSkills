# bioSkills

A collection of skills that guide AI coding agents (Claude Code, Codex, Gemini) through common bioinformatics tasks using Biopython and Bioconductor.

## Purpose

These skills help AI agents write correct, idiomatic code for bioinformatics workflows. Whether you're an undergrad learning computational biology, a biohacker running analyses, or a PhD processing large-scale data, these skills enable your AI assistant to help effectively.

## Installation

### Dependencies

**Biopython (Python):**
```bash
pip install biopython
```

**Bioconductor (R):** *(only needed for R-based skills)*
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
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
git clone https://github.com/your-username/bioSkills.git
cd bioSkills

# Install globally
./install-codex.sh

# Or install to a specific project
./install-codex.sh --project /path/to/your/project
```

### Gemini CLI

```bash
git clone https://github.com/your-username/bioSkills.git
cd bioSkills

# Install globally
./install-gemini.sh

# Or install to a specific project
./install-gemini.sh --project /path/to/your/project
```

Codex and Gemini installers convert to the Agent Skills standard (`examples/` -> `scripts/`, `usage-guide.md` -> `references/`).

### Other Agents

*Installation guides for Copilot and Cursor coming soon.*

## Available Skills

### sequence-io/
Reading, writing, and converting biological sequence files.

| Skill | Description |
|-------|-------------|
| read-sequences | Parse FASTA, FASTQ, GenBank, ABI, SFF, PDB and 40+ formats; index(), index_db(), low-level parsers |
| write-sequences | Write SeqRecord objects to sequence files |
| format-conversion | Convert between file formats |
| compressed-files | Handle gzip/bzip2/BGZF compressed files (BGZF supports indexing) |
| fastq-quality | Analyze quality scores, filter by quality, convert encodings |
| filter-sequences | Filter sequences by length, ID, GC content, patterns |
| batch-processing | Process multiple files, merge, split, batch convert |
| sequence-statistics | Calculate N50, length/GC distributions, summary reports |
| paired-end-fastq | Handle R1/R2 pairs, interleave/deinterleave, sync filtering |

### sequence-manipulation/
Working with sequence data programmatically.

| Skill | Description |
|-------|-------------|
| seq-objects | Create and modify Seq, MutableSeq, SeqRecord objects |
| transcription-translation | DNA to RNA to protein, codon tables, ORF finding |
| reverse-complement | Reverse complement, palindrome detection, primer design |
| sequence-slicing | Slice, extract, concatenate sequences |
| motif-search | Find patterns with regex, PWM/PSSM, parse motif files (JASPAR, MEME) |
| sequence-properties | GC content, GC skew, molecular weight, Tm, protein analysis |
| codon-usage | CAI, RSCU, codon optimization for expression |

## Usage

Once skills are deployed, ask your agent naturally:

```
"Parse my FASTA file and show sequence lengths"
"Read the GenBank file and extract CDS features"
"Filter FASTQ reads by quality score"
"Translate this coding sequence to protein"
"Calculate the CAI for E. coli expression"
"Find all occurrences of GAATTC in my sequence"
```

The agent will use the skill patterns to generate correct Biopython/Bioconductor code.

## Project Structure

```
bioSkills/
├── README.md              # This file
├── CLAUDE.md              # Project guidance for Claude Code
├── install-claude.sh      # Installation script for Claude Code
├── install-codex.sh       # Installation script for Codex CLI
├── install-gemini.sh      # Installation script for Gemini CLI
├── sequence-io/           # Sequence file operations (9 skills)
│   ├── read-sequences/
│   ├── write-sequences/
│   ├── format-conversion/
│   ├── compressed-files/
│   ├── fastq-quality/
│   ├── filter-sequences/
│   ├── batch-processing/
│   ├── sequence-statistics/
│   └── paired-end-fastq/
├── sequence-manipulation/ # Sequence manipulation (7 skills)
│   ├── seq-objects/
│   ├── transcription-translation/
│   ├── reverse-complement/
│   ├── sequence-slicing/
│   ├── motif-search/
│   ├── sequence-properties/
│   └── codon-usage/
│       ├── SKILL.md       # Agent instructions (YAML frontmatter + markdown)
│       ├── usage-guide.md # Human documentation
│       └── examples/      # Sample scripts
└── planning.md            # Roadmap for additional skill categories
```

## Contributing

Skills should be:
- Discrete but not overly granular
- Biopython-first, Bioconductor when R's statistical capabilities are needed
- Compatible with multiple AI agents
- Include examples and test cases
