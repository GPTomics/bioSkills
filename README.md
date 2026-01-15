# bioSkills

A collection of skills that guide AI coding agents (Claude Code, Copilot, Codex) through common bioinformatics tasks using Biopython and Bioconductor.

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

Use the install script to deploy skills:

```bash
# Clone the repository
git clone https://github.com/your-username/bioSkills.git
cd bioSkills

# Install globally (available in all projects)
./install.sh

# Or install to a specific project
./install.sh --project /path/to/your/project

# List available skills
./install.sh --list
```

Skills are auto-invoked by Claude Code when your request matches their description.

### Other Agents

*Installation guides for Copilot, Cursor, and Codex coming soon.*

## Available Skills

### sequence-io/
Reading, writing, and converting biological sequence files.

| Skill | Description |
|-------|-------------|
| read-sequences | Parse FASTA, FASTQ, GenBank, and 40+ formats |

*More skills coming soon*

## Usage

Once skills are deployed, ask your agent naturally:

```
"Parse my FASTA file and show sequence lengths"
"Read the GenBank file and extract CDS features"
"Filter FASTQ reads by quality score"
```

The agent will use the skill patterns to generate correct Biopython/Bioconductor code.

## Project Structure

```
bioSkills/
├── README.md              # This file
├── install.sh             # Installation script
├── planning.md            # Skills roadmap
├── sequence-io/           # Sequence file operations
│   └── read-sequences/
│       ├── SKILL.md       # Agent instructions (with YAML frontmatter)
│       ├── usage-guide.md # Human documentation
│       ├── tests.json     # Validation tests
│       └── examples/      # Sample files and scripts
└── ...                    # More skill categories
```

## Contributing

Skills should be:
- Discrete but not overly granular
- Biopython-first, Bioconductor when R's statistical capabilities are needed
- Compatible with multiple AI agents
- Include examples and test cases
