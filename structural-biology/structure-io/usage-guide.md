# Structure I/O - Usage Guide

## Overview

This skill handles reading, writing, and downloading protein structure files. It supports PDB, mmCIF, MMTF, and BinaryCIF formats, with tools for downloading structures from RCSB PDB and extracting metadata.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Download structure 1ABC from PDB"
- "Parse this PDB file and show the chains"
- "Convert this mmCIF file to PDB format"

## Example Prompts

### Reading Structures
> "Read 1abc.pdb and tell me how many atoms it has"

> "Parse this mmCIF file and list the chains"

> "Load this structure and show the resolution"

### Downloading Structures
> "Download PDB 4HHB"

> "Get the biological assembly for 1ABC"

> "Download these 10 PDB structures"

### Writing Structures
> "Save only chain A to a new PDB file"

> "Convert this PDB to mmCIF format"

> "Write the structure without hydrogens"

### Extracting Metadata
> "What is the resolution of this structure?"

> "Show me all header information"

> "Extract the experimental method from this mmCIF"

## What the Agent Will Do

1. Identify the file format from extension or content
2. Use appropriate parser (PDBParser, MMCIFParser, etc.)
3. Download from RCSB if PDB ID provided
4. Extract requested information
5. Write output with optional filtering using Select class

## Supported Formats

| Format | Read | Write | Use Case |
|--------|------|-------|----------|
| PDB | Yes | Yes | Legacy, widely supported |
| mmCIF | Yes | Yes | Modern standard, full metadata |
| MMTF | Yes | No | Compact binary |
| BinaryCIF | Yes | No | RCSB recommended binary |

## Tips

- **Use mmCIF for new work** - PDB format is frozen and has limitations
- **MMCIF2Dict for metadata** - Use for detailed structure information
- **Select class for filtering** - Write subsets without modifying structure
- **QUIET=True** - Suppress parser warnings for cleaner output
- **Biological assemblies** - Use `assembly_num` parameter for functional units
