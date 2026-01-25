# Tree I/O - Usage Guide

## Overview

This skill handles reading, writing, and converting phylogenetic tree files. It supports all major tree formats including Newick, Nexus, PhyloXML, and NeXML.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Read this Newick tree file and show me the taxa"
- "Convert my Nexus tree to Newick format"
- "Parse this PhyloXML file and list the bootstrap values"

## Example Prompts

### Reading Trees
> "Read the tree from tree.nwk and show its structure"

> "Parse all trees from my MrBayes output file"

> "Load this Newick string and print an ASCII diagram"

### Writing Trees
> "Save this tree as a PhyloXML file with metadata"

> "Export my modified tree to Newick format"

### Format Conversion
> "Convert my Nexus file from MrBayes to Newick"

> "Transform this Newick tree to PhyloXML to add annotations"

### Tree Inspection
> "How many taxa are in this tree?"

> "Show me all the branch lengths in this tree"

> "Is this tree bifurcating?"

## What the Agent Will Do

1. Identify the tree format from file extension or content
2. Parse the tree file using appropriate reader
3. Extract requested information (taxa, branch lengths, etc.)
4. Convert formats if requested
5. Write output in the specified format

## Supported Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| Newick | .nwk, .tre | Universal exchange format |
| Nexus | .nex | MrBayes, PAUP output |
| PhyloXML | .xml | Annotated trees with metadata |
| NeXML | .nexml | Modern XML trees |

## Tips

- **Newick** is the most portable format - use for sharing trees
- **PhyloXML** preserves the most metadata - use for annotated trees
- **Nexus** files from MrBayes often contain multiple posterior trees
- Use `Phylo.parse()` for multiple trees, `Phylo.read()` for single trees
- Use `format(tree, 'newick')` to get tree as string without writing to file
- Branch lengths of `None` mean unspecified, not zero
