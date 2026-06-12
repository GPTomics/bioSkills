# Tree I/O - Usage Guide

## Overview

This skill reads, writes, and converts phylogenetic tree files, and chooses a parser that preserves the annotations an analysis depends on. The central idea is that a tree file is a lossy serialization: the topology almost always survives a conversion, but branch supports, posterior probabilities, HPD intervals, and per-branch rates do not, and which ones survive is decided by the parsing tool, not by the format string. The most consequential decision is therefore tool selection: reading a BEAST or MrBayes tree with Bio.Phylo or `ape::read.nexus` silently discards the credible intervals, while treeio (R) and DendroPy (Python) keep them.

A second, subtler trap is that a bare number in plain Newick has no fixed meaning: `(A,B)95` could be a bootstrap, a posterior, or a clade name, and IQ-TREE writes two values (`SH-aLRT/UFBoot`) into the same slot, so a careless parser truncates or mislabels support.

## Prerequisites

- Python: `pip install biopython` (general I/O), `pip install dendropy` (annotation-preserving reads).
- R alternative for annotated trees: `treeio` and `ape` (Bioconductor/CRAN).
- Conceptual: know what wrote the file. BEAST/MrBayes output is annotated Nexus; the annotations ride inside `[&...]` comments and only purpose-built readers parse them.

## Quick Start

Tell your AI agent what you want to do:
- "Read this Newick tree file and show me the taxa"
- "Convert my Nexus tree to Newick but keep the posterior probabilities"
- "Parse all trees from my MrBayes posterior file"
- "My BEAST MCC tree lost its HPD bars after conversion - how do I recover them?"
- "Are the numbers on my internal nodes bootstraps or clade names?"

## Example Prompts

### Reading and inspecting
> "Read the tree from tree.nwk and show its structure and branch lengths"
> "Parse all trees from my MrBayes .trees posterior and tell me how many there are"

### Preserving annotations
> "I need plain Newick for a downstream tool but my tree is a BEAST MCC with HPD intervals - extract the node ages and credible intervals first"
> "Read this BEAST tree and list every annotation it carries before I decide what to drop"

### Format conversion
> "Convert my Nexus file from MrBayes to Newick and tell me exactly what is lost"
> "Transform this Newick tree to phyloXML so I can add typed taxonomy"

### Disambiguating support
> "The internal labels on my IQ-TREE output are like 87.5/98 - split them into SH-aLRT and UFBoot"

## What the Agent Will Do

1. Identify the format and, more importantly, the software that wrote the file.
2. Choose a parser: Bio.Phylo/ape for plain topology and standard bootstraps; treeio/DendroPy when `[&...]` annotations matter.
3. For annotated trees, extract the posteriors/HPDs/rates into a side table before any down-conversion.
4. Convert only among formats of sufficient capability, and state what each conversion drops.
5. Verify support landed in the right slot (`.confidence` vs `.name`) and that tip labels survive a round-trip against any metadata table.

## Tips

- The tool, not the format, decides whether `[&...]` metadata survives. Never bridge a BEAST tree through Bio.Phylo or `ape::read.nexus` when the HPDs matter.
- Run a round-trip test (read, write, read) and diff the annotations you care about, not just the topology; the topology surviving gives false confidence.
- A `.trees` file is the full posterior; a `.tree`/`.mcc` is the single annotated summary. Use `Phylo.parse`/`TreeList` for the former, and do not confuse them.
- Sanitize tip names to `[A-Za-z0-9_.]` and remember the Newick underscore-space convention, which silently desyncs labels from a metadata join key.
- A support number is meaningless without its provenance. Preserve which test produced it (bootstrap, UFBoot, SH-aLRT, posterior), not just the value.

## Related Skills

- tree-manipulation - rooting, pruning, and collapsing where rooted/unrooted and polytomy choices bite
- tree-visualization - ggtree and ETE consume the annotations preserved here
- bayesian-inference - produces the BEAST/MrBayes annotated trees whose metadata must survive
- divergence-dating - produces MCC trees with HPD intervals on node ages
- sequence-io/read-sequences - taxon-name sanitization shares the whitespace and non-ASCII traps
