---
name: bio-phylo-tree-io
description: Read, write, and convert phylogenetic tree files with Biopython Bio.Phylo, and choose an annotation-preserving parser (treeio, DendroPy) when metadata matters. Covers why a tree file is a lossy serialization, why format conversion silently drops BEAST/MrBayes node annotations (posteriors, HPD intervals, rates), the Newick support-vs-label ambiguity that mislabels bootstrap values, and the Nexus TRANSLATE and rooted/unrooted traps. Use when parsing Newick, Nexus, NHX, phyloXML, or NeXML, converting between formats, handling posterior tree sets, or moving annotated BEAST trees without losing the credible intervals. Routes annotation-critical reads to DendroPy or treeio and orthology/alignment context to sibling skills.
tool_type: python
primary_tool: Bio.Phylo
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+. Annotation-preserving alternatives: DendroPy 5+ (Python), treeio 1.26+ / ape 5.8+ (R).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(module.function)` to check signatures
- R: `packageVersion('treeio')` then `?read.beast` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Bio.Phylo stores a Newick `[&...]` bracket as opaque `.comment` text and does NOT parse BEAST key-values; DendroPy `extract_comment_metadata=True` and treeio `read.beast` do.

# Tree I/O -- A Tree File Is a Lossy Serialization

**"Read and convert my tree files"** -> Parse a tree into an in-memory object and re-serialize it, knowing which annotations each format and parser preserves.
- Python: `Phylo.read('tree.nwk', 'newick')`, `Phylo.convert(...)` (Bio.Phylo)
- Annotation-critical: `dendropy.Tree.get(..., extract_comment_metadata=True)` or treeio `read.beast()`

Scope: reading, writing, converting, and inspecting tree files, and selecting a parser that keeps the annotations the analysis needs. Rooting, pruning, collapsing -> tree-manipulation. Plotting and mapping annotations onto branches -> tree-visualization. Producing BEAST/MrBayes annotated trees -> bayesian-inference, divergence-dating. Taxon-name sanitization shares the whitespace traps in sequence-io/read-sequences.

## The Single Most Important Modern Insight

A tree file is a lossy serialization of a richer in-memory object. The biologist cares about the topology plus its annotations -- branch supports, posterior probabilities, 95% HPD intervals on node heights, per-branch rates, divergence dates, taxon metadata -- and formats differ enormously in which of these they can hold, while parsers differ in which they actually read back. Three load-bearing facts:

1. **Conversion is a silent data-destroying operation.** Reading a BEAST MCC tree and writing plain Newick produces a topologically identical tree that plots fine, but the HPD intervals, clade posteriors, and per-branch rates are gone and unrecoverable without re-running a multi-day MCMC. The loss is invisible until a reviewer asks where the credible intervals went.
2. **The tool, not the format string, decides whether `[&...]` metadata survives.** In Python the naive default (Bio.Phylo) drops BEAST key-values; in R the naive default (`ape::read.nexus`) drops them; the tools built to preserve them are DendroPy (`extract_comment_metadata=True`) and treeio (`read.beast`). Route annotated trees through those.
3. **In plain Newick a bare number has no fixed meaning.** In `(A,B)95:0.3` the `95` could be a bootstrap, a posterior, an internal clade name, or a second branch length. Only the tool that wrote the file knows; a parser that guesses wrong turns supports into names silently. IQ-TREE overloads the slot further, writing `SH-aLRT/UFBoot` (e.g. `87.5/98`), which a single-value parser truncates or chokes on.

## Tool Taxonomy

| Tool (lang) | BEAST `[&...]` | NHX | phyloXML richness | When |
|-------------|----------------|-----|-------------------|------|
| treeio (R) | YES, structured (`read.beast`/`read.mrbayes`/`read.iqtree`) | YES | via tidytree/ggtree | the default whenever annotations matter; feeds ggtree |
| DendroPy (Py) | YES, structured (`.annotations`) | YES | no | Python work needing metadata, posterior sets, tree distances, conversion with annotations |
| ETE3/ETE4 (Py) | partial (custom features) | YES, native | no | reconciliation, NHX round-trips, programmatic node features |
| Bio.Phylo (Py) | NO key-value parse (opaque `.comment`) | no | YES, richest | general pipelines, conversion among its 5 formats, phyloXML annotation |
| ape (R) | NO (drops `[&...]`) | no | no | fast topology/branch-length analysis; pair with treeio for annotated files |

References: Bio.Phylo Talevich 2012; DendroPy Sukumaran 2010; ETE3 Huerta-Cepas 2016; ape Paradis 2019; treeio Wang 2020.

One-line decision rule: if the file came from BEAST, MrBayes, or RevBayes or carries `[&...]`/`&&NHX` that matters, read it with treeio (R) or DendroPy (Py); otherwise Bio.Phylo (Py) or ape (R) is fine. Never route a BEAST MCC tree through `Bio.Phylo` or `ape::read.nexus` when the HPDs are needed.

## Format Capability (What Each Format Can Hold)

| Capability | Newick | NHX | Nexus (BEAST-annotated) | phyloXML | NeXML |
|------------|--------|-----|-------------------------|----------|-------|
| Topology + branch lengths | yes | yes | yes | yes | yes |
| One support value | ambiguous slot | `B=` tag | comment key | typed `<confidence>` | typed meta |
| Multiple supports per node | no (single slot) | tags | comment keys | yes (n elements) | yes |
| Posterior / HPD intervals | no | no | yes (`_95%_HPD={}`) | via property | via meta |
| Per-branch rates / dates | no | tags | yes | yes | yes |
| Taxonomy (NCBI id/rank) | no | `S=`/`T=` | no | yes, typed | yes |
| Schema-validated | no | no | loose | yes (XSD) | yes (XSD) |

Newick/NHX are compact and grep-able; phyloXML/NeXML are verbose XML (often 5-20x larger) but typed and validatable -- prefer XML for archiving/exchange where machine-checkable semantics matter, Newick/Nexus for pipeline interchange. BEAST/FigTree metadata rides inside a Nexus (or Newick) `[&...]` comment, so "it is a Nexus file" says nothing about whether annotations survive -- only the parser does.

## Reading, Writing, and Converting (Bio.Phylo)

**Goal:** Move trees between formats and inspect them without assuming a single tree or losing annotations.

**Approach:** Use `Phylo.read` for exactly one tree and `Phylo.parse` for many (posterior sets); use `Phylo.convert` only among formats of equal or greater capability; check `.confidence` vs `.name` to confirm support was read into the right slot.

```python
from Bio import Phylo

tree = Phylo.read('tree.nwk', 'newick')          # exactly one tree; raises if 0 or >1
posterior = list(Phylo.parse('run.trees', 'nexus'))   # many trees: posterior/bootstrap set
Phylo.write(tree, 'tree.xml', 'phyloxml')        # phyloXML is Bio.Phylo's richest format

Phylo.convert('tree.nex', 'nexus', 'tree.nwk', 'newick')   # WARNING: Newick cannot hold [&...]; annotations dropped

for clade in tree.get_nonterminals():
    print(clade.confidence, clade.name)          # confirm the support landed in .confidence, not .name
```

Supported format strings: `newick`, `nexus`, `phyloxml`, `nexml`, `cdao`. Colors and branch widths persist only in phyloXML.

## Preserve BEAST/MrBayes Annotations Before Down-Converting

**Goal:** Keep posteriors, HPD intervals, and rates when a downstream tool wants plain Newick.

**Approach:** Read with an annotation-aware parser, extract the numbers into a side table that travels with the analysis, and only then write a stripped topology -- never down-convert first.

```python
import dendropy

tree = dendropy.Tree.get(path='mcc.tree', schema='nexus', extract_comment_metadata=True)
for node in tree:
    if node.annotations.get_value('posterior') is not None:
        post = node.annotations.get_value('posterior')
        hpd = node.annotations.get_value('height_95%_HPD')   # raw BEAST key; treeio typically exposes it as height_0.95_HPD (exact name varies by source program and version -- introspect the columns)
        # persist post/hpd to a side table keyed by the clade before any conversion
tree.write(path='topology.nwk', schema='newick', suppress_annotations=True)   # intentional, after extraction
```

In R the equivalent is treeio `read.beast('mcc.tree')` then `get.data()` / `as_tibble()`, feeding ggtree (tree-visualization); `write.beast()` re-serializes with annotations intact.

## Per-Method Failure Modes

### BEAST/MrBayes MCC to Plain Newick Erases the Credible Intervals
**Trigger:** `Phylo.convert`, `ape::read.nexus` + `write.tree`, or any "just give me the topology" step on an annotated tree.
**Mechanism:** The HPDs, posteriors, and rates live only in the `[&...]` comments, which plain Newick cannot hold and stripping parsers discard.
**Symptom:** The output plots fine but the credible intervals are gone, irrecoverable without re-running the MCMC.
**Fix:** Read with treeio `read.beast` or DendroPy `extract_comment_metadata=True`; extract the numbers to a side table; keep the original `.tree` as the source of truth.

### Support Value Read as a Node Name (or Truncated)
**Trigger:** Parsing a tree whose internal-node slot holds a bootstrap, a posterior, a clade name, or IQ-TREE's `SH-aLRT/UFBoot` dual value.
**Mechanism:** The Newick grammar gives one slot for all of these; the parser must be told which it is, and a single-value reader truncates the `/`-delimited dual support.
**Symptom:** Supports appear as `.name` strings, or only one of two IQ-TREE values survives, or the parse errors on `/`.
**Fix:** Know what wrote the file; in Bio.Phylo inspect `.confidence` vs `.name`; in treeio use `read.iqtree`/`read.raxml`, which split dual support correctly.

### Whitespace, Underscore, or Non-ASCII Taxon Names
**Trigger:** Tip names with spaces, parentheses, commas, or accented characters; reliance on the Newick underscore-space convention.
**Mechanism:** Naive CLI tools split unquoted spaces, and underscore-to-space auto-conversion silently desyncs tip labels from a metadata join key.
**Symptom:** Downstream tools error or a metadata merge matches nothing.
**Fix:** Sanitize to `[A-Za-z0-9_.]`, single-quote when spaces are unavoidable, and round-trip-test the labels against the metadata table before any join.

### Nexus TRANSLATE-Table Desync and Rooted/Unrooted Confusion
**Trigger:** Hand-editing or merging Nexus tree blocks; assuming topology shape implies rootedness.
**Mechanism:** The integer-to-name TRANSLATE map can decouple from the tree and silently relabel tips; Newick does not flag rootedness (a basal trifurcation conventionally signals unrooted, but tools disagree), while Nexus carries an explicit `[&R]`/`[&U]`.
**Symptom:** Tips are mislabeled after a merge, or a rooting-sensitive analysis runs on the wrong assumption without erroring.
**Fix:** Parse with a translate-aware reader (treeio/DendroPy/ape apply it); verify tip-label sets match across merged trees; set rootedness explicitly rather than trusting topology shape.

## Quantitative and Practical Notes

| Item | Guidance | Why |
|------|----------|-----|
| Support-value scales | bootstrap/UFBoot in [0,100], posterior in [0,1], SH-aLRT in [0,100] | a number is meaningless without knowing which test produced it; preserve provenance, not just the value |
| Multi-tree files | use `Phylo.parse` / DendroPy `TreeList` / treeio `read.beast`; check object length | a single-tree reader on a `.trees` posterior returns only the first or errors |
| Round-trip test | read -> write -> read and diff the annotations, not just the topology | topology almost always survives and gives false confidence |
| Posterior set vs MCC | `.trees` is the full posterior; `.tree`/`.mcc` is the single annotated summary | `read.beast` on a full posterior is huge; usually the MCC is wanted |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| HPD bars missing after conversion | converted a BEAST tree to Newick | extract annotations with treeio/DendroPy first |
| `Phylo.read` raises on a `.trees` file | multiple trees in the file | use `Phylo.parse` and iterate |
| Bootstrap values show up as taxon names | node-label slot read as `.name` | set/inspect confidence parsing; use a software-specific reader |
| Metadata join matches nothing | underscore/space relabeling of tips | sanitize and round-trip-test labels before joining |
| Parser errors on `[` | strict parser chokes on FigTree comment | strip comments only after extracting needed metadata |

## References

Cock PJA, Antao T, Chang JT, et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422-1423.
Talevich E, Invergo BM, Cock PJA, Chapman BA. 2012. Bio.Phylo: a unified toolkit for processing, analyzing and visualizing phylogenetic trees in Biopython. *BMC Bioinformatics* 13:209.
Sukumaran J, Holder MT. 2010. DendroPy: a Python library for phylogenetic computing. *Bioinformatics* 26(12):1569-1571.
Huerta-Cepas J, Serra F, Bork P. 2016. ETE 3: reconstruction, analysis, and visualization of phylogenomic data. *Molecular Biology and Evolution* 33(6):1635-1638.
Paradis E, Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35(3):526-528.
Wang L-G, Lam TT-Y, Xu S, et al. 2020. Treeio: an R package for phylogenetic tree input and output with richly annotated and associated data. *Molecular Biology and Evolution* 37(2):599-603.
Maddison DR, Swofford DL, Maddison WP. 1997. NEXUS: an extensible file format for systematic information. *Systematic Biology* 46(4):590-621.
Han MV, Zmasek CM. 2009. phyloXML: XML for evolutionary biology and comparative genomics. *BMC Bioinformatics* 10:356.
Vos RA, Balhoff JP, Caravas JA, et al. 2012. NeXML: rich, extensible, and verifiable representation of comparative data and metadata. *Systematic Biology* 61(4):675-689.

## Related Skills

- tree-manipulation - rooting, pruning, and collapsing where rooted/unrooted and polytomy choices bite
- tree-visualization - ggtree and ETE consume the annotations preserved here
- bayesian-inference - produces the BEAST/MrBayes annotated trees whose metadata must survive
- divergence-dating - produces MCC trees with HPD intervals on node ages
- sequence-io/read-sequences - taxon-name sanitization shares the whitespace and non-ASCII traps
