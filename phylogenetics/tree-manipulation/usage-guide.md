# Tree Manipulation - Usage Guide

## Overview

This skill edits the structure of an existing phylogenetic tree: rooting, pruning/subsetting, extracting clades or induced subtrees, collapsing branches into polytomies, resolving polytomies, and ladderizing. The idea that matters most is that rooting is not editing; it is a separate statistical inference and the highest-error decision in the whole tree. Most inference returns an unrooted tree, so placing the root is what creates every ancestor/descendant claim, every "basal" lineage, and every character polarity, and the deep nodes near an outgroup-defined root are the least trustworthy part of the figure. The other operations are comparatively mechanical but each carries a silent trap: pruning can corrupt every patristic distance, collapsing by support produces uncertainty (soft) polytomies that are not radiations, and ladderizing can manufacture an apparent evolutionary trend that the data do not contain.

## Prerequisites

```bash
pip install biopython          # Bio.Phylo (primary)
pip install dendropy ete3       # annotation-aware / NHX alternatives
# R: install.packages(c('ape','phangorn','phytools'))
# CLI (optional): Newick Utilities (nw_reroot, nw_prune); MAD, MinVar-Rooting (FastRoot.py), RootDigger for outgroup-free / likelihood rooting
```

Conceptual prerequisites: a tree object (from tree-io), an understanding that branch support is not root confidence, and a priori knowledge of which taxa are outgroups when outgroup rooting is intended.

## Quick Start

Tell your AI agent what you want to do:
- "Root this tree using OutA and OutB as the outgroup and check ingroup monophyly"
- "Root at the midpoint, but warn me if this tree is not clock-like"
- "Prune everything except the primates and keep the branch lengths correct"
- "Collapse branches with bootstrap below 70 into polytomies"

## Example Prompts

### Rooting
> "Root this tree with OutA and OutB as outgroups, but tell me if they are not monophyletic"

> "I have no outgroup and the tree is deep with rate variation: root it with MAD and MinVar and tell me if they agree"

> "Is this tree rooted, and how reliable is the root placement?"

### Pruning and subsetting
> "Remove all the bacterial taxa and confirm the patristic distances are unchanged"

> "Give me the tree restricted to just these eight taxa, even though they are not a clade"

> "Extract the mammal clade as a self-contained subtree"

### Collapsing and resolving
> "Collapse all branches with UFBoot2 below 95 into polytomies"

> "Collapse the zero-length branches but leave the real short branches alone"

> "I need a fully binary tree for a downstream method: resolve the polytomies and warn me about the bias"

### Ladderizing
> "Ladderize the tree for the figure, and confirm this changes no biology"

## What the Agent Will Do

1. Read the tree from a file or string (routing annotation-critical reads to tree-io).
2. For rooting: pick a method from the data (close monophyletic outgroup preferred; MAD/MinVar when outgroup-free; non-reversible likelihood when an alignment and a confidence value are needed), verify ingroup monophyly, and flag a distant/lonely outgroup as an LBA risk.
3. For pruning: drop tips while suppressing degree-2 nodes and summing their branch lengths, and verify a known pairwise distance is unchanged.
4. For collapsing: collapse below a stated support cutoff, label the result a soft (uncertainty) polytomy, and refuse to describe it as a radiation.
5. Verify tree integrity (rootedness, bifurcation, taxon set) after each edit and save or return the result.

## Tips

- Justify the root separately from branch support: bootstrap 95 on a clade says nothing about whether the root is correct.
- Prefer multiple, closely related, monophyletic outgroups over a single distant one, which is a long branch that maximizes long-branch-attraction exposure and can root inside the ingroup.
- Avoid distant outgroups: they sit on long branches that attract ingroup long branches and root the tree inside the ingroup.
- Use midpoint rooting only for shallow, clock-like data (e.g. a single virus over months); for deep trees with rate variation use MAD or MinVar and prefer their agreement.
- When MAD and MinVar disagree, treat the root as poorly determined and hedge every "basal" claim.
- When pruning, confirm distances are preserved (ape `drop.tip` and Bio.Phylo `prune` are correct by default; ete3 `prune` needs `preserve_branch_length=True`).
- Do not re-interpret retained bootstrap values as support for a pruned subset; re-run inference on the subset if subset support is needed.
- Collapsing by support makes soft polytomies (uncertainty), not hard ones (simultaneous radiation); never conflate them.
- Ladderize only for legibility; rotation changes no topology, branch lengths, or biology.

## Related Skills

- tree-io - reading and writing the trees these edits consume and produce
- tree-visualization - ladderize as a perception choice and mapping support onto branches
- divergence-dating - clock-based and Bayesian rooting co-estimated with node ages
- modern-tree-inference - produces the unrooted ML tree and support values these edits act on
