---
name: bio-phylo-tree-manipulation
description: Modify phylogenetic tree structure using Biopython Bio.Phylo. Use when rooting trees with outgroups, midpoint, or MAD methods, pruning taxa, collapsing clades, ladderizing branches, or extracting subtrees. Includes rooting method decision guidance.
tool_type: python
primary_tool: Bio.Phylo
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Tree Manipulation

**"Root and prune my phylogenetic tree"** → Modify tree topology by rooting with outgroups or midpoint, pruning taxa, collapsing low-support clades, ladderizing branches, or extracting subtrees.
- Python: `tree.root_with_outgroup()`, `tree.prune()`, `tree.ladderize()` (Bio.Phylo)

Modify phylogenetic tree structure: rooting, pruning, ladderizing, and subtree extraction.

## Rooting Method Decision

| Method | Mechanism | When to Use | Pitfalls |
|--------|-----------|-------------|----------|
| Outgroup | External taxon(a) known to be sister to ingroup | **Default method** when a suitable outgroup exists | Too-distant outgroup causes LBA; single outgroup has ~50% misplacement risk |
| Midpoint | Root halfway between most divergent taxa | Quick exploratory; when outgroup is unknown; viral phylogenetics | Assumes clocklike evolution; fails with rate variation |
| MAD | Minimizes ancestor branch deviation across all tip pairs | More robust than midpoint; good complement to outgroup | Still assumes approximately clocklike rates |
| RootDigger | Non-reversible substitution model (CLI tool) | When outgroup is unavailable and clock assumption is untenable | Computationally intensive |

**Best practice:** Use outgroup rooting with multiple sister-group taxa as the primary method. Validate with MAD. If results conflict, investigate rate variation. A single outgroup taxon gives only ~50% accuracy for root placement, so always prefer multiple outgroups that form a monophyletic group.

**Common mistake:** Using a very distant outgroup (e.g., a bacterial sequence to root a eukaryotic tree) introduces a long branch that attracts other long branches, distorting topology near the root.

## Required Import

```python
from Bio import Phylo
from io import StringIO
```

## Rooting Trees

### Root with Outgroup (Preferred)

**Goal:** Root the tree using known sister-group taxa as outgroup.

**Approach:** Verify outgroup taxa are monophyletic, then reroot. Prefer multiple outgroup taxa over a single taxon to improve root placement accuracy.

```python
tree = Phylo.read('tree.nwk', 'newick')

# Root with single taxon (acceptable but less reliable)
tree.root_with_outgroup({'name': 'Outgroup'})

# Root with multiple taxa (preferred -- must be monophyletic)
outgroup = [{'name': 'TaxonA'}, {'name': 'TaxonB'}]
if tree.is_monophyletic(outgroup):
    tree.root_with_outgroup(*outgroup)
else:
    print('Outgroup is not monophyletic -- check taxon selection')
```

### Root at Midpoint

Appropriate when no suitable outgroup exists and evolution is approximately clocklike (e.g., viral phylogenetics). Not reliable when substitution rates vary substantially across lineages.

```python
tree = Phylo.read('tree.nwk', 'newick')
tree.root_at_midpoint()
```

### MAD Rooting (CLI Alternative)

Minimal Ancestor Deviation rooting is more robust to rate variation than midpoint. Install via `pip install mad` or use the standalone tool. Useful as a validation of outgroup rooting results.

### Check Rooting Status

```python
# 2 children at root = rooted; 3+ = unrooted
root = tree.root
print(f'Root has {len(root.clades)} children')
print(f'Is bifurcating: {tree.is_bifurcating()}')
```

## Ladderizing

Sort clades for consistent visual presentation.

```python
tree = Phylo.read('tree.nwk', 'newick')

# Larger clades at bottom
tree.ladderize()

# Larger clades at top
tree.ladderize(reverse=True)

Phylo.write(tree, 'ladderized.nwk', 'newick')
```

## Pruning Trees

### Remove Specific Taxa

```python
tree = Phylo.read('tree.nwk', 'newick')

# Find and remove a taxon
target = tree.find_any(name='TaxonToRemove')
if target:
    tree.prune(target)

# Remove multiple taxa
for name in ['TaxonA', 'TaxonB', 'TaxonC']:
    target = tree.find_any(name=name)
    if target:
        tree.prune(target)
```

### Keep Only Specified Taxa

```python
tree = Phylo.read('tree.nwk', 'newick')
keep_taxa = {'Human', 'Chimp', 'Gorilla'}

terminals = tree.get_terminals()
for term in terminals:
    if term.name not in keep_taxa:
        tree.prune(term)
```

## Collapsing Clades

Collapse branches below a support threshold into polytomies. The threshold depends on the support measure used: >= 70 for standard bootstrap, >= 95 for UFBoot2 (these are not equivalent scales).

```python
tree = Phylo.read('tree.nwk', 'newick')

# Collapse single clade
target = tree.find_any(name='SomeInternalNode')
if target:
    tree.collapse(target)

# Collapse by branch length
tree.collapse_all(lambda c: c.branch_length and c.branch_length < 0.01)

# Collapse poorly-supported nodes
# Use 70 for standard bootstrap; use 95 for UFBoot2 values
tree.collapse_all(lambda c: c.confidence is not None and c.confidence < 70)
```

## Extracting Subtrees

### Get Clade as Subtree

```python
tree = Phylo.read('tree.nwk', 'newick')

# Find common ancestor of taxa
clade = tree.common_ancestor({'name': 'Human'}, {'name': 'Chimp'})

# The clade itself can be treated as a subtree
Phylo.draw_ascii(clade)

# Get all terminals in this clade
subtree_taxa = [t.name for t in clade.get_terminals()]
print(f'Subtree contains: {subtree_taxa}')
```

### Extract Subtree by Common Ancestor

```python
tree = Phylo.read('tree.nwk', 'newick')

# Find MRCA (Most Recent Common Ancestor)
taxa = [{'name': 'Human'}, {'name': 'Chimp'}, {'name': 'Gorilla'}]
mrca = tree.common_ancestor(*taxa)
print(f'MRCA branch length: {mrca.branch_length}')
```

## Tree Traversal

```python
tree = Phylo.read('tree.nwk', 'newick')

# Iterate all clades (preorder by default)
for clade in tree.find_clades():
    print(clade.name, clade.branch_length)

# Level-order traversal (breadth-first)
for clade in tree.find_clades(order='level'):
    print(clade.name)

# Postorder traversal
for clade in tree.find_clades(order='postorder'):
    print(clade.name)

# Only terminal nodes
for term in tree.get_terminals():
    print(term.name)

# Only internal nodes
for internal in tree.get_nonterminals():
    print(internal)
```

## Finding Clades

```python
tree = Phylo.read('tree.nwk', 'newick')

# Find by name
clade = tree.find_any(name='Human')

# Find all matching criteria
matches = tree.find_clades(branch_length=lambda x: x and x > 0.5)
for m in matches:
    print(f'{m.name}: {m.branch_length}')

# Find by terminal status
terminals = list(tree.find_clades(terminal=True))
internals = list(tree.find_clades(terminal=False))
```

## Getting Path Between Nodes

```python
tree = Phylo.read('tree.nwk', 'newick')

# Path from root to a node
target = tree.find_any(name='Human')
path = tree.get_path(target)
print(f'Path from root to Human: {len(path)} nodes')
for clade in path:
    print(f'  {clade.name}: {clade.branch_length}')

# Trace path between any two nodes
human = tree.find_any(name='Human')
mouse = tree.find_any(name='Mouse')
trace = tree.trace(human, mouse)
print(f'Path Human to Mouse: {len(trace)} nodes')
```

## Checking Tree Properties

```python
tree = Phylo.read('tree.nwk', 'newick')

# Check if monophyletic
taxa = [tree.find_any(name='Human'), tree.find_any(name='Chimp')]
taxa = [t for t in taxa if t is not None]
print(f'Is monophyletic: {tree.is_monophyletic(taxa)}')

# Check if bifurcating
print(f'Is bifurcating: {tree.is_bifurcating()}')

# Check if preterminal (parent of only terminals)
for clade in tree.get_nonterminals():
    print(f'{clade}: is_preterminal={clade.is_preterminal()}')
```

## Modifying Branch Lengths

```python
tree = Phylo.read('tree.nwk', 'newick')

# Set missing branch lengths
for clade in tree.find_clades():
    if clade.branch_length is None:
        clade.branch_length = 0.0

# Scale all branch lengths
scale_factor = 100  # Convert to percent divergence
for clade in tree.find_clades():
    if clade.branch_length:
        clade.branch_length *= scale_factor

# Remove branch lengths (convert to cladogram)
for clade in tree.find_clades():
    clade.branch_length = None
```

## Renaming Taxa

```python
tree = Phylo.read('tree.nwk', 'newick')

# Rename individual taxon
target = tree.find_any(name='OldName')
if target:
    target.name = 'NewName'

# Batch rename from mapping
name_map = {'Hsap': 'Human', 'Ptro': 'Chimp', 'Mmus': 'Mouse'}
for term in tree.get_terminals():
    if term.name in name_map:
        term.name = name_map[term.name]

Phylo.write(tree, 'renamed.nwk', 'newick')
```

## Counting Nodes

```python
tree = Phylo.read('tree.nwk', 'newick')

n_terminals = len(tree.get_terminals())
n_internals = len(tree.get_nonterminals())
n_total = tree.count_terminals() + len(tree.get_nonterminals())

print(f'Terminals: {n_terminals}')
print(f'Internal nodes: {n_internals}')
print(f'Total nodes: {n_total}')
```

## Tree Depths

```python
tree = Phylo.read('tree.nwk', 'newick')

# Get depths from root
depths = tree.depths()
for clade, depth in depths.items():
    if clade.is_terminal():
        print(f'{clade.name}: depth={depth:.3f}')

# Get maximum depth (tree height)
max_depth = max(depths.values())
print(f'Tree height: {max_depth:.3f}')
```

## Splitting Clades

```python
tree = Phylo.read('tree.nwk', 'newick')

# Split a terminal into multiple children
target = tree.find_any(name='TaxonA')
if target and target.is_terminal():
    target.split(n=2, branch_length=0.05)  # Creates 2 children

# Split with specific branch lengths
target.split(branch_length=[0.1, 0.2, 0.3])  # Creates 3 children
```

## Generating Random Trees

```python
from Bio.Phylo.BaseTree import Tree

# Generate random bifurcating tree
taxa = ['Human', 'Chimp', 'Gorilla', 'Mouse', 'Rat']
random_tree = Tree.randomized(taxa)
Phylo.draw_ascii(random_tree)

# With branch lengths
random_tree = Tree.randomized(taxa, branch_length=1.0)
```

## Quick Reference: Tree Methods

| Method | Description |
|--------|-------------|
| `root_with_outgroup()` | Reroot using outgroup |
| `root_at_midpoint()` | Reroot at midpoint |
| `ladderize()` | Sort branches by size |
| `prune()` | Remove a clade |
| `collapse()` | Collapse a clade into polytomy |
| `collapse_all()` | Collapse all matching clades |
| `split()` | Split clade into children |
| `trace()` | Get path between two clades |
| `Tree.randomized()` | Generate random tree |
| `common_ancestor()` | Find MRCA of taxa |
| `find_any()` | Find first matching clade |
| `find_clades()` | Find all matching clades |
| `get_path()` | Get path from root to clade |
| `depths()` | Get depth of all clades |
| `is_monophyletic()` | Check if taxa form clade |
| `is_bifurcating()` | Check if tree is binary |

## Related Skills

- tree-io - Read and write tree files
- tree-visualization - Draw modified trees
- distance-calculations - Build trees from alignments
- modern-tree-inference - ML tree inference (rooting before/after inference)
- species-trees - Coalescent methods and concordance factors
