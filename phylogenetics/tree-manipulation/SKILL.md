---
name: bio-phylo-tree-manipulation
description: Edit phylogenetic tree structure with Biopython Bio.Phylo, and treat rooting as a separate statistical inference rather than a display choice. Covers why most inference returns an unrooted tree so placing the root creates every ancestor/descendant and basal claim; why a distant or lonely outgroup misroots inside the ingroup via long-branch attraction; the outgroup/midpoint/MAD/MinVar/non-reversible-likelihood rooting tradeoffs; why pruning must suppress degree-2 nodes and sum their branch lengths or all patristic distances silently corrupt; and why collapsing by support makes SOFT (uncertainty) polytomies, not HARD (radiation) ones. Use when rooting, re-rooting, pruning or subsetting taxa, extracting a clade or induced subtree, collapsing low-support branches, resolving polytomies, or ladderizing. Routes clock-based rooting to divergence-dating, inference to modern-tree-inference, and reading/plotting to tree-io and tree-visualization.
tool_type: mixed
primary_tool: Bio.Phylo
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+. Alternatives: ape 5.8+ / phangorn / phytools (R), DendroPy 5+ and ete3 (Python), Newick Utilities 1.6+, MAD and RootDigger as standalone CLIs.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(module.function)` to check signatures
- R: `packageVersion('ape')` then `?drop.tip` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

ape has NO midpoint function of its own (use phangorn::midpoint or phytools::midpoint.root); ete3 prune drops branch lengths unless preserve_branch_length=True; ape di2multi filters by branch LENGTH (tol), not support.

# Tree Manipulation -- Rooting Is an Inference, Not a Cosmetic Operation

**"Root and prune my tree"** -> Edit an existing tree object, treating the root placement as a separate statistical inference and every structural edit as potentially distance-corrupting.
- Python: `tree.root_with_outgroup(...)`, `tree.root_at_midpoint()`, `tree.prune(...)`, `tree.collapse_all(...)` (Bio.Phylo)
- R: `root()`, `phangorn::midpoint()`, `drop.tip()`, `di2multi()` (ape/phangorn); CLI: `nw_reroot`, `nw_prune` (Newick Utilities); outgroup-free rooters: `mad`, `FastRoot.py`, `rootdigger`

Scope: editing an existing tree -- rooting, pruning/subsetting, extracting clades or induced subtrees, collapsing branches into polytomies, resolving polytomies, ladderizing. Reading/writing/converting tree files -> tree-io. Plotting and mapping annotations -> tree-visualization. Inferring the tree in the first place -> modern-tree-inference. Clock-based or Bayesian rooting co-estimated with dates -> divergence-dating.

## The Single Most Important Modern Insight

Rooting is a separate statistical inference layered on top of the topology, and it is the highest-error decision in the entire tree. Nearly every inference method (ML under a time-reversible model, neighbor-joining, most Bayesian runs) returns an UNROOTED tree: a reversible model is blind to the direction of time, so the likelihood is identical wherever the root sits. The unrooted tree states which taxa are neighbors; it is silent about who is ancestor and who is descendant. Placing the root is what converts a neighbor-graph into an evolutionary narrative, and it CREATES every "X is basal", every "the common ancestor had trait T", every character-polarity, every divergence ordering. Three load-bearing facts:

1. **Rooting needs its own justification and its own uncertainty statement, distinct from branch support.** Bootstrap 95 on an ingroup clade says nothing about whether the root is in the right place. The relationships can be robust while the root -- the thing every downstream story hangs on -- is the least reliable part of the figure.
2. **The deep nodes near an outgroup-defined root are the LEAST trustworthy.** A distant outgroup sits on a long branch, and long branches attract each other (long-branch attraction): the outgroup pulls the fastest ingroup taxon toward the root and makes it look artifactually "basal". Accuracy degrades monotonically as outgroup distance grows.
3. **The choice of rooting METHOD is a modeling choice with assumptions** (clock? outgroup monophyly? non-reversible model?), and the wrong method silently fabricates basal lineages. Pruning, collapsing, and ladderizing are comparatively mechanical -- but each has its own silent traps (distance corruption, soft/hard polytomy conflation, perceived trends).

## Rooting Method Decision

| Method | Core assumption | Needs | Robust to rate variation | Gives root uncertainty | Best when | Avoid when |
|--------|-----------------|-------|--------------------------|------------------------|-----------|------------|
| Outgroup (multiple, close) | outgroup truly outside ingroup, not LBA-misplaced | a priori outgroup taxa + sequences | moderate (depends on outgroup branch) | weakly (via ingroup-monophyly check) | closely related, balanced, monophyletic outgroups exist | only a distant/lonely outgroup available |
| Outgroup (single) | same, but one long branch | one outgroup taxon | poor (max LBA exposure) | no | a close single sister exists, nothing better | only a distant/lonely outgroup (one long branch, max LBA exposure) |
| Midpoint | strict molecular clock | tree with branch lengths | no | no | shallow, clock-like, intraspecific/viral data | deep trees, heterotachy, any long branch |
| MinVar | clock deviations are random/unbiased | tree with branch lengths | better than midpoint | no | outgroup-free, noise looks like unbiased clock scatter | strong lineage-specific rate shifts |
| MAD | minimal relative ancestor deviation | tree with branch lengths | yes (tolerates heterotachy) | yes (root ambiguity index) | outgroup-free deep/prokaryotic trees; default no-outgroup choice | strongly structured rate variation |
| Non-reversible likelihood (RootDigger / IQ-TREE) | non-reversible model carries directional signal | tree + ALIGNMENT; enough data | yes (model-based) | yes (per-branch likelihood confidence) | the alignment is available and a root confidence is needed | little data; weak signal; compute-limited |
| Relaxed clock (BEAST) -> divergence-dating | explicit clock + tree prior | alignment + calibrations/tip dates; MCMC | yes | yes (posterior over roots) | time-scaled / dated / phylodynamic analyses | a quick structural edit, no dating intended |

Rule: never report a root from a single method without a sanity check. If a close, monophyletic outgroup is available, use it. If not, run MAD and MinVar and prefer agreement; disagreement means the root is poorly determined and all "basal" claims must be hedged. With an alignment and a need for a confidence value, use RootDigger (Bettisworth and Stamatakis 2021) or non-reversible IQ-TREE. MAD = Tria et al. 2017; MinVar = Mai et al. 2017; outgroup-free CLIs and Newick Utilities = Junier and Zdobnov 2010.

## Tool Taxonomy

| Tool (lang) | Rooting | Pruning / clade | Collapse / resolve | When |
|-------------|---------|-----------------|--------------------|------|
| Bio.Phylo (Py) | `root_with_outgroup`, `root_at_midpoint` | `prune`, `common_ancestor` | `collapse_all` | general Python pipelines; the default here |
| ape / phangorn / phytools (R) | `root`, `phangorn::midpoint`, `phytools::midpoint.root` | `drop.tip`, `extract.clade` | `di2multi`, `multi2di` | R workflows; `drop.tip` sums suppressed branch lengths correctly |
| DendroPy (Py) | `reroot_at_edge`, `reroot_at_midpoint` | `retain_taxa_with_labels(suppress_unifurcations=True)` | edge collapse, `resolve_polytomies` | metadata-aware edits, posterior tree sets |
| ete3 (Py) | `set_outgroup`, `set_outgroup(get_midpoint_outgroup())` | `prune([...], preserve_branch_length=True)`, `detach` | `delete`, `resolve_polytomy` | NHX features; MUST pass preserve_branch_length |
| Newick Utilities (CLI) | `nw_reroot` | `nw_prune`, `nw_clade` | `nw_ed`, `nw_condense` | streaming/Unix-pipeline edits on many trees |
| MAD / MinVar / RootDigger (CLI) | `mad`, `FastRoot.py`, `rootdigger` | -- | -- | outgroup-free or likelihood rooting when no trustworthy outgroup exists |

## Root with an Outgroup (Multiple, Monophyletic)

**Goal:** Root the tree using known sister-group taxa, preferring multiple close outgroups and verifying ingroup monophyly first.

**Approach:** Confirm the outgroup taxa form a monophyletic group, root on the branch separating them from the ingroup, then check the ingroup is recovered as monophyletic -- if not, the rooting is suspect.

```python
from Bio import Phylo

tree = Phylo.read('tree.nwk', 'newick')
outgroup = [{'name': 'OutA'}, {'name': 'OutB'}]      # multiple close outgroups beat a single long branch

if tree.is_monophyletic([tree.find_any(name='OutA'), tree.find_any(name='OutB')]):
    tree.root_with_outgroup(*outgroup)               # root_with_outgroup, NOT root_with_midpoint
else:
    print('outgroup not monophyletic: root placement is unreliable, re-check taxon choice')
```

## Root at the Midpoint (Clock-Limited Fallback)

**Goal:** Root an outgroup-free tree where evolution is approximately clock-like (shallow/viral data).

**Approach:** Place the root at the midpoint of the longest tip-to-tip path; trust it only when no single long branch can hijack that path. For deep trees with rate variation prefer MAD/MinVar.

```python
tree = Phylo.read('tree.nwk', 'newick')
tree.root_at_midpoint()                              # assumes a clock; a long branch slides the root onto the fast lineage
# Outgroup-free and clock-relaxed (deep trees): standalone CLIs run on the Newick file
# MAD:    mad tree.nwk           -> tree.nwk.rooted   (per-branch root ambiguity index)
# MinVar: FastRoot.py -i tree.nwk -m MV -o rooted.nwk
```

## Prune Taxa With Branch-Length Preservation

**Goal:** Remove tips and keep every surviving patristic distance unchanged.

**Approach:** Drop the tip and SUPPRESS the resulting degree-2 ("knee") node, ADDING its branch length to the child so path lengths are conserved. Bio.Phylo `prune` and ape `drop.tip` do this by default; ete3 `prune` needs `preserve_branch_length=True`.

```python
tree = Phylo.read('tree.nwk', 'newick')
keep = {'Human', 'Chimp', 'Mouse'}

for term in list(tree.get_terminals()):
    if term.name not in keep:
        tree.prune(term)                             # collapses the degree-2 parent and sums branch lengths
# ape (R):   drop.tip(phy, c('X','Y'))               # sums suppressed branch lengths by default
# ete3 (Py): tree.prune(list(keep), preserve_branch_length=True)   # the flag is mandatory, else distances shrink
```

Non-monophyletic targets cannot be "extracted as a clade" -- there is no node whose descendants are exactly those taxa. Prune to the taxon set to get the induced subtree instead; `common_ancestor` of non-monophyletic taxa returns an MRCA whose clade contains EXTRA taxa.

## Collapse Low-Support Branches Into SOFT Polytomies

**Goal:** Replace poorly-supported resolved nodes with multifurcations that honestly say "unresolved".

**Approach:** Collapse any internal branch whose support is below a stated cutoff. The result is a SOFT (uncertainty) polytomy, never a HARD (simultaneous-radiation) one -- label it as such.

```python
tree = Phylo.read('tree.nwk', 'newick')              # support parsed into clade.confidence

tree.collapse_all(lambda c: c.confidence is not None and c.confidence < 70)   # 70 for std bootstrap; use 95 for UFBoot2
tree.collapse_all(lambda c: c.branch_length is not None and c.branch_length < 1e-8)   # collapse genuinely-zero branches
# ape (R): di2multi filters by LENGTH (tol), not support -> zero out low-support branch lengths first, THEN di2multi(phy)
```

Resolving the inverse direction (`multi2di` / `resolve_polytomy`) invents an arbitrary order with zero-length branches; analyzing one random resolution treats an arbitrary choice as fact. If a binary tree is required, integrate over many random resolutions and treat the inserted zero-length branches as "no information", not instantaneous divergence.

## Per-Method Failure Modes

### Distant Outgroup Roots Inside the Ingroup
**Trigger:** A single distant outgroup, or many outgroups all far from the ingroup, used to root.
**Mechanism:** The long outgroup branch attracts the fastest ingroup taxon (LBA), pulling the root into the ingroup; far outgroups attach at essentially random positions (DeSalle et al. 2023).
**Symptom:** Ingroup not recovered as monophyletic; a fast taxon appears "basal"; the root jumps as the outgroup set changes.
**Fix:** Use multiple CLOSELY-related, monophyletic, roughly-equidistant outgroups; verify ingroup monophyly; treat a lonely distant outgroup as a red flag.

### Single Distant Outgroup -- Long-Branch Misrooting
**Trigger:** Rooting on one distant outgroup taxon.
**Mechanism:** One unbroken long branch maximizes LBA exposure, with no way to subdivide it or check ingroup monophyly, so the root is drawn toward other long branches.
**Symptom:** The root lands inside the ingroup or on a spurious deep branch; deep nodes near the root are unstable across analyses.
**Fix:** Add multiple closer, monophyletic outgroups to subdivide the long branch and enable a monophyly check, or use an outgroup-free method (MAD/MinVar).

### Midpoint Misroots Under Rate Variation
**Trigger:** Midpoint rooting a deep tree with heterotachy or any long branch.
**Mechanism:** Midpoint assumes a clock; a long branch hijacks the longest path and slides the root onto the fast lineage.
**Symptom:** A rate-elevated taxon appears earliest-diverging; the root sits on a suspiciously long branch.
**Fix:** Use MAD or MinVar (clock-relaxed) or a close outgroup; report root uncertainty and hedge "basal".

### Soft vs Hard Polytomy Conflation
**Trigger:** Collapsing branches below a support threshold, then describing the multifurcation as a radiation.
**Mechanism:** Threshold-collapse encodes UNCERTAINTY (soft polytomy = "we cannot resolve the order"); a hard polytomy is a biological claim of simultaneous divergence -- different meanings.
**Symptom:** A paper claims simultaneous divergence from what is just unresolved data.
**Fix:** Label threshold-collapsed nodes as unresolved/soft; never read them as a biological radiation.

### Pruning Leaves a Spurious Node or Drops Distances
**Trigger:** A tool that leaves the degree-2 node in place, or removes it without summing branch lengths (notably ete3 `prune` without `preserve_branch_length=True`).
**Mechanism:** The suppressed knee's branch length is not added to its child, so every path through that lineage shortens.
**Symptom:** Patristic distances shrink; subsequent midpoint/MAD/MinVar rooting on the pruned tree is now wrong.
**Fix:** Use Bio.Phylo `prune` / ape `drop.tip` (correct by default) or pass ete3 `preserve_branch_length=True`; verify a known pairwise distance is unchanged.

### Randomly Resolving a Polytomy Biases Downstream
**Trigger:** `multi2di` / `resolve_polytomy` to satisfy a binary-tree requirement, then analyzing the single tree.
**Mechanism:** An arbitrary order with zero-length branches is invented; the downstream result depends on a topology the data never supported.
**Symptom:** Results that change under a different random resolution.
**Fix:** Integrate over many resolutions and summarize; treat zero-length branches as "no info", not instantaneous divergence.

### Ladderizing Manufactures an Apparent Trend
**Trigger:** Ladderizing a figure so one lineage sits visually at the top/bottom.
**Mechanism:** Rotation about internal nodes changes no topology, branch lengths, or bipartitions, but readers unconsciously read the staircase as an early-to-late sequence.
**Symptom:** Reviewers infer a basal-to-derived trend that does not exist.
**Fix:** Ladderize only for legibility; state that rotation changes no biology. -> tree-visualization.

## Quantitative Thresholds

| Quantity | Value | Source / rationale |
|----------|-------|--------------------|
| Midpoint root recovery, single-outgroup source data | ~54% (barely a coin flip) | Hess and Russo 2007 |
| Midpoint root recovery, multiple-outgroup source data | ~82% (inconsistent) to ~94% (consistent) | Hess and Russo 2007 |
| MAD accuracy on benchmarks | >~70%, beating midpoint | Tria et al. 2017 |
| MinVar vs midpoint | matched or beat midpoint in all simulated conditions | Mai et al. 2017 |
| Bootstrap collapse cutoff | <50% (near-uninformative floor) or <70% (reliability boundary); state which | common practice |
| UFBoot2 collapse cutoff | <95% (a different scale from bootstrap -- not 70) | Hoang 2018 |
| Bayesian posterior collapse cutoff | <0.95 | common practice |
| Zero-length collapse tolerance | ~1e-8 (or machine epsilon) | removes genuinely-zero branches without deleting real short ones |
| Outgroup distance | accuracy degrades monotonically with distance | prefer the closest credible outgroup |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `ape::midpoint` function-not-found | midpoint is not in base ape | use `phangorn::midpoint` or `phytools::midpoint.root` |
| Pruned tree has wrong distances | degree-2 node not suppressed / not summed | Bio.Phylo `prune` / ape `drop.tip` default, or ete3 `preserve_branch_length=True` |
| "Extract the clade of X,Y,Z" returns extra taxa or None | X,Y,Z are not monophyletic | prune to the taxon set (induced subtree) instead of extracting a clade |
| `di2multi(tol)` did not drop low-support nodes | `tol` filters branch LENGTH, not support | zero out low-support branch lengths first, then `di2multi` |
| Bootstrap re-read as support for a pruned subtree | support was computed on the full taxon set | re-run inference on the subset for valid support |
| MAD and MinVar disagree on the root | the root is genuinely poorly determined | treat as uncertain; seek a close outgroup or RootDigger confidence; hedge all ancestral claims |
| `root_with_midpoint` AttributeError | no such method | the methods are `root_at_midpoint()` and `root_with_outgroup()` |

## References

Tria FDK, Landan G, Dagan T. 2017. Phylogenetic rooting using minimal ancestor deviation. *Nature Ecology & Evolution* 1:0193.
Mai U, Sayyari E, Mirarab S. 2017. Minimum variance rooting of phylogenetic trees and implications for species tree reconstruction. *PLOS ONE* 12(8):e0182238.
Bettisworth B, Stamatakis A. 2021. Root Digger: a root placement program for phylogenetic trees. *BMC Bioinformatics* 22:225.
Hess PN, De Moraes Russo CA. 2007. An empirical test of the midpoint rooting method. *Biological Journal of the Linnean Society* 92(4):669-674.
DeSalle R, Narechania A, Tessler M. 2023. Multiple outgroups can cause random rooting in phylogenomics. *Molecular Phylogenetics and Evolution* 184:107806.
Junier T, Zdobnov EM. 2010. The Newick utilities: high-throughput phylogenetic tree processing in the Unix shell. *Bioinformatics* 26(13):1669-1670.
Huerta-Cepas J, Serra F, Bork P. 2016. ETE 3: reconstruction, analysis, and visualization of phylogenomic data. *Molecular Biology and Evolution* 33(6):1635-1638.
Paradis E, Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35(3):526-528.

## Related Skills

- tree-io - reading and writing the trees these edits consume and produce, without dropping annotations
- tree-visualization - ladderize is a perception choice; mapping support onto branches
- divergence-dating - clock-based and Bayesian rooting co-estimated with node ages, not bolted on by midpoint
- modern-tree-inference - produces the unrooted ML tree and the support values these edits act on
