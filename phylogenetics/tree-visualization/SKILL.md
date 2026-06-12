---
name: bio-phylo-tree-visualization
description: Draw and export phylogenetic trees with Bio.Phylo plus matplotlib, and route rich figures to ggtree, ETE4, or iTOL. Covers why a tree figure is an argument not a neutral picture, the cladogram-vs-phylogram-vs-chronogram choice that hides or reveals rate and time, how ladderization manufactures a false arrow of progress, why an unlabeled support number always flatters the result (bootstrap vs posterior vs SH-aLRT vs UFBoot are different scales), why Bio.Phylo silently drops BEAST HPD bars so annotated Bayesian trees must go through treeio plus ggtree, and the tip-count and raster-vs-vector thresholds for legible publication figures. Use when drawing a tree, choosing a layout, coloring branches, showing or labeling support, exporting vector figures, or deciding a drawing tool. Routes annotation-preserving reads to tree-io, and rooting and ladderizing to tree-manipulation.
tool_type: python
primary_tool: Bio.Phylo
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, matplotlib 3.8+. Rich-figure alternatives: ggtree 3.12+ / treeio 1.28+ / ggtreeExtra 1.14+ (Bioconductor), ete4 4.x (Python), iTOL v6 (web), FigTree (desktop).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Phylo.draw)` to check signatures
- R: `packageVersion('ggtree')` then `?ggtree` to verify layout and geom arguments

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Bio.Phylo `draw()` renders rectangular layouts only and parses no BEAST/MrBayes annotation block; circular/fan/unrooted layouts and HPD bars require ggtree, ETE4, or iTOL. The ETE3->ETE4 migration changed APIs; introspect rather than assume ETE3 signatures.

# Tree Visualization -- A Figure Is an Argument

**"Draw and export my tree figure"** -> Render a topology under a chosen layout and branch-length encoding, decide which support to show and how to label it, and export in a format that survives print.
- Python quick look: `Phylo.draw(tree, axes=ax)` (Bio.Phylo + matplotlib)
- Publication with metadata, rich support, or BEAST HPD bars: ggtree + treeio (R)

Scope: drawing, styling, annotating, and exporting tree figures, and choosing the drawing tool. Reading/converting files and preserving `[&...]` annotations before plotting -> tree-io. Rooting, ladderizing, pruning, collapsing low-support nodes -> tree-manipulation. ggplot2 grammar, themes, and `ggsave` underlying ggtree -> data-visualization/ggplot2-fundamentals. Composing a tree plus aligned panels into one figure -> data-visualization/multipanel-figures.

## The Single Most Important Modern Insight

A phylogenetic tree figure is a rhetorical and modeling artifact, not a neutral picture of the data. A Newick or Nexus file encodes a topology, optionally branch lengths, optionally one node label; it does NOT encode which child is drawn on top, whether the horizontal axis means evolutionary distance or just nesting depth, where the root sits, the aspect ratio, or which of several support measures a bare "95" refers to. Every one of those is supplied by the person drawing the tree, so the figure is the data filtered through a dozen interpretive choices the reader cannot see. Four argument-level decisions dominate, and the same file can be drawn to tell mutually contradictory stories that are all technically correct renderings of one topology:

1. **The geometry is a claim about what the lengths mean.** A cladogram asserts only nesting; a phylogram asserts nesting plus amount of change (substitutions per site); a chronogram asserts timing and demands a clock model plus calibrations. Drawing a phylogram as a cladogram erases rate variation, the very signal that flags long-branch attraction; drawing untrusted branch lengths as a phylogram fabricates a quantitative claim.
2. **Ladderization implies a directionality that is not there.** Sibling subclades are unordered, so any node rotation is the same tree, but the eye reads top-to-bottom ordering as a march of progress toward the bottom-most tip. Rotation can also bury non-monophyly: a paraphyletic group can be rotated to look contiguous, or a clean clade made to look scattered.
3. **Which support is shown, and whether it is labeled, changes the conclusion.** Support is not accuracy, and the measures are not interchangeable: a posterior of 0.95 is generally weaker evidence than a bootstrap of 95 for the same data, and UFBoot 95, SH-aLRT 80, and TBE live on different scales. A bare integer with no legend is the quietest lie because the reader assumes bootstrap and the ambiguity always flatters the result.
4. **The drawing tool is decided by the I/O layer, not the surrounding language.** If the tree carries BEAST/MrBayes HPD intervals and posteriors, only a tool whose data model can carry those annotations (treeio + ggtree) can plot them; Bio.Phylo flattens the tree to topology + length + one label and silently drops the uncertainty that was the result.

The skill's job is to make these choices conscious and DECLARE them in the caption, not to default into them. Default ladderization plus cladogram-ish rendering plus unlabeled nodes is the path of least resistance, and it is an argument the user did not intend to make.

## Tool Taxonomy

The central decision is which drawing tool to reach for, and it is coupled to the I/O decision: a tool that flattens the tree to topology plus one label cannot draw BEAST HPD bars or posteriors it cannot represent. The plotting path through treeio + ggtree is the same I/O choice made in tree-io.

| Tool (lang) | Citation | Layouts and annotation power | When to choose |
|-------------|----------|------------------------------|----------------|
| Bio.Phylo + matplotlib (Py) | Talevich 2012 | rectangular phylogram/cladogram only; one label per node via `label_func`/`branch_labels`; no circular/unrooted, no BEAST-block parse | a fast scripted look inside a Python pipeline; headless/CI rendering of simple trees; NOT publication figures with metadata or HPD bars |
| ggtree + treeio + ggtreeExtra (R) | Yu 2017; Wang 2020; Xu 2021 | rectangular, slanted, circular, fan, radial/unrooted; `%<+%` attaches a data.frame, `geom_tiplab`/`geom_nodelab`/`geom_cladelab`, `gheatmap`, `geom_facet` (Cartesian) or `geom_fruit` rings (circular); treeio reads BEAST/MrBayes/IQ-TREE so HPD bars and posterior plot directly | the publication standard whenever a figure integrates metadata, rich/dual support, HPD bars, or aligned panels; the only clean path for BEAST HPD bars and posterior |
| ETE4 (Py; ETE3 paper) | Huerta-Cepas 2016 | rectangular and circular; programmatic per-node `NodeStyle`/faces and layout functions; integrated NCBI taxonomy | styling thousands of nodes by rule, taxonomy-driven annotation, or rendering many trees headlessly in Python |
| iTOL v6 (web) | Letunic 2024 | rectangular, circular, unrooted, very large trees; tab-delimited dataset templates for color strips, heatmaps, bars, binary symbols, clade collapse; SVG/PDF/PNG/EPS export | large trees and template-based metadata with no code; collaborative web-driven polished figures |
| FigTree (desktop GUI) | Rambaut (software) | rectangular, polar, radial; reads BEAST NEXUS node annotations to display HPD bars and posterior on screen; manual rotate/collapse/color | interactive inspection of a single BEAST/MrBayes tree before scripting the final figure; not reproducible, use for exploration not the pipeline |

Decision rule: large tree or no-code web annotation -> iTOL v6. Publication figure with metadata, dual support, HPD bars, or aligned heatmaps -> ggtree + treeio (+ ggtreeExtra for circular rings). Programmatic styling of thousands of nodes or NCBI-taxonomy annotation -> ETE4. Fast topology look from a Python script -> Bio.Phylo. Interactive rooting/inspection of one BEAST tree -> FigTree, then move the reproducible figure to ggtree. The recurring mistake is defaulting to Bio.Phylo because the pipeline is Python, then discovering it cannot draw a circular layout, place a heatmap, or show the BEAST HPD bars that are the whole point.

## What Each Layout Reveals and Hides

| Layout | Reveals | Hides or distorts | Good for | Bad for |
|--------|---------|-------------------|----------|---------|
| Rectangular phylogram | branch lengths read directly off a linear axis; honest distance comparison | tip labels overplot past a few hundred tips; uses vertical space | the default when length matters and tip count is moderate (<~150) | very large trees |
| Slanted / triangular | compactness; quick topology read | the diagonal implies a ladder/progression; imprecise length reading | quick schematic topology views | precise branch-length comparison |
| Circular / fan | hundreds to thousands of tips in one panel; metadata rings | radial distortion -- the same length subtends a larger arc near the rim, so distance is perceptually compressed near the root; near-root structure cramped | big trees where the message is broad structure plus aligned metadata | quantitative branch-length comparison |
| Unrooted / radial | honestly shows NO assumed root or time direction; overall shape and long-branch outliers | no time direction; clade membership harder to trace; clusters over-read as clades | exploratory views, "no root committed" | directional stories; any "basal/early-diverging" claim |
| Chronogram (time-tree) | node timing on a time axis | requires a clock model and calibrations the reader must trust; without HPD bars timing looks falsely precise | dated BEAST/treePL/MCMCtree trees where timing is the message | any tree lacking a clock; never draw without age uncertainty bars |

Circular is still a ROOTED tree bent into a ring (center = root = past); unrooted/radial explicitly refuses a root. Presenting an unrooted radial tree and narrating "X is basal" is a contradiction.

## Bio.Phylo + matplotlib Recipes

Quick text and ASCII inspection, no figure needed:

```python
from Bio import Phylo

tree = Phylo.read('tree.nwk', 'newick')
print(tree)                  # indented text summary
Phylo.draw_ascii(tree)       # ASCII-art diagram, useful in a terminal or log
```

Draw to a vector file (always pass an axes and `do_show=False` for headless/scripted use):

```python
from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read('tree.nwk', 'newick')
tree.ladderize()             # legibility only; ordering carries NO phylogenetic meaning -- say so in the caption

fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(tree, axes=ax, do_show=False)
ax.set_title('Phylogenetic tree (phylogram, branch length = subs/site)')
fig.savefig('tree.pdf', bbox_inches='tight')   # vector: text and lines stay sharp at any size
plt.close(fig)
```

Label tips, and show support with its measure named (never a bare integer):

```python
def tip_only(clade):
    return clade.name if clade.is_terminal() else ''

def support_label(clade):
    # the measure MUST be stated in the legend/caption; here values are bootstrap percentages
    if not clade.is_terminal() and clade.confidence is not None:
        return f'{clade.confidence:.0f}'
    return ''

fig, ax = plt.subplots(figsize=(12, 10))
Phylo.draw(tree, axes=ax, do_show=False, label_func=tip_only, branch_labels=support_label)
ax.set_title('Bootstrap support shown at internal nodes')
fig.savefig('supported_tree.svg', bbox_inches='tight')
plt.close(fig)
```

Color branches by group (convert to phyloXML for native color support):

```python
from Bio.Phylo.PhyloXML import BranchColor

xtree = tree.as_phyloxml()                     # phyloXML carries branch color through draw()
for clade in xtree.find_clades():
    if clade.name and clade.name.startswith('Homo'):
        clade.color = BranchColor.from_name('red')

fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(xtree, axes=ax, do_show=False)
fig.savefig('colored_tree.pdf', bbox_inches='tight')
plt.close(fig)
```

Scale the panel to tip count so labels stay legible, and drop the axis frame:

```python
n_tips = len(tree.get_terminals())
height = max(8, n_tips * 0.25)                  # ~0.25 in/tip keeps ~6-8 pt labels from colliding

fig, ax = plt.subplots(figsize=(10, height))
Phylo.draw(tree, axes=ax, do_show=False)
ax.axis('off')
fig.savefig('scaled_tree.pdf', bbox_inches='tight')
plt.close(fig)
```

For circular/fan/unrooted layouts, metadata heatmaps, dual support, or BEAST HPD bars, Bio.Phylo cannot help -- route to ggtree + treeio (R), ETE4, or iTOL.

## Per-Method Failure Modes

### Drawing Meaningless Branch Lengths as a Phylogram
**Trigger:** Rendering a constraint tree, a supertree, or a tree whose lengths are non-comparable with the horizontal axis proportional to length.
**Mechanism:** The phylogram geometry asserts a quantitative claim about evolutionary distance that the lengths do not support.
**Symptom:** A topology-only or arbitrary-length tree appears to make precise distance statements.
**Fix:** Draw as a cladogram and SAY so (`branch.length='none'` in ggtree); or, if lengths are real, keep them and state the unit (subs/site vs time) plus a scale bar.

### Unlabeled or Mislabeled Support
**Trigger:** A node shows "98" with no legend, or two measures are printed without saying which is which.
**Mechanism:** The reader defaults to assuming bootstrap, but it may be a posterior (much weaker for the same number), an SH-aLRT (cutoff 80), a UFBoot (cutoff 95, not the BP-70 scale), or a TBE.
**Symptom:** A weakly resolved bush reads as a confident comb because the displayed number is over-read.
**Fix:** Always state the measure(s) and their order (e.g. "SH-aLRT/UFBoot" at each node); collapse nodes below threshold into polytomies rather than drawing fake resolution, since bootstrap, posterior, SH-aLRT, and UFBoot sit on different scales.

### Tip-Label Overplotting on Large Trees
**Trigger:** A rectangular layout past a few hundred tips.
**Mechanism:** Horizontal labels collide into an unreadable black band; authors then shrink the font to illegibility or silently drop labels.
**Symptom:** Tip labels are unreadable or missing.
**Fix:** Rotate labels, switch to circular/fan with radial labels, collapse uninformative clades, annotate by colored strips/rings (ggtreeExtra, iTOL datasets) instead of text, or move to iTOL which is built for large trees.

### Non-Monophyly Hidden by Ladderization or Rotation
**Trigger:** Rotating nodes to make a paraphyletic or polyphyletic group look contiguous (or a clade look scattered).
**Mechanism:** Node rotation is information-free with respect to the tree, but the eye reads contiguity as a clade and ordering as direction.
**Symptom:** A group that is not monophyletic appears unified, or a directional narrative is implied.
**Fix:** Color by group and let the topology speak; never narrate ordering as meaning; state "tips ladderized for legibility, ordering carries no phylogenetic meaning."

### Chronogram Without Age Uncertainty
**Trigger:** Drawing a BEAST/MrBayes time-tree with point-estimate node ages and no HPD bars.
**Mechanism:** The 95% HPD intervals on node heights ARE the result; omitting them asserts false precision.
**Symptom:** Node ages look known to the day; a reviewer asks where the credible intervals went.
**Fix:** Draw the HPD bars via treeio `read.beast()` -> ggtree `geom_range('height_0.95_HPD')` (introspect the column name; it varies by source program and treeio version), or enable node bars in FigTree; never draw an annotated Bayesian tree with Bio.Phylo, which drops the annotations silently.

### Raster Export for Publication
**Trigger:** Saving a tree as a 150-dpi PNG for a paper.
**Mechanism:** Raster pixelates thin branches and small tip labels at print size.
**Symptom:** Fuzzy labels and pixelated branches; journal rejection or blurry print.
**Fix:** Export SVG/PDF/EPS (vector) so text and lines stay sharp at any scale; only if forced to rasterize, do so at final size with >=600 dpi for line art.

## Quantitative Thresholds

These are operational defaults for a standard portrait panel with ~6-8 pt labels; adapt to font, page, and journal specs.

| Quantity | Threshold | Rationale / source |
|----------|-----------|--------------------|
| Tips, rectangular horizontal labels comfortable | <=~50 | labels do not collide at legible font |
| Tips, start rotating/shrinking labels | ~50-150 | label height x tip count approaches panel height |
| Tips, switch to circular/fan or collapse clades | ~150-500 | rectangular labels collide; radial labels recover space |
| Tips, per-tip text impractical; use strips/rings | >~500-1000 | annotate by color/groups (iTOL, ggtreeExtra); iTOL is built for large trees |
| Publication export | vector (SVG/PDF/EPS) | text and lines stay sharp at any scale; the default |
| Forced raster, line-art/tree | >=600 dpi (many journals 600-1200) | ~300 dpi is the photo floor but too coarse for thin branches and small labels |
| Scale bar on any phylogram | mandatory (subs/site or a time axis) | without it the reader cannot recover true distances under aspect-ratio distortion |
| Bootstrap strong | >=95 (>=70 moderate) | Hillis-Bull heuristic; state the measure (Felsenstein 1985) |
| UFBoot strong | >=95 (not the BP-70 scale) | recalibrated; IQ-TREE joint rule SH-aLRT>=80 AND UFBoot>=95 |
| SH-aLRT strong | >=80 | IQ-TREE recommendation |
| Posterior probability strong | >=0.95, but WEAKER than bootstrap 95 | PP systematically higher for the same data; do not equate |

Lock the branch-length scale; do not let the figure engine non-uniformly stretch a phylogram to fill a fixed panel, which distorts the lengths the figure exists to communicate.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| BEAST HPD bars absent from a Python figure | drew an annotated tree with Bio.Phylo | route through treeio `read.beast` + ggtree `geom_range` |
| Figure not saving / blank | `do_show=True` opens a window instead of writing | pass `do_show=False`, then `fig.savefig(...)` |
| Branch colors not appearing | plain Newick tree has no color slot | convert with `tree.as_phyloxml()` and set `clade.color` |
| Labels overlap into a black band | too many tips for rectangular layout | increase panel height, rotate labels, or switch to circular/iTOL |
| "Basal" claim on a radial tree | narrated an unrooted layout as if rooted | root explicitly (tree-manipulation) and show the root before any directional claim |
| Support number misread | bare integer with no measure stated | label the measure(s) and order; collapse sub-threshold nodes to polytomies |
| `Phylo.draw` has no circular option | Bio.Phylo is rectangular-only | use ggtree `layout='circular'`, ETE4, or iTOL |

## References

Talevich E, Invergo BM, Cock PJA, Chapman BA. 2012. Bio.Phylo: a unified toolkit for processing, analyzing and visualizing phylogenetic trees in Biopython. *BMC Bioinformatics* 13:209.
Yu G, Smith DK, Zhu H, Guan Y, Lam TT-Y. 2017. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. *Methods in Ecology and Evolution* 8(1):28-36.
Yu G. 2020. Using ggtree to visualize data on tree-like structures. *Current Protocols in Bioinformatics* 69(1):e96.
Wang L-G, Lam TT-Y, Xu S, Dai Z, Zhou L, Feng T, Guo P, Dunn CW, Jones BR, Bradley T, Zhu H, Guan Y, Jiang Y, Yu G. 2020. treeio: an R package for phylogenetic tree input and output with richly annotated and associated data. *Molecular Biology and Evolution* 37(2):599-603.
Xu S, Dai Z, Guo P, Fu X, Liu S, Zhou L, Tang W, Feng T, Chen M, Zhan L, Wu T, Hu E, Jiang Y, Bo X, Yu G. 2021. ggtreeExtra: compact visualization of richly annotated phylogenetic data. *Molecular Biology and Evolution* 38(9):4039-4042.
Huerta-Cepas J, Serra F, Bork P. 2016. ETE 3: reconstruction, analysis, and visualization of phylogenomic data. *Molecular Biology and Evolution* 33(6):1635-1638.
Letunic I, Bork P. 2024. Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool. *Nucleic Acids Research* 52(W1):W78-W82.
Felsenstein J. 1985. Confidence limits on phylogenies: an approach using the bootstrap. *Evolution* 39(4):783-791.

## Related Skills

- tree-io - parsing and preserving BEAST/MrBayes/IQ-TREE annotations so they survive into the figure; the tightest coupling
- tree-manipulation - rooting, ladderizing, pruning, and collapsing low-support nodes before drawing
- data-visualization/ggplot2-fundamentals - the grammar, themes, and ggsave vector export underlying ggtree
- data-visualization/multipanel-figures - composing a tree plus aligned metadata panels into one figure
