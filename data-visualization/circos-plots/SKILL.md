---
name: bio-data-visualization-circos-plots
description: Build circular genome visualizations using circlize (R), pyCirclize (Python), or Circos (Perl CLI) with ideogram tracks, multi-data tracks (scatter, histogram, heatmap), chord/link arcs for interactions, and explicit circos.clear() between plots. Covers when circular is appropriate vs when Cartesian wins (Cleveland-McGill 1984), karyograms, and chromosome adjacency in chord diagrams. Use when adjacency on the circle conveys meaning — chromosome-level overview, structural variants, Hi-C interactions, cross-genome comparisons.
tool_type: mixed
primary_tool: circlize
---

## Version Compatibility

Reference examples tested with: circlize 0.4.16+ (R), pyCirclize 1.4+ (Python), Circos 0.69-9 (Perl CLI), ComplexHeatmap 2.18+ (uses circlize for color mapping).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Circular Genome Plots (Circos)

**"Make a circos plot"** -> Render genome chromosomes around a circle with stacked tracks (histogram, scatter, heatmap) and arcs/chords showing interactions. Krzywinski 2009 *Genome Res* 19:1639 introduced the genre for genome-scale comparative views. The single decision that matters: **does the circular layout convey meaning that Cartesian cannot?**

- R: `circlize::circos.initializeWithIdeogram` + `circos.genomicTrack*` (Gu 2014)
- Python: `pyCirclize.Gcircle`
- CLI: Circos (Perl); config-driven; most flexible but steepest learning

## The Single Most Important Modern Insight -- Circular Plots Often Hide What Cartesian Reveals

Cleveland-McGill 1984 *J Am Stat Assoc* 79:531 effectiveness rankings establish that **position-on-common-scale (Cartesian) is the most accurate visual channel**; circular position requires mental "unwrapping" and impairs precise value comparison. Heer-Bostock 2010 *CHI* replicated the ranking in modern crowd studies. Use circular only when adjacency on the circle conveys meaning that linear cannot.

Use circular ONLY when:
- **Chromosome adjacency matters** (whole-genome SVs, Hi-C contacts where genome circularity is the geometry)
- **Pairwise interactions between many entities** (chord diagrams; chromosome translocations)
- **Aesthetic / overview** infographic for cover figure

Do NOT use circular for:
- Comparing values across categories (Cartesian bar/dot wins)
- Time series (linear axis wins)
- Anything where precise value reading matters

The circos plot is a beautiful but dangerous default. The most-cited published critique is the genre being applied where it adds no information.

## circlize (R) — Modern Default

**Goal:** Render a multi-track circos plot with ideograms, gene-density histogram, variant-density heatmap, and inter-chromosomal SV links.

**Approach:** Initialize with chromosome ideograms via `circos.initializeWithIdeogram`; add tracks with `circos.genomicTrack` + appropriate panel function; add links with `circos.link`; **always call `circos.clear()` after the plot completes.**

```r
library(circlize)

# 1. Initialize with hg38 ideograms
pdf('circos.pdf', width = 8, height = 8)
circos.par(start.degree = 90,             # 12 o'clock start
           gap.degree = c(rep(1, 23), 5)) # bigger gap before chr1 for visual break
# hg38 specifically benefits from explicit chromosome.index to skip unmapped contigs
circos.initializeWithIdeogram(species = 'hg38',
                               chromosome.index = paste0('chr', c(1:22, 'X', 'Y')),
                               plotType = c('axis', 'labels', 'ideogram'))

# 2. Gene-density histogram (outermost data track)
circos.genomicDensity(gene_bed, col = '#0072B2', track.height = 0.08)

# 3. Variant-density heatmap
circos.genomicHeatmap(variant_bed,
                       col = colorRamp2(c(0, 100), c('white', '#D55E00')),
                       heatmap_height = 0.08, side = 'inside')

# 4. CNV scatter
circos.genomicTrack(cnv_bed, ylim = c(-2, 2),
                     panel.fun = function(region, value, ...) {
                         circos.genomicPoints(region, value,
                                              col = ifelse(value > 0.3, '#D55E00',
                                                           ifelse(value < -0.3, '#0072B2', 'grey60')),
                                              pch = 16, cex = 0.4)
                     },
                     track.height = 0.1)

# 5. Inter-chromosomal SV links
for (i in seq_len(nrow(sv_df))) {
    circos.link(sv_df$chr1[i], c(sv_df$start1[i], sv_df$end1[i]),
                sv_df$chr2[i], c(sv_df$start2[i], sv_df$end2[i]),
                col = '#888888', lwd = 0.4)
}

# 6. CRITICAL -- clear global state
circos.clear()
dev.off()
```

## The `circos.clear()` Trap

`circos.par()` settings (start.degree, gap.degree, canvas.xlim, canvas.ylim, clock.wise, circle.margin) are GLOBAL state. After a plot completes, those settings persist into the next plot.

Forgetting `circos.clear()` produces:
- Next `circos.par()` calls silently fail to take effect (warning, easily missed in loops)
- Re-initialization may error or render at wrong angles
- Loop-rendered figures inherit state from the previous iteration

**Always call `circos.clear()` after every plot.** Make it the last line of the plotting block alongside `dev.off()`.

## pyCirclize (Python)

```python
from pycirclize import Circos
import matplotlib.pyplot as plt

sectors = {'chr1': 248956422, 'chr2': 242193529, ...}
circos = Circos(sectors, space=2)                       # space = degree gap between sectors

for sector in circos.sectors:
    sector.text(sector.name, r=110, size=8)
    # outer ideogram
    sector.axis(r_lim=(95, 100), fc='lightgrey')
    # data track
    track = sector.add_track((75, 90))
    track.bar(positions, heights, width=bin_size, color='#0072B2')

# Inter-sector links (chord diagram)
circos.link(('chr1', 1e8, 1.1e8), ('chr5', 2e8, 2.1e8),
            color='#888888', alpha=0.5)

fig = circos.plotfig()
fig.savefig('circos_py.pdf', bbox_inches='tight')
```

pyCirclize is a younger package than circlize but actively developed (Shimoyama 2024+). API more Pythonic than circlize-via-rpy2.

## Circos CLI (Perl) — Most Powerful, Steepest Curve

```bash
# config: circos.conf with karyotype, ideogram, plots, links sections
circos -conf circos.conf -outputfile output.png
```

Circos (Krzywinski 2009) is the original; supports unlimited tracks and arbitrary geometries via configuration. For publication-grade complex figures the Perl tool remains the most powerful. For Python/R workflows, circlize/pyCirclize are more accessible.

## Decision Tree by Use Case

| Use case | Recommended | Why |
|----------|-------------|-----|
| Whole-genome CNV summary | circlize/pyCirclize | Standard genre |
| SV link diagram | Chord arcs in circos | Inter-chromosomal adjacency |
| Hi-C contact summary at chromosome level | circos heatmap track | Adjacency matters |
| Per-sample mutation overview | Circular karyogram | Aesthetic; comparable to OncoPrint |
| Cohort-wide gene expression comparison | NOT circular | Use heatmap (Cartesian wins) |
| Time-series of any kind | NOT circular | Use line plot |
| Pathway diagram | NOT circular | Use Cytoscape |

## Ideogram + Karyogram Without Circos

For per-chromosome data display where circularity is not required, `karyoploteR` (Gel 2017 *Bioinformatics* 33:3088) renders linear ideograms with stacked data tracks — often the better choice for CNV per-chromosome views.

```r
library(karyoploteR)
kp <- plotKaryotype(genome = 'hg38', chromosomes = c('chr1', 'chr7', 'chr17'))
kpAddBaseNumbers(kp)
kpLines(kp, data = cnv_gr, y = cnv_gr$log2)
kpAddCytobandLabels(kp)
```

See `copy-number/cnv-visualization` for karyoploteR in depth.

## Per-Method Failure Modes

### circos.clear() forgotten in a loop

**Trigger:** Plotting multiple circos figures in a `for` loop without `circos.clear()` between.

**Mechanism:** circos.par settings (gap.degree, start.degree, clock.wise) persist across plots.

**Symptom:** Plots 2..N inherit state from plot 1; gap sizes, rotation differ unexpectedly.

**Fix:** End every plot block with `circos.clear()`. Make it a hygiene rule.

### Using circular when Cartesian would be better

**Trigger:** "Circos plot of gene expression across 20 conditions."

**Mechanism:** Circular impairs value comparison (Cleveland-McGill 1984; Heer-Bostock 2010).

**Symptom:** Reviewer or coauthor says "I can't tell which condition is highest."

**Fix:** Use clustered heatmap. Reserve circos for genome-adjacency or chord-diagram use cases.

### Too many links produce a black blob

**Trigger:** Plotting 10000+ chord links between chromosomes.

**Mechanism:** Overlap saturates the center; no individual link visible.

**Symptom:** Center of circos is uniformly dark.

**Fix:** Filter to top-confidence links; OR color-bin by interaction strength with alpha; OR aggregate to chromosome-level summary then link.

### Sector ordering arbitrary

**Trigger:** Default sector order is input order.

**Mechanism:** circlize / pyCirclize do not auto-order chromosomes 1..22, X, Y.

**Symptom:** Chromosomes appear in genome-build-file order.

**Fix:** Explicit `chromosome.index = c(paste0('chr', 1:22), 'chrX', 'chrY')`.

### Wrong species ideogram

**Trigger:** `species = 'hg19'` when data is hg38-coordinate.

**Mechanism:** circlize fetches cytoband data per species; mismatch renders correct ideogram but wrong banding for the data.

**Symptom:** Cytoband boundaries don't match published references.

**Fix:** Match `species` to data coordinate system. For non-standard genomes, supply custom cytoband file. For `species = 'hg38'` specifically, always pass `chromosome.index = paste0('chr', c(1:22, 'X', 'Y'))` to skip unmapped contigs (jokergoo/circlize issue #46).

### Ideogram covers data track

**Trigger:** Default ideogram track height too large; data track squeezed.

**Mechanism:** circos.initializeWithIdeogram uses ~5% of radius; left-over for data.

**Symptom:** Data values invisible because track is too narrow.

**Fix:** Reduce `cytoband.height` in initialization; OR use `plotType = c('axis', 'labels')` to omit ideogram entirely.

### Chromosome label collisions for small chromosomes

**Trigger:** Default label position; small chromosomes (chr21, chr22, chrY) have overlapping labels.

**Mechanism:** Labels drawn at sector midpoints regardless of sector width.

**Symptom:** Labels overlap.

**Fix:** `circos.par(gap.degree = c(rep(1, 22), 10, 10, 10))` for larger gaps before small chromosomes; OR reduce label font size.

## Reconciliation

| Pattern | Cause | Action |
|---------|-------|--------|
| circlize and pyCirclize differ in default rotation | Different start-angle convention | Set `start.degree=90` (R) / equivalent (Python) explicitly |
| Cytoband colors don't match UCSC | Different species cytoband source | Verify species; for custom genomes supply band file |
| Inter-sector links arc the "long way around" | Default arc direction | Some chord packages support `direction = 'short'` |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Max sectors readable | ~30 | Visualization practical |
| Max links before blob | ~5000 | Practical; depends on alpha |
| Cytoband.height default | 0.05 of radius | circlize default |
| When circular adds value | Adjacency-meaningful only | Cleveland-McGill 1984 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Subsequent plots use wrong rotation | `circos.clear()` forgotten | Always end with `circos.clear()` |
| Chromosomes out of order | Default = input order | Explicit chromosome.index |
| Cytoband mismatch | Wrong species | Match species to data coords |
| Center of circos black | Too many links | Filter or aggregate |
| Reviewer asks "why circular?" | Cartesian would have been clearer | Migrate to heatmap unless adjacency matters |
| Small-chromosome label overlap | Default label position | Larger gap.degree before small sectors |

## References

- Cleveland WS, McGill R. 1984. Graphical perception: theory, experimentation, and application to the development of graphical methods. *J Am Stat Assoc* 79(387):531-554.
- Gel B, Serra E. 2017. karyoploteR: an R/Bioconductor package to plot customizable genomes. *Bioinformatics* 33(19):3088-3090.
- Gu Z, Gu L, Eils R, Schlesner M, Brors B. 2014. circlize implements and enhances circular visualization in R. *Bioinformatics* 30(19):2811-2812.
- Heer J, Bostock M. 2010. Crowdsourcing graphical perception: using Mechanical Turk to assess visualization design. *Proc CHI* 203-212.
- Krzywinski M, Schein J, Birol I, et al. 2009. Circos: an information aesthetic for comparative genomics. *Genome Res* 19(9):1639-1645.
- Shimoyama Y. 2024. pyCirclize: Circular visualization in Python. *GitHub* https://github.com/moshi4/pyCirclize

## Related Skills

- copy-number/cnv-visualization - karyoploteR linear alternative for CNV
- variant-calling/structural-variant-calling - SV data for circos links
- hi-c-analysis/hic-visualization - Hi-C contact data circular display
- data-visualization/genome-tracks - Linear track alternative
- data-visualization/color-palettes - Sector and link palettes
