---
name: bio-copy-number-cnv-visualization
description: Visualize copy number profiles, segments, allele-specific tracks, and cohort patterns from CNVkit, GATK, ASCAT, FACETS, Sequenza, and other callers. Covers genome-wide and per-chromosome log2 scatter plots, B-allele-frequency/minor-allele-fraction tracks, ideograms, cohort heatmaps, circos views, and caller-native plots. Use when creating publication CNV figures, choosing which plot answers a given question, diagnosing a wrong diploid baseline visually, displaying loss of heterozygosity, or deciding what depth-only plots cannot reveal.
tool_type: mixed
primary_tool: matplotlib
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, pandas 2.2+, numpy 1.26+, seaborn 0.13+, CNVkit 0.9.10+, GATK 4.5+; R 4.3+ with ggplot2 3.5+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show matplotlib pandas` then `help(function)` for signatures
- R: `packageVersion('ggplot2')` then `?function_name`
- CLI: `cnvkit.py version`, `gatk --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example rather than retrying.

# CNV Visualization

**"Plot my copy number profile"** -> A CNV figure is an argument, not a picture. The plot type, the y-axis quantity, and where the diploid baseline sits all determine what the reader can conclude. The single most important rule: a depth-only log2 plot cannot show loss of heterozygosity, cannot show tumor purity, and silently misleads if the diploid baseline is centered on a non-diploid mode.

- CLI: `cnvkit.py scatter` / `diagram` / `heatmap`; `gatk PlotModeledSegments`
- Python: `matplotlib` for custom genome-wide and allele-specific tracks
- R: `ggplot2`, `karyoploteR` for publication ideograms

## Plot Selection — What Each View Reveals and Hides

| Plot | Answers | Reveals | Cannot show |
|------|---------|---------|-------------|
| Genome-wide log2 scatter + segments | Where are the gains/losses? | Focal vs broad events, noise level | LOH, purity, allele-specific state |
| Per-chromosome scatter | Is this focal event real and where are its boundaries? | Breakpoints, bin support, weight | Absolute CN without purity |
| BAF / minor-allele-fraction track | Is there allelic imbalance / LOH? | CN-neutral LOH, mirrored imbalance | Total copy number alone |
| Combined log2 + BAF (two-panel) | What is the allele-specific state? | Gains vs CN-LOH vs balanced | — (this is the complete view) |
| Cohort heatmap | What is recurrent across samples? | Shared arm/focal events | Per-sample breakpoint detail |
| Ideogram / diagram | Where do events sit relative to cytobands/genes? | Gene-level context | Quantitative amplitude |
| Circos | Genome-wide CNV + SV breakpoints together | CNV-SV co-localization | Fine amplitude detail |
| Caller-native (GATK/ASCAT/FACETS) | Did the caller fit correctly? | Model fit, segment confidence | — (diagnostic, not publication) |

The decision rule: if the biological question involves LOH, allele-specific gain, or whole-genome doubling, a log2-only plot is insufficient — pair it with a BAF track.

## CNVkit Built-in Plots

```bash
cnvkit.py scatter sample.cnr -s sample.cns -o scatter.png            # genome-wide
cnvkit.py scatter sample.cnr -s sample.cns -c chr17 -o chr17.png      # one chromosome
cnvkit.py scatter sample.cnr -s sample.cns -v sample.vcf.gz -o baf.png  # with BAF panel
cnvkit.py diagram sample.cnr -s sample.cns -o diagram.pdf             # ideogram
cnvkit.py heatmap cohort/*.cns -d -o cohort_heatmap.pdf              # cohort, desaturated
```

Passing `-v` with a VCF adds a B-allele-frequency panel — use it whenever LOH matters.

## Genome-Wide log2 Profile with Segments

**Goal:** Render a publication genome-wide CNV profile with colored segments.

**Approach:** Map per-bin log2 to cumulative genomic coordinates, plot bins as faint points, overlay segment medians as colored horizontal lines, mark chromosome boundaries.

```python
import pandas as pd
import matplotlib.pyplot as plt

def plot_genome_profile(cnr_file, cns_file, output=None, gain=0.3, loss=-0.3):
    '''Genome-wide log2 scatter with segment overlay.'''
    cnr = pd.read_csv(cnr_file, sep='\t')
    cns = pd.read_csv(cns_file, sep='\t')
    chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

    offsets, cum = {}, 0
    for c in chroms:
        sub = cnr[cnr['chromosome'] == c]
        if sub.empty:
            continue
        offsets[c] = cum
        cum += sub['end'].max()

    cnr = cnr[cnr['chromosome'].isin(offsets)].copy()
    cnr['x'] = cnr.apply(lambda r: offsets[r['chromosome']] + r['start'], axis=1)

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.scatter(cnr['x'], cnr['log2'], s=1, c='0.7', alpha=0.5, rasterized=True)
    for _, seg in cns.iterrows():
        if seg['chromosome'] not in offsets:
            continue
        x0 = offsets[seg['chromosome']] + seg['start']
        x1 = offsets[seg['chromosome']] + seg['end']
        color = 'red' if seg['log2'] > gain else 'blue' if seg['log2'] < loss else '0.3'
        ax.hlines(seg['log2'], x0, x1, colors=color, linewidth=2.5)
    for c, x in offsets.items():
        ax.axvline(x, color='0.9', linewidth=0.5)
    ax.axhline(0, color='black', linewidth=0.6)
    ax.set_ylim(-2, 2)
    ax.set_ylabel('log2 copy ratio')
    ax.set_xlabel('genomic position')
    fig.tight_layout()
    if output:
        fig.savefig(output, dpi=200)
    return fig, ax
```

## Combined log2 + B-Allele-Frequency Panel

**Goal:** Show total copy number and allelic imbalance together so CN-neutral LOH and allele-specific gains are visible.

**Approach:** Stack two axes — log2 on top, BAF below. Mirror BAF about 0.5 so allelic imbalance reads as deviation from the center line.

```python
import numpy as np

def plot_log2_baf(cnr_file, baf_df, output=None):
    '''Two-panel plot: log2 copy ratio above, B-allele frequency below.
    baf_df: columns chromosome, position, baf (germline-het sites only).'''
    cnr = pd.read_csv(cnr_file, sep='\t')
    fig, (ax_cn, ax_baf) = plt.subplots(2, 1, figsize=(16, 6), sharex=True)

    ax_cn.scatter(range(len(cnr)), cnr['log2'], s=1, c='0.6', alpha=0.5)
    ax_cn.axhline(0, color='black', linewidth=0.6)
    ax_cn.set_ylabel('log2 ratio')
    ax_cn.set_ylim(-2, 2)

    # Plot BAF and its mirror; a tight band at 0.5 = balanced, split bands = imbalance/LOH
    ax_baf.scatter(range(len(baf_df)), baf_df['baf'], s=2, c='0.4', alpha=0.5)
    ax_baf.scatter(range(len(baf_df)), 1 - baf_df['baf'], s=2, c='0.4', alpha=0.5)
    ax_baf.axhline(0.5, color='black', linewidth=0.6)
    ax_baf.set_ylabel('B-allele frequency')
    ax_baf.set_ylim(0, 1)
    ax_baf.set_xlabel('het SNP index')
    fig.tight_layout()
    if output:
        fig.savefig(output, dpi=200)
    return fig
```

A copy-neutral LOH region shows log2 ~ 0 but BAF splitting away from 0.5 — invisible on any log2-only plot.

## Cohort Heatmap

**Goal:** Show recurrent CNV patterns across a cohort.

**Approach:** Resample every sample's segments onto a common genomic bin grid, stack into a samples-by-bins matrix, render with a diverging colormap centered at zero.

```python
import seaborn as sns

def plot_cohort_heatmap(cns_files, bin_size=1_000_000, output=None):
    '''Recurrent-CNV heatmap across a cohort on a uniform bin grid.'''
    chroms = [f'chr{i}' for i in range(1, 23)]
    columns = []
    for c in chroms:
        columns += [(c, b) for b in range(0, 250_000_000, bin_size)]
    matrix = {}
    for f in cns_files:
        name = f.split('/')[-1].replace('.cns', '')
        cns = pd.read_csv(f, sep='\t')
        row = {}
        for c, b in columns:
            hits = cns[(cns['chromosome'] == c) &
                       (cns['start'] < b + bin_size) & (cns['end'] > b)]
            row[(c, b)] = hits['log2'].mean() if not hits.empty else 0.0
        matrix[name] = row
    df = pd.DataFrame(matrix).T
    fig, ax = plt.subplots(figsize=(14, max(4, 0.3 * len(cns_files))))
    sns.heatmap(df, cmap='RdBu_r', center=0, vmin=-1.5, vmax=1.5,
                xticklabels=False, ax=ax)
    ax.set_xlabel('genomic bin')
    ax.set_ylabel('sample')
    if output:
        fig.savefig(output, dpi=200, bbox_inches='tight')
    return fig
```

## Caller-Native Diagnostic Plots

```bash
# GATK: denoised ratios + modeled segments with allelic info
gatk PlotModeledSegments --denoised-copy-ratios tumor.denoisedCR.tsv \
    --allelic-counts tumor.hets.tsv --segments tumor.modelFinal.seg \
    --sequence-dictionary reference.dict --output-prefix tumor -O plots/
```

ASCAT (`ascat.runAscat` ASPCF and sunrise plots), Sequenza (`sequenza.results` chromosome view and the cellularity/ploidy contour), and FACETS (`plotSample`) emit diagnostic plots — always inspect these to confirm the purity/ploidy fit before trusting downstream calls. They are diagnostic, not publication, figures.

## Failure Modes

### The diploid-baseline centering trap

**Trigger:** Plotting log2 from a hyper-aneuploid or whole-genome-doubled tumor with the y-axis centered on the data median or mode.

**Mechanism:** If most of the genome is at tetraploid baseline, centering on the mode places 4 copies at log2 0. Every true diploid region then appears deleted and the plot tells the opposite story.

**Symptom:** A genome-wide pattern of "loss" (or "gain") inconsistent with the BAF track or with the caller's ploidy estimate.

**Fix:** Anchor the y-axis baseline to the caller's ploidy estimate, not the data mode. Always show a BAF track alongside; if BAF says balanced where log2 says deleted, the centering is wrong.

### log2 axis presented as if it were absolute copy number

**Trigger:** Labeling a log2 y-axis "copy number" or comparing log2 amplitudes across samples of different purity.

**Mechanism:** log2 ratio compresses with decreasing purity — a true CN=4 amplification at 40% purity has roughly half the log2 amplitude of the same event at 80% purity.

**Symptom:** A real amplification looks weaker than a passenger gain in a purer sample; cross-sample amplitude comparisons are meaningless.

**Fix:** For cross-sample or absolute claims, plot integer copy number from a purity-corrected caller, not raw log2. Label log2 axes "log2 copy ratio".

### Cohort heatmap binning erases focal events

**Trigger:** Large uniform bins (e.g. 1-3 Mb) in a cohort heatmap.

**Mechanism:** A focal amplification (e.g. a few hundred kb at MYC or ERBB2) is averaged with flanking neutral sequence and disappears.

**Symptom:** Known recurrent focal drivers absent from the heatmap; only arm-level events visible.

**Fix:** Use a bin size matched to the question — Mb bins for arm-level surveys, gene-centric or GISTIC peak regions for focal drivers. Consider a separate gene-level panel.

## Quantitative Thresholds

| Choice | Value | Rationale |
|--------|-------|-----------|
| log2 y-axis range | -2 to 2 | Covers homozygous loss to ~8-copy gain; clip extreme amplicons separately |
| Gain/loss plot coloring | log2 > 0.3 / < -0.3 | Visual convention (no single primary citation); not a calling threshold (see cnvkit-analysis) |
| Cohort heatmap bin (arm-level) | ~1 Mb | Balances resolution and matrix size |
| Cohort heatmap colormap center | 0 | Diverging map must be zero-centered or gains/losses are not comparable |
| Rasterize scatter points | yes, for > ~50k bins | Keeps vector PDFs openable; segments stay vector |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Whole genome looks lost or gained | y-axis centered on a non-diploid mode | Anchor baseline to caller ploidy; add a BAF track |
| LOH region not visible | log2-only plot | Add a BAF / minor-allele-fraction panel |
| Focal driver missing from heatmap | Bins too large | Use gene-level or GISTIC-peak bins |
| Giant unopenable PDF | 100k+ vector scatter points | `rasterized=True` on the scatter |
| Chromosomes out of order / overlapping | String-sorted contig names | Explicit chromosome order list |
| Amplitudes incomparable across samples | Plotting raw log2 across mixed purity | Plot purity-corrected integer CN |

## References

- Talevich E et al 2016. CNVkit: genome-wide copy number detection from targeted DNA sequencing. PLoS Comput Biol 12:e1004873
- Van Loo P et al 2010. Allele-specific copy number analysis of tumors. PNAS 107:16910 (BAF interpretation)
- Gel B, Serra E 2017. karyoploteR: an R/Bioconductor package to plot customizable genomes. Bioinformatics 33:3088

## Related Skills

- copy-number/cnvkit-analysis - Generates the .cnr/.cns inputs and built-in plots
- copy-number/gatk-cnv - GATK denoised ratios and modeled-segment plots
- copy-number/allele-specific-copy-number - Source of BAF/MAF tracks and ploidy estimates
- copy-number/recurrent-cnv - Cohort-level recurrence underlying heatmaps
- data-visualization/ggplot2-fundamentals - General publication-figure grammar
- data-visualization/circos-plots - Circular genome layouts for CNV + SV
