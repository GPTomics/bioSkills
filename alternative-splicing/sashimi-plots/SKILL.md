---
name: bio-sashimi-plots
description: Creates sashimi-style plots showing RNA-seq read coverage and splice junction counts using ggsashimi (general-purpose, condition-grouped overlays), rmats2sashimiplot (rMATS-output-aware), MAJIQ-VOILA (LSV posteriors interactive HTML), leafviz (leafcutter clusters Shiny), Jutils (tool-agnostic heatmaps and sashimi for rMATS/leafcutter/SUPPA2/MAJIQ output), or pyGenomeTracks (multi-track publication figures). Tool choice depends on the upstream differential-splicing tool's output format and the publication vs interactive use case. Use when visualizing specific splicing events, validating differential splicing calls, or producing publication-quality figures.
tool_type: python
primary_tool: ggsashimi
---

## Version Compatibility

Reference examples tested with: ggsashimi 1.1+, rmats2sashimiplot 3.0+, MAJIQ 3.0+, leafcutter 0.2.9+, pyGenomeTracks 3.8+, ggplot2 3.5+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Sashimi Plot Visualization

Visualize RNA-seq coverage tracks with splice junction arcs labeled by read count. Sashimi plots originated with MISO (Katz 2010 *Nat Methods*); modern tools differ in input handling, group aggregation logic, and customization. Tool choice is not interchangeable — some tools work only with specific upstream output formats.

## Tool Selection Matrix

| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| ggsashimi | Publication-quality grouped overlays from any BAM | BAMs + region | `--overlay` aggregates samples within a group; clean PDFs | No native rMATS/MAJIQ integration; need to extract coords manually |
| rmats2sashimiplot | One-line plot from rMATS output | rMATS event file + BAMs | No manual coord extraction | rMATS-specific; doesn't handle leafcutter or MAJIQ |
| MAJIQ-VOILA | Interactive LSV browsing with posterior PSI distributions | MAJIQ build + psi/deltapsi | Splice-graph topology; LSV-aware; posterior violins | Static figures; non-academic license |
| leafviz | Cluster-level interactive browsing with NMD annotation | leafcutter differential output | Filter table + sashimi-like plots; NMD-aware | leafcutter-specific |
| Jutils | Unified output across rMATS, leafcutter, SUPPA2, MAJIQ | Tool-specific differential output | Heatmaps, Venn, sashimi tool-agnostically | Output less polished than ggsashimi |
| pyGenomeTracks | Multi-track publication figures (RNA-seq + ChIP/ATAC) | BigWig + BED + GTF | Combine RNA with chromatin tracks | Not splicing-specific; configure tracks manually |
| IGV (interactive) | Quick ad-hoc inspection | BAM + region | Scrollable, instant | Not for publication figures |
| MISO sashimi | Historical | MISO output | Original sashimi format | MISO unmaintained; no longer recommended |

## Decision Tree by Goal

| Goal | Recommended tool |
|------|-------------------|
| Validate a specific rMATS hit | rmats2sashimiplot (one-line) or ggsashimi (custom) |
| Validate a leafcutter cluster | leafviz (interactive) or ggsashimi with cluster coordinates |
| Validate a MAJIQ LSV (complex topology) | MAJIQ-VOILA (only tool that shows full LSV graph) |
| Publication-quality two-condition comparison | ggsashimi `-O 3 -A mean_j` for grouped overlay |
| Multi-track figure (RNA-seq + H3K4me3 + ATAC) | pyGenomeTracks |
| Quick ad-hoc browsing during development | IGV sashimi |
| Tool-agnostic batch heatmap of significant events | Jutils |
| Interactive cohort-level filtering of leafcutter results | leafviz Shiny |

## ggsashimi for Publication Overlays

**Goal:** Generate publication-quality sashimi plot for a region with samples grouped by condition and per-sample tracks aggregated.

**Approach:** Define samples + groups + colors in a TSV (no header), then call ggsashimi with coordinates, GTF, and visual flags.

```python
import subprocess
import pandas as pd

groups = pd.DataFrame({
    'bam': ['ctrl1.bam', 'ctrl2.bam', 'ctrl3.bam', 'trt1.bam', 'trt2.bam', 'trt3.bam'],
    'group': ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment'],
    'color': ['#1f77b4'] * 3 + ['#ff7f0e'] * 3
})
groups.to_csv('sashimi_groups.tsv', sep='\t', index=False, header=False)

subprocess.run([
    'ggsashimi.py',
    '-b', 'sashimi_groups.tsv',
    '-c', 'chr17:43094000-43125000',
    '-o', 'BRCA1_sashimi',
    '-M', '10',
    '--alpha', '0.25',
    '--height', '3',
    '--width', '10',
    '--shrink',
    '--fix-y-scale',
    '--ann-height', '4',
    '-g', 'gencode_v45.gtf',
    '--base-size', '14',
    '-O', '3',
    '-A', 'mean_j',
    '-F', 'pdf'
], check=True)
```

Key ggsashimi flags (Garrido-Martin 2018 *PLoS Comput Biol*):
- `--overlay 3` (or `-O 3`): aggregate multiple samples within a group into a single overlay track with summary statistics — its signature feature
- `-A mean_j`: junction aggregation method (`mean`, `median`, `mean_j` accounts for sample-wise normalization); use `mean_j` for biological replicates
- `--shrink`: rescale long introns (>2x flanking exons) for compact display
- `--fix-y-scale`: identical y-axis across groups (essential for visual comparison)
- `--alpha 0.25`: transparency for per-sample coverage in overlay mode
- `-M 10`: minimum junction reads to display (lower = noisier; 5-10 typical; raise to 20+ for crowded plots)
- `--ann-height`: gene annotation track height
- `-F pdf`: output format (pdf, png, svg, eps)

## Batch Plotting from rMATS Hits

**Goal:** Auto-generate sashimi plots for all significant rMATS differential events.

**Approach:** Parse SE.MATS.JC.txt, expand coordinates to flanking exons + 500nt context, iterate ggsashimi.

```python
import subprocess
import pandas as pd
from pathlib import Path

diff = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')
sig = diff[(diff['FDR'] < 0.05) & (diff['IncLevelDifference'].abs() > 0.10)]

Path('sashimi_plots').mkdir(exist_ok=True)
for idx, ev in sig.head(25).iterrows():
    region = f'{ev["chr"]}:{ev["upstreamES"] - 500}-{ev["downstreamEE"] + 500}'
    safe_name = f'{ev["geneSymbol"]}_{ev["chr"]}_{ev["upstreamES"]}'
    subprocess.run([
        'ggsashimi.py',
        '-b', 'sashimi_groups.tsv',
        '-c', region,
        '-o', f'sashimi_plots/{safe_name}',
        '-M', '5',
        '--shrink',
        '--fix-y-scale',
        '-O', '3',
        '-A', 'mean_j',
        '-g', 'annotation.gtf',
        '-F', 'pdf'
    ], check=True)
```

For MXE events, plot from upstreamES of exon 1 to downstreamEE of exon 2 to show both alternative exons in the same figure.

## rmats2sashimiplot

**Goal:** Plot directly from rMATS event coordinates without manual region calculation.

**Approach:** Pass rMATS event file + BAM lists + event type; rmats2sashimiplot extracts coordinates and produces per-event PDFs.

```bash
rmats2sashimiplot \
    --b1 ctrl1.bam,ctrl2.bam,ctrl3.bam \
    --b2 trt1.bam,trt2.bam,trt3.bam \
    -t SE \
    -e rmats_output/SE.MATS.JC.txt \
    --l1 Control \
    --l2 Treatment \
    -o sashimi_rmats \
    --exon_s 1 \
    --intron_s 5 \
    --color '#1f77b4,#ff7f0e' \
    --group-info group_def.txt
```

`--exon_s 1 --intron_s 5` shrinks intron-to-exon visual ratio 5:1 (introns drawn 1/5 their actual length). The `--group-info` flag (newer versions) allows custom replicate groupings.

## MAJIQ-VOILA Interactive HTML

**Goal:** Browse LSV posterior PSI distributions interactively with splice-graph topology.

**Approach:** Run `voila` on MAJIQ output to generate self-contained HTML.

```bash
# MAJIQ V3 (June 2025+) uses Zarr-format splicegraph (V2's .sql is deprecated)
voila view -p 5000 -j 8 build/splicegraph.zarr psi_output/sample.psi.voila -o voila_psi_html

voila view -p 5000 -j 8 build/splicegraph.zarr deltapsi_output/group1_group2.deltapsi.voila -o voila_dpsi_html
```

VOILA shows:
- Complete LSV graphs (single source / single target nodes)
- Per-junction posterior PSI violin plots
- ΔPSI distributions across all conditions
- Confidence by junction within an LSV

**The only tool that visualizes complex multi-junction LSVs intuitively.** For events that don't fit canonical SE/A5SS/A3SS, VOILA is the visualization of choice.

## leafviz Shiny App

**Goal:** Browse leafcutter clusters with intron-level effects, sashimi-like plots, and NMD annotation.

**Approach:** Prepare leafviz input from leafcutter differential output, then launch Shiny.

```bash
prepare_results.R \
    -o leafviz \
    -m groups.txt \
    leafcutter_perind_numers.counts.gz \
    ds_results_cluster_significance.txt \
    ds_results_effect_sizes.txt \
    annotation_codes
```

```r
library(leafviz)
run_leafviz('leafviz.RData')
```

Standalone alternative: `jackhump/leafviz` GitHub repo for the lightweight installable subset. Useful for cohort-level interactive filtering.

## Jutils for Tool-Agnostic Output

**Goal:** Visualize differential splicing output uniformly across rMATS, leafcutter, SUPPA2, and MAJIQ.

**Approach:** Convert tool output to Jutils' standard format, then plot.

```bash
jutils convert -t rmats -i SE.MATS.JC.txt -o rmats_jutils.tsv
jutils heatmap -i rmats_jutils.tsv -o heatmap.pdf --top 50
jutils sashimi -i rmats_jutils.tsv -b sashimi_groups.tsv -g annotation.gtf -o sashimi_jutils/
jutils venn -i rmats_jutils.tsv leafcutter_jutils.tsv -o overlap_venn.pdf
```

(Yang 2021 *Bioinformatics*) Useful when comparing multiple tools' outputs across publications or doing meta-analysis.

## pyGenomeTracks for Multi-Track Figures

**Goal:** Combine splicing with chromatin or coverage tracks for publication figures.

**Approach:** Define tracks in an INI file (genes, BAM, BigWig, BED), then run `pyGenomeTracks --tracks tracks.ini --region ... -o figure.pdf`.

```ini
[gene_models]
file = annotation.gtf
height = 3
title = GENCODE v45
fontsize = 10
file_type = gtf

[ctrl_coverage]
file = ctrl_merged.bw
title = Control
color = #1f77b4
height = 3
file_type = bigwig

[trt_coverage]
file = trt_merged.bw
title = Treatment
color = #ff7f0e
height = 3
file_type = bigwig

[junctions]
file = junctions.bedpe
title = Junctions
height = 2
file_type = links
links_type = arcs
```

The `junctions.bedpe` file must be in **BEDPE format** (6 columns: chr1 start1 end1 chr2 start2 end2 [+ optional score]). Convert from regtools .bed12 junctions:

```bash
# Convert regtools junctions BED12 to BEDPE for pyGenomeTracks.
# regtools BED12 column 11 is blockSizes (anchor_left, anchor_right);
# column 12 is blockStarts (0, intron_length + anchor_left).
# Intron start = chromStart + anchor_left = $2 + a[1]
# Intron end   = chromStart + blockStarts[2] = $2 + b[2]
awk 'BEGIN{OFS="\t"} {split($11,a,","); split($12,b,","); s=$2+a[1]; e=$2+b[2]; print $1, s, s+1, $1, e-1, e, $5}' \
    regtools_junctions.bed > junctions.bedpe
```

```bash
pyGenomeTracks --tracks tracks.ini --region chr17:43094000-43125000 -o figure.pdf
```

## Reading Sashimi Plots (Interpretation Guide)

| Visual element | What it represents |
|----------------|--------------------|
| Filled coverage track | Read coverage at each genomic position (depth-normalized in `-A` mode) |
| Arc / curve between exons | Junction-spanning reads; arc connects donor to acceptor |
| Number on arc | Count of junction-spanning reads (raw, not normalized, unless `-A` set) |
| Arc thickness | Often proportional to read count (tool-dependent) |
| Gene model below | Exons (boxes) and introns (lines) from GTF |
| Multiple parallel tracks | Per-sample (default) or per-group (with `-O`) |

**Junction count interpretation:** the number on an arc is the absolute count of reads whose CIGAR string contained an `N` operation matching that intron coordinate. Higher = more usage. Compare counts on inclusion vs skipping arcs to estimate PSI visually.

**Color convention:** by convention, control = blue (`#1f77b4`), treatment = orange (`#ff7f0e`); always document. Use ColorBrewer or matplotlib defaults for >2 groups.

## Per-Tool Failure Modes

### ggsashimi: Off-Strand Junction Artifacts

**Trigger:** Stranded RNA-seq library plotted without strand specification.

**Mechanism:** ggsashimi reads BAM strand from CIGAR + flag; without strand info, antisense junctions appear as artifacts.

**Symptom:** Implausible junctions in regions with overlapping antisense genes; "noise" arcs at unexpected locations.

**Fix:** Use `-S RF` (reverse-forward) for Illumina TruSeq stranded; verify with RSeQC `infer_experiment.py`. Alternatively, pre-filter BAM by strand with `samtools view -f 16` / `-F 16`.

### rmats2sashimiplot: Wrong Coordinate Convention

**Trigger:** Older versions or non-default rMATS output.

**Mechanism:** rmats2sashimiplot expects 1-based coordinates from rMATS' .MATS.JC.txt; rMATS outputs 0-based half-open in some columns.

**Symptom:** Plot region shifted by 1 nt; arcs misaligned with gene model.

**Fix:** Verify rmats2sashimiplot version matches rMATS-turbo output convention; use ggsashimi for cleaner control.

### MAJIQ-VOILA: Browser Memory

**Trigger:** Loading large VOILA HTML in browser (cohort with hundreds of LSVs).

**Mechanism:** VOILA HTML embeds all LSV data; large cohorts produce >100 MB HTMLs.

**Symptom:** Browser unresponsive on opening; "page unresponsive" warnings.

**Fix:** Filter LSVs in MAJIQ before voila step (`--changing-pvalue-threshold 0.95` and `--changing-between-group-dpsi-threshold 0.2`); split into per-gene HTMLs.

### leafviz: Annotation Codes Mismatch

**Trigger:** Using leafviz with annotation_codes from different GENCODE version than leafcutter clusters.

**Mechanism:** annotation_codes encodes intron-to-event-class mapping per GTF version.

**Symptom:** Many clusters show as "unannotated" despite being in canonical GTF.

**Fix:** Generate annotation_codes from the same GENCODE version used in differential analysis.

## Customization Reference

| Visual goal | ggsashimi flag |
|-------------|-----------------|
| Reduce intron whitespace | `--shrink` |
| Identical y-axis across groups | `--fix-y-scale` |
| Per-group overlay aggregation | `-O 3 -A mean_j` |
| Larger figure | `--width 12 --height 4` |
| Bigger fonts | `--base-size 16` |
| Vector output | `-F pdf` or `-F svg` |
| Custom palette | Edit colors in groups TSV |
| Filter junction noise | `-M 10` (raise to 20+) |
| Transparency | `--alpha 0.25` |
| GTF feature filter | `--gtf-filter protein_coding` |

## Best Practices

| Tip | Rationale |
|-----|-----------|
| Use `--shrink` for genes with large introns | Keeps exons visible (TTN, brain genes with multi-kb introns) |
| `--fix-y-scale` for cross-group comparisons | Otherwise auto-rescaling visually exaggerates differences |
| Aggregate replicates with `-O 3 -A mean_j` | Reduces clutter; per-sample variance still shown via alpha |
| Limit to 3-4 groups per figure | More becomes hard to read |
| Include 200-500 nt flanking exons | Show full splicing context |
| For MXE events, plot both alternative exons | Otherwise only half of the event is visible |
| Check accessibility colors | Use ColorBrewer-safe palettes for color-blind readers |
| Always include a legend | Sashimi figures without legends are uninformative for non-experts |
| Specify output format explicitly | PDF for publication; PNG for slides; SVG for editing |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `ggsashimi: 'samtools' not found` | samtools not in PATH | Install via conda; `which samtools` to verify |
| `ggsashimi: empty plot` | Region has no reads or wrong chromosome name | Check BAM with `samtools view sample.bam chr1:100-200`; chrom name match (chr1 vs 1) |
| `rmats2sashimiplot: KeyError 'IJC_SAMPLE_1'` | Old rmats2sashimiplot with new rMATS output | Update both to matching versions |
| `voila: out of memory` | Large LSV cohort | Filter by deltapsi threshold before voila |
| `pyGenomeTracks: ini parse error` | Missing closing bracket or invalid track type | Validate INI syntax; check `pyGenomeTracks --listTracks` for supported types |
| `leafviz: missing exon file` | annotation_codes path wrong | Re-run `prepare_results.R` with correct paths |

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No junctions shown | Default `-M 10` too strict | Lower to `-M 3` or `-M 5` |
| Plot too crowded | Many samples without aggregation | Use `-O 3` to overlay groups |
| Annotation missing or wrong gene | GTF lacks gene_name attribute or wrong build | Verify GTF version vs BAM reference; switch to `--gtf-filter protein_coding` |
| Memory issues on large regions | >100 kb regions with many samples | Plot smaller windows or pre-extract reads with samtools view |
| Y-axis dominated by one peak | Outlier sample | Use `-A mean_j` to aggregate; or filter outlier |

## Related Skills

- differential-splicing - Identify events to plot; sashimi plots are validation
- splicing-quantification - Context for PSI values; sashimi provides visual confirmation
- data-visualization/genome-browser-tracks - Multi-track figure design (pyGenomeTracks)
- data-visualization/ggplot2-fundamentals - ggsashimi customization (extends ggplot2)
- data-visualization/color-palettes - Accessible color choices
- data-visualization/specialized-omics-plots - Volcano, heatmap complement to sashimi

## References

- Katz et al 2010 *Nat Methods* - MISO sashimi plot original
- Garrido-Martin et al 2018 *PLoS Comput Biol* - ggsashimi
- Yang et al 2021 *Bioinformatics* - Jutils
- Vaquero-Garcia et al 2016 *eLife* - MAJIQ / VOILA
- Li et al 2018 *Nat Genet* - leafcutter / leafviz
- Ramirez et al 2018 *Nat Commun* - pyGenomeTracks
