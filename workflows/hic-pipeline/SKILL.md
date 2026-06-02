---
name: bio-workflows-hic-pipeline
description: End-to-end Hi-C analysis workflow from FASTQ to compartments, TADs, and loops, with the decision of WHICH features the sequencing depth can support. Covers pairtools read-pair processing and library QC, cooler matrices, ICE balancing and distance-decay expected, A/B compartments, TAD boundaries, loop calling, and the routing of HiChIP/PLAC-seq/Capture Hi-C to protein-directed loop callers. Use when processing Hi-C data end to end, deciding a resolution for a given depth, or choosing between bulk-Hi-C and protein-directed loop calling.
tool_type: mixed
primary_tool: cooler
workflow: true
depends_on:
  - hi-c-analysis/contact-pairs
  - hi-c-analysis/hic-data-io
  - hi-c-analysis/matrix-operations
  - hi-c-analysis/compartment-analysis
  - hi-c-analysis/tad-detection
  - hi-c-analysis/loop-calling
  - hi-c-analysis/hic-visualization
  - hi-c-analysis/hic-differential
  - hi-c-analysis/hichip-plac-loops
qc_checkpoints:
  - after_pairs: "Long-range cis (>=20kb) fraction, not just %valid; trans is genome-size-dependent"
  - after_balance: "balance=True returns finite weights; masked bins are NaN by design"
  - after_analysis: "Eigenvector sign phased by GC; feature scale matches the resolution"
---

## Version Compatibility

Reference examples tested with: BWA-MEM2 2.2.1+, cooler 0.10+, cooltools 0.7+, bioframe 0.7+, matplotlib 3.8+, pairtools 1.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Hi-C Pipeline

**"Analyze my Hi-C data from FASTQ to 3D genome features"** -> Process read pairs and judge library quality, build and balance a cooler, then call ONLY the features the depth can support: compartments are cheap, TADs need moderate depth, de-novo loops need billions of contacts.

Complete workflow for Hi-C chromosome conformation capture analysis.

## The Decision That Frames the Whole Pipeline -- depth dictates the feature

Contacts scale with the SQUARE of the bin count, so the affordable resolution is set by depth, not ambition. Read `hi-c-analysis/matrix-operations` for the budget; the rule of thumb is ~1000 contacts/bin. Compartments (100kb-1Mb) come from almost any library; TAD boundaries (10-40kb) need a moderate map; de-novo loop calling (5-10kb dots) needed ~5 billion contacts in Rao 2014. On a shallow map, do NOT de-novo call loops -- run aggregate peak analysis (APA) on known CTCF/cohesin anchors instead (`hi-c-analysis/loop-calling`).

Protein-directed assays branch here: HiChIP, PLAC-seq, and Capture Hi-C have non-uniform peak-anchored coverage, so generic dots/HiCCUPS use the wrong null. Route them to `hi-c-analysis/hichip-plac-loops` (FitHiChIP/MAPS/CHiCAGO), NOT to step 6 below.

## Workflow Overview

```
Hi-C FASTQ files
    |
    v
[1. Alignment & Pairs] --> bwa-mem2 -SP5M + pairtools (parse/sort/dedup/split)
    |                       QC: long-range cis fraction is the one-number readout
    v
[2. Matrix Generation] --> cooler cload + zoomify (sum RAW, re-balance per resolution)
    |
    v
[3. Balancing] --------> ICE (cooler balance); REQUIRED before any analysis
    |
    v
[4. Compartments 100kb] -> eigs_cis, sign-phased by GC (E1 is a choice, not an output)
    |
    v
[5. TADs 10kb] ---------> insulation score across a window sweep (boundaries, not domains)
    |
    v
[6. Loops 10kb] --------> cooltools dots IF deep; else APA on known anchors
    |
    v
Hi-C features (compartments / boundaries / loops)
```

## Step 1: Alignment and Pair Processing

**Goal:** Turn raw Hi-C FASTQ into a deduplicated, classified `.pairs` list and judge whether the library worked.

**Approach:** Align the two mates independently with `bwa-mem2 -SP5M` (proper pairing would destroy long-range contacts), then parse, sort, deduplicate, and split with pairtools, reading the long-range cis fraction as the go/no-go.

```bash
# Align Hi-C reads (each end separately, then combine)
bwa-mem2 mem -SP5M -t 16 reference.fa reads_R1.fastq.gz | \
    pairtools parse --min-mapq 40 --walks-policy 5unique \
    --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 \
    --chroms-path reference.genome | \
    pairtools sort --nproc 16 --tmpdir ./tmp | \
    pairtools dedup --nproc-in 8 --nproc-out 8 \
    --mark-dups --output-stats stats.txt | \
    pairtools split --nproc-in 8 --output-pairs sample.pairs.gz

# Alternative: align both ends
bwa-mem2 mem -SP5M -t 16 reference.fa \
    reads_R1.fastq.gz reads_R2.fastq.gz | \
    pairtools parse --min-mapq 40 --walks-policy 5unique \
    --max-inter-align-gap 30 --chroms-path reference.genome | \
    pairtools sort | \
    pairtools dedup --mark-dups --output-stats stats.txt | \
    pairtools split --output-pairs sample.pairs.gz
```

**QC Checkpoint:** read `pairtools stats` as the go/no-go. The one-number readout is the LONG-RANGE cis fraction (>=20kb), not bare %valid: short-range cis is inflated by dangling ends and self-circles. Trans fraction is a noise floor but its acceptable value is genome-size-dependent (a human <10% threshold is meaningless for a microbe). High duplicate rate = low library complexity (not rescuable by sequencing deeper). See `hi-c-analysis/contact-pairs` for the orientation-balance QC and the Micro-C/Arima variants.

## Step 2: Generate Contact Matrix

```bash
# Create cooler file at multiple resolutions
cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    reference.genome:1000 \
    sample.pairs.gz \
    sample.1000.cool

# Multi-resolution (mcool)
cooler zoomify sample.1000.cool \
    -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000 \
    -o sample.mcool
```

## Step 3: Normalization (ICE Balancing)

```python
import cooler
import cooltools

# Load matrix
clr = cooler.Cooler('sample.mcool::/resolutions/10000')

# Balance (ICE normalization)
cooler.balance_cooler(clr, store=True, mad_max=5)

# Or via command line
# cooler balance sample.mcool::/resolutions/10000
```

## Step 4: Compartment Analysis

**Goal:** Assign each genomic bin to the active (A) or inactive (B) compartment with a non-arbitrary sign.

**Approach:** At a coarse 100kb resolution, compute the cis eigenvector and orient it with a GC phasing track so positive E1 is the active compartment (the sign is arbitrary without it).

```python
import cooler
import cooltools
import bioframe
import numpy as np

# Compartments are coarse-scale: 100kb, balanced matrix
clr = cooler.Cooler('sample.mcool::/resolutions/100000')

# Phasing track is NOT optional: the eigenvector sign is arbitrary. A GC track
# (matching the cooler binning exactly) orients positive E1 to the active (A) compartment.
view_df = bioframe.make_viewframe(clr.chromsizes)
gc = bioframe.frac_gc(clr.bins()[:][['chrom', 'start', 'end']], bioframe.load_fasta('reference.fa'))

eig_values, eig_vectors = cooltools.eigs_cis(clr, gc, view_df=view_df, n_eigs=3)
compartments = eig_vectors[['chrom', 'start', 'end', 'E1']].copy()
compartments['compartment'] = np.where(compartments['E1'] > 0, 'A', 'B')
compartments.to_csv('compartments.tsv', sep='\t', index=False)
```

## Step 5: TAD Detection

**Goal:** Locate domain boundaries at the sub-Mb scale.

**Approach:** Compute the insulation score across a window sweep at 10kb and read the `is_boundary`/`boundary_strength` columns the function returns directly (report boundaries, not a fixed domain partition).

```python
import cooltools

# Load matrix at TAD resolution
clr = cooler.Cooler('sample.mcool::/resolutions/10000')

# Insulation across a window sweep; the function already returns boundary columns
# (is_boundary_<W>, boundary_strength_<W>) -- there is no separate find_boundaries call.
ins = cooltools.insulation(clr, window_bp=[100000, 200000, 500000])

# Boundaries at the 200kb window; keep the continuous strength (comparable across samples)
boundaries = ins[ins['is_boundary_200000']][['chrom', 'start', 'end', 'boundary_strength_200000']]
boundaries.to_csv('tad_boundaries.tsv', sep='\t', index=False)

# Alternative: use HiCExplorer
# hicFindTADs -m sample.cool --outPrefix tads --correctForMultipleTesting fdr
```

## Step 6: Loop Calling

**Goal:** Detect focal CTCF/enhancer-promoter contacts, but only when the map is deep enough to support de-novo calling.

**Approach:** Compute a distance-matched expected, then run `cooltools dots` on a deep map; on a shallow library, skip de-novo calling and run APA on known anchors instead.

```python
import cooltools

# Load high-resolution matrix
clr = cooler.Cooler('sample.mcool::/resolutions/10000')

# De-novo dot calling is only honest on a DEEP map (Rao 2014 used ~5B contacts).
# On a shallow library, skip this and run APA on known anchors (see loop-calling).
expected = cooltools.expected_cis(clr)
loops = cooltools.dots(clr, expected, max_loci_separation=2000000, nproc=4)
loops.to_csv('loops.tsv', sep='\t', index=False)

# Alternative caller (template matching): chromosight
# chromosight detect --pattern loops sample.mcool::/resolutions/10000 loops
# For HiChIP/PLAC-seq/Capture Hi-C do NOT use dots -> hi-c-analysis/hichip-plac-loops
```

## Step 7: Visualization

```python
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cooltools.lib.plotting   # registers 'fall'; needs matplotlib < 3.9 with cooltools 0.7.x (else use 'afmhot_r')

# A square balanced map on a log scale; importing cooltools.lib.plotting registers 'fall'.
# Show O/E with a symmetric diverging cmap to see compartments/loops (see hic-visualization).
mat = clr.matrix(balance=True).fetch('chr1:50000000-60000000')
fig, ax = plt.subplots(figsize=(8, 8))
ax.matshow(mat, norm=LogNorm(vmax=mat[mat > 0].max() * 0.5), cmap='fall')
plt.savefig('hic_matrix.pdf')

# Triangle/track-stacked browser views: use pyGenomeTracks or HiCExplorer hicPlotTADs
# (data-visualization/genome-tracks), not a hand-rolled rotation.
```

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

THREADS=16
REF="reference.fa"
GENOME="reference.genome"
R1="sample_R1.fastq.gz"
R2="sample_R2.fastq.gz"
OUTDIR="hic_results"

mkdir -p ${OUTDIR}/{pairs,cool,analysis}

# Step 1: Alignment and pairs
echo "=== Alignment ==="
bwa-mem2 mem -SP5M -t ${THREADS} ${REF} ${R1} ${R2} | \
    pairtools parse --min-mapq 40 --walks-policy 5unique \
    --chroms-path ${GENOME} | \
    pairtools sort --nproc ${THREADS} --tmpdir ./tmp | \
    pairtools dedup --mark-dups --output-stats ${OUTDIR}/pairs/stats.txt | \
    pairtools split --output-pairs ${OUTDIR}/pairs/sample.pairs.gz

# Step 2: Generate matrix
echo "=== Matrix Generation ==="
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
    ${GENOME}:1000 ${OUTDIR}/pairs/sample.pairs.gz ${OUTDIR}/cool/sample.1000.cool

cooler zoomify ${OUTDIR}/cool/sample.1000.cool \
    -r 1000,5000,10000,25000,50000,100000,500000 \
    -o ${OUTDIR}/cool/sample.mcool

# Step 3: Balance
echo "=== Balancing ==="
for res in 10000 25000 100000; do
    cooler balance ${OUTDIR}/cool/sample.mcool::/resolutions/${res}
done

echo "=== Pipeline Complete ==="
echo "Run Python script for compartments, TADs, and loops"
```

## Python Analysis Script

```python
import cooler
import cooltools
import bioframe
import os

outdir = 'hic_results/analysis'
os.makedirs(outdir, exist_ok=True)

# Compartments (100kb) -- pass a GC phasing track (Step 4) so the sign is meaningful;
# eigs_cis without phasing returns an arbitrary-sign eigenvector.
print('Compartments...')
clr = cooler.Cooler('hic_results/cool/sample.mcool::/resolutions/100000')
gc = bioframe.frac_gc(clr.bins()[:][['chrom', 'start', 'end']], bioframe.load_fasta('reference.fa'))
eig_values, eig_vectors = cooltools.eigs_cis(clr, gc, n_eigs=3)
eig_vectors.to_csv(f'{outdir}/compartments.tsv', sep='\t', index=False)

# TADs (10kb)
print('TADs...')
clr = cooler.Cooler('hic_results/cool/sample.mcool::/resolutions/10000')
insulation = cooltools.insulation(clr, window_bp=[100000, 200000])
insulation.to_csv(f'{outdir}/insulation.tsv', sep='\t')

# Loops (10kb)
print('Loops...')
expected = cooltools.expected_cis(clr)
loops = cooltools.dots(clr, expected, nproc=4)
loops.to_csv(f'{outdir}/loops.tsv', sep='\t')

print(f'Results saved to {outdir}/')
```

## Related Skills

- hi-c-analysis/contact-pairs - Read-pair processing and the library-QC decision
- hi-c-analysis/hic-data-io - Cooler file operations and format conversion
- hi-c-analysis/matrix-operations - ICE balancing, expected/P(s), and the resolution-vs-depth budget
- hi-c-analysis/compartment-analysis - Sign-phased A/B compartments and saddle strength
- hi-c-analysis/tad-detection - Insulation-score boundaries across a window sweep
- hi-c-analysis/loop-calling - Dot calling and APA validation
- hi-c-analysis/hic-visualization - Normalization-aware contact-map plotting
- hi-c-analysis/hic-differential - Scale-matched comparison between conditions
- hi-c-analysis/hichip-plac-loops - Protein-directed loops (HiChIP/PLAC-seq/Capture Hi-C)
