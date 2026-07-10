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
# Pass BOTH mates to ONE bwa-mem2 call. -SP5M: -S/-P make bwa treat the mates as single-end and
# skip proper-pair rescue (so long-range contacts survive), while both sides are still emitted for
# pairtools to form the pair; -5 reports the 5'-most alignment of split reads, -M flags secondaries.
bwa-mem2 mem -SP5M -t 16 reference.fa reads_R1.fastq.gz reads_R2.fastq.gz | \
    pairtools parse --min-mapq 40 --walks-policy 5unique \
    --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 \
    --chroms-path reference.genome | \
    pairtools sort --nproc 16 --tmpdir ./tmp | \
    pairtools dedup --nproc-in 8 --nproc-out 8 \
    --mark-dups --output-stats stats.txt | \
    pairtools split --nproc-in 8 --output-pairs sample.pairs.gz
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

**Goal:** ICE-balance the matrix so every bin has equal marginal coverage, without letting empty/artifact bins corrupt the result.

**Approach:** Mask low-coverage and blacklist bins BEFORE balancing, then balance per resolution. ICE assumes equal visibility per bin, so an unmasked empty or repeat/blacklist bin is iteratively up-weighted into a bright stripe artifact; `mad_max` filters bins whose coverage is `mad_max` MADs below the median, and a blacklist/`--blacklist` (or pre-masking bad bins) removes known artifacts. Balancing is REQUIRED before any analysis, but it does NOT make two maps comparable across conditions — that needs depth-matching + distance-stratified normalization (hi-c-analysis/hic-differential).

```python
import cooler
import cooltools

# Mask before balancing: mad_max drops low-coverage bins that would otherwise become stripes.
# Balance EVERY resolution the downstream steps analyze -- weights are resolution-specific, and an
# unbalanced cooler has no 'weight' column, so cooltools (eigs_cis, insulation, dots) fails on it.
# Steps 4-6 below use 100kb (compartments) and 10kb (loops and insulation).
for res in (10000, 25000, 100000):
    clr = cooler.Cooler(f'sample.mcool::/resolutions/{res}')
    cooler.balance_cooler(clr, store=True, mad_max=5, ignore_diags=2, min_nnz=10)   # masked bins are NaN by design

# CLI equivalent (add --blacklist regions.bed to remove known-artifact bins first):
# for res in 10000 25000 100000; do cooler balance --mad-max 5 --ignore-diags 2 --min-nnz 10 sample.mcool::/resolutions/${res}; done
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
# Masked bins have E1 = NaN; NaN > 0 is False, so guard or a bare np.where mislabels them all 'B'.
compartments['compartment'] = np.where(compartments['E1'].isna(), None, np.where(compartments['E1'] > 0, 'A', 'B'))
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

# Boundaries at the 200kb window; keep the continuous strength (comparable across samples).
# is_boundary is NaN for bad/low-mappability bins; fillna(False) before masking or pandas raises.
boundaries = ins[ins['is_boundary_200000'].fillna(False).astype(bool)][['chrom', 'start', 'end', 'boundary_strength_200000']]
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
import cooltools.lib.plotting   # registers the 'fall' cmap; if unavailable in your stack, use 'afmhot_r'

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

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Long-range contacts missing / map looks like short-range only | Mates aligned as a proper pair instead of independently | Align with `bwa-mem2 mem -SP5M` (each end separately, no proper-pair rescue) |
| Library "passed" %valid but is unusable | Judged on bare %valid; short-range cis is inflated by dangling ends/self-circles | Read the long-range cis (>=20kb) fraction as the go/no-go |
| Bright stripes/plaid artifacts after balancing | Empty/repeat/blacklist bins not masked before ICE | Mask with `mad_max`/`--blacklist`/`min_nnz` before balancing |
| A/B compartments flipped between samples | Eigenvector sign is arbitrary without phasing | Phase E1 by a GC track matching the cooler binning exactly |
| De-novo loops look sparse/noisy | Called `dots` on a shallow map | Only de-novo call on deep maps (~billions of contacts, Rao 2014); else APA on known anchors |
| Cross-condition differences dominated by depth | Compared balanced maps directly | Depth-match + distance-stratified normalize first (hi-c-analysis/hic-differential) |
| HiChIP/PLAC "loops" full of false positives | Generic dots/HiCCUPS null on peak-anchored coverage | Route to FitHiChIP/MAPS/CHiCAGO (hi-c-analysis/hichip-plac-loops) |

## References

- Rao SSP, Huntley MH, Durand NC, et al (2014) A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. *Cell* 159:1665-1680. DOI 10.1016/j.cell.2014.11.021. (depth-vs-resolution; ~5B contacts for kilobase loops.)
- Imakaev M, Fudenberg G, McCord RP, et al (2012) Iterative correction of Hi-C data reveals hallmarks of chromosome organization. *Nature Methods* 9:999-1003. DOI 10.1038/nmeth.2148. (ICE balancing.)
- Abdennur N, Mirny LA (2020) Cooler: scalable storage for Hi-C data and other genomically labeled arrays. *Bioinformatics* 36:311-316. DOI 10.1093/bioinformatics/btz540.
- Open2C, Abdennur N, Fudenberg G, et al (2024) Cooltools: enabling high-resolution Hi-C analysis in Python. *PLOS Computational Biology* 20:e1012067. DOI 10.1371/journal.pcbi.1012067.
- Open2C, Abdennur N, Fudenberg G, et al (2024) Pairtools: from sequencing data to chromosome contacts. *PLOS Computational Biology* 20:e1012164. DOI 10.1371/journal.pcbi.1012164.

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
