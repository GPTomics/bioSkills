---
name: bio-atac-seq-nucleosome-positioning
description: Map nucleosome center positions, occupancy, and fuzziness from ATAC-seq fragment-size patterns using NucleoATAC, ATACseqQC, DANPOS3, or scprinter. Use when characterizing nucleosome organization at promoters and enhancers, calling +1/-1 nucleosomes flanking NFRs, generating V-plots for chromatin structure visualization, or comparing nucleosome positioning between conditions.
tool_type: mixed
primary_tool: NucleoATAC
---

## Version Compatibility

Reference examples tested with: NucleoATAC 0.3.4+, ATACseqQC 1.26+, DANPOS 3.1+, samtools 1.19+, pysam 0.22+, pyBigWig 0.3+, BSgenome.Hsapiens.UCSC.hg38 1.4+, TxDb.Hsapiens.UCSC.hg38.knownGene 3.18+.

NucleoATAC is unmaintained since 2018 but remains the canonical ATAC-specific nucleosome caller; ATACseqQC, DANPOS3, and scprinter are actively developed alternatives. Verify versions before use:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Nucleosome Positioning

**"Where are the nucleosomes in my ATAC-seq data?"** -> Use fragment-size classes (Tn5 cuts twice through naked DNA generating short fragments; once on each side of a single nucleosome generating ~147+linker fragments) to call nucleosome centers, occupancy scores, and the spacing pattern around regulatory elements.

- CLI: `nucleoatac run --bed regions.bed --bam sample.bam --fasta genome.fa`
- R: `ATACseqQC::splitGAlignmentsByCut()` -> fragment classes; `factorFootprints()` -> per-TF flanking nuc analysis
- CLI: `danpos3 dpos sample.bam` (alternative; supports MNase, ATAC, DNase)
- Python: `scprinter` for multi-scale nucleosome inference

## Nucleosome Physics for ATAC

A nucleosome wraps ~147 bp DNA in 1.65 turns. Adjacent nucleosomes are separated by 20-50 bp linker; mean **nucleosome repeat length (NRL)** is species-dependent:

| Cell type / organism | NRL | Notes |
|----------------------|-----|-------|
| Yeast S. cerevisiae | 165 bp | Tightly packed; less linker |
| Drosophila S2 | 175-185 bp | |
| Mouse ES cells | 188-196 bp | |
| Human HEK293 / K562 | 196-200 bp | Standard somatic |
| Human cortical neurons | 211 bp | Longer linker |
| Sperm chromatin | 240-250 bp | Tight packaging via protamines |
| Active gene bodies | -10 bp shorter than genome avg | Active transcription disrupts |

NRL determines fragment-size peak positions. ATAC mono-nucleosome peak is at NRL (NOT 147 bp -- that's the protected length; ATAC fragments span the full nucleosome+linker). Di-nuc is at 2x NRL minus a small overlap.

## Fragment-Size Classes (Buenrostro 2013, refined)

| Class | Fragment range | Origin | Use |
|-------|---------------|--------|-----|
| Sub-nucleosomal / NFR | < 100 bp | Two Tn5 cuts in naked accessible DNA | TF binding, footprinting |
| Mono-nucleosomal | 180-247 bp | Tn5 cuts on each side of one nucleosome | Nucleosome positioning |
| Di-nucleosomal | 315-473 bp | Tn5 cuts span two nucleosomes | Phasing, NRL estimation |
| Tri-nucleosomal | 558-615 bp | Three nucleosomes | Heterochromatin / phasing |
| > 700 bp | Rare | Often artefact (chimeric); discard | -- |

Mono-nucleosome window 180-247 bp is the Buenrostro 2013 convention; ATACseqQC uses 180-250. Adjust for the target organism's NRL.

## V-Plot Interpretation

V-plots (fragment-size vs position) are diagnostic. X-axis is position relative to a feature (TSS, motif center); Y-axis is fragment size. Aggregate density forms characteristic patterns:

| Pattern | Visual | Meaning |
|---------|--------|---------|
| V (apex at center, low size at center, increasing flanks) | Classic V | TF or NFR at center, flanking nucleosomes |
| W (two V's flanking center) | W-shape | NFR at center plus +1 / -1 nucleosomes |
| Inverted V (peak at center) | Mountain | Fragment fully enclosed at feature; e.g. nucleosome-bound TF | 
| Flat band at 200 bp | Horizontal line | Constitutive nucleosome (no positioning relative to feature) |
| 10.4 bp helical phasing on V | Sub-peaks at 50, 60, 70, 80 bp size | Tn5 helical preference visible; high-quality library |

V-plots are the primary diagnostic for whether nucleosome-positioning analysis will succeed. Flat-band patterns mean no positioning information; classic V/W patterns mean positioning is recoverable.

## Algorithmic Taxonomy

| Tool | Method | Resolution | Strength | Fails when |
|------|--------|-----------|----------|------------|
| NucleoATAC | Cross-correlation with idealized V-plot template; per-base occupancy + nucleosome calls | Single-bp | ATAC-specific; provides occupancy + fuzziness | Unmaintained since 2018; pegs Python 2/3.6; struggles on chromatin without clear NRL |
| ATACseqQC | Fragment-size split + Tn5-shifted GAlignments + V-plot from BAM | Region-level | R/Bioconductor; integrates with TxDb / motif analysis | No per-base nucleosome calls; visualization-focused |
| DANPOS3 | Smoothing + peak call on cleavage signal; tested on MNase, ATAC, DNase | ~50 bp | Robust differential mode (`dpeak`); MNase legacy; broadly maintained | Designed for MNase-Seq; ATAC adaptation needs careful parameter tuning |
| scprinter | CNN multi-scale; resolves co-occurring TF + nucleosome footprints | Single-bp | Modern; single-cell aware; multi-scale | Newer; benchmarks evolving; GPU recommended |
| custom (pysam V-plot) | Fragment counting + 2D density | Region-level | Maximally flexible; reproducible | Requires manual calling logic; slow |

Methodology evolves; verify against current Schep 2015 (NucleoATAC), Chen 2013 (DANPOS), Bao 2024 (scprinter) before locking pipelines.

## +1 Nucleosome Calling

The +1 nucleosome (first nucleosome downstream of TSS, immediately bordering the NFR) is the most-studied positioning feature. Its position relative to TSS determines transcription initiation kinetics.

**Canonical +1 position:** +50 to +60 bp from TSS in metazoa; -100 to -120 bp from TATA in yeast; varies by gene type (Pol II vs Pol III, housekeeping vs developmental).

**Calling strategy:**

**Goal:** Identify each gene's +1 nucleosome, the first nucleosome downstream of the TSS that flanks the NFR.

**Approach:** Build gene-body intervals slopped around TSSs, run NucleoATAC over them to call per-base nucleosome positions, then pick the most-downstream-of-TSS nucleosome per gene.

```bash
# 1. Define gene-body intervals
bedtools slop -i genes.bed -g chrom.sizes -l 200 -r 1000 > gene_bodies.bed

# 2. Run NucleoATAC
nucleoatac run --bed gene_bodies.bed --bam sample.dedup.bam --fasta genome.fa \
    --out tss_nuc/ --cores 8

# 3. The first nucleosome downstream of each TSS in nucpos.bed is +1
```

A failure to detect a clear +1 peak in aggregate V-plot suggests TSS annotation is wrong or library is over-transposed.

## Per-Tool Failure Modes

### NucleoATAC -- Region size and depth dependence

**Trigger:** Short region BED (< 1 kb per region); shallow library (< 25M nuclear reads).

**Mechanism:** NucleoATAC fits an idealized V-plot template per region. Short regions provide too few fragments for stable correlation; shallow data provides noisy templates.

**Symptom:** No nucleosome calls in shallow regions; "occupancy" track is flat at zero.

**Fix:** Use regions >= 500 bp; merge adjacent peaks via bedtools to ensure region size; require >= 30M nuclear reads.

### NucleoATAC -- Maintenance status

**Trigger:** Installing NucleoATAC in 2025+.

**Mechanism:** Last release 2018; pegs Python 3.6 in some installs; depends on outdated NumPy API.

**Fix:** Use a dedicated conda env (`conda create -n nucleoatac python=3.7 numpy=1.18 scipy=1.5 pysam`); accept it works but is no longer updated. Consider scprinter or DANPOS3 alternatives for new projects.

### ATACseqQC factorFootprints -- Asymmetric nucleosome flanks

**Trigger:** Pioneer-factor binding sites where one face is on a nucleosome.

**Mechanism:** factorFootprints assumes symmetric flanking nucleosomes. Pioneer TFs (FOXA1, GATA) only have nucleosome on one side -> asymmetric output.

**Symptom:** Single shoulder in flanking signal; unbalanced V-plot.

**Fix:** Treat asymmetry as biological signal, not artefact. For pioneers, use stranded analysis.

### DANPOS dpos with default parameters -- ATAC mismatch

**Trigger:** Running `danpos3 dpos` with MNase defaults on ATAC.

**Mechanism:** DANPOS3's smoothing window and peak-calling defaults are tuned for MNase signal (smoother coverage). ATAC's sharper signal requires `--smooth_width 80 --width 145` or similar; otherwise calls are over-smoothed.

**Fix:** Use ATAC-tuned parameters. See DANPOS docs for ATAC-specific recipe; or use NucleoATAC instead.

### Mono-nucleosome filter window mis-set

**Trigger:** Using strict 147 bp filter for mono-nuc fraction; using 100-180 bp instead of 180-247.

**Mechanism:** Mono-nuc fragments are 180-247 bp because they span the nucleosome AND a linker. Filtering tighter excludes the legitimate signal.

**Symptom:** Mono-nuc count is much lower than expected (< 30% of NFR count).

**Fix:** Use Buenrostro 2013 windows: NFR < 100, mono 180-247, di 315-473.

## Decision Tree by Goal

| Goal | Recommended workflow |
|------|---------------------|
| Per-base nucleosome occupancy track | NucleoATAC (with caveat about maintenance); or scprinter |
| V-plot at TSS or motif center | ATACseqQC vPlot |
| Differential nucleosome positioning between conditions | DANPOS3 dpeak |
| +1 nucleosome calling at all genes | NucleoATAC + post-process to first nuc downstream of TSS |
| Single-cell nucleosome positioning | scprinter |
| Quick fragment-size QC plot | ATACseqQC fragSizeDist |
| NRL estimation | Custom Fourier / autocorrelation on fragment-end coverage |
| Nucleosome-aware peak calling | MACS3 hmmratac (peak-calling skill) |

## Estimating NRL from Fragment-Size Distribution

**Goal:** Estimate the nucleosome repeat length from ATAC fragment-size periodicity.

**Approach:** Collect proper-pair fragment lengths from the BAM, build a histogram, find density peaks via scipy find_peaks, and read off the mono-nucleosome peak position within the 150-250 bp window.

```python
import numpy as np, pysam
from scipy.signal import find_peaks

bam = pysam.AlignmentFile('sample.bam', 'rb')
frag_lengths = [abs(r.template_length) for r in bam.fetch()
                if r.is_proper_pair and r.is_read1 and 0 < abs(r.template_length) < 1500]
hist, edges = np.histogram(frag_lengths, bins=300, range=(0, 1500))
peaks, _ = find_peaks(hist, distance=50, prominence=hist.max() * 0.05)
peak_positions = edges[peaks] + (edges[1] - edges[0]) / 2
# Mono peak should be ~NRL; di peak ~2*NRL
mono = peak_positions[(peak_positions > 150) & (peak_positions < 250)][0]
print(f'Estimated NRL: {mono:.0f} bp')
```

NRL inferred this way is approximate; for precision use autocorrelation on cumulative cleavage coverage instead.

## V-Plot in Python

**Goal:** Build a fragment-size-by-position density plot to diagnose nucleosome positioning around a feature.

**Approach:** Iterate proper-pair fragments in a flank window around each feature center, accumulate counts into a (fragment_size x position) grid, and render the 2D density.

```python
import numpy as np, pysam, matplotlib.pyplot as plt

def vplot(bam_path, regions_bed, max_size=600, flank=1000):
    bam = pysam.AlignmentFile(bam_path, 'rb')
    grid = np.zeros((max_size, 2 * flank))
    for line in open(regions_bed):
        chrom, start, *_ = line.strip().split('\t')
        center = int(start)
        for r in bam.fetch(chrom, max(0, center - flank), center + flank):
            if not r.is_proper_pair or not r.is_read1: continue
            size = abs(r.template_length)
            if size <= 0 or size >= max_size: continue
            frag_center = r.reference_start + size // 2
            x = frag_center - center + flank
            if 0 <= x < 2 * flank:
                grid[size, x] += 1
    return grid

g = vplot('sample.bam', 'tss.bed')
plt.imshow(g, aspect='auto', origin='lower', cmap='magma',
           extent=[-1000, 1000, 0, 600])
plt.xlabel('Distance from feature (bp)')
plt.ylabel('Fragment size (bp)')
plt.savefig('vplot.png', dpi=200, bbox_inches='tight')
```

V-plot quality is the most useful diagnostic before nucleosome calling. Classic V at TSS = positioning info recoverable; flat band = not.

## Differential Nucleosome Positioning (DANPOS3 dpeak)

```bash
# Compare control vs treatment nucleosome positions
danpos3 dpeak \
    -b condition2.bam:condition1.bam \
    -c control.bam \
    -o danpos_diff/ \
    --paired 1 \
    --width 145 --smooth_width 80
```

DANPOS reports four event types: shifted nucleosomes, gained, lost, fuzziness change. ENCODE has no official threshold; require >= 30 bp shift and FDR < 0.05 for nucleosome shift calls.

**Full ATAC-tuned DANPOS3 recipe:**
```bash
# --width 145: nucleosome footprint width
# --smooth_width 80: ATAC-specific (MNase default 60 too narrow)
# -jd 145: min distance between adjacent nuc calls (single-dash short flag)
# --pheight 0.95: peak height fraction for calling
# --frsz 200: fragment size used (mono-nuc)
danpos3 dpos sample.bam \
    --paired 1 \
    --width 145 \
    --smooth_width 80 \
    -jd 145 \
    --pheight 0.95 \
    --frsz 200 \
    --out danpos_out/
```

Verify exact flags with `danpos3 dpos --help`; DANPOS3 documentation has been spotty and flag names can drift across releases.

Adapted from DANPOS3 docs for ATAC; `--smooth_width 80` widens the smoothing kernel to match ATAC's sharper signal vs MNase's broader cleavage. `-jd 145` (single-dash short, alternative `--distance 145`) enforces nucleosome spacing >= 145 bp (one nucleosome footprint).

## Histone Variant Detection from Fragment Size

**Trigger:** Suspected H2A.Z- or H3.3-containing nucleosomes; differential nucleosome composition between conditions.

**Mechanism:** H2A.Z replacement of H2A produces nucleosomes with weaker DNA-histone interaction (lower thermal stability); fragment-size shifts ~10-15 bp shorter than canonical H2A nucleosomes (Voong 2017 supplementary). H3.3 replacement is more subtle but H3.3-H2A.Z double-variant nucleosomes are particularly destabilized at active promoters.

**Detection:** Aggregate fragment-size distribution at H2A.Z ChIP-seq peaks vs H3K4me3-only peaks; the H2A.Z population shows mean fragment size ~10 bp shorter. ATAC alone CANNOT definitively call H2A.Z; H2A.Z ChIP-seq is needed for ground truth. ATAC fragment-size analysis is a hypothesis generator.

```python
# Per-region fragment-size mean as H2A.Z indicator
def region_frag_size(bam, region):
    sizes = [abs(r.template_length) for r in bam.fetch(*region)
             if r.is_proper_pair and r.is_read1 and 100 < abs(r.template_length) < 300]
    return np.mean(sizes) if sizes else np.nan

# Compare H2A.Z-positive vs H2A.Z-negative TSSs
```

## Long-Read Single-Molecule Chromatin (Fiber-seq, NanoNOMe)

**Alternative to short-read ATAC for nucleosome positioning:**

| Method | Tech | Resolution | Strength |
|--------|------|-----------|----------|
| Fiber-seq (Stergachis 2020) | PacBio HiFi + DNA methylation footprinting | Per-molecule single-bp | Reads continuous chromatin fiber up to 20 kb; resolves haplotype-specific positioning |
| NanoNOMe (Lee 2020) | Nanopore + GpC methyltransferase | Per-molecule single-bp | Same single-molecule but cheaper than PacBio |
| MOSE (Dong 2024) | Nanopore + methylation | Improved Nanopore | Newer; benchmarks emerging |

Fiber-seq can detect nucleosome occupancy directly per single chromatin molecule (no aggregation needed). Resolves cell-cycle-dependent and stochastic positioning that bulk ATAC averages out. Preferred for fine-structure analysis of regulatory elements.

For most labs, short-read ATAC + NucleoATAC remains primary; Fiber-seq is special-purpose when single-molecule resolution is essential.

## Nucleosome Fuzziness

Fuzziness measures how sharply positioned a nucleosome is across cells. Defined as the standard deviation of per-cell nucleosome center positions.

| Fuzziness range | Interpretation |
|-----------------|----------------|
| < 20 bp | Sharply positioned (rare in metazoa; common at +1 in yeast) |
| 20-50 bp | Standard well-positioned |
| 50-100 bp | Fuzzy; constitutive but non-stable |
| > 100 bp | Effectively unpositioned |

NucleoATAC reports per-nucleosome fuzziness as `nuc_size`; values around the NRL are expected.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `nucleoatac run` ImportError on numpy | Python 3.6 incompatibility | Use dedicated conda env with pinned versions |
| Empty .nucpos.bed output | Region BED too short or library too shallow | Verify region size >= 500 bp; depth >= 30M |
| V-plot shows horizontal band, no V | No positioning info; library over-transposed or wrong feature center | Check feature BED; verify TSS positions are correct |
| Mono-nuc count very low | Wrong fragment-size window (used 100-180 instead of 180-247) | Use Buenrostro windows |
| factorFootprints asymmetric | Pioneer TF; this is biological | Treat as signal, not artefact |
| DANPOS calls many shifts | MNase parameters used on ATAC | Tune `--smooth_width 80 --width 145` for ATAC |
| splitGAlignmentsByCut error in ATACseqQC | BAM is single-end | Mono-nuc analysis requires paired-end |
| +1 nucleosome not visible at TSS aggregate | TSS list mixes coding + non-coding strands; or wrong genome build | Restrict to protein-coding TSSs in matched build |

## References

- Schep AN et al 2015 Genome Res 25:1757 (NucleoATAC)
- Chen K et al 2013 Genome Res 23:341 (DANPOS)
- Buenrostro JD et al 2013 Nat Methods 10:1213 (ATAC fragment-size classes)
- Ou J et al 2018 BMC Genomics 19:169 (ATACseqQC)
- Bao Y et al 2024 bioRxiv (scprinter)
- Mavrich TN et al 2008 Nature 453:358 (+1 nucleosome positioning)
- Voong LN et al 2017 Cell 169:1422 (high-resolution nucleosome mapping)
- Teif VB et al 2012 Nat Struct Mol Biol 19:1185 (NRL variation across cell types)

## Related Skills

- atac-seq/atac-qc - Fragment-size periodicity QC
- atac-seq/atac-peak-calling - Nucleosome-aware MACS3 hmmratac
- atac-seq/footprinting - Per-TF flanking nucleosome analysis
- atac-seq/single-cell-atac - scprinter for sc nucleosome positioning
- chip-seq/peak-annotation - Annotate nucleosome positions to genes
- alignment-files/bam-statistics - Insert-size statistics upstream
