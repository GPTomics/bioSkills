---
name: bio-genome-intervals-proximity-operations
description: Performs proximity operations on genomic intervals with bedtools (closest, window, flank, slop) and pybedtools - nearest-feature queries with signed/strand-aware distance, fixed-radius window searches, strand-aware promoter construction, and interval extension. Covers the closest -d/-D a/b/ref/-t/-k/-io/-iu/-id flags, the -D ref strand sign-flip, silent chromosome-end clipping in slop/flank, -t all tie double-counting, and the critical distinction between a geometry answer (nearest TSS) and a biology answer (which gene an element regulates). Use when assigning peaks or variants to genes, defining promoters from a gene model, building distance-to-TSS distributions, finding features within a window, or extending intervals - and when deciding whether nearest-gene is a fair prior (GWAS locus) or a trap (distal enhancer).
tool_type: mixed
primary_tool: bedtools
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, pybedtools 0.10+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `bedtools --version` then `bedtools <subcommand> --help` to confirm flags
- Python: `pip show pybedtools` then `help(pybedtools.BedTool.closest)` to check signatures

`flank` and `slop` REQUIRE a chrom-sizes (`genome.txt`, two columns: chrom<TAB>length) file via `-g`; `closest` requires both inputs coordinate-sorted (`sort -k1,1 -k2,2n`). If code throws an error, introspect the installed tool and adapt rather than retrying.

# Proximity Operations

**"Which gene is nearest to each peak, and is that the gene it regulates?"** -> Compute interval geometry (nearest feature, signed distance, window membership, strand-aware promoters) with bedtools, then decide honestly whether geometry answers the biological question.
- CLI: `bedtools closest -D b -t first -a peaks.bed -b genes.bed`, `bedtools window -w 50000`, `bedtools slop -s -l 2000 -r 200 -g genome.txt`
- Python: `peaks.closest(genes.sort(), D='b', t='first')`, `peaks.window(genes, w=50000)`, `tss.slop(g='genome.txt', s=True, l=2000, r=200)` (pybedtools)

## The Single Most Important Modern Insight -- closest Answers a GEOMETRY Question Misread as a BIOLOGY Question

`bedtools closest` answers "what is the nearest annotated TSS?" - a coordinate fact. The user almost always wants "which gene does this element regulate?" - a biology claim. For distal regulatory elements these disagree the **majority of the time**. In the CRISPRi-FlowFISH gold standard (Fulco 2019 *Nat Genet* 51:1664), assigning each tested distal element to the **closest expressed gene gave only ~47% precision and ~37% recall** - the nearest gene was the wrong target most of the time, and the method missed nearly two-thirds of real links. Enhancers routinely skip intervening genes: the canonical case is the obesity-associated *FTO* intron regulating **IRX3 ~500 kb away**, not *FTO* (Smemo 2014 *Nature* 507:371). Do the bedtools arithmetic flawlessly here, then route real enhancer->gene linking to activity/contact/QTL methods (ABC: Fulco 2019, Nasser 2021; PCHi-C; eQTL-coloc) at atac-seq/enhancer-gene-linking - never present "nearest gene" as a regulatory call for a distal element.

The deeper twist - **two regimes, opposite advice, identical command:**
- **Enhancer -> target (closest is a TRAP).** Distal ATAC/H3K27ac peaks, enhancer GWAS variants: nearest gene is wrong most of the time. Use as a candidate generator, validate with ABC/PCHi-C/eQTL.
- **GWAS locus -> gene (closest is a fair PRIOR).** For a fine-mapped, colocalized credible-set SNP, the nearest **protein-coding** gene is right ~50-65% of the time - a strong, hard-to-beat baseline for which gene a locus implicates. Route to causal-genomics for the rigorous version, but nearest-coding-gene is a defensible first pass.

Conflating the two regimes is the real error. The discriminator: is the question the *target of an enhancer* (distrust nearest) or *the gene under a GWAS peak* (nearest is a fine first pass)?

## Operation Taxonomy

| Operation | What it computes | Strand-aware? | Needs genome file? |
|-----------|------------------|---------------|--------------------|
| closest | For each A, the nearest B (+ optional signed distance) | optional (`-s`/`-S`, `-D a`/`-D b`) | no |
| window | For each A, all B within +-W bp (fuzzy intersect) | optional (`-sw`/`-sm`/`-Sm`) | no |
| slop | Grow each interval by N bp, keeping it one feature | optional (`-s`) | yes (`-g`) |
| flank | Emit the regions BESIDE each interval, dropping the body | optional (`-s`) | yes (`-g`) |

`closest`/`window` are queries (A vs B); `slop`/`flank` are transforms (A only, + genome file). The slop-vs-flank distinction trips people: `slop -b 1000` makes a peak 2 kb wider (one feature); `flank -b 1000` returns only the left/right neighboring 1 kb regions and discards the peak itself (two features). `window -w 0` is approximately `intersect`.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Nearest gene to a promoter-proximal mark (H3K4me3, Pol II, CAGE) | `closest -D b -io -t first` | the peak really is at the gene it marks; closest is honest here |
| Distal enhancer / ATAC peak -> which gene? | `closest`/`window` as candidates, then -> atac-seq/enhancer-gene-linking | nearest is wrong the majority of the time (ABC/PCHi-C/eQTL link it) |
| GWAS credible-set SNP -> implicated gene | `closest` to nearest **protein-coding** gene, then -> causal-genomics/colocalization-analysis | nearest-coding-gene is a ~50-65% prior; a fair first pass |
| All candidate genes near an element | `window -w 50000` (or TAD-scale) | honest "candidate set", not a single call |
| Build promoters from a gene model | collapse to TSS, then `slop -s -l UP -r DOWN -g` | a promoter is an imposed definition, strand-aware, from the TSS |
| Distance-to-TSS distribution | `closest -D b -d` then plot signed distance | a distribution beats a binary "promoter vs distal" threshold |
| Upstream-only / downstream-only nearest | `closest -D b -iu` / `-id` | direction must be strand-relative (`-D b`), never `-D ref` |
| Peak-set GO enrichment from proximity | -> GREAT/rGREAT (regulatory-domain model) | avoids the `-t all` double-counting and distal mis-assignment |
| Regions flanking a feature (splice/boundary context) | `flank -s -b N -g` | the regions outside the feature, strand-aware |
| Peaks not yet called | -> chip-seq/peak-calling, atac-seq/atac-peak-calling | this skill operates on existing intervals |

## closest - Nearest Feature with Signed, Strand-Aware Distance

Default: for each A, report the single nearest B; **on ties, report ALL tied B** (`-t all` is the default - the double-counting trap below). Both inputs must be sorted. When A's chromosome has no B feature, bedtools prints `none` for B columns and **`-1`** for distance - filter this sentinel before any numeric summary.

```bash
# Nearest gene, signed distance by the GENE's strand, ignore overlaps, one row per peak
bedtools sort -i peaks.bed > peaks.sorted.bed
bedtools sort -i genes.bed  > genes.sorted.bed
bedtools closest -a peaks.sorted.bed -b genes.sorted.bed -D b -io -t first > nearest.bed
#                                                         ^^^^ sign by gene strand (biology, not coordinates)
#                                                              ^^^ closest non-overlapping gene
#                                                                  ^^^^^^^^ resolve ties deterministically (document this)

# k=3 nearest with unsigned distance (k>1 intentionally multiplies rows)
bedtools closest -a peaks.sorted.bed -b genes.sorted.bed -k 3 -d > top3.bed
```

```python
import pybedtools

peaks = pybedtools.BedTool('peaks.bed').sort()
genes = pybedtools.BedTool('genes.bed').sort()
near = peaks.closest(genes, D='b', io=True, t='first')      # -D b -io -t first
near = near.filter(lambda x: int(x.fields[-1]) != -1)        # drop the no-feature sentinel
near.saveas('nearest.bed')
```

Key flags: `-d` unsigned distance (overlaps = 0); `-D ref` signed by coordinate only (strand-agnostic - see Failure Modes); `-D a`/`-D b` signed by A's / B's strand; `-t all|first|last`; `-k N` k-nearest; `-io` ignore overlapping B; `-iu`/`-id` ignore upstream/downstream (require `-D`); `-fu`/`-fd` first upstream/downstream; `-s`/`-S` same/opposite strand; `-N` require different names; `-mdb each|all` and `-names`/`-filenames` for multiple `-b` files.

## window - Features Within a Search Radius

`window` reports all B within a window around each A (default 1000 bp each side). Use it for the honest "candidate genes near this element" framing.

```bash
# All genes within 50 kb of each peak, counted per peak
bedtools window -a peaks.bed -b genes.bed -w 50000 -c > peak_gene_counts.bed
```

Flags: `-w N` symmetric (default 1000); `-l N`/`-r N` asymmetric (coordinate left/right); `-sw` define `-l`/`-r` BY STRAND; `-sm`/`-Sm` keep only same/opposite-strand B; `-u` boolean (A once if any B); `-c` count of B per A; `-v` A with no B in window. `-sw` controls *where the window is*; `-sm`/`-Sm` control *which B count* - distinct concerns.

## slop / flank - Extend or Find Adjacent Regions (genome file REQUIRED)

Both need `-g genome.txt` precisely so they can clip at chromosome boundaries - extension past coordinate 0 or past chrom length is **silently truncated** (start floored at 0, end capped). Flags: `-b N` both sides; `-l N`/`-r N` per side (coordinate unless `-s`); `-s` strand-aware (on a `-`-strand feature `-l` adds to the END, so `-l` always means "upstream of the feature"); `-pct` treat N as a fraction of feature length; `-header` echo input header.

### Build Strand-Aware Promoters from a Gene Model

**Goal:** Produce a promoter BED (TSS -2000 / +200 bp, strand-aware) that is correct for both strands - the right way to define "promoter", which is a choice imposed on a TSS, not an annotated feature.

**Approach:** Collapse genes to their TSS first (start for `+`, end-1 for `-`), THEN `slop -s` so "upstream" tracks strand. Running `slop -b 2000` on a gene BODY is the wrong promoter (it grows the whole gene, ignores strand).

```bash
# 1) TSS BED from a BED6 gene model (strand-aware single base)
awk -v OFS='\t' '{ if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6 }' genes.bed > tss.bed

# 2) Promoter = TSS -2000 / +200, strand-aware (-l is always the upstream side under -s)
bedtools slop -i tss.bed -g genome.txt -s -l 2000 -r 200 > promoters.bed
```

```python
import pybedtools

UP = 2000   # bp upstream of TSS; common core-promoter convention, NOT a fact -- report it and tune per assay
DOWN = 200  # bp downstream of TSS; asymmetric on purpose (+1 nucleosome / 5'UTR sit downstream)

genes = pybedtools.BedTool('genes.bed')
tss = genes.each(lambda f: pybedtools.create_interval_from_list([f[0], str(f.start) if f.strand == '+' else str(f.end - 1), str(f.start + 1) if f.strand == '+' else str(f.end), f.name, f.score, f.strand])).saveas()
promoters = tss.slop(g='genome.txt', s=True, l=UP, r=DOWN).saveas('promoters.bed')
```

`flank` shares the flag vocabulary but emits the regions BESIDE each feature and drops the original (two intervals per input, used for splice/boundary context):

```bash
bedtools flank -i exons.bed -g genome.txt -s -b 1000 > exon_flanks.bed   # 1 kb each side, strand-aware
```

## Per-Method Failure Modes

### -D ref silently mis-signs minus-strand genes
**Trigger:** using `closest -D ref` and interpreting the sign as upstream/downstream. **Mechanism:** `-D ref` signs by genomic coordinate only (lower = negative); for a `-`-strand gene the TSS is at the HIGHER coordinate, so "upstream" runs to higher coordinates and the coordinate sign is inverted relative to biology. **Symptom:** half the genes (the `-`-strand ones) are folded the wrong way; symmetric QC (TSS-enrichment plot) still looks fine, but any "enhancers preferentially upstream" claim washes out or inverts. **Fix:** use `-D b` (sign by the gene's strand) for any upstream/downstream biology; reserve `-D ref` for pure left/right genomic distance.

### slop on a gene body is not a promoter
**Trigger:** `slop -b 2000` (or `-l 2000 -r 0` without `-s`) on a gene-body BED, called "the promoter". **Mechanism:** it grows the window around the whole gene, not the TSS, and without `-s` adds the "upstream" side to the wrong (3') end on `-`-strand genes. **Symptom:** a 100 kb gene becomes a 104 kb "promoter"; every `-`-strand promoter is shifted into the gene body. **Fix:** collapse to TSS first, then `slop -s -l UP -r DOWN`.

### slop/flank clip silently at chromosome ends
**Trigger:** fixed-width windows near contig starts / telomeres. **Mechanism:** slop/flank truncate at 0 and chrom length with no warning. **Symptom:** a TSS 800 bp from a contig start yields a 1000-bp (not 2000-bp) upstream window - quietly asymmetric, biasing per-window normalization (reads/kb, motif density); flank can drop a region entirely, breaking a 2:1 feature->flank assumption. **Fix:** after slop verify `end-start == requested width`; after flank verify the per-feature flank count; treat chrom-end features as edge cases.

### -t all double-counts ties into inflated enrichment
**Trigger:** letting default `-t all` rows flow into a per-gene tally, `wc -l` peak count, or GO/hypergeometric enrichment. **Mechanism:** a peak equidistant to two TSSs emits two rows; ties concentrate NON-randomly at bidirectional (head-to-head) promoters and gene-dense regions. **Symptom:** association counts inflated exactly where biology is most interesting; broken independence inflates significance. **Fix:** `-t first` (deterministic but arbitrary - document it) OR `-t all` then aggregate counting distinct PEAKS not rows; for enrichment prefer GREAT/rGREAT, whose regulatory-domain model exists to avoid this artifact.

### closest on unsorted input
**Trigger:** `closest` on a BED that was filtered/edited and not re-sorted. **Mechanism:** closest assumes coordinate-sorted input. **Symptom:** wrong nearest feature or an error. **Fix:** `bedtools sort` (or `.sort()` in pybedtools) both A and B first.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Promoter TSS -2000 / +200 bp (strand-aware) | common convention | a CHOICE, not a fact; asymmetric because core-promoter elements sit upstream and the +1 nucleosome / 5'UTR downstream. Report it; "% promoter-proximal" is sensitive to it |
| GREAT basal domain 5 kb up / 1 kb down, extension <=1 Mb | McLean 2010 *Nat Biotechnol* 28:495 | the principled "proximity++": asymmetric basal domain + extension to the neighbor, a far better proximity heuristic than raw closest |
| ChIPseeker default promoter +-3 kb | tool default | shows the convention spans an order of magnitude (+-500 bp to +-10 kb across tools) |
| ABC candidate window 5 Mb | Fulco 2019 | activity-by-contact scores all elements within 5 Mb of a gene's promoter - "distal" is tens of kb to megabases |
| Nearest gene precision/recall ~47% / ~37% (enhancers) | Fulco 2019 CRISPRi-FlowFISH | the empirical ceiling on nearest-gene for distal-enhancer targeting |
| Nearest protein-coding gene ~50-65% right (GWAS loci) | fine-mapping/coloc literature | the GWAS-regime baseline; strong, hard to beat, but imperfect |
| Distal flag |dist| > ~50-100 kb | field convention | beyond promoter scale; flag for ABC/PCHi-C/eQTL validation rather than trusting nearest |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Nearest-gene call wrong for an enhancer | geometry != regulation for distal elements | treat as candidate; route to atac-seq/enhancer-gene-linking (ABC/PCHi-C/eQTL) |
| Upstream/downstream asymmetry washes out or inverts | `-D ref` mis-signs `-`-strand genes | use `-D b` |
| Promoter window includes the whole gene | `slop -b` on a gene body, not the TSS | collapse to TSS, then `slop -s -l UP -r DOWN` |
| Per-gene counts inflated near bidirectional promoters | `-t all` rows counted as peaks | `-t first` or aggregate by distinct peak; or use GREAT |
| Asymmetric "fixed-width" windows near contig ends | silent slop/flank clipping | verify `end-start`; treat chrom-end features as edge cases |
| `none` / `-1` rows poison a mean distance | no B feature on that chromosome | filter the `-1` sentinel before summarizing |
| Wrong nearest feature, or closest errors | unsorted input | `bedtools sort` both A and B |
| Empty output | `chr1` vs `1` naming mismatch between A, B, genome file | harmonize chromosome naming across all files |

## References

- Quinlan AR, Hall IM. 2010. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26:841-842.
- Dale RK, Pedersen BS, Quinlan AR. 2011. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. *Bioinformatics* 27:3423-3424.
- Fulco CP, Nasser J, Jones TR, et al. 2019. Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations. *Nat Genet* 51:1664-1669.
- Nasser J, Bergman DT, Fulco CP, et al. 2021. Genome-wide enhancer maps link risk variants to disease genes. *Nature* 593:238-243.
- Smemo S, Tena JJ, Kim KH, et al. 2014. Obesity-associated variants within FTO form long-range functional connections with IRX3. *Nature* 507:371-375.
- McLean CY, Bristor D, Hiller M, et al. 2010. GREAT improves functional interpretation of cis-regulatory regions. *Nat Biotechnol* 28:495-501.

## Related Skills

- bed-file-basics - BED coordinate systems and the sort/conversion this skill depends on
- gtf-gff-handling - Extract TSS and gene models from GTF/GFF for promoter construction
- interval-arithmetic - intersect/merge/subtract; window -w 0 is approximately intersect
- chip-seq/peak-annotation - Assigns peaks to genes via the same closest-TSS logic and caveats
- atac-seq/enhancer-gene-linking - The real enhancer->gene science (ABC, contact, peak-gene correlation) this skill routes distal calls to
- atac-seq/footprinting - Uses strand-aware windows over motif/TSS sites
- data-visualization/genome-tracks - Render the promoter/proximity intervals built here
