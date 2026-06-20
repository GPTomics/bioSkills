---
name: bio-ribo-seq-ribosome-periodicity
description: Validate Ribo-seq library quality by measuring 3-nucleotide periodicity and calibrating read-length-specific P-site offsets. Use when checking whether footprints capture genuine translation, determining P-site offsets for downstream ORF/TE/stalling analysis, or deciding which read lengths to keep.
tool_type: mixed
primary_tool: riboWaltz
---

## Version Compatibility

Reference examples tested with: riboWaltz 2.0+, plastid 0.6+, numpy 1.26+, scipy 1.12+, pysam 0.22+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Ribosome Periodicity and P-site Calibration

**"Check if my Ribo-seq data shows triplet periodicity and get my P-site offsets"** -> Confirm footprints carry codon phase (the signature of genuine elongating ribosomes) and compute the read-length-specific offset from the read end to the P-site codon, the prerequisite for every codon-resolution analysis.
- R: `riboWaltz` for P-site offset calibration and per-length periodicity (the de-facto standard)
- Python: `plastid` `metagene` + `psite` CLI scripts as the alternative path

## Why periodicity is the QC gate

An elongating ribosome advances exactly one codon (3 nt) per translocation, so P-site-assigned footprints over a CDS pile up in one reading frame (frame 0 >> frames +1/+2). This sub-codon comb is what distinguishes Ribo-seq from RNA-seq, and a library without it cannot support frame-based ORF calling, TE, or dwell-time work regardless of read depth. Contamination (rRNA/tRNA), degraded RNA, and over-digestion all give phase-free, RNA-seq-like coverage.

## P-site geometry (what is being calibrated)

The ribosome has three tRNA sites: A (aminoacyl, decodes the incoming codon), P (peptidyl, holds the nascent chain), E (exit). Codon position is reported at the P-site, which sits at a fixed OFFSET inside the ~28-30 nt footprint. The offset is the distance from the mapped read end to the first nucleotide of the P-site codon. A-site offset = P-site + 3; E-site = P-site - 3. The A-site is the relevant site for tRNA/decoding effects (see ribosome-stalling).

## The decisions that shape periodicity QC

### P-site offset method

| Method | Map from | Best when | Caveat |
|--------|----------|-----------|--------|
| 5'-end + offset | 5' end | sharp 5' ends, classic RNase I libraries (~+12 for 28 nt) | breaks if the 5' end is ragged/variably trimmed |
| 3'-end + offset | 3' end | variable 5' trimming, sharper 3' end; standard for bacteria/MNase | offset still length-dependent; verify per length |
| auto (riboWaltz) | 5' or 3', chosen per length | default; let the data pick the more consistent end | reports both, decides per read-length population |
| center (plastid) | both ends | very noisy ends | loses sub-codon sharpness |

The canonical ~+12 nt offset for 28-29 nt mammalian footprints is a STARTING expectation, not a constant. Offsets are read-length-specific and dataset-specific and must be calibrated empirically; a fixed lookup table silently misassigns codons.

### Tool choice

| Tool | Language | Role |
|------|----------|------|
| riboWaltz | R | offset calibration + per-length frame % (primary) |
| plastid | Python | `metagene`/`psite` CLI offsets + count vectors |
| Ribo-seQC | R | one-shot HTML QC report (P-site, region, periodicity) |
| ribotricer | Python | phase-score check that is robust to P-site shift |

## Calibrate offsets and frame with riboWaltz

**Goal:** Determine the per-read-length P-site offset and the frame-0 fraction that together certify the library.

**Approach:** Convert BAMs to riboWaltz tables, filter lengths by periodicity, compute offsets with `extremity="auto"`, then read off frame percentages per length.

```r
library(riboWaltz)

annotation <- create_annotation(gtfpath = "annotation.gtf")
reads <- bamtolist(bamfolder = "bams", annotation = annotation)

# Keep only read lengths with strong frame-0 enrichment
# periodicity_threshold is a frame-0 percentage (here 50%); tune per dataset
reads <- length_filter(reads, length_filter_mode = "periodicity",
                       periodicity_threshold = 50)

# extremity="auto" picks the 5' or 3' end giving the most consistent per-length offset;
# the corrected offset refines the temporary one to the local maximum (occupancy correction)
offsets <- psite(reads, flanking = 6, extremity = "auto")
reads_psite <- psite_info(reads, offsets)

# Frame-0 fraction per read length is the primary, defensible periodicity metric
frames_by_length <- frame_psite_length(reads_psite, annotation,
                                       sample = names(reads)[1])
```

## Read off the periodicity metrics

**Goal:** Decide pass/fail and which lengths to retain.

**Approach:** Use the frame-0 fraction as the headline number and the metaheatmap/metaprofile as visual confirmation.

```r
# Pooled frame distribution and the start/stop metaprofile
frames <- frame_psite(reads_psite, annotation, sample = names(reads)[1])
metaprofile_psite(reads_psite, annotation, sample = names(reads)[1],
                  utr5l = 25, cdsl = 40, utr3l = 25)
metaheatmap_psite(reads_psite, annotation, sample = names(reads)[1])
```

Frame-0 fraction rule of thumb: good libraries put roughly >60-70% of in-CDS P-sites in frame 0 (vs the 33% null); ~45-60% is marginal; near-uniform 33/33/33 is uninterpretable at codon level. Report per read length, not just pooled.

If NO read length clears the periodicity threshold (length_filter returns empty), the library is RNA-seq-like and supports only gene-level counting, not codon-resolution ORF/TE/stalling analysis; that is the verdict, not a reason to keep lowering the threshold. Lower it only to inspect the best-available length, not to rescue an aperiodic library.

A bimodal length distribution is expected, not an error: alongside the ~28-30 nt footprint there is a ~21 nt population from ribosomes with an open (empty) A-site (Lareau 2014). Inspect the ~21 nt class per length rather than discarding it as contamination; its phase and offset differ from the long footprints and it carries elongation-state information.

## Alternative: plastid offsets via the CLI

**Goal:** Get per-length offsets without R, using plastid's verified workflow.

**Approach:** Build a start-codon ROI with `metagene generate`, then run the `psite` script, which writes an offsets table and per-length profile plots.

```bash
# CLI-first: there is NO top-level plastid.metagene_analysis() function
metagene generate cds_start --landmark cds_start --annotation_files annotation.gtf
psite cds_start_rois.txt psite_out --min_length 26 --max_length 34 \
    --require_upstream --count_files riboseq.bam
```

```python
# Apply the calibrated offsets in Python
from plastid import BAMGenomeArray, VariableFivePrimeMapFactory, GTF2_TranscriptAssembler

ga = BAMGenomeArray('riboseq.bam')
ga.set_mapping(VariableFivePrimeMapFactory.from_file(open('psite_out_p_offsets.txt')))
transcripts = list(GTF2_TranscriptAssembler('annotation.gtf'))
# Per-transcript P-site counts: vec = transcript.get_counts(ga)
```

## Compute a body-coverage periodicity score

**Goal:** Quantify periodicity strength from the CDS body, not the initiation peak.

**Approach:** Build per-nucleotide P-site coverage along the CDS, trim the start/stop peaks, then take the frame-0 fraction or the spectral power at period 3.

```python
import numpy as np

def body_frame_fraction(psite_coverage, trim_start=45, trim_stop=15):
    '''Frame-0 fraction over CDS-body P-site coverage.

    The start (initiation) and stop (termination) peaks dwarf the body and carry
    their own phase, so they are trimmed (trim in nt; ~15 codons start, ~5 codons stop).
    '''
    body = psite_coverage[trim_start:len(psite_coverage) - trim_stop]
    frames = [body[f::3].sum() for f in range(3)]
    total = sum(frames)
    return frames[0] / total if total else 0.0
```

Running an FFT on the start-codon metagene is the wrong signal: that profile is dominated by a single initiation peak, not sustained codon phase. The spectral test must run on uniform CDS-body P-site coverage; otherwise report the frame-0 fraction directly.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `ImportError: cannot import name 'metagene_analysis'` | No such function exists in plastid | Use the `metagene generate` + `psite` CLI, or riboWaltz |
| Periodicity "score" always ~0 or meaningless | FFT run on the start-codon metagene, or frames never populated | Score CDS-body P-site coverage; use frame_psite_length |
| Offset works for one length, breaks others | A single hardcoded offset (e.g. 12) applied to all lengths | Calibrate per read length; A-site = P-site + 3 |
| Strong "periodicity" that is just the start peak | Start/stop codon peaks not trimmed | Trim ~15 codons at start, ~5 at stop before scoring |
| Bacterial library looks aperiodic | MNase data with ragged 5' ends mapped 5'-anchored | Anchor on the 3' end; expect weaker periodicity than RNase I |
| Short/long read lengths dilute the signal | Phase-free length tails kept in the analysis | length_filter mode "periodicity"; analyze per length |

## Related Skills

- riboseq-preprocessing - Produce the aligned BAM and inspect the read-length distribution
- orf-detection - Consumes the per-length P-site offsets to call translated ORFs
- translation-efficiency - Needs correct P-site positioning for CDS footprint counts
- ribosome-stalling - Uses the calibrated A-site offset for codon occupancy

## References

- Lauria F, Tebaldi T, Bernabò P, Groen EJN, Gillingwater TH, Viero G. 2018. riboWaltz: Optimization of ribosome P-site positioning in ribosome profiling data. PLoS Comput Biol 14(8):e1006169. doi:10.1371/journal.pcbi.1006169
- Dunn JG, Weissman JS. 2016. Plastid: nucleotide-resolution analysis of next-generation sequencing and genomics data. BMC Genomics 17(1):958. doi:10.1186/s12864-016-3278-x
- Calviello L, Sydow D, Harnett D, Ohler U. 2019. Ribo-seQC: comprehensive analysis of cytoplasmic and organellar ribosome profiling data. bioRxiv 601468. doi:10.1101/601468
- Ingolia NT, Ghaemmaghami S, Newman JRS, Weissman JS. 2009. Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling. Science 324(5924):218-223. doi:10.1126/science.1168978
- Lareau LF, Hite DH, Hogan GJ, Brown PO. 2014. Distinct stages of the translation elongation cycle revealed by sequencing ribosome-protected mRNA fragments. eLife 3:e01257. doi:10.7554/eLife.01257
