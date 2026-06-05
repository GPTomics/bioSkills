---
name: bio-cfdna-preprocessing
description: Decides how to preprocess plasma cfDNA sequencing data so the recoverable signal survives - library-prep-aware fragment expectations (dsDNA vs ssDNA/adaptase prep), UMI/duplex consensus with fgbio (ExtractUmisFromBam, GroupReadsByUmi --strategy paired for duplex, CallMolecularConsensusReads vs CallDuplexConsensusReads, FilterConsensusReads min-reads "total s1 s2"), the align->group->consensus->RE-align ordering, and the cfDNA dedup trap where naive coordinate dedup collapses nucleosome-coincident independent molecules. Covers when single-strand consensus suffices vs when duplex is mandatory, the singleton/sensitivity tax at low input, and reading the insert-size histogram as a pre-analytical QC instrument. Use when processing plasma cfDNA reads before fragmentomics, ctDNA mutation calling, or tumor-fraction estimation.
tool_type: mixed
primary_tool: fgbio
---

## Version Compatibility

Reference examples tested with: bwa 0.7.17+, fgbio 2.1+, numpy 1.26+, pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: fgbio flag semantics drift across major versions - in `CallDuplexConsensusReads` the `--min-reads` is a permissive PRE-filter (fgbio issue #1009), and in `FilterConsensusReads` the `--min-reads` (`-M`) takes up to three values "total strand1 strand2" where if values two and three differ the more stringent value must come first. Confirm both against the installed `fgbio --help` before scripting.

# cfDNA Preprocessing

**"Preprocess my plasma cfDNA reads"** -> Convert reads to consensus molecules and an analysis-ready BAM without destroying the fragment-size signal or collapsing independent molecules.
- CLI: `fgbio ExtractUmisFromBam` -> `bwa mem -Y` -> `fgbio GroupReadsByUmi` -> `fgbio Call*ConsensusReads` -> RE-align -> `fgbio FilterConsensusReads`
- Python: shell the fgbio/bwa chain with `subprocess`; read insert sizes with `pysam` for QC

## The Single Most Important Modern Insight -- Pre-Analytics and Library Prep Set the Ceiling; Consensus Is Error Suppression, Not Just Dedup

What a cfDNA assay can ever see is fixed before any bioinformatics runs. cfDNA is not randomly sheared - it is a nucleosome footprint with a mononucleosome mode at ~167 bp and a ctDNA-enriched short tail (134-144 bp), so fragment length is structured biological signal, not noise to be normalized away. Two upstream choices determine what survives: (1) the blood draw and plasma prep (a delayed draw dumps leukocyte gDNA into the denominator and irreversibly dilutes tumor fraction - no algorithm recovers it; defer that QC to liquid-biopsy/analytical-validation), and (2) the library chemistry (dsDNA ligation polishes native ends and discards the sub-100 bp population; ssDNA/adaptase prep recovers short and damaged molecules). Reverse engineering a fragment the prep threw away, or a molecule a bad draw never delivered, is impossible.

The second reframe: UMI/duplex consensus is error suppression, not merely duplicate removal. Grouping reads by UMI and majority-voting a consensus erases PCR/sequencing errors that are not shared across a family - but a single-strand consensus votes unanimously for damage (C->T deamination, G->T 8-oxoG) that was on the template before amplification. Only DUPLEX consensus, requiring the same base on both independently-copied strands, removes those lesions. Single-strand consensus cannot. Choosing simplex vs duplex is choosing an error floor, and it trades against molecular recovery at low input (the singleton tax below).

## Library-Prep and Consensus Landscape

| Choice | What it does | Recoverable distribution / error floor |
|--------|--------------|----------------------------------------|
| dsDNA ligation prep (NEB/KAPA-style) | needs duplex substrate; end-repair/A-tail polishes native ends | clean ~167 bp mode; sub-100 bp tail and native-end signal LOST |
| ssDNA prep, SRSLY/Kircher-Meyer lineage | denatures, ligates single strands; retains native ends + sawtooth | recovers short/nicked/damaged + sub-nucleosomal; prep for fragmentomics |
| ssDNA prep, adaptase/tail-based (Swift/Accel-1S) | single-strand but adaptase chemistry shifts apparent size | ~10 bp short of the canonical mode; sawtooth blunted (a chemistry signature, not a bug) |
| Single-strand UMI consensus (CallMolecularConsensusReads) | majority-vote one source strand's reads | removes PCR/sequencer error (~1e-4 to 1e-5); CANNOT remove deamination/oxidation damage |
| DUPLEX consensus (CallDuplexConsensusReads) | combine both-strand single-strand consensuses | error must occur identically on both strands to survive (~<1e-7); removes deamination/oxidation |

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Deep targeted panel with single-strand UMIs | adjacency group -> CallMolecularConsensusReads -> FilterConsensusReads | simplex consensus suppresses PCR/sequencer error to ~1e-4; the panel-VAF workhorse |
| Need VAF below ~0.1% / MRD-grade specificity | duplex prep -> `--strategy paired` -> CallDuplexConsensusReads -> filter `2 1 1` | duplex removes damage artifacts; reaches the ~1e-7 floor single-strand cannot |
| sWGS / ULP-WGS for tumor fraction | minimal processing: trim -> align -> light dedup; NO consensus | TF from copy number needs even coverage, not error suppression; defer to tumor-fraction-estimation |
| Degraded / low-input / FFPE-adjacent sample | ssDNA prep (SRSLY) to recover short+damaged molecules | dsDNA prep discards exactly the molecules a degraded sample has left |
| Fragmentomics / end-motif readout | ssDNA prep (native ends), NO in-silico size selection upstream | dsDNA fills jagged ends; size-selecting conditions on length and biases every feature |
| Picogram input, detection (not genotyping) | call permissively (`--min-reads 1`), accept singletons | requiring duplicate observation discards the only evidence for a low-VAF variant |
| No UMIs at all, quantitative readout | do NOT coordinate-dedup; document the bias (see Failure Modes) | nucleosome-positioned ends make naive dedup delete real molecules |

Methodology evolves: confirm current fgbio flag semantics and prep-vendor size behavior against live docs before committing a pipeline.

## The fgbio Consensus Pipeline

**Goal:** Turn UMI-tagged raw reads into error-suppressed consensus molecules with correct coordinates.

**Approach:** Extract UMIs into the `RX` tag, align (consensus needs coordinates to group), group by UMI + approximate position, call consensus (which emits UNMAPPED reads because the consensus sequence differs from any input read), RE-align the consensus, then apply the real quality gate with `FilterConsensusReads`. The two-pass alignment is mandatory.

```bash
# 1. Extract inline UMIs into RX. Read structure tokens: M=UMI, S=skip/stem, T=template, +=all remaining.
#    6M11S+T per end = 6 bp UMI, 11 bp stem (the S that bleeds into T if omitted), rest = insert.
fgbio ExtractUmisFromBam --input raw.unmapped.bam --output with_umis.bam \
    --read-structure 6M11S+T 6M11S+T --single-tag RX

# 2. Align. -Y soft-clips supplementaries so tag-bearing short-fragment sequence is not dropped.
bwa mem -t 8 -Y reference.fa with_umis.bam | samtools sort -o aligned.bam -
samtools index aligned.bam

# 3. Group by UMI. adjacency = simplex default; paired = MANDATORY for duplex (reconstructs strand pairing).
fgbio GroupReadsByUmi --input aligned.bam --output grouped.bam \
    --strategy paired --edits 1          # use --strategy adjacency for single-strand UMIs

# 4a. SIMPLEX: single-strand consensus, --min-reads takes ONE value.
fgbio CallMolecularConsensusReads --input grouped.bam --output consensus.unmapped.bam --min-reads 1

# 4b. DUPLEX: --min-reads here is a permissive PRE-filter (fgbio #1009) - call low, filter later.
fgbio CallDuplexConsensusReads --input grouped.bam --output consensus.unmapped.bam --min-reads 1

# 5. RE-align: consensus reads are emitted UNMAPPED by design. Re-map, then ZipperBams
#    transfers the consensus/UMI tags from the unmapped BAM onto the new alignments.
samtools fastq consensus.unmapped.bam | bwa mem -t 8 -Y -p reference.fa - \
  | fgbio ZipperBams --unmapped consensus.unmapped.bam --ref reference.fa \
  | samtools sort -o consensus.bam -

# 6. The REAL quality gate. --min-reads "total strand1 strand2"; "2 1 1" = true duplex (both strands seen).
#    If values two and three differ, the more stringent must come first (e.g. "6 3 0", not "0 3").
fgbio FilterConsensusReads --input consensus.bam --output filtered.bam --ref reference.fa \
    --min-reads 2 1 1 --max-read-error-rate 0.025 --max-base-error-rate 0.1 \
    --min-base-quality 40 --reverse-per-base-tags
```

Key flags: `--strategy paired` is non-negotiable for duplex (`adjacency` cannot reconstruct A/B strand pairing). `--min-reads 1 1 0` on the filter accepts single-strand consensus too (one strand may be absent); `2 1 1` requires both strands = true duplex. `--reverse-per-base-tags` makes per-base depth/error tags read in genomic orientation after alignment.

## Reading the Insert-Size Histogram (QC)

**Goal:** Use the fragment-length distribution as a pre-analytical instrument before trusting any downstream number.

**Approach:** Tabulate proper-pair template lengths from the BAM, locate the mode, and compare the shape against the expected nucleosome footprint - the mode, the sawtooth, and the long-fragment fraction each diagnose a specific failure.

```python
import pysam
import numpy as np

def insert_size_qc(bam_path, max_size=600):
    '''Summarize cfDNA fragment lengths as a QC readout.'''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    sizes = [r.template_length for r in bam.fetch()
             if r.is_proper_pair and not r.is_secondary and 0 < r.template_length <= max_size]
    bam.close()
    sizes = np.array(sizes)
    long_frac = np.mean(sizes > 250)  # excess >250 bp signals gDNA/leukocyte-lysis contamination
    return {'n': len(sizes), 'mode_bp': int(np.bincount(sizes).argmax()),
            'median_bp': float(np.median(sizes)), 'frac_over_250bp': float(long_frac)}
```

A mode at ~167 bp with a ~10.4 bp sawtooth below it is healthy. A mode drifting up plus excess mass >250 bp is gDNA contamination. A mode ~10 bp low is the adaptase chemistry signature, not contamination. A ~120-130 bp spike is adapter dimer (trim adapters before reading the histogram).

## Per-Method Failure Modes

### The dedup trap (naive coordinate dedup on no-UMI cfDNA)
**Trigger:** running Picard MarkDuplicates / `samtools markdup` on cfDNA without UMIs. **Mechanism:** cfDNA ends pile up non-randomly at nucleosome/linker boundaries, so many independent molecules genuinely share the same ~167 bp start+end; coordinate dedup declares them PCR duplicates and keeps one. **Symptom:** deflated unique-molecule counts and depressed VAF for true low-frequency variants - worst exactly where coverage and nucleosome positioning are strongest. **Fix:** use UMIs and group by family; if none, do not dedup by position alone for quantitative readouts and document the bias. This is the strongest single argument for putting UMIs on a cfDNA assay.

### Strict filtering on the consensus caller instead of FilterConsensusReads
**Trigger:** setting `CallDuplexConsensusReads --min-reads` high. **Mechanism:** it is a PRE-filter (fgbio #1009) that discards data before the real filter sees it. **Symptom:** lower molecular recovery than expected with no specificity gain. **Fix:** call permissively (`--min-reads 1`), apply strictness in `FilterConsensusReads`.

### adjacency grouping on duplex data
**Trigger:** `--strategy adjacency` for a duplex library. **Mechanism:** adjacency cannot link the two strands of one duplex (which carry the UMI in opposite order). **Symptom:** duplex consensus finds no two-strand molecules; output collapses to simplex. **Fix:** `--strategy paired`.

### Blanket MAPQ filtering kills short ctDNA fragments
**Trigger:** a global `MAPQ >= 30/60` filter. **Mechanism:** a 40-80 bp insert has fewer anchoring bases and lands in the low-MAPQ tail even when correctly placed. **Symptom:** the short, ctDNA-enriched fragments are preferentially deleted - the molecules size-selection was meant to keep. **Fix:** filter MAPQ with the size distribution in mind; use a gentler threshold for fragmentomics.

### Reflexive `--min-reads >= 2` at low input (the singleton tax)
**Trigger:** requiring duplicate observation per family on a picogram-input library. **Mechanism:** a large fraction of unique cfDNA molecules are sequenced once (singletons); requiring two reads discards genuine, often the only, evidence for a variant. **Symptom:** sensitivity loss disguised as quality. **Fix:** for detection favor recovery (call permissively); reserve `2 1 1` for genotyping/MRD where error suppression dominates.

### Fragmentomics on a size-selected library
**Trigger:** computing end-motif/ratio/VAF features after in-vitro (Pippin) or in-silico size selection. **Mechanism:** selection conditions on length, so length-derived features are biased by construction. **Symptom:** distorted end-motif spectra and fragment-ratio features that do not reproduce. **Fix:** size-select for detection sensitivity only; never report length-derived features from a length-selected library; keep the selection step in metadata.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Mononucleosome mode 167 bp; ~10.4 bp sub-mode periodicity | Snyder 2016 *Cell* 164:57 | core particle (~147 bp) + linker; periodicity = helical pitch of nucleosome-bound DNA, a QC sanity check |
| TF/CTCF-footprint short fragments 35-80 bp | Snyder 2016 *Cell* 164:57 | sub-nucleosomal protection; only surfaces with ssDNA prep |
| ctDNA principal length 134-144 bp vs 167 bp germline | Underhill 2016 *PLoS Genet* 12:e1006162 | tumor fragments shorter; BRAF V600E mutant 132-145 bp vs WT 165 bp |
| ctDNA-enrichment selection windows 90-150 bp (and 250-320 bp) | Mouliere 2018 *Sci Transl Med* 10:eaat4921 | selecting the short window enriches mutant allele fraction at no extra sequencing cost |
| ssDNA prep mitochondrial / microbial cfDNA enrichment ~10.7x / ~71.3x | Burnham 2016 *Sci Rep* 6:27859 | dsDNA prep is blind to the ultrashort fraction these reside in |
| iDES error-suppression gain ~3x (barcode) x ~3x (in-silico) ~= ~15x | Newman 2016 *Nat Biotechnol* 34:547 | family-size and background polishing are complementary, not redundant |
| Duplex error floor <1 per 1e7 nt | Kennedy 2014 *Nat Protoc* 9:2586 | both-strand concordance requirement |
| FilterConsensusReads defaults: max-read-error 0.025, max-base-error 0.1, max-no-calls 0.2 | fgbio docs | per-read vs per-base vs fraction/count switch (<1.0 = fraction, >=1.0 = count) |
| gDNA-contamination signature: excess mass >180-250 bp / ladder distortion | Snyder 2016 (biology); pre-analytical QC | leukocyte-lysis gDNA is long; a "too clean, too long" library is contaminated, not pristine |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Consensus BAM has garbage coordinates | skipped re-alignment (consensus emitted unmapped) | RE-align consensus reads before FilterConsensusReads |
| UMI bases corrupt mapping / consensus | omitted the `S` skip in the read structure | include the stem: e.g. `6M11S+T` not `6M+T` |
| Duplex run finds no two-strand molecules | grouped with `--strategy adjacency` | regroup with `--strategy paired` |
| FilterConsensusReads rejects valid strand spec | `--min-reads 0 3` (less stringent first) | put the more stringent value first: `3 0` (or `6 3 0`) |
| Deflated VAF / low unique-molecule count | coordinate dedup on no-UMI cfDNA | use UMI families; never naive-dedup quantitative cfDNA |
| Short ctDNA fragments missing after filtering | blanket high MAPQ filter | lower MAPQ threshold; account for short-fragment mapping |
| Apparent ~10 bp size shift "needs correcting" | adaptase (Swift/Accel-1S) chemistry signature | record the prep; do not correct a real chemistry effect |

## References

- Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. 2016. Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. *Cell* 164:57-68. -- 167 bp mode, ~10.4 bp periodicity, 35-80 bp TF/CTCF footprints; the nucleosome-positioning basis of the dedup trap.
- Underhill HR, Kitzman JO, Hellwig S, Welker NC, Daza R, Baker DN, Gligorich KM, Rostomily RC, Bronner MP, Shendure J. 2016. Fragment length of circulating tumor DNA. *PLoS Genet* 12:e1006162. -- ctDNA 134-144 bp vs 167 bp germline.
- Mouliere F, Chandrananda D, Piskorz AM, Moore EK, Morris J, Ahlborn LB, Mair R, Goranova T, Marass F, Heider K, et al. 2018. Enhanced detection of circulating tumor DNA by fragment size analysis. *Sci Transl Med* 10:eaat4921. -- 90-150 bp and 250-320 bp ctDNA-enrichment windows.
- Burnham P, Kim MS, Agbor-Enoh S, Luikart H, Valantine HA, Khush KK, De Vlaminck I. 2016. Single-stranded DNA library preparation uncovers the origin and diversity of ultrashort cell-free DNA in plasma. *Sci Rep* 6:27859. -- ssDNA prep recovers the ultrashort mitochondrial/microbial fraction dsDNA prep discards.
- Newman AM, Lovejoy AF, Klass DM, Kurtz DM, Chabon JJ, Scherer F, Stehr H, Liu CL, Bratman SV, Say C, et al. 2016. Integrated digital error suppression for improved detection of circulating tumor DNA. *Nat Biotechnol* 34:547-555. -- barcode + in-silico error suppression are complementary (~3x each).
- Kennedy SR, Schmitt MW, Fox EJ, Kohrn BF, Salk JJ, Ahn EH, Prindle MJ, Kuong KJ, Shen JC, Risques RA, Loeb LA. 2014. Detecting ultralow-frequency mutations by Duplex Sequencing. *Nat Protoc* 9:2586-2606. -- both-strand concordance and the <1e-7 error floor that single-strand consensus cannot reach.

## Related Skills

- analytical-validation - the LoD/molecule-counting framework input quality feeds
- fragment-analysis - fragmentomics consumes the preprocessed fragment ends
- ctdna-mutation-detection - consensus reads feed low-VAF calling
- tumor-fraction-estimation - sWGS minimal-processing path
- alignment-files/duplicate-handling - general dedup vs the cfDNA UMI caveat
- read-qc/quality-reports - upstream read QC
