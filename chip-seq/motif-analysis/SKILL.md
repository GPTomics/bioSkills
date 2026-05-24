---
name: bio-chipseq-motif-analysis
description: Discovers de novo motifs and tests known motif enrichment in ChIP-seq, ATAC-seq, or other peak sequences using HOMER, MEME-ChIP (STREME, CentriMo, TOMTOM, FIMO), monaLisa, and AME. Handles background selection (GC-matched, dinucleotide-shuffled, Markov order-2, peak-flanks), motif databases (JASPAR 2024 CORE PWMs, JASPAR 2026 deep-learning collection, HOCOMOCO v12, HOMER built-in), centrally-enriched motif testing, and differential motif analysis. Use when identifying TF binding motifs in peaks, testing for known TF enrichment, scanning for motif instances, comparing motif content between conditions, or interpreting motifs from deep learning models.
tool_type: cli
primary_tool: HOMER
---

## Version Compatibility

Reference examples tested with: HOMER 4.11+, MEME suite 5.5+ (STREME replaces DREME from 5.4+), monaLisa 1.10+, JASPAR 2024 CORE, HOCOMOCO v12, BioPython 1.83+, bedtools 2.31+.

DREME was removed from MEME suite 5.4+; use STREME instead. Some tutorials still reference DREME — verify the installed version via `meme --version`. JASPAR 2026 (released late 2025) integrates 1259 BPNet ChIP models in a Deep Learning collection; the CORE collection remains the standard PWM source.

# Motif Analysis on ChIP-seq Peaks

**"Find enriched DNA binding motifs in my ChIP-seq peaks"** -> Discover de novo motif patterns and test for known TF motif enrichment in peak sequences, with appropriate background to control for compositional and positional biases.

- CLI (HOMER, fast): `findMotifsGenome.pl peaks.bed hg38 outdir/ -size 200 -p 8`
- CLI (MEME-ChIP, comprehensive): `meme-chip -db JASPAR.meme peaks.fa`
- R (regression-based, selective enrichment): `monaLisa::calcBinnedMotifEnrR(seqs, bins, pwms)`
- CLI (deep-learning-derived motifs): TF-MoDISco on BPNet attribution scores (see chip-deep-learning)

Motif discovery is sensitive to background choice and peak quality. Hyper-ChIPable artifacts at rRNA / housekeeping loci often produce false-positive motifs (GC-rich or A-T-rich biases of those regions). Filter peaks against blacklists and inspect peak distribution before running motif discovery.

## Tool Taxonomy

| Tool | Discovery type | Background handling | Strength | Fails when |
|------|----------------|---------------------|----------|------------|
| **HOMER findMotifsGenome.pl** | De novo + known | GC-matched genomic regions (auto) | Fast (multi-core); integrated vertebrate/insect/plant DBs; one-command full report | Background can include unmasked repeats producing motif artifacts; `-size given` slow; auto background may include peaks themselves |
| **MEME-ChIP** | De novo (STREME, MEME) + central enrichment (CentriMo) + DB comparison (TOMTOM) + scanning (FIMO) | Markov order-2 from input; shuffled (preserves dinucleotide) | Comprehensive single command; rigorous statistics; HTML report | Slower; sequences must be 100-500 bp; central enrichment requires summit-centered peaks |
| **STREME** (MEME 5.4+) | De novo (replaced DREME) | Markov order-2 | Bailey 2021 benchmark: more accurate than DREME/HOMER/MEME/Peak-motifs; handles 3-30 bp; scales to 100k+ sequences | Memory-hungry for very long sequences (>1 kb) |
| **MEME** (classical) | De novo (long, gapped) | Markov | Long motifs; gapped motifs | Slow (no parallel); replaced by STREME for short motifs |
| **DREME** | De novo (short) | Shuffled | Historical; small fast | Removed from MEME 5.4+; use STREME |
| **monaLisa** (Stadler lab) | Binned enrichment regression | Native (binned scoring) | Modern; regression-based; selectivity (TF-specific in differential peaks) | R-only; less integrated with browsers |
| **AME** (MEME suite) | Known motif differential | Matched background set required | Designed for two-set comparison (e.g., peaks vs. control regions) | Requires user-provided background set |
| **CentriMo** | Known motif central enrichment | Auto from input | Tests positional enrichment relative to peak center | Requires summit-centered peaks (200-500 bp) |
| **FIMO** | Motif scanning | Markov model | Genome-wide scanning at user-set p-value | Many false positives at p ≤ 1e-4; tighten to 1e-5 for whole-genome |
| **HOMER scanMotifGenomeWide.pl** | Motif scanning | None | Genome-wide scanning at fixed score threshold | Less calibrated than FIMO; HOMER's PWM format |
| **RSAT peak-motifs** | De novo + known | k-mer comparison | Web-server; multi-tool ensemble | Web limits; less reproducible from CLI |
| **TF-MoDISco** | DL attribution-based | Implicit in model | Motifs from BPNet/chromBPNet attribution scores; captures soft motif syntax | Requires trained DL model; see chip-deep-learning |

## Background Selection — The Biggest Source of Error

Motif enrichment p-values are conditional on the background distribution. Wrong background produces wrong motifs.

| Background | What it preserves | Use case | Limitation |
|------------|-------------------|----------|------------|
| **GC-matched genomic regions** | Mononucleotide composition; chromatin context | TF motifs; avoid GC-bias artifacts | Doesn't preserve dinucleotide (CpG, TpA) |
| **Dinucleotide-shuffled** | CpG and TpA frequencies | Short motifs; avoiding repeat-derived artifacts | Doesn't capture genomic position context |
| **Markov order-2 (trinucleotide)** | Trinucleotide context | Compositional control; STREME/MEME default | Doesn't capture chromatin context |
| **Peak-flanking sequences** (±500 bp upstream/downstream of peak) | Local genomic context | When peak GC differs from genome | May contain shared regulatory motifs if peaks cluster |
| **Repeat-masked input** | Sequence with repeats replaced by N | Avoid TE-derived motif artifacts | Loses motifs in evolved-from-repeat regulatory elements |
| **Input control peaks** | Open-chromatin / artifact regions | TF discrimination from generic chromatin | Hard to obtain; controversial |
| **Differential set** (AME) | Treatment-condition-specific peaks vs ctrl peaks | Differential motif enrichment | Requires a control peak set |

**Practical default:** STREME / MEME-ChIP with Markov order-2 (built-in default); HOMER with `-mask` flag (mask repeats); always inspect for repeat-derived motifs (e.g., Alu-derived AluY consensus, LINE motifs).

## Window Around Summit Matters

Motif enrichment improves dramatically when sequences are summit-centered:

| Window | When |
|--------|------|
| ±100 bp (200 bp total) | Sharp TFs (CTCF, p53); summit reliably reflects motif position |
| ±150-250 bp (300-500 bp) | Most TFs and sharp histones; balance of motif coverage and noise |
| Full peak width (`-size given`) | Variable-width broad marks; computationally expensive |
| Whole gene body (>1 kb) | Almost always wrong; dilutes motif signal |

**For ChIP-seq broad histone marks (H3K27me3, H3K9me3) motif analysis is generally NOT informative** — these marks reflect Polycomb / heterochromatin domains without sequence-specific binding. Motif analysis applies to TFs and to histone marks deposited by sequence-specific cofactors (H3K27ac partial, since BRD4 reads acetyl).

## HOMER Workflow

```bash
# De novo + known motif discovery, repeat-masked, GC-matched background
findMotifsGenome.pl peaks.narrowPeak hg38 homer_out/ \
    -size 200 \
    -mask \
    -p 8

# With user-supplied background (e.g., control peaks or random genomic)
findMotifsGenome.pl peaks.narrowPeak hg38 homer_out/ \
    -size 200 -mask -p 8 \
    -bg background_peaks.bed

# Known motifs only (skip de novo; faster)
findMotifsGenome.pl peaks.narrowPeak hg38 homer_known_only/ \
    -size 200 -mask -nomotif

# Differential motif analysis: peaks gained in condition A vs condition B
findMotifsGenome.pl gained_in_A.bed hg38 differential_motifs/ \
    -size 200 -mask -bg gained_in_B.bed
```

HOMER output files:
- `homerResults.html` — de novo motifs ranked by significance
- `knownResults.html` — known motif enrichment
- `homerMotifs.all.motifs` — all de novo motifs (PWM format)
- `knownResults.txt` — tab-separated known motif stats

## MEME-ChIP Workflow

```bash
# Center peaks to ±100 bp around summit (column 10 in narrowPeak)
awk 'BEGIN{OFS="\t"} {summit = $2 + $10; print $1, summit - 100, summit + 100, $4, $5, $6}' \
    peaks.narrowPeak > peaks_centered.bed
bedtools getfasta -fi hg38.fa -bed peaks_centered.bed -fo peaks_centered.fa

# Full MEME-ChIP analysis
meme-chip \
    -oc meme_chip_out/ \
    -db JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
    -meme-nmotifs 5 \
    -streme-nmotifs 10 \
    -minw 6 -maxw 20 \
    peaks_centered.fa

# MEME-ChIP runs: STREME (replaces DREME) + MEME + CentriMo + TOMTOM + FIMO
# CentriMo tests known motifs for central enrichment — strongest signal of
# direct binding vs. tethered/indirect binding
```

## monaLisa Workflow (R, Regression-Based)

**Goal:** Test which TF motifs are enriched in specific bins of peaks (e.g., bins of differential log2FC, or bins of accessibility).

**Approach:** monaLisa builds a per-motif regression of peak signal on motif occurrence, controlling for GC content. Selective for TFs that discriminate between bins.

```r
library(monaLisa)
library(JASPAR2024)
library(Biostrings)

# Load peaks and split into bins (e.g., quintiles of log2FC)
peaks <- rtracklayer::import('peaks.bed')
peaks$log2FC <- ...  # from differential analysis
bins <- bin(peaks$log2FC, binmode = 'equalN', nElement = 200)

# Get sequences around peak centers
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, resize(peaks, width = 500, fix = 'center'))

# Load JASPAR PWMs
pwms <- getMatrixSet(JASPAR2024, list(species = 9606, collection = 'CORE'))

# Compute binned motif enrichment with GC control
res <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms, BPPARAM = MulticoreParam(8))

# Plot heatmap of motif enrichment vs bins
plotMotifHeatmaps(x = res, which.plots = c('log2enr', 'negLog10P'),
                   width = 1.8, maxEnr = 2, maxSig = 10)
```

## Per-Tool Failure Modes

### HOMER -- Background includes peaks themselves

**Trigger:** Running `findMotifsGenome.pl` without `-bg` on a large peak set covering >5% of genome.

**Mechanism:** HOMER auto-samples GC-matched genomic regions for background, which can overlap the peak set itself.

**Symptom:** Even known TF motif p-values are weak (>1e-3); de novo motifs less enriched than expected.

**Fix:** Supply explicit `-bg` background (e.g., random genomic intervals matching peak count and width) OR use MEME-ChIP with internal shuffled background.

### HOMER / MEME -- Repeat-derived false-positive motifs

**Trigger:** Running on unmasked peaks; peaks cover transposable elements (Alu, LINE, LTR).

**Mechanism:** TEs contain over-represented k-mers that motif algorithms mistake for biology.

**Symptom:** Top de novo motif matches AluY consensus (~280 bp), LINE/L1, or LTR families.

**Fix:** Use `-mask` (HOMER) or pre-mask peak sequences with RepeatMasker; verify TOMTOM matches against legitimate TF databases.

### STREME -- Memory failure on long sequences

**Trigger:** Running STREME on full-peak sequences (>1 kb each) with > 50k peaks.

**Mechanism:** STREME holds suffix structures in memory; long sequences explode RAM.

**Fix:** Resize peaks to ±100-250 bp summit-centered before STREME; or downsample peak count.

### CentriMo -- No central enrichment due to wrong centering

**Trigger:** Using full peak coordinates (BED `start, end`) without recentering on summit.

**Mechanism:** Peak start coordinate is the left edge, not the summit; motif may be enriched near the summit but appears unenriched relative to the start.

**Symptom:** Known TF motifs show no central enrichment despite being clearly enriched overall.

**Fix:** Recenter sequences on summit: `summit = start + summit_offset` (narrowPeak column 10) before extracting FASTA.

### FIMO -- Massive false positives at default p ≤ 1e-4

**Trigger:** Genome-wide scanning at FIMO default p-value threshold.

**Mechanism:** At p ≤ 1e-4, expect ~3M random matches in a 3 Gb genome; most are false positives.

**Symptom:** FIMO output has millions of motif "hits"; can't distinguish real binding.

**Fix:** Tighten to `--thresh 1e-5` or stricter for genome-wide scans; or restrict scan to peaks: `fimo --bgfile motif_bg motif.meme peaks.fa`.

### monaLisa -- GC bins not respected

**Trigger:** Running `calcBinnedMotifEnrR` without GC binning when peaks have systematic GC differences (e.g., promoters vs distal enhancers).

**Mechanism:** GC-rich motifs are over-enriched in GC-rich bins simply by chance.

**Symptom:** Top motifs are CpG-rich families (e.g., E2F, NRF1) regardless of biology.

**Fix:** Use `background = 'genome'` with GC-matched genomic background; or `background = 'otherBins'` with stratified GC.

### Motif analysis on hyper-ChIPable regions

**Trigger:** Running motifs on peaks dominated by hyper-ChIPable artifacts (rRNA, tRNA, housekeeping).

**Mechanism:** These regions have systematic compositional biases (GC-rich, A-T-rich, repeat-derived); motifs from artifact peaks reflect compositional biology, not TF binding.

**Symptom:** Top de novo motif matches no known TF; high-GC or low-complexity consensus.

**Fix:** Apply blacklist + custom hyper-ChIPable filter (top-1% input signal) before motif analysis. See chipseq-qc.

## Reconciliation: When Motif Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| HOMER finds motif X; MEME-ChIP misses | Different background; HOMER may have permissive background | Run MEME-ChIP with explicit GC-matched background; check |
| MEME finds long gapped motif; STREME doesn't | MEME captures variable-length structure; STREME limited to 30 bp | Both are correct; report MEME for long motifs |
| Top de novo motif doesn't match TOMTOM databases | Novel motif OR repeat artifact OR compositional artifact | Inspect peaks for repeats; check input control; could be genuine novel TF |
| Known motif enriched but no de novo recovery | Insufficient enrichment for de novo; or motif is degenerate | Trust known motif enrichment; de novo needs strong signal |
| Differential motif gained in treatment but TF expression unchanged | TF post-translational regulation (binding mode change without expression change) | Check ChIP signal at known TF target genes; not a contradiction |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| HOMER "configureHomer.pl genome not installed" | Genome not configured | `perl configureHomer.pl -install hg38` (one-time) |
| MEME "sequence too short" | Peaks < motif min width | Resize peaks to ≥ 200 bp |
| MEME-ChIP "out of memory" | Too many long sequences | Resize peaks to ±100-250 bp; downsample |
| No enriched motifs | Peak quality / hyper-ChIPable / wrong background | Check FRiP, filter blacklist, supply explicit background |
| Top motif is GC-rich consensus | GC-bias in peaks not matched by background | GC-matched background (HOMER `-bg` or MEME shuffled with order-2) |
| FIMO produces millions of hits | p-value threshold too loose | `--thresh 1e-5` for whole-genome; restrict to peaks for finer p |
| TOMTOM matches always say "no match" | Motif database species mismatch | Use vertebrates / insects / plants DB matching organism |

## References

- Heinz S et al 2010 Mol Cell 38:576 (HOMER)
- Bailey TL & Elkan C 1994 Proc ISMB (MEME)
- Bailey TL 2021 Bioinformatics 37:2834 (STREME)
- Machanick P & Bailey TL 2011 Bioinformatics 27:1696 (MEME-ChIP)
- Bailey TL et al 2015 Nucleic Acids Res 43:W39 (MEME suite update)
- Grant CE et al 2011 Bioinformatics 27:1017 (FIMO)
- Bailey TL & Machanick P 2012 Nucleic Acids Res 40:e128 (CentriMo)
- McLeay RC & Bailey TL 2010 BMC Bioinformatics 11:165 (AME)
- Machlab D et al 2022 Nucleic Acids Res 50:e49 (monaLisa)
- Castro-Mondragon JA et al 2022 Nucleic Acids Res 50:D165 (JASPAR 2022; CORE collection)
- Avsec Ž et al 2021 Nat Genet 53:354 (BPNet; soft motif syntax)
- Shrikumar A et al 2020 bioRxiv (TF-MoDISco)

## Related Skills

- chip-seq/peak-calling - Upstream peak calling; recenter on summit for motif input
- chip-seq/chipseq-qc - Filter hyper-ChIPable artifacts before motif discovery
- chip-seq/chip-deep-learning - BPNet/chromBPNet for sequence-attribution motif discovery (TF-MoDISco)
- chip-seq/peak-annotation - Annotate peaks before motif discovery to filter promoters vs enhancers
- atac-seq/motif-deviation - chromVAR per-cell motif activity (ATAC-specific)
- atac-seq/footprinting - TOBIAS footprint analysis (ATAC; complementary to motif enrichment)
- sequence-manipulation/motif-search - General sequence motif scanning
- genome-intervals/proximity-operations - bedtools getfasta to extract peak sequences
