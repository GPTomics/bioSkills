---
name: bio-epitranscriptomics-m6a-peak-calling
description: Calls m6A peaks from MeRIP-seq / m6A-seq paired IP-vs-input data using exomePeak2 (transcript-aware, GC-bias-corrected Poisson GLM), MeTPeak (HMM over sliding windows), MACS3/MACS2 with --nomodel --broad --keep-dup all (genome-wide broad alternative), and DRACH motif enrichment via HOMER or ggseqlogo as a sanity check (NOT a filter). Covers BED12 vs narrowPeak output, exonic vs intronic peak handling, multi-tool reconciliation (intersection vs union), the m6A-vs-m6Am ambiguity at 5'UTR peaks that antibody methods cannot resolve, and orthogonal validation (miCLIP/GLORI/m6A-SAC-seq/m6Anet). Use when calling peaks from paired IP/input genome BAMs, choosing exomePeak2 (transcript-aware default) vs MACS3 (broad genomic) vs MeTPeak (HMM-smoothed low-coverage), confirming DRACH enrichment as a sanity check on the peak set, reconciling differing peak sets across tools, validating MeRIP peaks against single-base methods, interpreting 5' peaks where m6Am contamination is possible, or recommending a consensus strategy.
tool_type: mixed
primary_tool: exomePeak2
---

## Version Compatibility

Reference examples tested with: exomePeak2 1.14+ (Bioconductor 3.18+), MeTPeak (GitHub commit SHA-pinned; no Bioconductor release), MACS3 3.0+, MACS2 2.2.9+, samtools 1.19+, GenomicFeatures 1.54+, BSgenome.Hsapiens.UCSC.hg38 1.4+, HOMER 4.11+, ggseqlogo 0.2+, rtracklayer 1.62+, GenomicRanges 1.54+, ChIPseeker 1.38+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('exomePeak2')` then `?exomePeak2` to verify parameters
- CLI: `macs3 callpeak --help`, `findMotifsGenome.pl` to confirm flags

If R throws `unused argument` or `argument is missing`, the exomePeak2 / MeTPeak API may have moved between releases; consult `?exomePeak2` and the installed package's NAMESPACE.

exomePeak2 differential is triggered by populating `bam_treated_ip` / `bam_treated_input` alongside the control-arm `bam_ip` / `bam_input` (see m6a-differential); there is NO `mode=` argument. MeTPeak is GitHub-only and unversioned — pin a commit SHA in any reproducible analysis; defaults are `WINDOW_WIDTH=50, SLIDING_STEP=50, FRAGMENT_LENGTH=100`. MACS3 supersedes MACS2 in development activity but both are widely used; default `--keep-dup` is 1 in BOTH and MUST be overridden to `all` for MeRIP. HOMER `findMotifsGenome.pl` is mature; the `-rna` flag is the RNA-mode entry.

# m6A Peak Calling

**"Find m6A sites in my MeRIP data"** -> Compare paired IP and input read distributions per transcript window, call windows with statistically significant IP enrichment as peaks, annotate against transcript features (5'UTR / CDS / 3'UTR / stop-codon), and confirm DRACH motif enrichment on the peak set relative to background as a sanity check that the antibody-IP worked. CRITICAL: DRACH is a sanity check on the WHOLE peak set, NOT a filter on individual peaks (filtering by DRACH drops 5-10% of real m6A sites and creates circular validation). Equally critical: peaks within ~50 nt of TSS may be PCIF1 m6Am, not METTL3 m6A — anti-m6A antibodies cross-react.

- R: `exomePeak2::exomePeak2()` -- transcript-aware GC-corrected GLM (field default)
- R: `MeTPeak::metpeak()` -- HMM-smoothed sliding-window alternative
- CLI: `macs3 callpeak --nomodel --keep-dup all --broad` -- broad genomic alternative
- CLI: `findMotifsGenome.pl peaks.bed -rna` (HOMER) -- DRACH enrichment sanity check
- R: `ggseqlogo::ggseqlogo()` -- peak-centre 5-mer sequence logo

## The Single Most Important Modern Insight -- MeRIP cannot distinguish m6A from m6Am near the 5' end, and DRACH is a sanity check not a filter

Anti-m6A antibodies (Synaptic Systems 202-003, Abcam ab151230, NEB EpiMark E1610, Cell Signaling 56593, Active Motif 61755) cross-react with m6Am — the cap-adjacent N6,2'-O-dimethyladenosine at the +1 nucleotide of capped mRNAs, installed by PCIF1 / CAPAM (Akichika 2019 *Science* 363:eaav0080; Boulias 2019 *Mol Cell* 75:631). METTL3 / METTL14 do NOT methylate the cap-adjacent position; PCIF1 does. Peaks within the first ~50 nt of a transcript are ambiguous between m6A and m6Am. METTL3-KO will REMOVE internal m6A peaks but LEAVE the cap m6Am peaks intact, which has caused multiple papers to mis-attribute METTL3-independent peaks to non-canonical writers when they are really PCIF1 m6Am. Linder 2015 *Nat Methods* 12:767 (miCLIP, the antibody-based single-base method) explicitly notes the cross-reactivity. Separately: while ~70% of mammalian METTL3-deposited m6A sites sit within the DRACH consensus (D=A/G/U, R=A/G, A=methylated, C=C, H=A/C/U), a non-trivial fraction are non-DRACH. Post-hoc filtering peaks by DRACH content drops real m6A peaks; DRACH should be reported as enrichment-relative-to-background (HOMER / MEME / ggseqlogo) on the PEAK SET as a whole, NOT as a per-peak filter. For unambiguous internal-m6A studies, restrict analysis to peaks past the first ~50 nt AND validate at high-stakes sites with an orthogonal method (miCLIP for single-base antibody validation; GLORI Liu 2023 *Nat Biotechnol* 41:355 for absolute stoichiometry; m6A-SAC-seq Hu 2022 *Nat Biotechnol* 40:1210; m6Anet for ONT direct-RNA confirmation).

## Algorithmic Taxonomy

| Tool / mode | Mechanism | Inputs | Output | Strength | Fails when |
|-------------|-----------|--------|--------|----------|------------|
| exomePeak2 (Liu 2022 *NAR Genom Bioinform* 4:lqac046) | Transcript-windowed Poisson GLM with on-the-fly GC-bias correction; supersedes exomePeak v1 | paired IP/input BAM + TxDb / GTF | BED12 + RDS + per-peak fold-change / FDR | Transcript-aware; GC-aware; integrates motif annotation; modern field default | Slow on very large datasets; argument signatures shift between Bioconductor minor releases — verify against `?exomePeak2` |
| MeTPeak (Cui 2016 *Bioinformatics* 32:i378) | HMM over sliding windows with Beta-binomial emission per window | paired IP/input BAM + GTF (or TxDb) | BED12 | HMM smooths spatial dependency; better at low coverage | GitHub-only; unversioned; default window/step is 50/50 not the small values some tutorials suggest |
| MACS3 / MACS2 broad mode (Zhang 2008 *Genome Biol* 9:R137) | Sliding-window negative-binomial test; `--broad` extends; `--nomodel` disables ChIP fragment-shift model | paired IP/input BAM | narrowPeak / broadPeak | Battle-tested ChIP-seq lineage; very fast | Not transcript-aware; misses GC-confounded peaks; default `--keep-dup 1` collapses MeRIP signal at high-coverage transcripts |
| MeRIPtools (R wrapper) | Bundles exomePeak / MeTPeak / motif / annotation steps | FASTQ -> peaks pipeline | full report | Reproducible end-to-end | Less flexibility than calling tools separately |
| MoAIMS | Mixture model alternative | paired IP/input BAM | peaks | Smaller user base; less benchmarked | Niche use |
| m6Aboost (R) | Boost peak-calling sensitivity by leveraging DRACH motif as a prior | paired IP/input BAM + motif file | refined peak set | Improves sensitivity in low-coverage regions | Builds DRACH into the calling — DON'T use as evidence for DRACH enrichment downstream (circular) |
| m6ACali (recent ML peak filter; verify current citation against the project repo) | ML-based false-peak filter trained on exomePeak2 + MACS2 outputs across many cell lines | called peak set + IP/input BAM | refined peak set | Modern QC layer; cuts antibody artifact peaks | Trained on specific antibody clones; verify it generalises to the antibody used |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| Standard mammalian MeRIP, 3+ replicates per arm | exomePeak2 + DRACH confirmation; reconcile with MeTPeak as a second opinion | MACS3 misses GC-confounded peaks; MeTPeak alone gives less GC awareness |
| Viral / kilobase-broad peaks | MACS3 `--broad --broad-cutoff 0.1 --keep-dup all` | exomePeak2 splits broad enrichment into many small per-window peaks |
| Low-coverage / scarce-sample MeRIP | MeTPeak HMM smooths; exomePeak2 as cross-check | MACS3 default loses sensitivity at low coverage |
| Need single-nucleotide resolution | NOT MeRIP -- switch to miCLIP, m6Anet, GLORI, SAC-seq | MeRIP windows are ~100-200 nt; cannot resolve to base |
| 5'UTR / cap-proximal peaks | Run normally BUT flag 5' peaks (within ~50 nt of TSS) as m6A-or-m6Am ambiguous; validate with PCIF1-KO if available | Antibody cross-reacts with m6Am; cannot assign without orthogonal data |
| Cross-tool reconciliation | Call with exomePeak2 + MeTPeak + MACS3 broad; intersect (NOT union) for high-confidence; report each separately | Union inflates false-positive rate; single-tool reports under-call ~30% of consensus peaks |
| Validation of high-stakes peak set | Cross-check against published m6A-Atlas / REPIC databases; orthogonal validation at top hits (miCLIP / GLORI / m6A-SAC-seq) | Single-method-single-study peaks have ~50% inter-study overlap (McIntyre 2020) |
| Wanting absolute stoichiometry | NOT MeRIP -- use GLORI (Liu 2023 *Nat Biotechnol*), SAC-seq (Hu 2022), eTAM-seq (Xiao 2023) | MeRIP IP fold-change is relative, not absolute |
| METTL3-KO validation experiment | Call peaks in WT and KO; peaks that DISAPPEAR in KO are m6A-dependent; peaks that REMAIN are antibody artifacts OR m6Am OR non-METTL3 modifications (METTL16, METTL5) | Calling only WT and assuming all peaks are m6A; many are not |
| Anti-m6A vs anti-m6Am specific analysis | For internal m6A only: exclude TSS-proximal peaks AND use GLORI / SAC-seq orthogonal validation; for m6Am specifically: m6Am-seq (Sun H 2021 *Nat Commun* 12:4778), m6ACE-seq, or PCIF1-KO subtraction | Antibody alone CANNOT distinguish — this is a chemistry problem, not a software problem |

Methodology evolves; before any high-stakes peak-calling analysis, web-search "exomePeak2 Bioconductor release notes" and "MeRIP peak caller benchmark 2024" for current consensus parameters.

## exomePeak2 Standard Workflow

**Goal:** Produce a transcript-aware set of m6A peaks with FDR and IP/input fold-change from paired IP/input genome BAM files, suitable as input to differential analysis, motif scanning, or downstream visualisation.

**Approach:** Build a TxDb from the matched GTF; pass paired IP/input BAM vectors with `bam_ip` and `bam_input`; let exomePeak2 handle GC-bias correction internally; export BED12 + RDS for downstream use; annotate against transcript features.

```r
library(exomePeak2)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb <- makeTxDbFromGFF('refs/annotation.gtf', format='gtf')

result <- exomePeak2(
    bam_ip       = c('aligned/IP_rep1.bam', 'aligned/IP_rep2.bam', 'aligned/IP_rep3.bam'),
    bam_input    = c('aligned/Input_rep1.bam', 'aligned/Input_rep2.bam', 'aligned/Input_rep3.bam'),
    txdb         = txdb,
    genome       = BSgenome.Hsapiens.UCSC.hg38,
    paired_end   = TRUE,
    library_type = 'unstranded',
    save_dir     = 'exomepeak2_output',
    experiment_name = 'm6a_run1'
)

# Inspect peaks: GRanges with metadata
peaks <- result
length(peaks)
head(as.data.frame(peaks))
```

`exomePeak2()` writes BED12 + RDS + per-peak fold-change / FDR to `save_dir/experiment_name/`. The `txdb` argument is the canonical interface; older tutorials use `gff=` (path) which is deprecated. `BSgenome` is required for GC correction; if unavailable, exomePeak2 falls back to a less-accurate GC-uncorrected mode.

## MeTPeak HMM-Smoothed Alternative

**Goal:** Call peaks with HMM smoothing across spatial windows, useful when coverage is low or when peak boundaries are ambiguous from per-window negative-binomial tests alone.

**Approach:** MeTPeak accepts IP and INPUT BAM vectors plus either a GTF file path (`GENE_ANNO_GTF=`) or a TxDb object (`TXDB=`); defaults are `WINDOW_WIDTH=50, SLIDING_STEP=50, FRAGMENT_LENGTH=100, MINIMAL_PEAK_LENGTH=FRAGMENT_LENGTH/2, MINIMAL_MAPQ=30`. The example below uses GTF input with explicit defaults shown for transparency; adjust only with rationale.

```r
library(MeTPeak)

metpeak(
    IP_BAM              = c('aligned/IP_rep1.bam', 'aligned/IP_rep2.bam', 'aligned/IP_rep3.bam'),
    INPUT_BAM           = c('aligned/Input_rep1.bam', 'aligned/Input_rep2.bam', 'aligned/Input_rep3.bam'),
    GENE_ANNO_GTF       = 'refs/annotation.gtf',
    OUTPUT_DIR          = 'metpeak_output',
    EXPERIMENT_NAME     = 'm6a_run1',
    WINDOW_WIDTH        = 50,
    SLIDING_STEP        = 50,
    FRAGMENT_LENGTH     = 100,
    MINIMAL_PEAK_LENGTH = 50,
    PEAK_CUTOFF_PVALUE  = 1e-5,
    PEAK_CUTOFF_FDR     = 0.05,
    FOLD_ENRICHMENT     = 1
)
```

MeTPeak accepts either `GENE_ANNO_GTF=` (a file path) or `TXDB=` (a TxDb object); the GTF path is the more common usage and is documented in the GitHub README. Output BED12 is at `OUTPUT_DIR/EXPERIMENT_NAME/peak.bed`. MeTPeak is GitHub-only (`compgenomics/MeTPeak`) and unversioned; pin a commit SHA in reproducible analyses.

## MACS3 Broad-Peak Alternative

**Goal:** Call broad MeRIP peaks across the genome using the ChIP-seq sliding-window negative-binomial framework; useful for viral genomes where peaks span kilobases, or as a cross-caller second opinion.

**Approach:** MACS3 `callpeak` with `--nomodel --extsize 150 --keep-dup all --broad`. The `--nomodel` flag disables the ChIP-seq fragment-shift model (designed for DNA); `--keep-dup all` is REQUIRED because MACS3 default deduplicates and collapses MeRIP signal at high-coverage transcripts.

```bash
macs3 callpeak \
    --treatment aligned/IP_rep1.bam aligned/IP_rep2.bam aligned/IP_rep3.bam \
    --control   aligned/Input_rep1.bam aligned/Input_rep2.bam aligned/Input_rep3.bam \
    --format BAMPE \
    --gsize hs \
    --nomodel \
    --extsize 150 \
    --keep-dup all \
    --broad \
    --broad-cutoff 0.1 \
    --qvalue 0.05 \
    --outdir macs3_output \
    --name m6a_run1
```

`--keep-dup all` is non-negotiable for MeRIP; the default `--keep-dup 1` (keep one read per position) destroys real signal. `--broad` and `--broad-cutoff 0.1` extend narrow peaks into broader regions, appropriate for MeRIP fragments. Output: narrowPeak + broadPeak + gappedPeak files in `macs3_output/`.

## DRACH Motif Confirmation (Sanity Check, NOT a Filter)

**Goal:** Confirm that the called peak set is enriched for the DRACH consensus motif relative to genomic background, validating antibody specificity. NEVER post-hoc remove individual non-DRACH peaks.

**Approach:** Run HOMER `findMotifsGenome.pl` with `-rna` mode on peak centres ±50 bp against a length-matched random-shuffled background; expect DRACH-like consensus in the top motifs with E-value < 1e-50 for a well-behaved MeRIP dataset.

```bash
findMotifsGenome.pl \
    exomepeak2_output/m6a_run1/peaks.bed \
    hg38 \
    motif_output \
    -rna \
    -size 100 \
    -len 5,6 \
    -p 8
```

```r
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggseqlogo)
library(rtracklayer)

peaks <- import('exomepeak2_output/m6a_run1/peaks.bed')

peak_centres <- resize(peaks, width=5, fix='center')
genome <- BSgenome.Hsapiens.UCSC.hg38
seqs <- as.character(getSeq(genome, peak_centres))

ggseqlogo(seqs, method='probability') +
    ggplot2::labs(title='Peak-centre 5-mer (DRACH consensus expected)')
```

If DRACH enrichment is NOT observed in the peak set, the IP failed OR the wrong antibody was used OR the assay actually captured a different modification (m6A is centred in coding regions / 3'UTR; m1A is centred at TSS — different metagene signature). Do NOT proceed to differential analysis until DRACH is confirmed on the peak set as a whole.

## Multi-Tool Consensus

**Goal:** Build a high-confidence peak set by intersecting multiple peak callers run on the same data; report both per-tool peak counts and the intersection.

**Approach:** Use `bedtools intersect` to combine exomePeak2, MeTPeak, and MACS3 broadPeak outputs; require at least 2-of-3 caller agreement for the high-confidence set.

```bash
bedtools intersect \
    -a exomepeak2_output/m6a_run1/peaks.bed \
    -b metpeak_output/m6a_run1/peak.bed macs3_output/m6a_run1_peaks.broadPeak \
    -wa \
    -u \
    -f 0.5 > consensus_at_least_2of3.bed

wc -l \
    exomepeak2_output/m6a_run1/peaks.bed \
    metpeak_output/m6a_run1/peak.bed \
    macs3_output/m6a_run1_peaks.broadPeak \
    consensus_at_least_2of3.bed
```

Intersection (consensus) is the conservative choice; union inflates the false-positive rate (every caller's idiosyncratic false positives propagate). For published analyses, report both per-tool count AND consensus count; the ratio is a useful caller-agreement metric.

## Per-Method Failure Modes

### Anti-m6A antibody cross-reacts with m6Am at 5' peaks

**Trigger:** Peak called at the very 5' end of a transcript (within ~50 nt of TSS), attributed to METTL3-deposited internal m6A.

**Mechanism:** Anti-m6A antibodies cross-react with m6Am (the cap-adjacent N6,2'-O-dimethyladenosine at the +1 nucleotide of capped mRNAs), installed by PCIF1 / CAPAM (Akichika 2019 *Science* 363:eaav0080). METTL3 / METTL14 do NOT methylate the cap-adjacent position. Linder 2015 *Nat Methods* 12:767 explicitly noted the cross-reactivity; Mauer 2017 *Nature* 541:371 deepened the m6Am story.

**Symptom:** 5' peaks remain in METTL3-KO cells; agent infers METTL3-independent m6A writers when the peaks are really PCIF1 m6Am.

**Fix:** Flag peaks within ~50 nt of TSS as m6A-or-m6Am ambiguous. For unambiguous internal-m6A studies, exclude TSS-proximal peaks; for m6Am-specific studies, use PCIF1-KO subtraction or m6Am-seq (Sun H 2021); for orthogonal stoichiometric resolution, use GLORI / SAC-seq.

### Peaks called at GC-extreme transcripts by MACS3 but not exomePeak2

**Trigger:** MACS3 broadPeak output has peaks at high-GC or low-GC transcripts that exomePeak2 does NOT call on the same BAMs.

**Mechanism:** exomePeak2 implements internal GC-bias correction; MACS3 does not. The IP step has GC bias because anti-m6A antibody pull-down efficiency varies with local GC content. MACS3 attributes GC-driven IP enrichment to true methylation; exomePeak2 corrects for it.

**Fix:** Trust exomePeak2 for transcript-aware analyses. Use MACS3 only as a second opinion or for broad / viral analyses. If a MACS3-only peak is biologically interesting, validate orthogonally.

### DRACH filtering applied post-hoc

**Trigger:** Filtering called MeRIP peaks to retain only those containing a DRACH motif within the peak window.

**Mechanism:** ~70% of mammalian METTL3-deposited m6A sites sit within DRACH (Linder 2015), but a non-trivial fraction are non-DRACH. MeRIP peaks span ~100-200 nt; the methylated A is somewhere within the window but the peak boundary is not the modification position. Filtering by "no DRACH in window" rejects real peaks where DRACH sits at the edge or where the m6A is non-DRACH. The filter also creates circular validation — the filtered set is enriched for DRACH by construction.

**Symptom:** Peak count drops 30-50% after DRACH filtering; reviewers note the filter is not a standard convention.

**Fix:** Report DRACH enrichment as a SANITY CHECK on the peak set (HOMER E-value < 1e-50 expected) — never as a per-peak filter. For single-nucleotide methods (miCLIP, m6Anet, GLORI), the DRACH constraint is built into the calling model and need not be re-applied.

### MACS3 default `--keep-dup 1` on MeRIP

**Trigger:** MACS3 invoked without `--keep-dup all`; default behaviour is to retain one read per position.

**Mechanism:** MeRIP libraries have NO UMI (standard protocol); positional duplicates in MeRIP are a mix of PCR duplicates and biological re-sampling of high-coverage transcripts. Default `--keep-dup 1` collapses the latter and destroys real signal at the most-abundant transcripts.

**Symptom:** MACS3 reports very few peaks at housekeeping genes (GAPDH, ACTB) despite obvious IP enrichment in IGV; peak count drops 5-20x compared to `--keep-dup all`.

**Fix:** Always pass `--keep-dup all` to MACS3 / MACS2 for MeRIP. This is the same logic as the do-NOT-dedup rule in merip-preprocessing.

### MeTPeak called with wrong argument name for annotation

**Trigger:** `metpeak(IP_BAM=..., INPUT_BAM=..., txdb=txdb_object)` (lowercase `txdb=` rather than `TXDB=`).

**Mechanism:** MeTPeak's annotation arguments are `GENE_ANNO_GTF=` (GTF file path) or `TXDB=` (uppercase, TxDb object). Lowercase `txdb=` is not recognised. exomePeak2 uses lowercase `txdb=` instead, which leads to easy confusion between the two skills' APIs.

**Fix:** Pass either `GENE_ANNO_GTF='annotation.gtf'` (string path) or `TXDB=txdb_object` (capitalised). Verify argument casing against the installed MeTPeak source.

### Single-tool peak set reported without cross-caller agreement

**Trigger:** Paper reports "18,234 m6A peaks" from exomePeak2 alone; downstream biology built on this single set.

**Mechanism:** Cross-tool concordance on MeRIP is empirically ~70% between exomePeak2 and MeTPeak on the same BAMs; single-tool reports include ~30% tool-specific calls. McIntyre 2020 *Sci Rep* 10:6590 documents this directly.

**Fix:** Run at least 2 callers (typically exomePeak2 + MeTPeak; or exomePeak2 + MACS3 broad); report the intersection as high-confidence; full per-tool counts as supplementary.

### Transcriptome alignment passed to peak caller

**Trigger:** A transcriptome BAM (reads aligned to a transcripts.fa file) passed to exomePeak2 or MACS3 instead of a genome BAM.

**Mechanism:** Peak callers operate on genome BAMs and use the GTF to project peaks back to transcript features. Transcriptome BAMs have reads in per-transcript coordinates; the GTF cannot resolve these back to genome coordinates without re-alignment.

**Fix:** Align to GENOME with STAR / HISAT2 for downstream MeRIP peak calling (see merip-preprocessing). Transcriptome BAMs are for m6anet-analysis only.

## Reconciliation: When Peak Callers Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| exomePeak2 calls a peak; MACS3 does not | GC-confounded; MACS3 calibration ignores GC | Trust exomePeak2 for transcript-aware analysis |
| MACS3 calls a broad peak; exomePeak2 calls 3-5 sub-peaks | Broad biology (e.g., long 3'UTR); exomePeak2 splits broad enrichment into per-window peaks | Merge exomePeak2 sub-peaks within 200 nt; report as broad |
| miCLIP single-nt site OUTSIDE MeRIP peak window | MeRIP window slightly shifted; OR low-coverage transcript | Inspect MeRIP coverage; consider miCLIP authoritative for single-nt identity |
| m6Anet high-confidence site at non-MeRIP-peak location | Non-DRACH site invisible to MeRIP antibody preferences; OR low MeRIP coverage; OR site is real but stoichiometry too low for MeRIP detection | Trust m6Anet for DRACH-context sites at adequate coverage |
| Two replicates concordant; third diverges | Failed IP in third replicate | Inspect IP enrichment QC in merip-preprocessing for the diverging replicate; consider exclusion |
| GLORI calls a high-stoichiometry site that MeRIP misses | MeRIP under-calls at low-expression transcripts; OR site outside common-core consensus | Trust GLORI for absolute stoichiometry; MeRIP is qualitative |
| exomePeak2 gives different peak counts across Bioconductor releases on same BAMs | Default-parameter or signature shifts between minor Bioconductor releases | Verify against `?exomePeak2`; pin Bioconductor + exomePeak2 version in reproducible analyses |
| MACS3 narrowPeak vs broadPeak give very different counts | `--broad` extends peaks; narrow mode requires strong per-window enrichment | For MeRIP, use `--broad`; broadPeak is the appropriate output format |
| 5' peaks dominate the peak set | High m6Am signal (PCIF1) OR TSS-proximal antibody binding artifact | Flag 5' peaks; restrict downstream to internal peaks for METTL3 biology |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| exomePeak2 default FDR | 0.05 | exomePeak2 default; standard in field |
| exomePeak2 fold-change minimum | log2(IP/input) > 0 (any enrichment) by default; stringent set: log2 > 1 | exomePeak2 default vs commonly-applied stringent threshold |
| MACS3 broad-cutoff for m6A | 0.1 | MACS2/3 default for `--broad` mode |
| MACS3 `--extsize` for MeRIP | 150 | MeRIP fragment length convention |
| MACS3 `--qvalue` | 0.05 | Convention (applies to narrow-peak output; broad-peak file uses `--broad-cutoff` instead) |
| DRACH enrichment E-value (HOMER) | < 1e-50 | Convention for "well-validated antibody dataset"; below this threshold the IP likely failed |
| Minimum coverage per peak window | 30 reads in BOTH IP and Input | Standard convention; below this, statistical calls are noisy |
| Peak reproducibility threshold | >=2 of N replicates | Stringency-vs-recall trade-off; report multiple thresholds |
| 5'UTR peak ambiguity zone | first 50 nt of transcript | Conservative; antibody cross-reactivity with m6Am peaks here |
| Minimum peak length | >=20 nt | Below this, suspect noise or peak-edge artefact |
| MeTPeak default p-value cutoff | 1e-5 | MeTPeak default |
| MeTPeak default FDR cutoff | 0.05 | MeTPeak default |
| MeTPeak `WINDOW_WIDTH` / `SLIDING_STEP` | 50 / 50 (defaults) | MeTPeak GitHub source; do not assume smaller values without rationale |
| MeTPeak `FRAGMENT_LENGTH` | 100 (default) | MeTPeak GitHub source; matches non-stranded MeRIP convention |
| Common-core m6A sites (cross-method) | ~6,000-15,000 in HEK293T | The intersection of MeRIP + miCLIP + GLORI + SAC-seq; peaks outside the common core need orthogonal validation |
| Cross-replicate peak overlap | ~80% within lab; ~30-60% between labs (median ~45%) | McIntyre 2020 *Sci Rep* 10:6590; bounds reproducibility expectations |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| exomePeak2 errors on TxDb chromosome mismatch | BAM uses `chr1`; GTF uses `1` (or reverse) | Verify with `samtools view -H` and `head genes.gtf`; reconcile with `seqlevelsStyle()` in R |
| exomePeak2 reports zero peaks | TxDb built from incorrect GTF; OR `paired_end` flag wrong; OR all peaks below FDR threshold | Verify TxDb covers same chromosomes as BAM; verify `paired_end=TRUE` for PE; relax `pvalue_cutoff` |
| MACS3 reports zero peaks | Default `--keep-dup 1` collapsed signal; OR `--nomodel` not set; OR `--gsize` wrong | Add `--keep-dup all --nomodel`; confirm `--gsize hs` (2.7e9) or `mm` |
| MeTPeak install fails | GitHub-only, requires devtools | `devtools::install_github('compgenomics/MeTPeak')` |
| MeTPeak `metpeak()` error on lowercase `txdb=` | MeTPeak uses uppercase `TXDB=` (TxDb object) or `GENE_ANNO_GTF=` (file path); exomePeak2 uses lowercase `txdb=` — easy to confuse | Use `GENE_ANNO_GTF='annotation.gtf'` or `TXDB=txdb_object` |
| HOMER DRACH motif not detected | Antibody failure OR wrong protocol OR insufficient peaks for motif detection | Re-inspect IP/input fingerprint in merip-preprocessing; verify peak count >>100 |
| Peak file is empty BED | exomePeak2 silently filtered all peaks | Check `pvalue_cutoff` / `fold_enrichment` arguments; lower thresholds |
| Cross-tool peak intersect tiny | Tool-default thresholds differ; OR fragment length parameters differ | Harmonise thresholds across callers; reconcile via IDR-equivalent |
| exomePeak2 takes >24h on whole-genome BAMs | Large BAM + many transcripts | Subset to expressed transcripts; OR use more cores via downstream parallelisation |
| `findMotifsGenome.pl` error on chromosome naming | hg38 vs HG38 vs Hg38 | Use lowercase `hg38` consistently |
| ggseqlogo throws "sequences must be equal length" | Mixed-length sequences passed | Resize peak ranges to fixed width: `resize(peaks, width=5, fix='center')` |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "How many peak callers were used?" | Two minimum (exomePeak2 + MeTPeak or + MACS3); intersection reported as high-confidence; per-tool counts as supplementary |
| "Was filtering by DRACH applied?" | No — DRACH reported as enrichment-on-the-peak-set (HOMER E-value); per-peak DRACH filtering drops 5-10% real m6A sites |
| "How were 5'UTR peaks handled?" | Flagged peaks within 50 nt of TSS as m6A-or-m6Am ambiguous; restricted internal-m6A analyses to peaks past the 5'UTR; cited Linder 2015 / Mauer 2017 cross-reactivity |
| "What's the FDR threshold?" | exomePeak2 default FDR 0.05; MACS3 q-value 0.05; peaks reported with both fold-change and FDR |
| "Was orthogonal validation done?" | High-stakes sites validated against published miCLIP / GLORI / SAC-seq / m6A-Atlas; cross-method overlap reported |
| "Why exomePeak2 over MeTPeak?" | exomePeak2 implements GC-bias correction; MeTPeak does not. For low-coverage datasets MeTPeak HMM smoothing helps; for typical datasets exomePeak2 is the field default |
| "Were failed IPs checked for?" | Replicate IPs inspected via plotFingerprint AND per-transcript IP/input ratio distribution in merip-preprocessing BEFORE peak calling |
| "What's the cross-replicate peak overlap?" | Reported per pair; expect ~80% within-lab per McIntyre 2020 |
| "Was the peak set intersected with m6A-Atlas?" | Yes — common-core overlap reported as a confidence anchor; novel peaks flagged for orthogonal validation |
| "Why weren't m6A-CLIP peaks called here?" | miCLIP / m6A-CLIP single-nucleotide methods live in `clip-seq/peak-calling`; this skill is for fragment-level MeRIP peak calling |

## References

- Dominissini D, Moshitch-Moshkovitz S, Schwartz S et al (2012) Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. *Nature* 485(7397):201-206. doi:10.1038/nature11112
- Meyer KD, Saletore Y, Zumbo P, Elemento O, Mason CE, Jaffrey SR (2012) Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. *Cell* 149(7):1635-1646. doi:10.1016/j.cell.2012.05.003
- Linder B, Grozhik AV, Olarerin-George AO, Meydan C, Mason CE, Jaffrey SR (2015) Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. *Nat Methods* 12(8):767-772. doi:10.1038/nmeth.3453
- Mauer J, Luo X, Blanjoie A et al (2017) Reversible methylation of m6Am in the 5' cap controls mRNA stability. *Nature* 541(7637):371-375. doi:10.1038/nature21022
- Liu J, Yue Y, Han D et al (2014) A METTL3-METTL14 complex mediates mammalian nuclear RNA N6-adenosine methylation. *Nat Chem Biol* 10(2):93-95. doi:10.1038/nchembio.1432
- Cui X, Meng J, Zhang S, Chen Y, Huang Y (2016) A novel algorithm for calling mRNA m6A peaks by modeling biological variances in MeRIP-seq data. *Bioinformatics* 32(12):i378-i385. doi:10.1093/bioinformatics/btw281
- Meng J, Lu Z, Liu H, Zhang L, Zhang S, Chen Y, Rao MK, Huang Y (2014) A protocol for RNA methylation differential analysis with MeRIP-Seq data and exomePeak R/Bioconductor package. *Methods* 69(3):274-281. doi:10.1016/j.ymeth.2014.06.008
- Liu J, Zhang Z, Meng J et al (2022) exomePeak2: a peak calling and differential analysis tool for MeRIP-Seq with bias awareness. *NAR Genom Bioinform* 4(3):lqac046. doi:10.1093/nargab/lqac046
- Zhang Y, Liu T, Meyer CA et al (2008) Model-based analysis of ChIP-Seq (MACS). *Genome Biol* 9(9):R137. doi:10.1186/gb-2008-9-9-r137
- Akichika S, Hirano S, Shichino Y et al (2019) Cap-specific terminal N6-methylation of RNA by an RNA polymerase II-associated methyltransferase. *Science* 363(6423):eaav0080. doi:10.1126/science.aav0080
- Boulias K, Toczydłowska-Socha D, Hawley BR et al (2019) Identification of the m6Am methyltransferase PCIF1 reveals the location and functions of m6Am in the transcriptome. *Mol Cell* 75(3):631-643.e8. doi:10.1016/j.molcel.2019.06.006
- Liu C, Sun H, Yi Y et al (2023) Absolute quantification of single-base m6A methylation in the mammalian transcriptome using GLORI. *Nat Biotechnol* 41(3):355-366. doi:10.1038/s41587-022-01487-9
- Hu L, Liu S, Peng Y et al (2022) m6A RNA modifications are measured at single-base resolution across the mammalian transcriptome. *Nat Biotechnol* 40(8):1210-1219. doi:10.1038/s41587-022-01243-z
- Garcia-Campos MA, Edelheit S, Toth U et al (2019) Deciphering the m6A code via antibody-independent quantitative profiling. *Cell* 178(3):731-747.e16. doi:10.1016/j.cell.2019.06.013
- McIntyre ABR, Gokhale NS, Cerchietti L, Jaffrey SR, Horner SM, Mason CE (2020) Limits in the detection of m6A changes using MeRIP/m6A-seq. *Sci Rep* 10(1):6590. doi:10.1038/s41598-020-63355-3
- Heinz S, Benner C, Spann N et al (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. *Mol Cell* 38(4):576-589. doi:10.1016/j.molcel.2010.05.004

## Related Skills

- merip-preprocessing - IP/input BAM preparation, library complexity QC, and IP enrichment QC upstream of peak calling
- m6a-differential - Compare peak sets between conditions; uses peaks called here as input
- m6anet-analysis - Orthogonal validation of MeRIP peaks via ONT direct-RNA single-nucleotide resolution
- modification-visualization - Metagene, browser-track, and peak-centred heatmap rendering of called peaks
- clip-seq/peak-calling - miCLIP / m6A-CLIP single-nucleotide validation methods (PureCLIP, PEKA, paraclu)
- clip-seq/clip-motif-analysis - Antibody-CLIP motif analysis context for cross-validation
- chip-seq/peak-calling - General sliding-window IP-vs-input peak calling (MACS3 design lineage)
- chip-seq/peak-annotation - Annotation of peaks against gene features (re-usable for m6A)
- read-alignment/star-alignment - Splice-aware STAR alignment defaults referenced by merip-preprocessing
- rna-quantification/featurecounts-counting - Peak count matrix construction (input to m6a-differential)
- pathway-analysis/go-enrichment - GO enrichment on gene lists derived from peak-bearing transcripts
- variant-calling/vcf-basics - Cross-reference m6A peaks against A-to-I edit sites (sometimes confounded)
- data-visualization/multipanel-figures - Combining metagene + heatmap + volcano for figures
