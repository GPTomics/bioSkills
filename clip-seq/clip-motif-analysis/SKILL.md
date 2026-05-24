---
name: bio-clip-seq-clip-motif-analysis
description: Discover RBP binding motifs from CLIP-seq peaks or single-nucleotide crosslink sites using HOMER, MEME/STREME, kpLogo, mCross (CL-position-registered motifs), PEKA (positional k-mer enrichment), RBPamp (affinity), and RNA Bind-n-Seq (RBNS) cross-validation. Use when characterizing RBP sequence specificity, registering motifs to crosslink positions, validating in vivo CLIP motifs against in vitro RBNS Kd, reconciling motif disagreements across tools, or correcting for the uracil crosslinking bias that contaminates raw CLIP motif logos.
tool_type: mixed
primary_tool: HOMER
---

## Version Compatibility

Reference examples tested with: HOMER 4.11+, MEME Suite 5.5+ (STREME, MEME-ChIP, FIMO), bedtools 2.31+, kpLogo 1.1+, mCross v1+, PEKA v1+, RBPamp 0.9+, ggseqlogo 0.1+, biopython 1.83+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws unexpected errors, introspect the installed tool and adapt the example to match the actual API rather than retrying.

# CLIP-seq Motif Analysis

**"Find enriched RNA motifs at my RBP binding sites"** -> Discover the in vivo sequence preference of an RNA-binding protein from CLIP-seq peaks or single-nucleotide crosslink sites. The fundamental confound is the uracil bias of UV254 crosslinking: U is the most-crosslinked base (>50% of CL events), so naive motif logos centered on CL positions are U-enriched even for non-U-binding RBPs. Modern tools (mCross, PEKA) register motifs relative to the CL position and correct for this bias; legacy tools (HOMER, MEME) need careful background selection.

- CLI (de novo, peak-based, HOMER RNA mode): `findMotifs.pl peaks.fa fasta motif_out -rna -len 5,6,7,8 -p 4`
- CLI (de novo, peak-based, MEME-ChIP / STREME): `streme --rna --oc streme_out -p peaks.fa -n background.fa --minw 5 --maxw 10`
- CLI (positional, single-nt CL-registered, mCross): `mCross -i crosslinks.bed -g genome.fa -k 7 -o mcross_out` (jointly models motif + CL position)
- CLI (positional k-mer, no input control needed, PEKA): `peka --peak_file_name peaks.bed --crosslinks_file_name crosslinks.bed --genome_file_name genome.fa --regions_file_name regions.bed --kmer_length 5 --percentile 30 --outpath peka_out` (the short flags `-i/-x/-g/-r/-k/-p` shown in earlier docs are not all stable; use the long forms or check `peka --help`)
- CLI (affinity-weighted, RBPamp): `rbpamp run -i peaks.fa -k 7 -o rbpamp_out` (joint affinity + motif model)
- CLI (positional logo, kpLogo): `kpLogo crosslinks_with_kmer_scores.txt -o kplogo_out` (position-specific significance)

ENCODE eCLIP motif convention: extract sequences from peaks (CLIPper + SMInput stringent), use shuffled or expression-matched background, run HOMER `-rna -len 5,6,7,8`. The Yeo lab eCLIP papers typically report the top 1-3 HOMER hits per RBP and validate with mCross for CL-registered position. For ASB and variant-effect work, mCross is mandatory.

## Algorithmic Taxonomy

| Tool | Input | Background | CL-position aware | Output | Strength | Fails when |
|------|-------|------------|-------------------|--------|----------|------------|
| HOMER findMotifs.pl | Peak FASTA | Auto-shuffled or user-provided | No | PWM + de novo + known scan | Mature, fast, RNA mode (U instead of T) | Background generation is heuristic; uracil bias inflates U-content motifs |
| MEME (MEME Suite) | Peak FASTA | Optional Markov model | No | PWM | Statistically rigorous; broadly cited | Slow on large peak sets (> 1000 sequences); single-motif default |
| STREME (MEME Suite) | Peak FASTA | Shuffled or user-provided | No | PWM | Fast successor to MEME-ChIP; tolerant of large sets | Same uracil bias issue as HOMER |
| MEME-ChIP | Peak FASTA | Markov order 1-3 | No | Multiple motifs | Wrapper for MEME + CentriMo + Tomtom comparison | Designed for ChIP; for CLIP use STREME directly |
| kpLogo | k-mer scores with positions | NA (significance from kpLogo internals) | Yes (positional) | Position-specific logo | Visualizes position-specific signal; complements PWM | Requires k-mer-score input file; not a peak-to-motif tool by itself |
| mCross (Feng 2019) | Single-nucleotide CL sites + genome | Generated internally | Yes (jointly models motif + CL position) | Registered PWM + CL offset | The standard for in vivo RBP motif registration; SRSF1 example showed clustered GGA half-sites | Requires high-quality single-nt CL sites (PureCLIP or CTK CITS) |
| PEKA (Mihai 2022) | Peaks + crosslink BED | Low-count crosslinks within same dataset | Yes (positional k-mer) | Position-enriched k-mers with significance | No external input required; cross-validates with mCross | Less granular than mCross's jointly-modeled PWM |
| RBPamp (Jens 2020) | Peak FASTA or CL sites | Joint affinity model | No (affinity-based) | Affinity-weighted motif + Kd estimate | Provides Kd; compatible with RBNS calibration | Slow; small community; needs validation |
| RCAS | Peak BED + GTF | Auto | No | Motif + region annotation | One-stop for motif + region distribution | Less granular than per-tool dedicated runs |
| MEME FIMO (motif scan) | PWM + FASTA | Markov bg | No | Match locations and scores | Scan known motifs across genome | Not a de novo tool; complements HOMER known scan |
| RBPNet (Jens & Gagneur 2024) | Raw sequence | Trained model | Yes (sequence-to-signal at single nt) | Predicted CL count distribution | Deep learning; predicts per-base affinity profile | See clip-seq/clip-deep-learning |

Methodology evolves; verify motif comparison databases (CISBP-RNA, ATtRACT, oRNAment, mCrossBase, RBPDB) for known motif validation. RBNS in vitro Kd values (Lambert 2014; Dominguez 2018 for 78 RBPs) are the most reliable orthogonal reference for in vivo CLIP motifs.

## Critical Choice: Peak-Based vs Crosslink-Site-Based Motif Discovery

**Peak-based (HOMER, MEME/STREME, RBPamp):** Extract sequences from peak BED, find enriched k-mers / PWMs. Produces classical motif logos. Fast and intuitive. Captures broad binding zones.

**Crosslink-site-based (mCross, PEKA, kpLogo):** Use single-nucleotide CL positions; analyze sequence in a window around each CL position. Produces motif logos REGISTERED to the CL site. Required for understanding the structural relationship between motif and crosslink.

The fundamental difference: peak-based logos show "what sequence is enriched near binding"; CL-registered logos show "what sequence the RBP contacts at the crosslink position." For most RBPs they agree on the core motif but diverge on flanking context.

## Uracil Crosslinking Bias

UV254 crosslinking is U-biased: ~60-80% of single-nt CL events occur at U residues (Konig 2010; Sugimoto 2017). Consequence: motif logos centered on raw CL positions are U-enriched at the center even for non-U-binding RBPs (e.g., PUM2 binds UGUANAUA but peak-center is shifted; HNRNPK binds C-rich tracts but centers can show spurious U).

**Mitigations:**
- Use peak boundaries, not single-nt CL, for HOMER/MEME (peak-level averaging dilutes the U bias).
- Use mCross or PEKA, which explicitly model the CL position offset from the motif center.
- Shuffle background should preserve mononucleotide composition (HOMER default does; MEME with `--markov-order 1` does).
- For PAR-CLIP, the bias is at U residues but T->C marks the EXACT CL position; mCross / PEKA register more cleanly than for iCLIP/eCLIP.

For naive HOMER/MEME output, examine the central column of the PWM: if it is U-skewed, suspect CL bias rather than true U preference. Cross-check against published RBNS Kd if available.

## Per-Tool Failure Modes

### HOMER -- Background mismatch inflates GC-biased motifs

**Trigger:** HOMER auto-generates a shuffled background that preserves only nucleotide frequency, not GC distribution per region; CLIP peaks come disproportionately from 3' UTRs which are AU-rich.

**Mechanism:** AU-rich foreground vs GC-shuffled background produces "AU-rich motif" as the top hit for any 3' UTR-binding RBP, regardless of biology.

**Symptom:** Top HOMER motif looks like the 3' UTR mean composition (AU-rich, no clear specificity); known motif scan (`-known`) finds the RBP's canonical motif at lower rank than the spurious "AU" motif.

**Fix:** Provide a GC-matched background: extract sequences from random 3' UTR regions of expressed transcripts matched by length distribution. Use `-bg matched_bg.fa` flag. Or use STREME with `-n matched_bg.fa`.

### MEME -- Slow on large peak sets

**Trigger:** MEME run on > 1000 peak sequences.

**Mechanism:** MEME is O(N^2) in sequence count; classical OOPS/ZOOPS/ANR models do not scale.

**Symptom:** Runtime > 12h; out-of-memory on large machines.

**Fix:** Use STREME (MEME Suite >= 5.0); it is the modern fast successor (sub-quadratic). Or restrict to top-1000 peaks by score.

### mCross -- Requires single-nt CL sites

**Trigger:** mCross run on a peak BED instead of crosslink BED.

**Mechanism:** mCross's model is `motif_PWM x CL_position_offset_PMF`. Without single-nt CL positions there is no offset axis; the model degenerates.

**Symptom:** mCross output looks identical to a peak-centered logo; the diagnostic CL-offset histogram is flat.

**Fix:** Provide CL sites from PureCLIP, CTK CITS, or PARalyzer. The Yeo lab and Zhang lab use the convention of PureCLIP CL sites as mCross input.

### PEKA -- Background from same dataset

**Trigger:** Concerned that PEKA's "low-count crosslinks within same dataset" background is biased.

**Mechanism:** PEKA does NOT need external input. It uses low-count crosslinks (below a quantile threshold) as the background and high-count crosslinks (above) as the foreground. This is by design (Mihai 2022) and cross-validated against mCross.

**Symptom:** False positive concern; user wants to verify.

**Fix:** Run BOTH PEKA and mCross; cross-validate. PEKA paper shows >90% concordance with mCross for ENCODE RBPs; treat them as orthogonal confirmation.

### RBPamp -- Slow convergence

**Trigger:** Large peak set (> 5000); RBPamp jointly optimizes affinity + motif iteratively.

**Mechanism:** Joint optimization is slow; some RBPs have multimodal binding (multiple motifs) that confuse single-PWM RBPamp.

**Symptom:** Runtime > 24h; output PWM disagrees with HOMER top motif.

**Fix:** Down-sample to top-2000 peaks; or run HOMER for initial PWM and use it as RBPamp seed.

### RBNS comparison -- in vitro vs in vivo divergence

**Trigger:** Comparing CLIP-derived motif to RBNS Kd-ranked motif (Dominguez 2018) and finding divergence (e.g., CLIP top motif rank 3 in RBNS; RBNS top motif weakly enriched in CLIP).

**Mechanism:** In vivo binding is shaped by structure accessibility (RNA secondary structure), cooperativity, and competing factors. RBNS measures intrinsic affinity in vitro. They diverge for:
- Structurally regulated RBPs (PUM2, RBFOX2 favor accessible loops)
- Cooperative binders (FUS, TDP-43, SR proteins)
- Co-factor-dependent (SRSF1 paralogs)

**Symptom:** CLIP motif logo plausible but rank-3 in RBNS; or RBNS top motif missing from CLIP peaks.

**Fix:** Accept the divergence as informative. Cite both. Use SHAPE-eCLIP or icSHAPE-MaP to test structural-accessibility hypothesis. For variant-effect studies, prefer RBNS-derived PWM (RBPamp Kd) over CLIP motif for absolute affinity prediction.

### Known motif scan -- FIMO threshold too lenient

**Trigger:** Used FIMO with default `-thresh 1e-4` and got matches everywhere.

**Mechanism:** Default threshold is too lenient for in vivo binding; CLIP background sequences contain many low-affinity matches.

**Symptom:** FIMO match rate > 50% of peaks AND > 20% of random background; no enrichment signal.

**Fix:** Use `-thresh 1e-5` or `-thresh 1e-6`; or use CentriMo for local enrichment around peak center (testing positional enrichment, not threshold).

## Decision Tree by Scenario

| Scenario | Tool | Why |
|----------|------|-----|
| Quick de novo motif from CLIPper peaks, ENCODE-comparable | HOMER `-rna -len 5,6,7,8` with GC-matched bg | Fast, mature, ENCODE convention |
| Motif registered to crosslink position (for variant effect) | mCross on PureCLIP single-nt sites | The 2019 Feng paper is the canonical reference |
| Confirm mCross result with independent method | PEKA on same peaks + crosslinks | Cross-validation; PEKA does not need input control |
| Affinity-weighted motif with Kd estimate | RBPamp | Provides quantitative Kd; integrates with RBNS comparison |
| Compare to in vitro RBNS Kd | RBPamp + Dominguez 2018 / Lambert 2014 tables | RBNS is the orthogonal in vitro standard |
| RBP with known motif - validation scan | FIMO -thresh 1e-5 against CISBP-RNA / ATtRACT PWM | Lower-confidence motifs need tight threshold |
| Multiple binding modes (SRSF1 GGA clusters) | mCross + manual inspection of CL-offset histogram | mCross reveals the multimodal structure |
| Position-specific logo (5' vs 3' of CL) | kpLogo with k-mer scores | Position-specific significance visualization |
| AU-rich 3' UTR binders (HuR, AUF1) - validate U bias is real | Compare HOMER motif center vs PUM2-style register | RBP-specific tradition |
| Allele-specific motif (variant effect) | mCross PWM + DeepRiPe/RBPNet for variant scoring | See clip-seq/clip-deep-learning |
| snoRNA / structural RNA - sequence motif less meaningful | Skip motif analysis; emphasize structural context (icSHAPE) | snoRNA binders read structure, not linear sequence |

## Reconciliation: When Motif Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| HOMER and STREME agree; mCross disagrees on flanking | mCross's CL-offset registers a flanking base that's identical in peak; HOMER sees both flanking bases as motif | Trust mCross for CL-registered position; HOMER is "what's nearby" |
| Top HOMER motif AU-rich; RBP knows to bind structure | GC-matched background not used | Re-run HOMER with `-bg matched_bg.fa` |
| HOMER finds canonical motif at rank 1; RBNS Kd ranks at rank 3 | In vitro vs in vivo specificity divergence | Both valid; cite RBNS for absolute affinity, CLIP for in vivo prevalence |
| mCross PWM has multiple CL-position peaks | Multiple binding modes (SRSF1 GGA half-sites; PTBP1 multimer) | Treat as biologically meaningful; report both |
| Skipper window-level motifs differ from CLIPper peak-level | Window size dilutes signal; peak boundaries sharpen | Both valid; report at the level of the downstream analysis |
| PAR-CLIP motifs cleaner than iCLIP for same RBP | PAR-CLIP T->C registers exactly; iCLIP truncates -1 from CL | Both correct; PAR-CLIP is the sharper single-nt method when 4SU labeling tolerated |
| HOMER finds different motif on top vs bottom strand | Strand information lost in upstream BED | Verify BED column 6 (strand) is preserved; HOMER respects strand |

**Operational rule for high-confidence motif reporting:** Report (a) HOMER top de novo motif with GC-matched background, (b) mCross PWM with CL-position offset histogram, (c) PEKA k-mer enrichment for orthogonal confirmation, (d) FIMO scan against published RBP PWM (CISBP-RNA / ATtRACT) for known-motif validation, (e) RBNS Kd ranking from Dominguez 2018 if available. Three independent methods agreeing on the core motif is the bar.

## Workflow: De Novo + Known + CL-Registered

**Goal:** Produce a publication-grade motif report combining peak-based de novo discovery, CL-position-registered motif (mCross), and known-motif validation.

**Approach:** Extract sequences with strand preservation, build a GC-matched 3' UTR background, run HOMER + STREME for de novo, mCross for CL-registered, FIMO scan against CISBP-RNA, and report all three views with information-content QC.

```bash
# Step 1: Extract peak sequences (use stringent peak set: log2 FC >= 3, -log10 p >= 3)
bedtools getfasta -fi genome.fa -bed peaks.stringent.bed -s -fo peaks.fa

# Step 2: GC-matched background (random regions from expressed transcripts)
# expressed.bed = transcripts with TPM >= 1 in the same cell type
shuffleBed -i peaks.stringent.bed -g chrom.sizes -incl expressed.bed -seed 42 > shuffled.bed
bedtools getfasta -fi genome.fa -bed shuffled.bed -s -fo background.fa

# Step 3: HOMER de novo + known
findMotifs.pl peaks.fa fasta homer_out \
    -rna -len 5,6,7,8 -p 8 \
    -fasta background.fa

# Step 4: STREME for cross-validation
streme --rna --oc streme_out -p peaks.fa -n background.fa --minw 5 --maxw 10

# Step 5: mCross requires single-nt CL sites (PureCLIP output)
# crosslinks.bed = PureCLIP -o sites.bed (single-nt CL positions, see clip-seq/crosslink-site-detection)
mCross \
    -i crosslinks.bed -g genome.fa -k 7 -n 5 -o mcross_out

# Step 6: PEKA orthogonal positional k-mer
peka --peak_file_name peaks.stringent.bed --crosslinks_file_name crosslinks.bed --genome_file_name genome.fa --regions_file_name regions.bed --kmer_length 5 --percentile 30 --outpath peka_out

# Step 7: Compare to RBNS in vitro Kd (Dominguez 2018 supplementary)
# Manual cross-reference to RBNS Kd-ranked top-5 motifs for the target RBP
```

## Quality Checks

```python
import numpy as np
import pandas as pd
from Bio import motifs

def parse_homer_motif(motif_file):
    '''Parse HOMER motif PWM and return information content'''
    pwm = []
    with open(motif_file) as f:
        for line in f:
            if line.startswith('>'):
                continue
            pwm.append([float(x) for x in line.strip().split()])
    pwm = np.array(pwm)
    # Information content per position (bits)
    epsilon = 1e-10
    ic = np.sum(pwm * np.log2(pwm / 0.25 + epsilon), axis=1)
    return {
        'length': pwm.shape[0],
        'mean_IC_per_position': np.mean(ic),
        'max_IC': np.max(ic),
        'min_IC': np.min(ic),
        'is_low_complexity': np.mean(ic) < 0.5
    }

# A high-quality CLIP motif has mean IC per position > 1.0 bit
# Mean IC < 0.5 suggests no real motif - check for background mismatch
```

| Metric | Good | Investigate |
|--------|------|-------------|
| Mean IC per position (HOMER PWM) | > 1.0 bit | < 0.5 = no real motif or background mismatch |
| Motif center U content | matches RBP biology | > 80% U in central position = U-CL bias suspect |
| FIMO scan rate in peaks vs background | peaks 3-10x bg | < 2x = motif too weak; > 50x = motif too narrow |
| mCross CL-offset histogram | single mode at offset N | flat = mCross unable to register |
| RBNS rank of CLIP top motif | top 3 | > 10 = strong in vivo / in vitro divergence; investigate |
| Replicate motif overlap | top motif identical in 2/2 reps | discordance = under-powered or one rep failed |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| HOMER top motif is "ATCC" or similar TruSeq-adapter-like | Adapter trim incomplete in preprocessing | Re-run cutadapt with correct adapter; verify trimming completed before peak calling |
| mCross "no crosslinks found" or output identical to peak-centered logo | Passed peak BED instead of single-nt crosslink-site BED | mCross requires PureCLIP or CTK CITS single-nt site BED; not a CLIPper peak BED. See clip-seq/crosslink-site-detection |
| mCross CL-offset histogram is flat | Strand information missing OR alignment soft-clipped 5' base | Verify BED column 6 (strand); confirm upstream `--alignEndsType EndToEnd` |
| All HOMER motifs are AU-rich, no specificity | Default shuffled background AU-rich for 3' UTR peaks | Provide GC-matched background from expressed 3' UTRs |
| mCross input "no crosslinks found" | Provided peak BED instead of single-nt CL BED | Pass PureCLIP / CTK CITS single-nt BED instead |
| MEME "out of memory" or > 12h runtime | Peak set too large | Switch to STREME; or down-sample to top-1000 peaks |
| FIMO match in nearly every peak | Threshold too lenient | Tighten to 1e-5 or 1e-6 |
| Motif width > 10 nt looks chimeric | HOMER reported overlapping co-motifs as one | Re-run with `-len 5,6,7,8` (force narrower); inspect rank 2-5 motifs |
| PWM info content < 0.5 bits | Background mismatch or peaks have no specificity | Verify peak set is stringent (log2 FC >= 3); regenerate background |
| Strand-specific motif different on +/- | BED column 6 missing or wrong | Verify strand encoding; re-extract with `bedtools getfasta -s` |
| Multiple distinct motifs per RBP | Biological (multimodal RBPs) OR contaminating multi-RBP CLIP | Validate with mCross multimode + compare to literature |

## References

- Heinz S et al 2010 Mol Cell 38:576 (HOMER motif discovery)
- Bailey TL et al 2015 Nucleic Acids Res 43:W39 (MEME Suite)
- Bailey TL et al 2021 Bioinformatics 37:2834 (STREME, MEME's fast successor)
- Wu X & Bartel DP 2017 Nucleic Acids Res 45:7331 (kpLogo positional logo)
- Feng H et al 2019 Mol Cell 74:1189 (mCross, jointly modeling motif + CL position)
- Kuret K, Amalietti AG, Jones DM, Capitanchik C, Ule J 2022 Genome Biol 23:160 (PEKA, positional k-mer no input)
- Jens M et al 2018 Cell Reports 24:2940 (RBPamp, affinity-weighted motif)
- Lambert N et al 2014 Mol Cell 54:887 (RNA Bind-n-Seq)
- Dominguez D et al 2018 Mol Cell 70:854 (78-RBP RBNS atlas)
- Sugimoto Y et al 2015 Nature 519:491 (iCLIP U-crosslink bias)
- Van Nostrand EL et al 2020 Nature 583:711 (ENCODE 150 RBP eCLIP + motifs)
- Konig J et al 2010 Nat Struct Mol Biol 17:909 (iCLIP, U bias origin)

## Related Skills

- clip-seq/clip-peak-calling - Upstream peak calls feed motif discovery
- clip-seq/crosslink-site-detection - Single-nt CL sites required for mCross / PEKA
- clip-seq/binding-site-annotation - Region-level context for motifs
- clip-seq/clip-deep-learning - Sequence-to-binding deep models (RBPNet, DeepRiPe)
- clip-seq/ago-clip-mirna-targets - miRNA seed motifs from CLEAR-CLIP chimeras
- chip-seq/motif-analysis - DNA-protein motif analogue
- pathway-analysis/go-enrichment - Functional context of motif-bearing genes
