---
name: bio-clip-seq-clip-qc
description: Comprehensive quality control for CLIP-seq libraries (eCLIP, iCLIP, iCLIP2, PAR-CLIP) covering library complexity (preseq), FRiP, IDR replicate reproducibility, read-distribution metagene, SMInput vs IgG control rationale, rRNA / snoRNA contamination, fragment-length distribution, and ENCODE-compliance thresholds. Use when assessing whether a CLIP library passed, deciding lenient vs stringent peak thresholds, comparing replicates with IDR rescue and self-consistency ratios, or distinguishing failed IP from over-amplified library.
tool_type: mixed
primary_tool: preseq
---

## Version Compatibility

Reference examples tested with: preseq 3.2+, picard 3.1+, samtools 1.19+, bedtools 2.31+, deeptools 3.5+, idr 2.0.4+, MultiQC 1.21+, RSeQC 5.0+, pysam 0.22+, fastp 0.23+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed binary and adapt the example to match the actual CLI rather than retrying.

# CLIP-seq Quality Control

**"Did my CLIP library pass?"** -> Assess preprocessing retention, alignment rate, library complexity, replicate reproducibility (IDR), fraction reads in peaks (FRiP), read-distribution metagene, rRNA/snoRNA contamination, fragment-length distribution, and SMInput vs IP enrichment. ENCODE eCLIP compliance is the canonical bar: >= 1M unique fragments per replicate, IDR rescue and self-consistency ratios both < 2, FRiP >= 0.005 (narrow-binding), library complexity rising linearly with depth on preseq lc_extrap. A library can fail at any of these stages, and the failure mode determines whether the data is salvageable.

- CLI (library complexity, primary QC): `preseq lc_extrap -B -P aligned.bam -o complexity.txt`
- CLI (FRiP after peak calling): `bedtools intersect -c -s -a peaks.bed -b dedup.bam | awk '{s+=$NF} END{print s}'` then divide by total reads
- CLI (IDR, ENCODE convention): `idr --samples rep1.sorted.bed rep2.sorted.bed --input-file-type bed --rank 5 --output-file idr.out --idr-threshold 0.05 --plot`
- CLI (read distribution / metagene): `RSeQC read_distribution.py -i dedup.bam -r gencode.v38.bed` + `geneBody_coverage.py -i dedup.bam -r housekeeping.bed -o gb`
- CLI (rRNA contamination check): `samtools idxstats dedup.bam | awk '$1 ~ /rRNA|45S|18S|28S/ { sum+=$3 } END {print sum}'`
- CLI (consolidated report): `multiqc <run_dir>` aggregates FastQC + cutadapt + STAR + umi_tools + preseq + samtools stats

The ENCODE eCLIP standards (encodeproject.org/eclip) define: >= 2 biological replicates with >= 1M unique fragments each (or saturated peak detection); IDR rescue and self-consistency ratios both < 2; narrow-binding RBPs FRiP >= 0.005. CLIP libraries should have 40-70% PCR duplication BY DESIGN - the IP enriches a small molecule pool, so high pre-dedup duplication is normal; low duplication suggests failed IP.

## QC Stage Hierarchy

CLIP QC progresses through five gates; failure at an earlier gate makes later gates meaningless.

| Gate | Metric | Tool | ENCODE threshold | Failure interpretation |
|------|--------|------|------------------|------------------------|
| 1. Preprocessing retention | % reads retained after UMI + adapter trim | cutadapt log | >= 70% | Adapter pattern wrong; degraded RNA |
| 2. Alignment rate | % reads aligned to genome (unique) | STAR Log.final.out | >= 60% for eCLIP, 70% for iCLIP/PAR-CLIP | Wrong genome; rRNA pre-map missing |
| 3. Library complexity | Predicted unique fragments at sequenced depth | preseq lc_extrap | >= 1M unique | Over-amplified or under-input library |
| 4. IP enrichment | log2(IP/SMInput) at expected sites; FRiP | bedtools + idr | FRiP >= 0.005; log2 >= 3 at top peaks | Failed antibody / antibody not IP-grade |
| 5. Reproducibility | IDR rescue + self-consistency ratios | idr | both < 2 | Biological variation too high; or low complexity |

A library failing at Gate 3 (complexity) cannot be rescued analytically; gates 4-5 fail downstream of complexity by construction.

## Library Complexity with preseq

**Goal:** Determine whether the CLIP library captured enough independent molecules to support genome-wide peak calling (ENCODE: >= 1M unique fragments per replicate).

**Approach:** Run `preseq lc_extrap` on the PRE-dedup BAM (preseq counts PCR duplicates to extrapolate); also compute picard ESTIMATED_LIBRARY_SIZE at sequenced depth. Flag any library predicted to plateau below 3M unique fragments at infinite depth.

```bash
# After alignment, BEFORE UMI dedup (preseq counts PCR duplicates)
preseq lc_extrap \
    -B -P \
    -o sample_complexity.txt \
    sample_aligned.bam

# Output columns:
# TOTAL_READS  EXPECTED_DISTINCT  LOWER_0.95CI  UPPER_0.95CI
# At 100M reads, EXPECTED_DISTINCT:
#   >= 10M  = excellent complexity
#   3-10M   = acceptable; restrict to high-expression transcripts
#   < 3M    = library failed; cannot rescue analytically

# picard direct estimate at current depth
picard EstimateLibraryComplexity \
    I=sample_aligned.bam \
    O=picard_complexity.txt
# ESTIMATED_LIBRARY_SIZE > 5M = healthy CLIP library
```

A linear plateau on preseq's curve at low depth indicates over-amplification; a curve still climbing at sequenced depth means more sequencing would yield more unique fragments.

## FRiP (Fraction Reads in Peaks)

FRiP measures how much of the IP signal falls into the called peak set. ENCODE eCLIP narrow-binding RBP minimum: FRiP >= 0.005. Atypical-binding RBPs (rare-transcript binders like TROVE2 on Y RNAs) are exempt.

```bash
# Reads in peaks (use stringent peaks: log2 FC >= 3, -log10 p >= 3)
reads_in_peaks=$(bedtools intersect -c -s -a peaks.stringent.bed -b dedup.bam | awk '{s+=$NF} END {print s}')
total_reads=$(samtools view -c -F 4 dedup.bam)
frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)
echo "FRiP: $frip"

# Per-region FRiP breakdown
for region in three_utr exon intron; do
    rip=$(bedtools intersect -c -s -a peaks_${region}.bed -b dedup.bam | awk '{s+=$NF} END {print s}')
    echo "${region}: $(echo "scale=4; $rip / $total_reads" | bc)"
done
```

| RBP class | Expected FRiP (ENCODE eCLIP) |
|-----------|------------------------------|
| Splicing factors (PTBP1, U2AF2) | 0.01 - 0.10 |
| 3' UTR mRNA stability (HuR, PUM2) | 0.02 - 0.20 |
| Translation factors (EIF3J) | 0.01 - 0.05 |
| Repeat binders (MATR3) | 0.05 - 0.30 (high; concentrated in repeats) |
| Mitochondrial (FASTKD2) | 0.05 - 0.40 (very high; chrM is small) |
| snoRNA binders (DKC1) | 0.10 - 0.50 (high; snoRNA is rare) |
| Failed IP (any RBP) | < 0.005 |

## IDR for CLIP Reproducibility

IDR (Li et al 2011) measures peak-rank reproducibility across replicates. ENCODE eCLIP convention applies IDR identically to ChIP-seq, using CLIPper + SMInput log2 FC + -log10 p as the ranking signal. The CLIP-specific consideration: rank by signalValue (log2 FC) or p-value, NOT by score column (CLIPper score is sparse and tied).

```bash
# Sort each replicate's peaks by signal
sort -k5,5gr rep1.compressed.bed > rep1.sorted
sort -k5,5gr rep2.compressed.bed > rep2.sorted

# True-replicates IDR (threshold 0.05)
idr --samples rep1.sorted rep2.sorted \
    --input-file-type bed --rank 5 \
    --output-file idr_true.out \
    --idr-threshold 0.05 \
    --plot --log-output-file idr.log

# Pseudo-replicates from each individual replicate (split BAM in half)
samtools view -b -h -s 1.5 rep1.dedup.bam > rep1.psr1.bam   # seed 1, fraction 0.5
samtools view -b -h -s 2.5 rep1.dedup.bam > rep1.psr2.bam   # seed 2 (different)
# Re-run peak calling on each pseudoreplicate, then IDR at threshold 0.10
```

**ENCODE consistency rules for eCLIP:**
- Nt = peaks passing IDR on true replicates
- Nself = peaks passing IDR on pseudo-replicates of each rep
- Library passes if: max(Nt, Nself) / min(Nt, Nself) <= 2
- If both ratios > 2: library rejected

## SMInput vs IgG Control: Which?

The eCLIP design uses a size-matched input (SMInput) from the SAME lysate, treated identically (UV, IP buffer, RNase, ligation, IP, but with NO antibody addition - just bead-only control or a non-specific control IP). This is fundamentally different from IgG controls and from RNA-seq.

| Control | What it measures | Pros | Cons |
|---------|------------------|------|------|
| SMInput | Background from non-specific binding + ligation/RT/gel biases at same size | Captures all CLIP-specific biases; ENCODE standard | Requires same-day prep; cannot use a previous IgG library |
| IgG-IP | Non-specific antibody binding to the same RBP-naive lysate | Direct nonspecificity measure | Yields very low (~3-10x less reads); high PCR dup; hard to normalize |
| Empty beads (mock) | Bead surface non-specificity only | Cleanest baseline | Misses real CLIP background (IP buffer + ligation step bias) |
| RNA-seq (matched cell type) | Transcript abundance | Easy to obtain | Misses CLIP-specific biases entirely; not a real CLIP control |
| Total RNA / nuclear RNA | Cellular RNA distribution | Easy | Same as RNA-seq above |

**Consensus (ENCODE / Hentze / Yeo):** SMInput. The bead-only IgG/mock alternatives are ill-suited for quantification because their library yields are 5-10x lower, dominated by PCR duplicates, and produce sparse read-density tracks. RNA-seq cannot replace SMInput because it does not capture the non-specific binding that occurs during IP.

```bash
# SMInput preparation - same as IP except no antibody added during incubation
# Same UV dose, same lysate aliquot, same RNase, same library prep
# Critical: same SDS-PAGE size cut from membrane as the IP

# Check SMInput vs IP enrichment. A ratio of whole-library totals only measures relative sequencing
# depth; normalize the in-peak fraction in each library instead:
#   log2( (IP_in_peaks/IP_total) / (SMI_in_peaks/SMI_total) )
# > 0.5 indicates the IP concentrates reads into peaks beyond SMInput (good)
# ~ 0 indicates failed IP (SMInput == IP)
```

## Read Distribution Metagene

```bash
# RSeQC read_distribution.py shows fractional read placement
read_distribution.py -i dedup.bam -r gencode.v38.bed > read_dist.txt

# Output reports:
#   CDS_Exons, 5'UTR_Exons, 3'UTR_Exons, Introns, TES_down_10kb, TES_down_1kb,
#   TSS_up_10kb, TSS_up_1kb, Intergenic_region
# CLIP-seq biology-specific patterns:
#   Splicing factor: > 60% Introns + 5' UTR exons (containing 5' splice sites)
#   3' UTR factor: > 50% 3'UTR_Exons
#   m6A reader: 3'UTR_Exons + Stop_codon region
#   Failed IP: matches RNA-seq distribution (no enrichment)

# GeneBody coverage for 5' vs 3' bias
geneBody_coverage.py \
    -i dedup.bam \
    -r housekeeping.bed \
    -o sample_gb
# Flat curve = no positional bias (normal for most RBPs)
# 3' end bias = polyA-dependent enrichment (suspect)
# 5' end bias = nascent / TSS-proximal (only sensible for some RBPs)
```

## Fragment-Length Distribution (Paired-End)

```bash
# Paired-end fragment-length distribution
samtools view -f 2 dedup.bam | awk '$9 > 0 && $9 < 500 {print $9}' | sort -n | uniq -c > fragment_lengths.txt

# CLIP normal range: 20-75 nt insert (rises from short trimmed reads to ~75 nt)
# Wider range (20-200) seen in high-RNase / long-fragment protocols
# Narrow peak at one length (e.g., all reads at 30 nt) = over-trimmed
# Bimodal at 30 nt and 150 nt = library preparation artifact
```

## Pre-Map rRNA / snoRNA Contamination Check

```bash
# rRNA reads as fraction of total
total=$(samtools view -c -F 4 dedup.bam)
rRNA=$(samtools idxstats dedup.bam | awk '$1 ~ /rRNA|45S|18S|28S|5_8S/ { sum+=$3 } END {print sum}')
rrna_frac=$(echo "scale=4; $rRNA / $total" | bc)
echo "rRNA fraction: $rrna_frac"

# eCLIP without pre-map: rRNA 5-30% normal
# After pre-map: < 2% expected
# > 30% indicates rRNA dominance - IP captured ribosomes preferentially
# Some RBPs (RPL/RPS) expected to bind ribosomes; verify against RBP biology
```

## Antibody Validation Sanity Check

The antibody is the single most common point of failure in CLIP. Even commercial "IP-grade" antibodies fail at 30-50% rates in practice. Cross-check IP enrichment against expected biology:

```python
import pandas as pd

# Load peak file
peaks = pd.read_csv('peaks.stringent.bed', sep='\t', header=None,
                    names=['chr','start','end','name','log2fc','strand'])

# Filter top 100 peaks
top_peaks = peaks.nlargest(100, 'log2fc')

# Check enrichment in expected biology
# 1. If RBP is splicing factor: top peaks should be intronic or splice-site flanking
# 2. If RBP is HuR: > 70% top peaks should be 3' UTR
# 3. If RBP is FASTKD2: top peaks should be chrM
# 4. If RBP is FUS: GUGGU motif should appear in top motif enrichment

# Quick check
chrom_dist = top_peaks['chr'].value_counts(normalize=True)
print(chrom_dist.head(10))
# Healthy RBP: chromosomes represented in proportion to expression
# Failed IP: top chromosome chrM (mt artifact) or chr21/22 (housekeeping bias) > 20%
```

GO term sanity: top-peak genes should enrich for the expected biology (e.g., HuR -> immune / inflammation / mRNA stability terms; PTBP1 -> RNA splicing terms). If the top GO term is "cellular metabolism" or "translation" for a splicing factor, the IP failed.

## Per-Stage Failure Modes

### Gate 1: Preprocessing retention < 70%

**Trigger:** cutadapt log shows > 30% reads filtered (too short or no adapter).

**Mechanism:** Adapter pattern wrong; or RNA was degraded before fragmentation (all reads adapter-only).

**Symptom:** Catastrophic loss at the adapter trim step.

**Fix:** Verify adapter sequence in library prep documentation; check library quality control (Bioanalyzer / TapeStation) for RIN value; re-prep if RIN < 7.

### Gate 2: Alignment rate < 60%

**Trigger:** STAR Log.final.out reports `Uniquely mapped reads %` < 60% for eCLIP, < 70% for iCLIP.

**Mechanism:** rRNA contamination dominating (no pre-map); wrong species genome; or chimeric library (contamination).

**Symptom:** Low alignment rate; samtools idxstats shows most reads as rRNA-aligned.

**Fix:** Run bowtie2 pre-map to rRNA + repeats index; verify genome species; check for sample swap with `samtools view -h sample.bam | grep '^@SQ'`.

### Gate 3: Library complexity < 1M unique

**Trigger:** preseq lc_extrap predicts plateau < 1M unique at sequenced depth; or picard reports ESTIMATED_LIBRARY_SIZE < 1M.

**Mechanism:** Over-amplification (> 25 PCR cycles); low input (< 5M cells); failed IP capturing only a few molecules.

**Symptom:** preseq curve flattens early; UMI clusters have median size > 8.

**Fix:** No analytic rescue. Re-prep with more input cells (10-20M for eCLIP), fewer PCR cycles (14-18 for eCLIP, 16-20 for iCLIP2). If forced to use this library, restrict downstream analysis to high-expression transcripts and acknowledge the limitation.

### Gate 4: FRiP < 0.005 for narrow-binding RBP

**Trigger:** FRiP value below ENCODE threshold for the RBP class.

**Mechanism:** IP failed to enrich specific binding sites; antibody cross-reactive or not IP-grade.

**Symptom:** Top peaks dominated by abundant transcripts (GAPDH, ACTB); GO enrichment generic; motif analysis returns AU-rich background.

**Fix:** Re-test antibody on a knockdown lysate (siRNA against the RBP) - if WB signal does not decrease, antibody is non-specific; switch antibody. ENCODE-validated antibodies are at encodeproject.org/biosamples.

### Gate 4 variant: IP/SMInput global log2 ~ 0

**Trigger:** Whole-genome log2(IP/SMInput) is near 0 instead of positive.

**Mechanism:** IP and SMInput are equivalent - the antibody is not enriching any RBP-bound RNA.

**Symptom:** Per-peak log2 FC scaled around 0; no obvious enrichment at known motif sites.

**Fix:** Same as above - antibody failure. Consider tagged-RBP system (Halo-CLIP, GoldCLIP, or knock-in endogenous tag).

### Gate 5: IDR fails - rescue/self-consistency > 2

**Trigger:** ENCODE IDR test reports rescue ratio > 2 OR self-consistency > 2.

**Mechanism:** Replicates inconsistent. Biological variation high; OR one replicate has lower library complexity; OR one replicate failed IP.

**Symptom:** Per-replicate peak counts differ > 2x; IDR plot shows the two replicates have very different peak rank distributions.

**Fix:** Run preseq and FRiP on each replicate independently; identify the failing one; re-sequence to higher depth or re-prep. ENCODE-style: down-sample both replicates to common depth, re-test IDR.

### High PCR duplication confused with low complexity

**Trigger:** Saw 60% PCR duplication and panicked.

**Mechanism:** CLIP libraries have 40-70% PCR duplication BY DESIGN - the IP enriches a small pool. UMI dedup recovers the unique molecules underneath.

**Symptom:** "60% duplication" headline. Actual unique count is what matters.

**Fix:** Confirm UMI dedup has been applied and unique fragment count >= 1M. The duplication rate metric is meaningless without UMI context for CLIP.

## Decision Tree by Failure

| Symptom | Likely failure | Action |
|---------|----------------|--------|
| > 50% reads "too short" in cutadapt | Adapter wrong or RNA degraded | Verify adapter; check RIN |
| Most reads align to rRNA | No pre-map; or IP captured ribosomes | bowtie2 pre-map; or accept if RBP is ribosomal |
| Unique frag count < 1M | Over-amplified or under-input | No rescue; re-prep |
| FRiP < 0.005 (narrow RBP) | Failed IP | Switch antibody |
| IP/SMInput global log2 ~ 0 | Antibody non-specific | Switch antibody or use tagged-RBP |
| IDR rescue > 2 | Replicate inconsistency | Identify failing replicate; re-sequence |
| Top GO terms generic | Failed IP or contamination | Check antibody on KD lysate WB |
| chrM peaks abundant for non-mt-RBP | Mitochondrial contamination | Filter chrM or accept if mt-RBP |
| Fragment length all 30 nt | Over-trimmed | Loosen quality trim from -q 20 to -q 6 |
| Strand-specific peaks lost | bedtools without -s | Add `-s` strand flag |

## Standard QC Report (MultiQC)

```bash
# Aggregate all QC tools into a single report
multiqc \
    fastqc_output/ \
    cutadapt_logs/ \
    star_logs/ \
    umi_tools_logs/ \
    preseq_output/ \
    samtools_stats/ \
    -o multiqc_report/
```

MultiQC reads logs from FastQC, cutadapt, STAR, umi_tools, preseq, samtools stats, picard, and others, and produces a single HTML report. For CLIP-seq, additionally run RSeQC `read_distribution.py` and document FRiP / IDR results manually.

## ENCODE-Style Comprehensive QC Table

| Metric | Threshold | Source |
|--------|-----------|--------|
| Biological replicates | >= 2 | ENCODE eCLIP |
| Read length | >= 50 nt (PE) | ENCODE |
| Unique fragments per replicate | >= 1M (or saturated) | ENCODE |
| Preprocessing retention | >= 70% | Practitioner consensus |
| Genome alignment rate (unique) | >= 60% eCLIP, 70% iCLIP | Practitioner consensus |
| Library complexity (preseq) | >= 10M expected unique at 100M | ENCODE |
| FRiP (narrow-binding RBP) | >= 0.005 | ENCODE |
| FRiP (atypical-binding RBP) | Exempt | ENCODE |
| log2(IP/SMInput) at stringent peaks | >= 3 | ENCODE |
| -log10 p at stringent peaks | >= 3 | ENCODE |
| IDR rescue ratio | < 2 | ENCODE |
| IDR self-consistency ratio | < 2 | ENCODE |
| rRNA fraction post pre-map | < 2% | Practitioner consensus |
| Read distribution match to RBP class | Y | Sanity check |
| Top motif consistent with literature | Y | Sanity check |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| preseq "all reads have duplicates" | Pre-deduplicated BAM passed to preseq | Re-run on PRE-dedup BAM |
| FRiP value > 1.0 | Peak set overlaps reads counted twice | Use unique-reads BAM; verify peak BED is non-redundant |
| IDR plot empty | Ranking column wrong; all peaks tied at same value | Rank by `--rank 5` (log2 FC); sort BED first |
| MultiQC missing modules | Logs not in expected directory structure | Verify log paths; rerun multiqc with explicit `-d` |
| Read distribution: > 50% intergenic for mRNA RBP | TxDb missing transcripts; rRNA contamination | Update GENCODE; pre-map rRNA |
| GeneBody coverage 3' biased | polyA-selected library; not CLIP | Verify library prep was random hexamer |
| Fragment-length distribution narrow at 30 nt | Adapter aggressively trimmed at 5' | Loosen trim parameters |
| Antibody KD WB shows no decrease | Antibody non-specific | Switch antibody; use ENCODE-validated |

## References

- Van Nostrand EL et al 2016 Nat Methods 13:508 (eCLIP, ENCODE QC standards)
- Li Q et al 2011 Ann Appl Stat 5:1752 (IDR framework)
- Landt SG et al 2012 Genome Res 22:1813 (ChIP/CLIP QC guidelines, IDR Nself rule)
- Daley T & Smith AD 2013 Nat Methods 10:325 (preseq library complexity)
- Wang L et al 2012 Bioinformatics 28:2184 (RSeQC)
- Ewels P et al 2016 Bioinformatics 32:3047 (MultiQC)
- Smith T et al 2017 Genome Res 27:491 (UMI-tools)
- ENCODE eCLIP Data Standards (encodeproject.org/eclip) - canonical thresholds
- Van Nostrand EL et al 2020 Nature 583:711 (ENCODE 150 RBP QC patterns)

## Related Skills

- clip-seq/clip-preprocessing - Gate 1 preprocessing logs
- clip-seq/clip-alignment - Gate 2 alignment metrics
- clip-seq/clip-peak-calling - Gates 4-5 peak QC + FRiP
- clip-seq/differential-clip - QC for differential analyses
- read-qc/quality-reports - General FastQC / MultiQC
- read-qc/contamination-screening - Cross-species contamination
- chip-seq/chipseq-qc - DNA-protein QC analogue
