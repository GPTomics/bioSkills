---
name: bio-chipseq-qc
description: Assesses ChIP-seq quality across antibody specificity, fragmentation, enrichment, replicate concordance, and library complexity. Computes FRiP, NSC/RSC (phantompeakqualtools), library complexity (NRF/PBC1/PBC2), deepTools plotFingerprint (JS distance, AUC, synthetic JS), ChIPQC, IDR with ENCODE Nself/Nt rules, and detects hyper-ChIPable artifacts. Use when validating an antibody, diagnosing failed peak calls, deciding whether to proceed with downstream analysis, grading against ENCODE thresholds, or auditing replicate concordance.
tool_type: mixed
primary_tool: deepTools
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: deepTools 3.5+, phantompeakqualtools 1.2.2+, ChIPQC 1.42+, IDR 2.0.4+, samtools 1.19+, bedtools 2.31+, pysam 0.22+, pybedtools 0.9+, MACS2 2.2.9+, MACS3 3.0.4+.

Verify versions before relying on numerical thresholds — phantompeakqualtools has known R-version compatibility issues with R ≥ 4.0 (use kundajelab fork or pin to R 3.6).

# ChIP-seq Quality Control

**"Should I trust this ChIP-seq experiment?"** -> Validate antibody, fragmentation, enrichment, replicate concordance, library complexity, and absence of hyper-ChIPable artifacts before committing to downstream peak calling and differential analysis.

- CLI: `Rscript run_spp.R -c=chip.bam -out=cc.txt` (NSC/RSC), `plotFingerprint -b chip.bam input.bam` (enrichment shape), `idr --samples rep1.np rep2.np` (replicate IDR)
- R: ChIPQC package (Carroll & Stark; computes the full ENCODE metric battery)
- Python: pysam + pybedtools for custom FRiP and library-complexity metrics

ChIP-seq fails for many independent reasons. The QC metrics below probe distinct failure modes — passing one metric does not rescue another. Antibody failure cannot be fixed by sequencing more.

## The Antibody Problem is the Real Problem

Every downstream metric is conditional on antibody specificity. "ChIP-grade" on a vendor datasheet is marketing, not validation. Run the cascade:

| Step | What | Why |
|------|------|-----|
| 1. Western blot | Expected MW + KO/KD negative | Confirms the antibody hits a band of the right size and loses signal in KO |
| 2. IP-Western | Pulls down the protein | Confirms IP recovery, not just recognition |
| 3. ChIP-qPCR | Known positive + known negative loci | First chromatin-context test; cheap |
| 4. ChIP-seq biological replicate | Two independent biological replicates | Reproducibility check |
| 5. KO/KD orthogonal | ChIP in KO/KD cells | Gold-standard: signal should drop to background |
| 6. Peptide array (histones) | Epicypher SNAP-ChIP or equivalent | Tests modification-state specificity |

Histone modification cross-reactivity is universal: H3K9me2 vs H3K9me3, H3K27me2 vs H3K27me3, and H3K4me1 vs H3K4me2 antibodies routinely show 10-30% cross-reactivity. Polyclonals vary lot-to-lot. CRISPR-knockout-validated lots from CST and Epicypher are the modern standard. Always record antibody catalog number + lot in methods.

## Fragment-Size Distribution is a Free Diagnostic

The fragment-size distribution from a properly prepared ChIP BAM is itself a quality readout:

| Distribution shape | Interpretation |
|--------------------|----------------|
| Sharp peak at ~50-100 bp (sub-nucleosomal) | Direct TF binding; expected for well-fragmented TF ChIP |
| Sharp peak at ~150 bp + secondary at ~300 bp | Mono- + di-nucleosomal; expected for histone ChIP |
| Bimodal at 150 + 300, no sub-nucleosomal | Histone-only signal; in TF ChIP, suggests trapping / hyper-ChIPable |
| Broad continuum 100-1000 bp | Over-sonication; biology lost; cannot be rescued |
| No peak structure, flat | Severe over-sonication or library prep failure |

```bash
# Quick diagnostic — count fragment sizes from properly-paired reads
samtools view -f 0x2 sample.bam | awk '{print $9}' | awk '$1>0' \
    | sort -n | uniq -c | awk '{print $2, $1}' > fragment_sizes.tsv
```

For CUT&Tag: 25-75 bp characteristic; fragments < 25 bp are Tn5 self-tagmentation noise (see cut-and-run-tag).

## QC Metric Battery with ENCODE Thresholds

| Metric | Tool | TF threshold | Histone threshold | Source / rationale |
|--------|------|--------------|-------------------|---------------------|
| **FRiP** (Fraction of Reads in Peaks) | bedtools / pysam / featureCounts | ≥ 0.01 minimum, > 0.05 ideal | ≥ 0.05, > 0.20 ideal; > 0.15 for H3K4me3 | Landt 2012; ENCODE flags experiments with FRiP < 1% |
| **NSC** (Normalized Strand Cross-correlation) | phantompeakqualtools | > 1.05 marginal, > 1.10 ideal | > 1.05 | Landt 2012; min = 1 (no enrichment); ratio of fragment-length CC to background |
| **RSC** (Relative Strand Cross-correlation) | phantompeakqualtools | > 0.8 marginal, > 1.0 ideal | > 0.8 | Landt 2012; ratio of (fragment - background) / (phantom - background) |
| **QualityTag** | phantompeakqualtools | ≥ 0 acceptable, 1-2 ideal | ≥ 0 | Composite based on NSC/RSC; -2 to 2 scale |
| **NRF** (Non-Redundant Fraction) | unique_pos / total | > 0.8 | > 0.8 | ENCODE; < 0.5 severe PCR bottleneck |
| **PBC1** (M1 / Mdistinct) | bedtools / pysam | > 0.8 | > 0.8 | ENCODE; fraction of singly-occupied positions |
| **PBC2** (M1 / M2) | bedtools / pysam | > 3 | > 3 | ENCODE; ratio of singletons to doubletons |
| **JS distance** (plotFingerprint) | deepTools | > 0.3 | > 0.05 (broad) to > 0.3 (narrow) | Distance between cumulative signal curves IP vs Input |
| **AUC** (plotFingerprint) | deepTools | < 0.6 | 0.6-0.9 | Input = ~0.5; lower AUC = more enrichment concentrated |
| **Synthetic JS** (plotFingerprint) | deepTools | Should ≈ measured JS | — | Sanity check vs simulated null |
| **Replicate Spearman correlation** | deepTools multiBamSummary / plotCorrelation | > 0.8 (true reps) | > 0.8 (true reps), > 0.6 (broad) | Replicates should correlate more than cross-condition |
| **Read count per replicate** | samtools flagstat | ≥ 20M unique mapped | 20M (narrow histone), 40-60M (broad histone) | ENCODE 2012 |

**Practical operational rule:** Compute the full battery. Reject any sample failing FRiP OR antibody validation OR fragment-size sanity check, regardless of other metrics. Failing one of NSC/RSC alone with strong FRiP can sometimes be rescued for narrow-peak biology; broad histones are more forgiving on NSC.

## Hyper-ChIPable Region Detection

Teytelman 2013 (PNAS): untagged GFP, no antibody, or non-existent targets all produce "binding" signal at highly-transcribed loci (rRNA, tRNA, histone gene clusters, snoRNA hosts, mtDNA, abundant housekeeping genes). ENCODE blacklist v2 (Amemiya 2019) catches repeat-driven artifacts but NOT these hyper-ChIPable transcribed regions.

**Detection:**

```bash
# Top 1% input signal as cell-type-specific custom blacklist
multiBigwigSummary BED-file -b input.bw -o input_signal.npz \
    --BED genes.bed --outRawCounts input_per_gene.tsv
awk 'NR > 1' input_per_gene.tsv | sort -k4,4nr | head -n $(($(wc -l < input_per_gene.tsv) / 100)) \
    > hyper_chipable.bed

# Intersect peaks against this list; flag peaks falling in hyper-ChIPable regions
bedtools intersect -a peaks.narrowPeak -b hyper_chipable.bed -u > suspicious_peaks.bed
```

**Disprove a suspicious peak:** Required for any claim at rRNA loci, tRNA clusters, HIST1/2 clusters, mitochondrial DNA:
1. Motif enrichment at peak (artifact has no enrichment)
2. KO/KD signal loss at peak (artifact persists)
3. Untagged-protein control ChIP shows no signal at this locus

Many "novel binding" claims at the rDNA repeat, mtDNA, and HIST1 cluster are spurious artifacts.

## Computing the Battery

### FRiP

```bash
total_reads=$(samtools view -c -F 260 chip.bam)
reads_in_peaks=$(bedtools intersect -a chip.bam -b peaks.narrowPeak -u | samtools view -c -)
frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)
```

### NSC / RSC / fragment length (phantompeakqualtools)

```bash
Rscript run_spp.R -c=chip.bam -savp=qc/chip_cc.pdf -out=qc/chip_cc.txt
# Output columns: filename | numReads | estFragLen | corr_estFragLen |
#                 phantomPeak | corr_phantomPeak | argmin_corr | min_corr |
#                 NSC | RSC | QualityTag
```

### Library complexity

```bash
# NRF
total=$(samtools view -c -F 260 chip.bam)
unique=$(samtools view -F 260 chip.bam | awk '{print $1, $3, $4}' | sort -u | wc -l)
nrf=$(echo "scale=4; $unique / $total" | bc)

# PBC1, PBC2 (singletons vs distinct positions vs doubletons)
samtools view -F 260 chip.bam | awk '{print $3":"$4}' | sort | uniq -c \
    | awk '{
        if ($1 == 1) m1++;
        if ($1 == 2) m2++;
        mdist++;
      } END {
        print "M1:", m1; print "M2:", m2; print "Mdistinct:", mdist;
        print "PBC1:", m1/mdist; print "PBC2:", m1/m2
      }'
```

### deepTools plotFingerprint

```bash
plotFingerprint \
    -b chip.bam input.bam \
    --labels ChIP Input \
    -o qc/fingerprint.pdf \
    --outRawCounts qc/fingerprint_counts.tab \
    --outQualityMetrics qc/fingerprint_qc.txt
# Inspect qc/fingerprint_qc.txt: AUC, JS distance, synthetic JS, X-intercept
# Good ChIP: AUC < 0.6 (TF), JS > 0.3 (TF); Input near diagonal (AUC ~ 0.5)
```

### Replicate Spearman correlation

```bash
multiBamSummary bins -b rep1.bam rep2.bam rep3.bam input.bam \
    --binSize 10000 -o results.npz
plotCorrelation -in results.npz --corMethod spearman \
    --whatToPlot heatmap --plotNumbers -o corr.pdf \
    --outFileCorMatrix corr_matrix.tab
# Replicates: > 0.8 (narrow), > 0.6 (broad)
# Cross-condition reps should correlate less than within-condition
```

### ChIPQC R package

```r
library(ChIPQC)
samples <- read.csv('samples.csv')
qc <- ChIPQC(samples, annotation = 'hg38')
ChIPQCreport(qc, reportFolder = 'ChIPQCreport')
# Generates the full ENCODE battery report per sample in one call
```

ChIPQC remains Bioconductor-maintained but mature; phantompeakqualtools is the canonical NSC/RSC source.

## IDR and Replicate Consistency Rules

For TFs: signal-ranked IDR with Nself/Nt consistency check; see chip-seq/peak-calling for full ENCODE workflow. Key thresholds:

- True replicate IDR threshold: 0.05
- Pseudoreplicate IDR threshold: 0.10 (per-rep self-consistency)
- **Nself/Nt rule:** `max(N1self, N2self) / min(N1self, N2self) ≤ 2` AND `max(Nt, max(Nself)) / min(Nt, min(Nself)) ≤ 2`. Failing both ratios rejects the library.

For histones: naive overlap with ≥ 40% reciprocal overlap (ENCODE default; commonly misquoted as 50%). IDR is too conservative for histone signal dynamic range.

## ENCODE 3 vs ENCODE 4 Thresholds (unchanged for most QC)

| Metric | ENCODE 3 | ENCODE 4 |
|--------|----------|----------|
| FRiP minimum | 1% | 1% (unchanged) |
| NSC threshold | > 1.05 | > 1.05 (unchanged) |
| RSC threshold | > 0.8 | > 0.8 (unchanged) |
| NRF threshold | > 0.8 | > 0.8 (unchanged) |
| Blacklist | v1 | v2 (Amemiya 2019) |
| Read depth (TF) | ≥ 20M unique mapped | ≥ 20M unchanged |
| Read depth (broad histone) | ≥ 40M | 40-60M recommended |

Most QC thresholds are stable across ENCODE versions; blacklist update is the main practical change.

## Per-Tool Failure Modes

### phantompeakqualtools / SPP -- R version incompatibility

**Trigger:** Running with R ≥ 4.0.

**Mechanism:** spp R package has unmaintained Boost / Rcpp dependencies; some shifts produce NaN cross-correlation values.

**Symptom:** NSC = NaN, RSC = NaN, or fragment length = 0 in output.

**Fix:** Pin to R 3.6 + spp 1.16 via conda env; OR use the kundajelab/phantompeakqualtools fork (current); OR substitute deepTools plotFingerprint for enrichment QC and `macs3 predictd` for fragment length.

### deepTools plotFingerprint -- Wrong baseline assumption for broad marks

**Trigger:** Interpreting JS distance with TF threshold (> 0.3) on broad histone mark.

**Mechanism:** Broad marks have less concentrated signal; JS distance is naturally lower (0.05-0.15 for H3K27me3) without indicating bad ChIP.

**Symptom:** Reports "failed JS distance" for high-quality broad-mark ChIP.

**Fix:** Use mark-specific thresholds: > 0.3 for TFs and sharp histones; > 0.05 for broad histones; check AUC instead (0.6-0.9 for broad; < 0.6 for TF/sharp).

### FRiP -- Computed before vs after blacklist filtering

**Trigger:** Calling FRiP from peak file pre- vs post-blacklist.

**Mechanism:** Hyper-ChIPable regions inflate "reads in peaks" because most reads at those loci are artifacts.

**Symptom:** FRiP looks great (>15%) but most of it is rRNA / mtDNA reads.

**Fix:** Apply blacklist + custom hyper-ChIPable filter BEFORE computing FRiP; or report both raw and filtered FRiP.

### NRF / PBC -- Computed after deduplication

**Trigger:** Running NRF on a MarkDuplicates-filtered BAM.

**Mechanism:** Library complexity metrics measure PCR redundancy; if duplicates are already removed, NRF = 1.0 by construction (uninformative).

**Symptom:** NRF reports 0.99-1.0; metric is meaningless.

**Fix:** Compute NRF / PBC1 / PBC2 on the PRE-deduplication BAM. ENCODE-compliant pipeline: filter -> MarkDuplicates (don't remove) -> compute NRF -> filter out duplicates -> call peaks.

### IDR -- Wrong rank column

**Trigger:** Sorting narrowPeak by signalValue (column 7) for IDR.

**Mechanism:** MACS signalValue scales with pile-up intensity which differs between libraries of different depth; rank correlation breaks.

**Symptom:** IDR returns 0 reproducible peaks despite good replicate Spearman correlation.

**Fix:** Sort by p-value (`-k8,8nr`), pass `--rank p.value` to IDR. ENCODE convention.

### ChIPQC -- Default annotation mismatch

**Trigger:** Using `annotation = 'hg19'` on hg38-aligned data.

**Mechanism:** ChIPQC computes feature-context enrichment from the specified annotation; mismatch silently corrupts enrichment metrics.

**Symptom:** Promoter / 5'UTR / 3'UTR enrichments look wrong; replicate report metrics drift.

**Fix:** Match `annotation` to the genome the BAMs were aligned to; for custom genomes pass a TxDb object explicitly.

## Reconciliation: When Metrics Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Good FRiP, bad NSC | High background but real enrichment | Acceptable for broad marks; for TFs, check phantompeakqualtools fragment length is reasonable |
| Good NSC, bad FRiP | Strong cross-correlation signal but few peaks pass q-value | Library shallow OR peak caller threshold too strict; try `-p 1e-2` |
| Good FRiP and NSC, bad replicate correlation | Real biology + replicate-specific batch effect | Check sample swap; check sequencing batch; consider PCA |
| Good Rep1, bad Rep2 | One replicate failed | Drop Rep2 + repeat; do NOT average metrics |
| All metrics fail | Antibody or fragmentation failure | Re-validate antibody (KO/KD); inspect fragment-size distribution; do not proceed |
| FRiP excellent at rRNA/mtDNA | Hyper-ChIPable artifact dominance | Build custom blacklist; recompute |

**Operational rule for proceeding with downstream analysis:** Require (1) antibody validated, (2) fragment-size distribution sane, (3) FRiP, NSC, RSC pass ENCODE thresholds, (4) Nself/Nt rule satisfied for TFs OR naive overlap concordance for histones, (5) hyper-ChIPable artifacts identified and either filtered or flagged.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Sequence chrM not found` (multi-tool) | chrM removed from BAM but kept in genome FASTA | Match chromosome naming convention; consistently include or exclude chrM |
| phantompeakqualtools hangs / OOM | Default tag chunking on deep libraries | Subsample to 15-25M reads (`samtools view -s`) before running |
| plotFingerprint blank or near-diagonal | Input control mislabeled as ChIP | Verify sample labels; AUC ~ 0.5 = Input-like signal |
| IDR runs but Nself ratio always > 2 | One pseudoreplicate dominates due to seed | Use sufficiently different seeds (`-s 1.5` and `-s 2.5`) |
| ChIPQC report missing peaks | samples.csv path columns wrong | Verify bamReads / Peaks paths; ChIPQC fails silently on missing files |
| Replicate Spearman > 0.95 | Technical (not biological) replicates | Treat as one sample; do not report as biological replicates |

## References

- Landt SG et al 2012 Genome Res 22:1813 (ENCODE/modENCODE QC guidelines, IDR Nself rule)
- Kharchenko PV et al 2008 Nat Biotechnol 26:1351 (SPP, NSC/RSC framework)
- Li Q et al 2011 Ann Appl Stat 5:1752 (IDR)
- Marinov GK et al 2014 G3 4:209 (large-scale ChIP-seq QC comparison)
- Teytelman L et al 2013 PNAS 110:18602 (hyper-ChIPable regions)
- Amemiya HM et al 2019 Sci Rep 9:9354 (ENCODE blacklist v2)
- Ramírez F et al 2016 Nucleic Acids Res 44:W160 (deepTools)
- Carroll TS et al 2014 Front Genet 5:75 (ChIPQC framework)
- Diaz A et al 2012 Stat Appl Genet Mol Biol 11:Article 9 (CHANCE; fingerprint-style QC)
- Park PJ 2009 Nat Rev Genet 10:669 (foundational review)
- Rothbart SB et al 2015 Sci Rep 5:8489 (peptide-array antibody specificity)

## Related Skills

- chip-seq/peak-calling - Use QC metrics to decide whether to proceed with peak calling
- chip-seq/cut-and-run-tag - CUT&RUN/CUT&Tag QC differs (spike-in % aligned, fragment-size signatures)
- chip-seq/spike-in-normalization - QC for spike-in carryover and Drosophila read depth
- chip-seq/differential-binding - Replicate concordance required before differential testing
- atac-seq/atac-qc - Parallel QC for ATAC-seq (no input control, different thresholds)
- alignment-files/bam-statistics - General BAM-level QC
- alignment-files/duplicate-handling - MarkDuplicates before NRF computation
