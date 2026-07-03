---
name: bio-epitranscriptomics-merip-preprocessing
description: Aligns and QCs methylated-RNA-immunoprecipitation (MeRIP / m6A-seq) IP and input libraries using STAR or HISAT2 splice-aware mapping, samtools sort/index, IP/input matched-pair tracking, antibody-lot metadata recording, replicate concordance via deepTools multiBamSummary + plotCorrelation, IP enrichment QC via plotFingerprint and per-transcript IP/input ratio distributions, library-complexity saturation curves via PreSeq, and the explicit do-NOT-deduplicate convention for standard non-UMI MeRIP. Use when preparing paired IP and input BAM files for exomePeak2 / MeTPeak / MACS3 peak calling, evaluating MeRIP replicate concordance and IP enrichment, deciding whether to deduplicate (standard MeRIP typically NOT), choosing genome-vs-transcriptome alignment for downstream peak vs m6Anet workflows, recording antibody clone and lot metadata for cross-batch reconciliation, detecting failed IPs via saturation curves and IP/input distribution shape, or generating IP-over-Input bigWig tracks for visualisation.
tool_type: cli
primary_tool: STAR
---

## Version Compatibility

Reference examples tested with: STAR 2.7.11+, HISAT2 2.2.1+, samtools 1.19+, deepTools 3.5+, PreSeq 3.2+, fastp 0.23+, Trim Galore 0.6.10+, Picard 3.1+, MultiQC 1.25+, bowtie2 2.5+, BWA-MEM2 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

STAR `--outSAMtype` accepts `BAM SortedByCoordinate` since 2.5.x; check the `Log.final.out` file for input/output statistics. deepTools `bamCompare --operation log2` is the modern syntax (older `--ratio log2` still works but is being phased out). PreSeq `c_curve` and `lc_extrap` have stable interfaces; HISAT2 reports unique vs multi-mapped in the summary log.

# MeRIP-seq Preprocessing

**"Get my MeRIP IP and input libraries ready for peak calling"** -> Trim adapters with MeRIP-appropriate defaults (do NOT trim UMIs unless the library is UMI-MeRIP — most are not), splice-aware-align IP and input to the GENOME (not transcriptome) with STAR / HISAT2, sort and index, evaluate replicate concordance and IP enrichment with deepTools, build a saturation curve per library so peak counts can be honestly compared across libraries, record antibody clone and lot metadata so cross-batch comparison is later auditable, and produce IP-over-Input log2 bigWig tracks for downstream visualisation. Crucially, do NOT deduplicate non-UMI MeRIP — see the failure-modes section.

- CLI: `STAR --runMode alignReads` -- splice-aware genome alignment, the field default
- CLI: `hisat2 -x index -1 R1.fq.gz -2 R2.fq.gz` -- graph-based alternative; lighter memory footprint
- CLI: `samtools sort -@ 8 -o sorted.bam in.bam && samtools index sorted.bam` -- post-alignment mechanics
- CLI: `deeptools multiBamSummary bins -b *.bam -o cov.npz` + `plotCorrelation` -- replicate Spearman
- CLI: `deeptools plotFingerprint -b IP.bam Input.bam -o fp.pdf` -- IP enrichment QC (ChIP-seq term; transfers cleanly)
- CLI: `preseq lc_extrap -B -o curve.txt sorted.bam` -- library complexity / saturation
- CLI: `deeptools bamCompare -b1 IP.bam -b2 Input.bam --operation log2 -o IP_over_Input.bw` -- downstream-ready coverage track

## The Single Most Important Modern Insight -- Peak counts are library-size-dependent; saturation curves are the only honest cross-library comparison

A MeRIP library sequenced to 20 million unique reads finds substantially fewer peaks than the same biology at 60 million reads. Per-sample peak counts reported without saturation curves (PreSeq `c_curve` / `lc_extrap`; Daley & Smith 2013 *Nat Methods* 10:325) are uninterpretable across studies and often across replicates within a study. Subsample BAMs to a common unique-read depth before peak calling for any cross-condition peak-count comparison, OR report peaks alongside the saturation curve. A corollary: **do NOT deduplicate standard MeRIP** — the typical MeRIP protocol (Synaptic Systems 202-003 / Abcam ab151230 / NEB EpiMark E1610 antibody pull-down on fragmented poly(A)-selected RNA) has NO unique molecular identifiers, and `picard MarkDuplicates` on such libraries collapses real biological replicates of high-coverage transcripts (the opposite of what dedup achieves in DNA ChIP-seq). UMI-MeRIP is the only exception — and most MeRIP libraries in print are NOT UMI. McIntyre et al. 2020 *Sci Rep* 10:6590 demonstrated that replicate-to-replicate peak overlap is ~80% within a single lab but drops to a median 45% between labs using nominally identical conditions; this irreducible technical noise constrains how strongly any single MeRIP study can support biological claims, and the preprocessing pipeline is where the variance is set.

## Algorithmic Taxonomy

| Tool / step | Mechanism | Output | Strength | Fails when |
|-------------|-----------|--------|----------|------------|
| STAR 2.7+ (Dobin 2013 *Bioinformatics* 29:15) | Two-pass splice-aware alignment with on-the-fly splice junction database | Sorted BAM + splice-junction TSV | Field default; multi-mapper retention configurable; STAR splice-junction DB | Memory-heavy (~30 GB human); slower than HISAT2 |
| HISAT2 2.2.1+ (Kim 2019 *Nat Biotechnol* 37:907) | Hierarchical graph FM-index; splice-aware | Sorted BAM | ~5x lighter memory than STAR; comparable accuracy | Less mature splice-junction handling for novel introns |
| BWA-MEM2 (Vasimuddin 2019 *IPDPS* 314) | DNA-style local alignment; NO splice awareness | Sorted BAM | Use ONLY for transcriptome-aligned MeRIP (rare) | Splits reads across exon junctions if used on genome |
| fastp 0.23+ (Chen 2018 *Bioinformatics* 34:i884) | Streaming adapter detection + quality trim | Trimmed FASTQ + JSON QC | Fast; JSON-readable QC output | UMI handling disabled by default; do NOT pass `--umi` for standard non-UMI MeRIP (the opposite of the failure direction in some other library types) |
| Trim Galore | Wrapper over cutadapt with paired-end auto-detect | Trimmed FASTQ | Conservative defaults; widely cited | Slower than fastp on large datasets |
| samtools sort / index | BAM coordinate sort + .bai index | Sorted BAM + index | Standard | None at default |
| Picard MarkDuplicates | Identifies PCR duplicates by 5' alignment start | Marked / removed BAM | Standard in DNA / ChIP | The dominant MeRIP convention is to SKIP dedup for non-UMI libraries (collapses real biology at high-coverage transcripts); a minority of pipelines dedup MeRIP — record the choice in metadata |
| deepTools multiBamSummary + plotCorrelation | Per-bin read counts; Spearman / Pearson matrix | Heatmap + clustering | Standard replicate-concordance plot | Bin size sensitive (use 10 kb for transcriptome-genome) |
| deepTools plotFingerprint (Diaz 2012 *Stat Appl Genet Mol Biol* 11:9) | Cumulative read-fraction vs cumulative-bin-fraction Lorenz curve | PDF + raw counts | Direct IP-vs-input enrichment QC; "good" IP has steep tail | A flat fingerprint = failed IP (or mock IgG) |
| deepTools bamCompare --operation log2 | Per-bin log2 (IP/input) bigWig | bigWig | Ready for downstream visualisation | Pseudocount choice matters at low-coverage bins |
| PreSeq c_curve / lc_extrap (Daley & Smith 2013 *Nat Methods* 10:325) | Capture-recapture; rational-function extrapolation | Curve TSV | The only honest library-complexity estimate | Requires uniquely-mapped reads to be reliable |
| MultiQC (Ewels 2016 *Bioinformatics* 32:3047) | Aggregator across tools | HTML report | Consolidates STAR + HISAT2 + samtools + deepTools + PreSeq into one report | None |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| Standard mammalian MeRIP, downstream exomePeak2 | STAR splice-aware -> genome BAM; do NOT deduplicate; build saturation curve | Transcriptome alignment breaks exomePeak2 (expects genome BAM + GTF); dedup collapses biology |
| Downstream m6Anet (ONT direct RNA) | Defer to `m6anet-analysis` -- alignment is to TRANSCRIPTOME with minimap2 `-ax map-ont -uf -k14` | Genome-aligned ONT input breaks m6Anet entirely (signal-level dataprep requires per-transcript coordinates) |
| Limited memory (<16 GB) | HISAT2 instead of STAR | STAR human genome index requires ~30 GB |
| UMI-MeRIP (rare) | Trim UMI to read header (umi_tools / fastp `--umi_loc`), align, then dedup with `umi_tools dedup` | Standard `picard MarkDuplicates` ignores UMI; effective dedup rate wrong |
| Cross-batch comparison (different antibody lots) | Record antibody clone + lot in sample-sheet metadata; include `batch` factor in downstream design | Pooling cross-batch counts without batch term inflates false-positive differential peaks |
| Spike-in normalisation (NEB EpiMark control oligos) | Align separately to spike-in reference; report IP-spike-in / Input-spike-in ratio per sample; use for absolute normalisation | Most users discard the NEB EpiMark Gluc / Cluc controls; they are the per-sample IP-efficiency QC anchor |
| Cross-library peak-count comparison | Subsample BAMs to common unique-read depth with `samtools view -s` BEFORE downstream peak calling, OR fit saturation curves and compare at common depth | Raw peak counts are sequencing-depth-dependent and not biologically interpretable |
| Suspect failed IP | Inspect deepTools `plotFingerprint` AND per-transcript IP/input ratio distribution; failed IP shows shallow Lorenz tail AND median IP/input ~1.0 | Trusting raw peak count alone — failed IPs still produce peaks |
| Viral / contamination-suspect samples | Build a combined host + viral index (and rRNA index) and check unmapped reads | Single-organism indexes hide systematic contamination |
| Aligning to transcriptome (rare; specific downstream tools) | BWA-MEM2 or bowtie2; defer to `read-alignment/` | STAR splice-aware on transcriptome causes spurious splice calls inside transcripts |

Methodology evolves; before any high-stakes preprocessing pipeline, web-search "STAR vs HISAT2 MeRIP 2024" and "MeRIP saturation curve preseq" for current consensus parameters.

## Adapter Trimming for MeRIP

**Goal:** Remove sequencing adapters and low-quality 3' ends WITHOUT removing biological signal (UMI-MeRIP must keep UMI sequences; standard MeRIP does not have UMIs and trimming should be minimal).

**Approach:** Use fastp or Trim Galore with adapter auto-detection; require minimum read length 25-30 nt (shorter reads multi-map and confound exomePeak2); for standard non-UMI MeRIP, do NOT pass `--umi` flags; preserve random-hexamer-priming artifacts ONLY if downstream pipeline expects them (most do not).

```bash
mkdir -p trimmed

for sample in IP_rep1 IP_rep2 IP_rep3 Input_rep1 Input_rep2 Input_rep3; do
    fastp \
        --in1 raw/${sample}_R1.fastq.gz \
        --in2 raw/${sample}_R2.fastq.gz \
        --out1 trimmed/${sample}_R1.fq.gz \
        --out2 trimmed/${sample}_R2.fq.gz \
        --html qc/${sample}_fastp.html \
        --json qc/${sample}_fastp.json \
        --length_required 25 \
        --detect_adapter_for_pe \
        --thread 8
done
```

For UMI-MeRIP (rare), insert `--umi --umi_loc read1 --umi_len 8` BEFORE the alignment step. Default fastp output preserves base quality information needed by downstream variant-aware tools; do NOT pass `--disable_quality_filtering` for MeRIP libraries.

## STAR Splice-Aware Alignment for IP and Input

**Goal:** Produce coordinate-sorted, indexed GENOME BAM files for each IP and input library with splice-junction-aware mapping, retaining a moderate number of multi-mappers for accurate per-window read counts at multi-isoform loci.

**Approach:** Build STAR genome index once with the matched GENCODE / Ensembl GTF used downstream; loop IP and input samples with identical parameters; retain up to 20 multi-mappers per read (MeRIP read counts at multi-isoform genes need this); request explicit BAM SortedByCoordinate; emit splice-junction tables for QC.

```bash
STAR \
    --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile annotation.gtf \
    --sjdbOverhang 100 \
    --runThreadN 12

mkdir -p aligned

for sample in IP_rep1 IP_rep2 IP_rep3 Input_rep1 Input_rep2 Input_rep3; do
    STAR \
        --runMode alignReads \
        --genomeDir star_index \
        --readFilesIn trimmed/${sample}_R1.fq.gz trimmed/${sample}_R2.fq.gz \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 20 \
        --outSAMattributes NH HI AS nM NM MD \
        --outFileNamePrefix aligned/${sample}_ \
        --runThreadN 12

    samtools index -@ 4 aligned/${sample}_Aligned.sortedByCoord.out.bam
done
```

`--outFilterMultimapNmax 20` is intentional: MeRIP at rRNA / snoRNA / pseudogene-rich loci needs multi-mapper retention. Reduce to 1 only if downstream analysis explicitly cannot tolerate multi-mappers. `--sjdbOverhang` should equal (read length - 1) but 100 is the common-enough default for 100-150 bp reads.

## HISAT2 Alternative for Memory-Constrained Environments

**Goal:** Achieve splice-aware alignment in ~5x less memory than STAR (12-16 GB suffices for human), with comparable accuracy for MeRIP applications.

**Approach:** Build HISAT2 graph index; align with `--dta` for downstream-transcript-assembly compatibility; pipe directly to samtools sort.

```bash
hisat2-build genome.fa hisat2_index/genome

for sample in IP_rep1 IP_rep2 IP_rep3 Input_rep1 Input_rep2 Input_rep3; do
    hisat2 \
        -x hisat2_index/genome \
        -1 trimmed/${sample}_R1.fq.gz \
        -2 trimmed/${sample}_R2.fq.gz \
        --dta \
        --summary-file qc/${sample}_hisat2.log \
        -p 12 | \
    samtools sort -@ 8 -o aligned/${sample}.sorted.bam -

    samtools index -@ 4 aligned/${sample}.sorted.bam
done
```

HISAT2 multi-mapper handling is governed by `-k`; the default reports the primary alignment only. For MeRIP, pass `-k 5` if multi-mapper-aware downstream counting is required.

## Per-Sample QC: flagstat and idxstats

```bash
mkdir -p qc

for bam in aligned/*sortedByCoord.out.bam aligned/*sorted.bam; do
    name=$(basename ${bam} .bam)
    samtools flagstat ${bam} > qc/${name}.flagstat
    samtools idxstats ${bam} > qc/${name}.idxstats
done
```

Inspect `flagstat` for properly-paired rate (>=85% indicates good pairing); inspect `idxstats` for unexpected chromosome-level read piles (rRNA bleed-through, mitochondrial domination — both are MeRIP red flags).

## Replicate Concordance via deepTools

**Goal:** Quantify how similar replicate IP libraries are to each other (and likewise input libraries) using a Spearman correlation matrix; flag a divergent replicate before it propagates into peak calling.

**Approach:** Compute genome-wide per-bin read counts at 10 kb resolution across all IP and input BAMs; convert to a clustered Spearman heatmap with deepTools `plotCorrelation`.

```bash
multiBamSummary bins \
    --bamfiles aligned/IP_rep1*.bam aligned/IP_rep2*.bam aligned/IP_rep3*.bam \
                aligned/Input_rep1*.bam aligned/Input_rep2*.bam aligned/Input_rep3*.bam \
    --binSize 10000 \
    --numberOfProcessors 8 \
    --outRawCounts qc/raw_bin_counts.tab \
    -o qc/cov.npz

plotCorrelation \
    --corData qc/cov.npz \
    --corMethod spearman \
    --skipZeros \
    --whatToPlot heatmap \
    --colorMap RdYlBu_r \
    --plotNumbers \
    -o qc/replicate_correlation.pdf
```

IP replicates within a condition should cluster (Spearman >= 0.85 typical); input replicates should cluster with each other; IP and input should NOT cluster together. A failed IP looks like input.

## IP Enrichment via plotFingerprint

**Goal:** Confirm IP libraries are enriched (a few transcripts have many reads) and input libraries are uniform (reads spread across transcripts); fail-fast on poor IP before peak calling.

**Approach:** deepTools `plotFingerprint` builds a cumulative Lorenz-style curve; a steep tail = signal concentrated in few regions (good IP); a diagonal = uniform coverage (input or failed IP). The framework is from ChIP-seq (Diaz 2012 *Stat Appl Genet Mol Biol* 11:9) and transfers cleanly to MeRIP.

```bash
plotFingerprint \
    --bamfiles aligned/IP_rep1*.bam aligned/IP_rep2*.bam aligned/IP_rep3*.bam \
                aligned/Input_rep1*.bam aligned/Input_rep2*.bam aligned/Input_rep3*.bam \
    --labels IP1 IP2 IP3 In1 In2 In3 \
    --numberOfProcessors 8 \
    --skipZeros \
    --outQualityMetrics qc/fingerprint_metrics.tab \
    -o qc/fingerprint.pdf
```

Good MeRIP IP: cumulative-fraction-of-reads vs cumulative-fraction-of-bins curve sits well below the diagonal in the right half (top-X% of bins capture >50% of reads). Input: near-diagonal. The `--outQualityMetrics` file reports JS distance and synthetic JS distance; the IP-vs-Input JS distance is a single-number IP-quality summary (higher = more concentrated signal).

## Library Complexity / Saturation Curves via PreSeq

**Goal:** Compute per-library complexity so peak counts can be honestly compared across libraries and conditions of different sequencing depth.

**Approach:** PreSeq `c_curve` (interpolation up to observed depth) and `lc_extrap` (extrapolation beyond observed) on the sorted BAM. Daley & Smith 2013 *Nat Methods* 10:325 capture-recapture model.

```bash
mkdir -p complexity

for bam in aligned/*.bam; do
    name=$(basename ${bam} .bam)

    preseq c_curve -B -o complexity/${name}_c_curve.txt ${bam}

    preseq lc_extrap -B -o complexity/${name}_lc_extrap.txt ${bam}
done
```

Inspect: the `lc_extrap` curve plots distinct molecules vs total reads; a plateau indicates saturation. For cross-condition peak-count comparison: pick a common depth (often 30M unique reads), subsample with `samtools view -s 0.<frac>` to that depth, THEN call peaks.

## IP-over-Input bigWig for Downstream Visualisation

**Goal:** Produce a per-bin log2 (IP / Input) coverage track per replicate, ready for downstream metagene / browser plots.

**Approach:** deepTools `bamCompare` with `--operation log2`; choose a sensible pseudocount to avoid divide-by-zero at low-coverage bins.

```bash
mkdir -p tracks

paste -d ' ' \
    <(printf '%s\n' IP_rep1 IP_rep2 IP_rep3) \
    <(printf '%s\n' Input_rep1 Input_rep2 Input_rep3) | \
while read ip input; do
    bamCompare \
        -b1 aligned/${ip}_Aligned.sortedByCoord.out.bam \
        -b2 aligned/${input}_Aligned.sortedByCoord.out.bam \
        --operation log2 \
        --pseudocount 1 \
        --binSize 25 \
        --normalizeUsing CPM \
        --numberOfProcessors 8 \
        -o tracks/${ip}_over_${input}.bw
done
```

`--pseudocount 1` prevents division-by-zero at zero-coverage bins; `--binSize 25` is fine-grained enough to preserve peak topology while keeping bigWig files reasonably sized.

## Per-Method Failure Modes

### Dedup applied to non-UMI MeRIP

**Trigger:** `picard MarkDuplicates REMOVE_DUPLICATES=true` invoked on a standard MeRIP BAM that has no UMI.

**Mechanism:** Standard MeRIP libraries have no unique molecular identifiers. PCR duplicates and biological re-sampling at high-coverage transcripts look identical at the alignment level. Dedup removes both, collapsing real coverage at the most-abundant transcripts to an artificially flat profile. This is the opposite of dedup's intent in DNA ChIP-seq.

**Symptom:** Coverage at housekeeping mRNAs (e.g., GAPDH, ACTB) drops 5-20x after dedup; downstream peak counts at highly-expressed transcripts collapse; volcano plot of differential peaks shows expression-driven false positives.

**Fix:** Skip dedup for standard non-UMI MeRIP. If the library is UMI-MeRIP, use `umi_tools dedup` (Smith 2017 *Genome Res* 27:491) which respects UMI rather than alignment position alone. Record dedup status in sample-sheet metadata.

### Transcriptome alignment for downstream peak calling

**Trigger:** STAR or bowtie2 aligned to transcriptome FASTA, then BAM passed to exomePeak2 / MeTPeak / MACS3.

**Mechanism:** exomePeak2 and MeTPeak expect a GENOME BAM plus GTF; they project peaks back to transcript features internally. A transcriptome BAM has reads in per-transcript coordinates which the GTF cannot resolve back to genome coordinates without re-alignment.

**Symptom:** exomePeak2 throws errors on TxDb-genome consistency; MeTPeak returns zero peaks; MACS3 calls peaks on transcript IDs as if they were chromosomes.

**Fix:** Align to GENOME with STAR / HISAT2 for downstream MeRIP peak calling. Transcriptome alignment is correct only for `m6anet-analysis` (ONT DRS) and rare quantification-only downstream tools.

### Failed IP indistinguishable from input

**Trigger:** A single replicate IP library has IP/input ratio distribution centred at 1.0 across all transcripts (no enrichment); fingerprint Lorenz curve sits at the diagonal.

**Mechanism:** Failed IP — antibody-RNA binding did not enrich m6A-containing fragments. Causes include antibody-batch defect, insufficient pulldown wash, RNA degradation during IP, or accidental mock IgG IP.

**Symptom:** plotFingerprint shows IP overlaying input on the Lorenz plot; per-transcript IP/input ratio histogram is centred at 1.0; downstream peak callers find few or no peaks AT THE FAILED REPLICATE while other replicates produce normal counts.

**Fix:** Identify the failed replicate via plotFingerprint AND IP/input ratio distribution BEFORE peak calling; exclude or re-do. Single failed IP in a 3-replicate design routinely produces "differential" peaks driven entirely by the failure.

### Antibody lot mismatch across samples

**Trigger:** A multi-condition MeRIP study uses Synaptic Systems 202-003 antibody lot A for the control IPs and lot B for the treatment IPs (because lot A ran out mid-study).

**Mechanism:** Anti-m6A polyclonals (Synaptic Systems 202-003, Abcam ab151230, NEB EpiMark E1610, Cell Signaling 56593, Active Motif 61755) have batch-to-batch variability in pulldown efficiency and m6A-vs-m6Am cross-reactivity. Pooling lot-A and lot-B counts in a downstream differential model attributes lot-effect to condition.

**Symptom:** "Differential" peaks at high abundance transcripts; effect sizes track antibody lot rather than condition; reanalysis with lot in the design matrix removes most differential peaks.

**Fix:** Record `antibody_clone` and `antibody_lot` per sample in metadata; include `lot` as a fixed effect in downstream differential analysis. Within a single study, ideally use a single lot for ALL replicates and ALL conditions.

### Peak counts compared across libraries of different depth

**Trigger:** "Condition A has 14,000 peaks; condition B has 22,000 peaks; condition B has more m6A."

**Mechanism:** Peak count is library-size-dependent. A library at 60M unique reads finds more peaks than 30M. Without rarefaction or saturation correction, peak-count comparisons across libraries are dominated by sequencing depth.

**Symptom:** Peak counts track total mapped reads more closely than they track biological condition; downstream "biological m6A change" claims do not survive rarefaction-to-common-depth.

**Fix:** Either rarefy all BAMs to common unique-read depth before peak calling, OR fit saturation curves with PreSeq `lc_extrap` and compare at matched depth, OR report peak count alongside the saturation curve.

### Random hexamer priming over-trim

**Trigger:** Aggressive 5' trimming of the first 6-12 nt to remove "random hexamer priming bias" applied to MeRIP libraries.

**Mechanism:** Random hexamer priming bias affects the 5' nucleotide composition of reads but does NOT degrade downstream peak-calling accuracy. Over-trimming removes biological signal and shortens reads enough to inflate multi-mapper fraction.

**Fix:** Standard adapter trimming with `--length_required 25` is sufficient; do not 5'-trim for hexamer bias unless downstream tooling explicitly requires unbiased 5' ends (most do not). The bias is a known artifact in the RNA-seq community and is robust to standard analytical pipelines.

## Reconciliation: When QC Signals Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| plotFingerprint diagonal but IP/input ratio shows enrichment | Mismatched chromosome naming (chr1 vs 1) between samples | Verify `samtools view -H bam | grep '@SQ'` matches across samples |
| Replicate Spearman 0.95 but plotFingerprint diverges | Replicates correlate in bulk but differ in IP enrichment depth | Check per-sample sequencing depth; reduce to common depth |
| Saturation curve plateaus early but peak count low | Library complexity exhausted (e.g., over-amplified PCR) | Inspect duplicate rate; re-prep library if possible |
| MultiQC reports input has higher mapping rate than IP | IP enriches non-canonical sequences (m6A on intronic RNA, mt-RNA) that map differently | Acceptable if STAR multi-mapper retention is on; verify with idxstats |
| Properly-paired rate < 60% | Insert size distribution off (RNA degradation; library prep failure) | Inspect `samtools view -f 0x2` count; re-prep if severe |
| HISAT2 reports many `discordant` pairs | Splice-junction not captured in index | Re-build with `--dta` and confirm GTF matches genome |
| plotFingerprint synthetic JS distance < 0.5 | Marginal IP enrichment; borderline failed | Inspect per-transcript IP/input ratio distribution; consider exclusion |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| Minimum read length after trimming | 25 nt | Below this, multi-mapping fraction inflates; downstream peak callers lose specificity |
| STAR `--outFilterMultimapNmax` for MeRIP | 20 | Retains multi-isoform mapping; tighten to 1 only when downstream cannot tolerate |
| `--sjdbOverhang` | read length - 1 | STAR convention; 100 is common for 100-150 bp reads |
| Properly-paired rate (samtools flagstat) | >=85% | Below indicates degraded RNA or library-prep failure |
| Replicate Spearman correlation (multiBamSummary 10 kb bins) | >=0.85 (IP-vs-IP within condition) | Below suggests one replicate is anomalous |
| plotFingerprint IP-vs-input JS distance | >=0.5 | Higher indicates better IP enrichment; <0.3 suggests failed IP |
| Saturation curve plateau depth | ~30-60M unique reads typical | Below this, peak calling under-samples; depth depends on cell type / antibody |
| Per-transcript IP/input ratio median | >1.5 (genome-wide median) | Lower suggests failed IP; conditions / cell lines vary |
| Minimum biological replicates | 3 (4-5 preferred) per condition | McIntyre 2020 *Sci Rep* 10:6590 — N=2 routinely under-powered |
| Dedup status for non-UMI MeRIP | OFF | Standard convention; UMI-MeRIP is the only exception |
| BAM sort order for downstream tools | Coordinate (SortedByCoordinate) | exomePeak2, MeTPeak, MACS3, deepTools all expect coordinate sort |
| `bamCompare --binSize` for downstream metagene | 25 | Fine enough to preserve peak topology; coarser only for whole-chromosome browser views |
| `bamCompare --pseudocount` | 1 | Prevents divide-by-zero at zero-coverage bins; larger values flatten signal |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| STAR runs out of memory on human genome | --genomeDir build needs ~30 GB RAM | Use HISAT2 (~12 GB) or run STAR on a high-memory node |
| `samtools index` fails with "is not coordinate sorted" | BAM is name-sorted or unsorted | Re-run `samtools sort` (not `sort -n`) |
| exomePeak2 errors on TxDb chromosome mismatch | BAM uses chr1, GTF uses 1 (or reverse) | Verify with `samtools view -H bam | head` and `head genes.gtf`; reconcile with `seqlevelsStyle()` in R or rename in shell |
| deepTools `bamCompare --ratio log2` deprecation warning | Newer deepTools uses `--operation log2` | Switch to `--operation log2` |
| PreSeq `lc_extrap` rejects with "low complexity" | Library too shallow OR genome too small (BAM under 1M unique reads) | Use `c_curve` only; or sequence deeper |
| MultiQC misses STAR Log.final.out | STAR output naming non-standard | Re-run with `--outFileNamePrefix` and rerun MultiQC; check `multiqc_config.yaml` search patterns |
| `picard MarkDuplicates` collapses all reads to 1 per position | Tiny BAM or single read pair per fragment | Verify BAM has many properly-paired reads; do NOT dedup non-UMI MeRIP regardless |
| Empty fingerprint output | All BAMs have identical bin coverage | Verify BAMs are different files; check `multiBamSummary --outRawCounts` |
| bigWig file size too large | Bin size too small at deep coverage | Increase `--binSize` from 25 to 50; bigWig is lossy at large bin sizes |
| Saturation curve never plateaus | Library deeply under-sampled | Sequence deeper OR accept curve does not plateau and report accordingly |
| `fastp --umi` errors on non-UMI library | UMI flag passed but library has no UMI | Drop `--umi` flag for standard non-UMI MeRIP |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Was deduplication applied?" | No — standard non-UMI MeRIP protocol; PCR duplicate vs biological resampling indistinguishable without UMI; dedup collapses real coverage at high-expression transcripts |
| "What is the IP enrichment QC?" | deepTools plotFingerprint reported per replicate; JS distance >=0.5 vs input |
| "Are the replicates concordant?" | Spearman correlation matrix reported via deepTools plotCorrelation on 10 kb bins; IP-IP within condition >=0.85 |
| "Saturation curve?" | PreSeq lc_extrap per library; libraries rarefied to common depth before downstream peak calling |
| "What antibody clone and lot?" | Recorded per sample in metadata; same lot for all replicates within study |
| "Why STAR instead of HISAT2?" | STAR splice-junction-DB-based vs HISAT2 graph-based; both valid for MeRIP; choice driven by memory budget |
| "How many biological replicates?" | N >=3 per condition (per McIntyre 2020); N=2 is under-powered for differential downstream |
| "Was alignment to genome or transcriptome?" | Genome (required for exomePeak2 / MeTPeak / MACS3 downstream); transcriptome alignment is for m6anet-analysis only |

## References

- Dobin A, Davis CA, Schlesinger F et al (2013) STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29(1):15-21. doi:10.1093/bioinformatics/bts635
- Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019) Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol* 37(8):907-915. doi:10.1038/s41587-019-0201-4
- Vasimuddin Md, Misra S, Li H, Aluru S (2019) Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. *IPDPS* 314-324. doi:10.1109/IPDPS.2019.00041
- Chen S, Zhou Y, Chen Y, Gu J (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884-i890. doi:10.1093/bioinformatics/bty560
- Ramírez F, Ryan DP, Grüning B et al (2016) deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Res* 44(W1):W160-W165. doi:10.1093/nar/gkw257
- Diaz A, Park K, Lim DA, Song JS (2012) Normalization, bias correction, and peak calling for ChIP-seq. *Stat Appl Genet Mol Biol* 11(3):Article 9. doi:10.1515/1544-6115.1750
- Daley T, Smith AD (2013) Predicting the molecular complexity of sequencing libraries. *Nat Methods* 10(4):325-327. doi:10.1038/nmeth.2375
- Smith T, Heger A, Sudbery I (2017) UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Res* 27(3):491-499. doi:10.1101/gr.209601.116
- McIntyre ABR, Gokhale NS, Cerchietti L, Jaffrey SR, Horner SM, Mason CE (2020) Limits in the detection of m6A changes using MeRIP/m6A-seq. *Sci Rep* 10(1):6590. doi:10.1038/s41598-020-63355-3
- Ewels P, Magnusson M, Lundin S, Käller M (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 32(19):3047-3048. doi:10.1093/bioinformatics/btw354
- Dominissini D, Moshitch-Moshkovitz S, Schwartz S et al (2012) Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. *Nature* 485(7397):201-206. doi:10.1038/nature11112
- Meyer KD, Saletore Y, Zumbo P, Elemento O, Mason CE, Jaffrey SR (2012) Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. *Cell* 149(7):1635-1646. doi:10.1016/j.cell.2012.05.003

## Related Skills

- m6a-peak-calling - Immediate downstream consumer of the IP/input BAM pairs
- m6a-differential - Downstream differential analysis on peak count matrices; design matrix relies on IP/input pairing recorded here
- m6anet-analysis - ONT DRS alternative; uses TRANSCRIPTOME alignment with minimap2, NOT the genome BAMs produced here
- modification-visualization - Uses the bigWig output of bamCompare for metagene plots and browser tracks
- read-qc/quality-reports - FastQC / MultiQC upstream of trimming
- read-alignment/star-alignment - General STAR splice-aware alignment patterns
- read-alignment/hisat2-alignment - HISAT2 graph-based alternative; general usage
- alignment-files/sam-bam-basics - General BAM mechanics, samtools fundamentals
- alignment-files/bam-statistics - flagstat / idxstats / per-chromosome counts
- alignment-files/duplicate-handling - General dedup philosophy (note: NOT applicable to non-UMI MeRIP)
- chip-seq/chipseq-qc - ChIP-seq IP QC concepts (FRiP, fingerprint, library complexity) that transfer directly
- chip-seq/peak-calling - General IP-vs-input peak-calling concepts
- clip-seq/clip-preprocessing - Antibody-RNA crosslink protocols (miCLIP / m6A-CLIP) overlap with MeRIP design
- rna-quantification/featurecounts-counting - Count matrix construction for downstream differential
- workflows/rnaseq-to-de - End-to-end pipeline orchestration patterns
