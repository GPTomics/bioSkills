---
name: bio-chipseq-spike-in-normalization
description: Normalizes ChIP-seq data using exogenous spike-in (ChIP-Rx with Drosophila chromatin per Orlando 2014 / Egan 2016; E. coli carryover for CUT&RUN/CUT&Tag). Distinguishes RRPM from Rx-Input scaling, integrates with DiffBind / DESeq2 / edgeR / csaw via sizeFactors and DiffBind library-size vectors, and applies the Patel et al 2024 *Nat Biotechnol* review's failure-mode framework to validate that normalization is correctly applied at the read level (not peak counts). Use when global signal shifts are expected (HDACi, BETi, EZH2i, dosage, target knockdown), when ChIPseqSpikeInFree detects post-hoc shifts, or when validating internal-control regions before publication.
tool_type: mixed
primary_tool: DiffBind
---

## Version Compatibility

Reference examples tested with: DiffBind 3.20+, DESeq2 1.42+, edgeR 4.0+, csaw 1.36+, ChIPseqSpikeInFree 1.6+, SpikChIP 1.0+, SpikeFlow (NAR Genom Bioinform 2024), samtools 1.19+, bowtie2 2.5+.

# ChIP-seq Spike-In Normalization

**"Account for global signal changes that defeat standard normalization"** -> Add exogenous reference chromatin (Drosophila for human/mouse ChIP-Rx; E. coli carryover for CUT&RUN/CUT&Tag) at fixed concentration BEFORE IP, derive scaling factors from spike-in read counts, and apply at the read or size-factor level (never to peak counts) to enable quantitative cross-condition comparison.

- CLI: align reads to combined target + spike genome; count spike reads via `samtools view -c`
- R (DiffBind integration): `dba.normalize(obj, spikein = TRUE)`
- R (DESeq2 / edgeR): `sizeFactors(dds) <- 1 / scale_factors` (note inverse)
- CLI (deepTools tracks): `bamCoverage --scaleFactor <derived>` (mutually exclusive with `--normalizeUsing`)
- Wrapper: SpikeFlow (Snakemake; 2024) automates end-to-end
- Post-hoc detection: ChIPseqSpikeInFree (when no spike-in was added)

The fundamental rule: spike-in scaling is applied at the READ level (via size factors or `--scaleFactor`), never multiplied into peak counts. ~25% of published spike-in ChIP papers violate this (Patel L et al 2024 *Nat Biotechnol* 42:1343).

## When Spike-In Is Required

| Experimental design | Spike-in needed? |
|---------------------|------------------|
| HDAC inhibitor -> global H3K27ac increase | Yes |
| BET inhibitor (JQ1, OTX015) -> global BRD4 / H3K27ac decrease | Yes |
| EZH2 inhibitor -> global H3K27me3 loss | Yes |
| DNMT inhibitor -> global 5mC loss; downstream histone mark shifts | Yes |
| Target factor knockdown / degron | Yes (or matched-input subtraction) |
| Cell-cycle synchronization / arrest | Yes |
| Dosage titration | Yes |
| Standard TF perturbation, local rebinding expected | No (reads-in-peaks RLE works) |
| Histone mark cross-cell-type comparison | Recommended |
| CUT&RUN/CUT&Tag standard | E. coli carryover (automatic); deliberate Drosophila for high-stakes |
| Replicate-only experiment, no condition comparison | No |

**Why this is necessary:** Standard normalization (RLE on reads-in-peaks, TMM on bins) assumes most regions don't change. When the perturbation IS the change-everything-globally biology, these methods force the median log2FC to zero, hiding the real effect.

## Spike-In Protocol Taxonomy

| Protocol | Spike organism | Added when | Notes |
|----------|----------------|------------|-------|
| **ChIP-Rx** (Orlando 2014) | Drosophila S2 nuclei | After lysis, before IP | 50k Drosophila nuclei per 5M target cells (Egan 2016) |
| **ChIP-Rx variant** (Bonhoure 2014) | Drosophila chromatin | After fragmentation, before IP | Different normalization layer |
| **CUT&RUN/Tag E. coli** | E. coli (carryover) | Automatic from bacterial pA-MNase/Tn5 | Free; variable across enzyme batches |
| **Heterologous spike-in** (Hu 2015) | Defined yeast / E. coli chromatin | Added at lysis | Less common; defined concentration |
| **xenoChIP** | Species swap (mouse cells + human chromatin spike) | Before IP | Niche; specific cancer xenograft contexts |

**The dominant standard for human/mouse ChIP is Drosophila (ChIP-Rx).** Drosophila is genetically distinct enough that mapping is unambiguous, and the genome size (~140 Mb) gives adequate read depth at small chromatin input.

## Scaling Factor Calculation

**RRPM (Orlando 2014):** reference-adjusted reads per million.

```
scale_factor_i = min(N_spike) / N_spike_i
```

Apply at the read level. The sample with the fewest spike reads gets scale_factor = 1; others get > 1 (compensating for less observed spike chromatin).

**Rx-Input (Fursova 2019):** additionally scales by input spike-in to correct IP efficiency variation.

```
RxInput_i = (N_spike_chip_i / N_total_chip_i) / (N_spike_input_i / N_total_input_i)
```

This is more rigorous when input controls are available; required for some inhibitor experiments where IP efficiency itself changes.

## Workflow: Drosophila ChIP-Rx Spike-In

**Goal:** Compute per-sample scaling factors from Drosophila spike-in reads and apply at the read level (not peak counts) to enable quantitative cross-condition ChIP-seq comparison.

**Approach:** Align reads to combined target + Drosophila genome, count spike reads at high mapq after deduplication, derive RRPM scaling factors (min/each), then apply via DESeq2 sizeFactors, DiffBind spike-in flag, or bamCoverage scaleFactor. Validate against internal-control regions (blacklist).

### Step 1: Alignment to combined genome

```bash
# Build combined index (target + Drosophila)
cat hg38.fa dm6.fa > hg38_dm6.fa
bowtie2-build hg38_dm6.fa hg38_dm6

# Align reads
bowtie2 -x hg38_dm6 -1 R1.fq -2 R2.fq -S aln.sam --very-sensitive --no-mixed
samtools view -bS aln.sam | samtools sort -o aln.bam
samtools index aln.bam
```

### Step 2: Filter, deduplicate, count spike reads

```bash
# Apply ENCODE filter (-F 1804 -q 30) BEFORE counting spike reads
samtools view -F 1804 -q 30 -b aln.bam > aln.filt.bam
samtools index aln.filt.bam

# Count Drosophila reads (NOT total reads)
DROSO_READS=$(samtools view -c aln.filt.bam chr2L chr2R chr3L chr3R chr4 chrX chrY)
echo "$SAMPLE: Drosophila reads = $DROSO_READS"

# Separate into target-only BAM for peak calling
samtools view -b aln.filt.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    > aln.filt.hg38.bam
samtools index aln.filt.hg38.bam
```

### Step 3: Compute scaling factors

```bash
# Per-sample Drosophila counts (assume saved in droso_counts.tsv)
# sample_id, droso_reads
# ctrl_1, 145000
# ctrl_2, 132000
# treat_1, 98000
# treat_2, 85000

awk 'BEGIN{min=1e10} NR>1{if($2<min) min=$2} END{print "min:", min}' droso_counts.tsv
# Use min as numerator: scale_factor_i = min / droso_reads_i
```

### Step 4: Apply scaling — three layers

**Layer 1: bigWig tracks**

```bash
SCALE=$(echo "scale=6; $MIN_DROSO / $SAMPLE_DROSO" | bc)
bamCoverage -b sample.bam -o sample.scaled.bw \
    --scaleFactor $SCALE --binSize 10 --extendReads 200
# DO NOT also pass --normalizeUsing; mutually exclusive
```

**Layer 2: DiffBind**

```r
library(DiffBind)
dba_obj <- dba(sampleSheet = 'samples.csv')   # spike-in BAM in sample sheet
dba_obj <- dba.count(dba_obj, summits = 250, bParallel = TRUE)

# Spike-in normalization
dba_obj <- dba.normalize(dba_obj, spikein = TRUE,
                          library = DBA_LIBSIZE_FULL,
                          normalize = DBA_NORM_LIB)

# Verify what was applied
dba.normalize(dba_obj, bRetrieve = TRUE)
```

**Layer 3: DESeq2 / edgeR direct**

```r
library(DESeq2)
# Read spike-in counts into a vector aligned with sample order
spike_reads <- c(ctrl_1 = 145000, ctrl_2 = 132000, treat_1 = 98000, treat_2 = 85000)
scale_factors <- min(spike_reads) / spike_reads

dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
# DESeq2 expects sizeFactors in INVERSE convention (sample with smallest factor gets largest sizeFactor)
sizeFactors(dds) <- 1 / scale_factors
dds <- DESeq(dds, fitType = 'parametric')
```

## Workflow: E. coli Spike-In (CUT&RUN/CUT&Tag Automatic)

E. coli DNA from bacterial pA-MNase/pA-Tn5 production is automatic spike-in carryover.

```bash
# Combined index
cat hg38.fa ecoli_k12.fa > hg38_ecoli.fa
bowtie2-build hg38_ecoli.fa hg38_ecoli

# Align as in ChIP-Rx; count E. coli reads
ECOLI_READS=$(samtools view -c aln.filt.bam ecoli_chr1)
TOTAL_READS=$(samtools view -c aln.filt.bam)
echo "E. coli fraction: $(echo "scale=4; $ECOLI_READS / $TOTAL_READS" | bc)"
# Target: 0.005-0.02 (0.5-2%); IgG: 0.02-0.05 (2-5%)

# Scale factor same as ChIP-Rx: min(ecoli) / per_sample_ecoli
# Apply at read or sizeFactors level
```

**E. coli carryover is variable** between enzyme production batches. For publication-grade cross-condition claims, supplement with deliberate Drosophila spike-in OR use a single enzyme lot across all experiments.

## ChIPseqSpikeInFree: Post-Hoc Detection

When no spike-in was added, ChIPseqSpikeInFree (Jin 2020) attempts post-hoc detection of global shifts by analyzing signal-distribution shape changes.

```r
library(ChIPseqSpikeInFree)

samples <- data.frame(
    ID = c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'),
    BAM = c('ctrl_1.bam', 'ctrl_2.bam', 'treat_1.bam', 'treat_2.bam'),
    ANTIBODY = rep('H3K27me3', 4),
    GROUP = c('Control', 'Control', 'Treatment', 'Treatment')
)

res <- ChIPseqSpikeInFree(bamFiles = samples$BAM, chromFile = 'hg38.chrom.sizes',
                          metaFile = 'metadata.txt', prefix = 'spikein_free_out')
# Output: per-sample scaling factor + global-shift detection
```

**Limitations:** Heuristic; not a substitute for true spike-in. Use as:
1. Sanity check when spike-in was forgotten
2. Initial diagnosis before deciding whether spike-in is needed in next experiment
3. NOT for publication-grade claims

## Internal-Control Sanity Check (Mandatory)

After applying spike-in scaling, internal-control regions should show NO signal change:

| Region type | Source | Expected behavior post-spike-in |
|-------------|--------|----------------------------------|
| ENCODE blacklist v2 | Amemiya 2019 | No change (artifact regions) |
| Constitutive housekeeping promoters | Eisenberg 2013 list (HK genes); U6 snRNA promoter | Minor change only |
| Custom hyper-ChIPable regions | Top-1% input signal | Stable signal at artifact regions |
| Untouched chromosome (e.g., chrY in cell types without expression) | Genome | No signal change |

```bash
# Compute mean signal at blacklist regions per condition; should be stable
bedtools multicov -bams ctrl_1.bam ctrl_2.bam treat_1.bam treat_2.bam \
    -bed hg38-blacklist.v2.bed > blacklist_signal.tsv
# Apply scaling factors to per-sample counts; verify no shift across conditions
```

If internal controls shift after scaling, the normalization is broken. Common causes (Patel L et al 2024 *Nat Biotechnol* 42:1343):
1. Scaling applied to peak counts instead of read counts
2. Spike-in reads not deduplicated before scaling
3. Spike-in genome not mapq-filtered (low-quality alignments inflated counts)
4. Spike-in saturated (>1M reads); titration not linear

## Per-Tool Failure Modes

### Scaling factor applied to peak counts instead of read counts

**Trigger:** Multiplying peak-by-sample count matrix entries by spike-in factor.

**Mechanism:** Peak counts already integrate over read counts; multiplying them double-corrects.

**Symptom:** Effect sizes 2-10× larger than expected biology; internal control regions also "shift" artifactually.

**Fix:** Apply via `sizeFactors(dds)` (DESeq2), `normFactors` (edgeR), or DiffBind's `dba.normalize(..., library=<numeric vector>, normalize=DBA_NORM_LIB)` to supply spike-in-derived library sizes, OR `--scaleFactor` (bamCoverage for tracks). Never multiply peak-level counts.

### Spike-in reads not deduplicated before scaling

**Trigger:** Counting all aligned reads to spike genome including duplicates.

**Mechanism:** PCR duplicates of spike-in reads vary independently of input chromatin amount.

**Symptom:** Scaling factors poorly correlated with library prep batch; high inter-replicate variability.

**Fix:** Deduplicate with MarkDuplicates; apply ENCODE filter `-F 1804 -q 30` before counting spike reads.

### Spike-in mapq filter too loose

**Trigger:** Counting all reads aligning to spike genome.

**Mechanism:** Low-mapq reads at low-complexity regions (E. coli rRNA, Drosophila satellite) are often misaligned from host genome.

**Fix:** Apply `-q 30` (high mapq) before counting spike reads.

### Inverse convention errors with DESeq2 / edgeR

**Trigger:** Passing `scale_factors` directly to `sizeFactors(dds)` without inversion.

**Mechanism:** DESeq2 / edgeR `sizeFactors` are MULTIPLIED to normalized counts; spike-in scaling factors should DIVIDE. Convention difference.

**Symptom:** Effect sizes inverted (treatment shifted in wrong direction).

**Fix:** `sizeFactors(dds) <- 1 / scale_factors` (inverse). Verify with internal-control sanity check.

### `--normalizeUsing` and `--scaleFactor` conflict in bamCoverage

**Trigger:** Passing both for spike-in scaled bigWig.

**Mechanism:** deepTools applies scaleFactor first, then normalizes; normalization undoes the scale.

**Fix:** Use ONE: `--scaleFactor` alone for spike-in; `--normalizeUsing` alone otherwise. Verify via `bamCoverage --help`.

### E. coli carryover inconsistent across enzyme batches

**Trigger:** Comparing CUT&Tag samples processed with different pA-Tn5 lots.

**Mechanism:** E. coli carryover varies between bacterial production batches; cross-batch comparison adds artificial variability.

**Fix:** Use single enzyme lot for cross-condition comparison; OR supplement E. coli with deliberate Drosophila spike-in.

### ChIPseqSpikeInFree applied as primary normalization

**Trigger:** No spike-in was added; ChIPseqSpikeInFree used for publication-grade scaling.

**Mechanism:** ChIPseqSpikeInFree infers global shift from signal-distribution shape; this is a heuristic, not a measurement.

**Fix:** Use only as diagnostic. For publication, re-do experiment with deliberate spike-in.

### Spike-in titration not verified linear

**Trigger:** Spike-in concentration too high (>5% of total reads) OR too low (<0.1%).

**Mechanism:** Outside linear range, scaling factor doesn't reflect actual ratio of input chromatin.

**Symptom:** Replicate-to-replicate scaling factor variability >2×.

**Fix:** Verify titration linearity by varying spike-in concentration on a single sample; only use spike-in counts in linear range (typically 0.5-5% of total reads).

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Spike-in scaled vs CPM give opposite signs | Global shift; CPM forced to median; spike-in revealed it | Spike-in is correct; CPM is fooled |
| Scaling factor varies wildly between reps | Spike-in saturated / not in linear range | Verify titration; subsample if needed |
| Internal-control signal shifts after scaling | Scaling applied wrong layer; reads not dedup'd; mapq too loose | Apply pre-test diagnostic; recompute |
| ChIPseqSpikeInFree predicts shift but spike-in says no | Both interpretations possible; trust spike-in when available | Spike-in measurement > distribution heuristic |
| DiffBind spike-in vs manual sizeFactors differ | DiffBind applies inverse convention internally | Verify via `dba.normalize(obj, bRetrieve=TRUE)` |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Spike-in BAM column missing in DiffBind sample sheet | `bamSpikeIn` (DiffBind 3.x) vs older `spikein` field | Use `spikein = TRUE` in `dba.normalize()` with appropriate column |
| Drosophila reads on chromosome X include host chrX | Combined genome chromosome naming collision | Prefix Drosophila chroms with `dm_` before combining |
| Scaling factors all close to 1 | Spike-in not added at fixed amount | Verify Egan 2016 protocol; titrate per cell count |
| Cross-condition results sign-flipped after scaling | Inverse convention bug | `sizeFactors(dds) <- 1 / scale_factors` |
| Blacklist signal shifts post-scaling | Normalization broken | Investigate per Patel L 2024 *Nat Biotechnol* 42:1343 failure modes |

## References

- Orlando DA et al 2014 Cell Rep 9:1163 (ChIP-Rx framework)
- Egan B et al 2016 PLoS One 11:e0166438 (ChIP-Rx protocol with cell counts)
- Bonhoure N et al 2014 Genome Res 24:1157 (alternative Drosophila spike-in)
- Fursova NA et al 2019 Mol Cell 74:1020 (Rx-Input scaling)
- Hu B et al 2015 PLoS One 10:e0145007 (heterologous spike-in)
- Jin H et al 2020 Bioinformatics 36:1270 (ChIPseqSpikeInFree)
- Blanco E et al 2021 NAR Genom Bioinform 3:lqab064 (SpikChIP)
- 2024 NAR Genom Bioinform 6:lqae118 (SpikeFlow)
- Patel L, Cao Y, Mendenhall EM, Benner C, Goren A 2024 Nat Biotechnol 42:1343 (spike-in normalization review; common failure modes; PMC12266361)
- Stark R & Brown G 2011 Bioconductor (DiffBind with spikein parameter)

## Related Skills

- chip-seq/peak-calling - Upstream peak calling
- chip-seq/chipseq-qc - Spike-in fraction QC
- chip-seq/differential-binding - Apply spike-in via DiffBind / DESeq2 / csaw
- chip-seq/cut-and-run-tag - E. coli spike-in carryover specifics
- chip-seq/super-enhancers - SE calling requires spike-in for cross-condition
- chip-seq/chipseq-visualization - Spike-in-scaled bigWig generation
- alignment-files/sam-bam-basics - Multi-genome alignment and chromosome filtering
- differential-expression/deseq2-basics - DESeq2 sizeFactors conventions
