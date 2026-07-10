---
name: bio-workflows-clip-pipeline
workflow: true
depends_on: [clip-seq/clip-preprocessing, clip-seq/clip-alignment, clip-seq/clip-qc, clip-seq/clip-peak-calling, clip-seq/crosslink-site-detection, clip-seq/binding-site-annotation, clip-seq/clip-motif-analysis, clip-seq/differential-clip]
qc_checkpoints: [preprocessing_retention, alignment_rate, library_complexity, frip, idr]
description: End-to-end CLIP-seq pipeline from FASTQ to ENCODE-compliant binding sites, single-nucleotide crosslink maps, annotation, motifs, and (optionally) differential binding. Use when running the full Yeo lab eCLIP / iCLIP / iCLIP2 / iCLIP3 / irCLIP / PAR-CLIP analysis with SMInput control, protocol-specific UMI extraction, ENCODE STAR parameters, CLIPper or Skipper peak calling with stringent log2 FC and -log10 p thresholds, IDR rescue and self-consistency QC, and downstream motif registration with mCross or PEKA.
tool_type: mixed
primary_tool: CLIPper
---

## Version Compatibility

Reference examples tested with: umi_tools 1.1.5+, cutadapt 4.6+, fastp 0.23+, STAR 2.7.11b+, samtools 1.19+, bedtools 2.31+, CLIPper 2.0+, Skipper (commit 2023.05+), PureCLIP 1.3.1+, HOMER 4.11+, ChIPseeker 1.40+, preseq 3.2+, picard 3.1+, idr 2.0.4+, MultiQC 1.21+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws unexpected errors, introspect the installed tool and adapt the example rather than retrying.

# CLIP-seq End-to-End Pipeline

**"Analyze my CLIP-seq data from raw FASTQ to ENCODE-compliant binding sites"** -> Orchestrate protocol-specific UMI extraction, 3'-only adapter trimming (preserving the R2 5' truncation = crosslink site -1), ENCODE STAR alignment, UMI-based deduplication, library complexity QC, peak calling against SMInput with stringent thresholds (log2 FC >= 3 AND -log10 p >= 3), single-nucleotide crosslink-site detection, ChIPseeker annotation with CLIP-appropriate `tssRegion`, motif discovery with GC-matched background and CL-position registration, and optional differential binding between conditions.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

A CLIP callset is decided at four seams, not inside the peak caller.

1. **The protocol variant is the master commitment made once at the top.** It fixes the UMI pattern, the STAR mismatch ceiling, and the crosslink signal (truncation vs PAR-CLIP T->C vs STAMP C->U edit). The CLIP Variant Selection table below is that decision; everything downstream inherits it.
2. **Every IP is normalized against its SMInput with ENCODE-stringent thresholds (log2 FC >= 3 AND -log10 p >= 3).** Peaks called without SMInput normalization are enrichment-uncontrolled and dominated by abundance.
3. **The R2 5' end IS the crosslink site (-1), so preprocessing and alignment must PRESERVE it.** Trim 3'-only and permissively (`-q 6`, never `-g` on R1), and align with STAR `--alignEndsType EndToEnd` — soft-clipping or aggressive 5' trimming destroys the truncation base and with it single-nucleotide resolution.
4. **40-70% PCR duplication is BY DESIGN; low duplication signals a FAILED IP, not a clean library.** The IP enriches a small molecule pool, so the real quality metric is the UNIQUE-fragment count after UMI dedup — never the raw duplication rate.

## Pipeline Overview

```
FASTQ + SMInput
  -> [clip-preprocessing]    UMI extract + 3' adapter trim (-q 6 -m 18) + two-pass for eCLIP
  -> [clip-alignment]        STAR ENCODE block (alignEndsType EndToEnd, mismatch 0.04 or 0.07 for PAR-CLIP) + UMI dedup
  -> [clip-qc]               preseq, FRiP, IDR rescue + self-consistency, read distribution
  -> [clip-peak-calling]     CLIPper + SMInput log2 norm (stringent: log2 FC >= 3, -log10 p >= 3) OR Skipper (substantially more sites)
  -> [crosslink-site-detection] PureCLIP or CTK CITS for single-nt CL positions
  -> [binding-site-annotation] ChIPseeker (tssRegion=c(-100,100), level=transcript) + RBP-Maps for splicing factors
  -> [clip-motif-analysis]   HOMER + mCross (registered) + RBNS Kd cross-check
  -> [differential-clip]     DEWSeq window-level NB with type:condition interaction (optional)
```

## CLIP Variant Selection

| Variant | When to use | UMI pattern | STAR mismatch ceiling | Detection signal |
|---------|-------------|-------------|----------------------|------------------|
| eCLIP (Van Nostrand 2016) | ENCODE comparability; SMInput available | 10 nt R1 | 0.04 | R2 5' truncation |
| iCLIP / iCLIP2 / iCLIP3 | Single-end; high motif specificity | NNNXXXXNN (3+4+2; demux first) | 0.04 | R1 5' truncation |
| irCLIP / FLASH | Non-radioactive; fast | Protocol-specific | 0.04 | Truncation |
| PAR-CLIP | Photoactivatable nucleoside (4SU); HEK293/K562 | 4 nt typical | 0.07 (raised for T->C) | T->C transitions |
| miCLIP / miCLIP2 | m6A modification | iCLIP-style | 0.04 | Truncation + C->T at m6A |
| STAMP / scSTAMP | Antibody-free; in vivo or single-cell | NA (no UV) | 0.04 (RNA-seq mode) | C->U editing (RBP-APOBEC1 fusion) |
| chimeric eCLIP / miR-eCLIP | Direct miRNA-target pairs | 10 nt R1 | 0.04 | Chimeric reads |

## Step 1: Quality Control of Raw FASTQ

```bash
# Initial QC
fastqc raw_R1.fq.gz raw_R2.fq.gz -o qc/raw/

# Inspect first 12 bases of 100 reads to verify UMI pattern matches the prep
zcat raw_R1.fq.gz | awk 'NR%4==2' | head -100 | cut -c1-12 | sort | uniq -c | sort -rn | head
# Random barcode positions show ~25% per base; library barcodes are fixed
```

## Step 2: Preprocessing (Protocol-Specific)

**Goal:** Convert raw CLIP FASTQ into UMI-deduplicated, alignment-ready FASTQ while preserving the R2 5' end (= crosslink site -1) that drives single-nucleotide resolution downstream.

**Approach:** Use the protocol-matched UMI pattern (10 nt eCLIP, NNNXXXXNN iCLIP, 4 nt PAR-CLIP), run `umi_tools extract` to move random barcodes to read names, then apply cutadapt with 3'-only adapter trimming at `-q 6 -m 18` (permissive 5' to protect the truncation base). eCLIP uses two-pass trimming to remove read-through inline adapters from R2 5' only; iCLIP and PAR-CLIP use single-pass.

```bash
# eCLIP: 10 nt UMI on R1; two-pass adapter trim for read-through
# See clip-seq/clip-preprocessing for protocol-specific patterns
umi_tools extract \
    --bc-pattern=NNNNNNNNNN \
    --stdin=raw_R1.fq.gz --read2-in=raw_R2.fq.gz \
    --stdout=R1.umi.fq.gz --read2-out=R2.umi.fq.gz \
    --log=qc/umi_extract.log

# Pass 1: 3' adapter on both reads
# -q 6 is intentionally permissive; aggressive trimming destroys R2 5' = CL site -1
cutadapt \
    -a AGATCGGAAGAGCACACGTCT \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --quality-base 33 -q 6 -m 18 \
    -j 8 \
    -o R1.p1.fq.gz -p R2.p1.fq.gz \
    R1.umi.fq.gz R2.umi.fq.gz \
    > qc/cutadapt_pass1.log 2>&1

# Pass 2: strip read-through 5' adapter from R2 only (NEVER -g on R1)
cutadapt \
    -G GATCGTCGGACTGTAGAACTCTGAAC \
    --quality-base 33 -q 6 -m 18 \
    -j 8 \
    -o R1.trim.fq.gz -p R2.trim.fq.gz \
    R1.p1.fq.gz R2.p1.fq.gz \
    >> qc/cutadapt_pass2.log 2>&1
```

For PAR-CLIP: same UMI extraction but downstream alignment raises `--outFilterMismatchNoverReadLmax` from 0.04 to 0.07 (the T->C signature would otherwise be filtered as sequencing error). See clip-seq/clip-preprocessing for full per-protocol guidance.

## Step 3: Alignment (ENCODE STAR Block)

```bash
# ENCODE eCLIP convention. Sacred: --alignEndsType EndToEnd (soft-clip would destroy truncation = CL site -1)
STAR --runMode alignReads \
    --runThreadN 16 \
    --genomeDir /path/to/STAR_hg38_index \
    --genomeLoad NoSharedMemory \
    --readFilesIn R1.trim.fq.gz R2.trim.fq.gz \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --outFilterMultimapNmax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterScoreMinOverLread 0.66 \
    --outFilterMatchNminOverLread 0.66 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outFileNamePrefix sample_

samtools index sample_Aligned.sortedByCoord.out.bam

# MAPQ >= 10 (255 = unique in STAR; lower = multi-mapper)
samtools view -b -q 10 sample_Aligned.sortedByCoord.out.bam > sample_q10.bam
samtools index sample_q10.bam

# UMI dedup. ENCODE convention: --method=unique
umi_tools dedup \
    --stdin=sample_q10.bam \
    --stdout=sample_dedup.bam \
    --method=unique \
    --paired \
    --log=qc/dedup.log
samtools index sample_dedup.bam
```

For PAR-CLIP: change `--outFilterMismatchNoverReadLmax 0.04` to `0.07`. For repeat-binding RBPs (MATR3, ZFP36, FUS at LINE-1, HNRNPK at SINEs): change `--outFilterMultimapNmax 1` to `100` and add `--outSAMmultNmax -1`, then run CLAM downstream for EM-based multi-mapper assignment. See clip-seq/clip-alignment for full guidance.

## Step 4: QC (Five Gates)

```bash
# Gate 1: preprocessing retention (cutadapt log, target >= 70%)
grep -E "passing filters|Pairs written" qc/cutadapt_pass1.log

# Gate 2: alignment rate (STAR Log.final.out, target >= 60% eCLIP, 70% iCLIP)
grep "Uniquely mapped reads %" sample_Log.final.out

# Gate 3: library complexity (preseq, target >= 1M unique at sequenced depth)
preseq lc_extrap -B -P sample_q10.bam -o qc/preseq.txt

# Gate 4: FRiP (after peak calling; target >= 0.005 narrow-binding RBP)
# Gate 5: IDR replicate reproducibility (after peak calling; target rescue and self-consistency < 2)

# Aggregate all QC into a single MultiQC report
multiqc qc/ -o qc/multiqc/
```

CLIP libraries have 40-70% PCR duplication BY DESIGN (the IP enriches a small molecule pool). Low duplication usually means failed IP, not a good library. The unique-fragment count after UMI dedup is the actual quality metric. See clip-seq/clip-qc for full five-gate diagnostic.

## Step 5: Peak Calling

```bash
# CLIPper (ENCODE canonical) + SMInput log2 normalization
clipper \
    -b sample_dedup.bam \
    -s GRCh38 \
    -o peaks/sample.clipper.bed \
    --FDR 0.05 \
    --save-pickle \
    --processors 8   # super-local p-values are hard-coded ON in current CLIPper; the --superlocal flag was removed

# ENCODE stringent: log2(IP/SMInput) >= 3 AND -log10 p >= 3
# (Yeo lab eclip-pipeline scripts implement the normalization; see clip-seq/clip-peak-calling)
python overlap_peakfi_with_bam_PE.py \
    peaks/sample.clipper.bed \
    sample_dedup.bam sminput_dedup.bam \
    sample_dedup.bam.readnum.txt sminput_dedup.bam.readnum.txt \
    peaks/sample.normed.bed

python compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.py \
    peaks/sample.normed.bed \
    peaks/sample.compressed.bed

# Stringent filter
awk 'BEGIN{FS=OFS="\t"} $5 >= 3 && $6 >= 3' peaks/sample.compressed.bed > peaks/sample.stringent.bed
```

For maximum sensitivity (substantially more sites than CLIPper for mRNA-binding RBPs), use the Skipper Snakemake workflow with the same SMInput control. Mandatory for FASTKD2 / mt-RBPs which CLIPper misses on chrM. See clip-seq/clip-peak-calling for the full caller taxonomy.

## Step 6: Single-Nucleotide Crosslink-Site Detection

```bash
# PureCLIP: HMM jointly modeling enrichment + truncation + CL motif.
# -iv learns HMM parameters on a CHROMOSOME SUBSET (semicolon-delimited) to cut memory/runtime
# (per PureCLIP docs); it is NOT a BED. To limit the callset to expressed regions, pre-filter the input BAM.
pureclip \
    -i sample_dedup.bam -bai sample_dedup.bam.bai \
    -g genome.fa \
    -ibam sminput_dedup.bam -ibai sminput_dedup.bam.bai \
    -o crosslinks/sample.sites.bed \
    -or crosslinks/sample.regions.bed \
    -nt 8 -dm 8 \
    -iv 'chr1;chr2;chr3;'
```

Single-nt CL sites feed mCross motif registration and allele-specific binding analyses. They are NOT a replacement for the broad peak list; complementary outputs. See clip-seq/crosslink-site-detection.

## Step 7: IDR Across Replicates

```bash
# Sort each replicate's compressed BED by signal (log2 FC, column 5)
sort -k5,5gr peaks/rep1.compressed.bed > peaks/rep1.sorted.bed
sort -k5,5gr peaks/rep2.compressed.bed > peaks/rep2.sorted.bed

# True replicates threshold 0.05
idr --samples peaks/rep1.sorted.bed peaks/rep2.sorted.bed \
    --input-file-type bed --rank 5 \
    --output-file qc/idr.true.out \
    --idr-threshold 0.05 \
    --plot --log-output-file qc/idr.log

# ENCODE rule: rescue + self-consistency ratios both < 2 to pass
# Pseudo-replicate IDR (split BAM in half) at threshold 0.10
```

## Step 8: Binding-Site Annotation

```r
# CLIP-appropriate ChIPseeker (tssRegion tight; level=transcript)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile('peaks/sample.stringent.bed')
anno <- annotatePeak(
    peaks,
    TxDb = txdb,
    level = 'transcript',
    tssRegion = c(-100, 100),
    genomicAnnotationPriority = c('Promoter','5UTR','3UTR','Exon','Intron','Downstream','Intergenic')
)
plotAnnoPie(anno)
```

Default ChIPseeker `tssRegion=c(-3000, 3000)` over-extends for CLIP (would label 30-50% peaks as "Promoter"). Splicing factors additionally need RBP-Maps (Yeo lab) for the 1400 nt cassette-exon regulatory metagene. See clip-seq/binding-site-annotation.

## Step 9: Motif Analysis (De Novo + CL-Registered)

```bash
# Extract peak sequences (strand-preserving)
bedtools getfasta -fi genome.fa -bed peaks/sample.stringent.bed -s -fo motifs/peaks.fa

# GC-matched 3' UTR background (NOT auto-shuffled, which biases to AU)
bedtools shuffle -i peaks/sample.stringent.bed -g chrom.sizes \
    -incl expressed_3utr.bed -seed 42 > motifs/background.bed
bedtools getfasta -fi genome.fa -bed motifs/background.bed -s -fo motifs/background.fa

# HOMER de novo
findMotifs.pl motifs/peaks.fa fasta motifs/homer \
    -rna -len 5,6,7,8 -p 8 -fasta motifs/background.fa

# mCross for CL-position-registered motif. mCross.pl takes a POSITIONAL FASTA of sequences
# pre-extracted/registered around the CL sites and an output stem (not a BED + genome + -i/-g/-k/-o):
#   bedtools slop -i crosslinks/sample.sites.bed -g genome.sizes -b 10 | bedtools getfasta -fi genome.fa -bed - -s > motifs/peakseqs.fa
mCross.pl motifs/peakseqs.fa motifs/mcross   # see clip-seq/clip-motif-analysis for options
```

UV254 crosslinking has a strong U bias (~60-80% CL events at U); naive logos centered on CL positions are U-enriched even for non-U-binding RBPs. mCross corrects this by registering motif relative to the CL offset. See clip-seq/clip-motif-analysis.

## Step 10: Differential Binding (Optional, Across Conditions)

```r
# DEWSeq window-level NB with the interaction-term design
# The interaction `~ type + condition + type:condition` tests whether IP/SMInput ratio shifts;
# naive `~ condition` confounds binding with expression changes.
library(DEWSeq)
counts <- read.table('counts/merged.tsv', sep='\t', header=TRUE, row.names=1)
colData <- data.frame(
    type = relevel(factor(c('ip','ip','ip','ip','sminput','sminput','sminput','sminput')), ref='sminput'),
    condition = relevel(factor(c('treat','treat','ctrl','ctrl','treat','treat','ctrl','ctrl')), ref='ctrl')
)
dds <- DESeqDataSetFromSlidingWindows(
    countData=counts, colData=colData,
    annotObj='annotation.txt',   # htseq-clip TAB annotation table (named columns), NOT a plain BED
    design = ~ type + condition + type:condition
)
dds <- DESeq(dds)
# with sminput/ctrl as the references, the interaction coefficient is typeip.conditiontreat
res <- results(dds, name='typeip.conditiontreat')
```

See clip-seq/differential-clip for full DEWSeq workflow and the htseq-clip preprocessing required upstream.

## Quality Checkpoints

| Step | Metric | ENCODE target |
|------|--------|---------------|
| Preprocessing | Retention after adapter trim | >= 70% |
| Alignment | Unique mapping rate | >= 60% (eCLIP); >= 70% (iCLIP) |
| Complexity | preseq predicted unique at 100M reads | >= 10M (good); >= 1M (minimum acceptable) |
| Peak calling | FRiP (narrow-binding RBP) | >= 0.005 |
| Peak calling | Stringent peaks log2(IP/SMI) | >= 3 |
| Peak calling | Stringent peaks -log10 p | >= 3 |
| IDR | Rescue ratio | < 2 |
| IDR | Self-consistency ratio | < 2 |
| Annotation | Top RBP-class match expectation | Y (HuR -> 3' UTR; PTBP1 -> intron; FASTKD2 -> chrM) |

## Per-Variant Adjustments

- **PAR-CLIP:** Raise STAR `--outFilterMismatchNoverReadLmax` from 0.04 to 0.07; downstream use PARalyzer or CTK CIMS substitution T->C
- **iCLIP / iCLIP2 multiplexed:** Demultiplex by inline library barcode (NNNXXXXNN) BEFORE umi_tools extract
- **HITS-CLIP:** Use deletion-tolerant aligner (BWA-aln); downstream CTK CIMS deletion mode
- **Repeat-binding RBPs:** STAR `--outFilterMultimapNmax 100 --outSAMmultNmax -1` + CLAM EM rescue
- **m6A profiling:** Switch to clip-seq/m6a-clip (miCLIP2 + m6Aboost or GLORI)
- **Antibody unavailable:** Switch to clip-seq/stamp-antibody-free (STAMP or TRIBE)
- **miRNA targets:** Switch to clip-seq/ago-clip-mirna-targets (chimeric eCLIP / miR-eCLIP)
- **Variant-effect prediction:** Use clip-seq/clip-deep-learning (RBPNet or RNAProt)

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Peaks everywhere, dominated by abundant transcripts | No SMInput normalization | Normalize IP against SMInput; keep log2 FC >= 3 AND -log10 p >= 3 |
| Single-nucleotide resolution lost | Soft-clipping or aggressive 5' trim destroyed the R2 truncation base | STAR `--alignEndsType EndToEnd`; 3'-only `-q 6` trim; never `-g` on R1 |
| PAR-CLIP T->C signal missing | Mismatch ceiling 0.04 filtered the transitions as error | Raise `--outFilterMismatchNoverReadLmax` to 0.07 |
| "Low-complexity" library discarded | Judged on raw duplication (40-70% is normal for CLIP) | Use the unique-fragment count after UMI dedup as the quality metric |
| Motif logo is all-U even for a non-U-binding RBP | Naive CL-centered logo + UV U-bias | mCross CL-registered motif + GC-matched (not shuffled) background |
| 30-50% of peaks labeled "Promoter" | Default `tssRegion=c(-3000,3000)` over-extends for CLIP | Tight `tssRegion=c(-100,100)`, `level='transcript'` |
| Differential binding confounded with expression | `~ condition` design | `~ type + condition + type:condition` interaction (DEWSeq) |

## References

- Van Nostrand EL, Pratt GA, Shishkin AA, et al (2016) Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP). *Nature Methods* 13:508-514. DOI 10.1038/nmeth.3810.
- Van Nostrand EL, Freese P, Pratt GA, et al (2020) A large-scale binding and functional map of human RNA-binding proteins. *Nature* 583:711-719. DOI 10.1038/s41586-020-2077-3. (ENCODE RBP; SMInput + IDR practice.)
- Krakau S, Richard H, Marsico A (2017) PureCLIP: capturing target-specific protein-RNA interaction footprints from single-nucleotide CLIP-seq data. *Genome Biology* 18:240. DOI 10.1186/s13059-017-1364-2.
- Li Q, Brown JB, Huang H, Bickel PJ (2011) Measuring reproducibility of high-throughput experiments. *Annals of Applied Statistics* 5:1752-1779. DOI 10.1214/11-AOAS466. (IDR.)

## Related Skills

- clip-seq/clip-preprocessing - UMI extraction and adapter trimming details
- clip-seq/clip-alignment - STAR ENCODE block + multi-mapper rescue
- clip-seq/clip-qc - Five-gate QC framework
- clip-seq/clip-peak-calling - CLIPper / Skipper / PureCLIP / CTK taxonomy
- clip-seq/crosslink-site-detection - Single-nt CL detection by chemistry
- clip-seq/binding-site-annotation - ChIPseeker + RBP-Maps
- clip-seq/clip-motif-analysis - HOMER + mCross + RBNS validation
- clip-seq/differential-clip - DEWSeq + Flipper for cross-condition
- clip-seq/m6a-clip - miCLIP2 / GLORI / DART for m6A modifications
- clip-seq/stamp-antibody-free - STAMP / TRIBE for antibody-free profiling
- clip-seq/ago-clip-mirna-targets - chimeric eCLIP for direct miRNA-target pairs
- clip-seq/clip-deep-learning - RBPNet / RNAProt for variant-effect prediction
- read-qc/quality-reports - FastQC / MultiQC upstream QC
- reporting/automated-qc-reports - MultiQC aggregates the per-tool QC into one report; gating stays in the five-gate framework, not MultiQC
- alternative-splicing/differential-splicing - Cassette exon tables for RBP-Maps
