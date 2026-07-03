---
name: bio-workflows-merip-pipeline
description: Orchestrates an end-to-end MeRIP-seq / m6A-seq analysis from raw FASTQ to differential m6A peak calls and metagene plots, chaining fastp adapter trimming, STAR splice-aware alignment (NO deduplication for non-UMI MeRIP), deepTools replicate-concordance + IP-enrichment QC, PreSeq saturation curves, exomePeak2 (transcript-aware Poisson GLM) peak calling, optional MACS3 broad-peak cross-check, DRACH motif confirmation as a sanity check (NOT a per-peak filter), exomePeak2 differential calling via the four-BAM-vector interface (bam_ip + bam_input control; bam_treated_ip + bam_treated_input treatment), ChIPseeker annotation, and the canonical Guitar metagene with stop-codon enrichment as the biological QC anchor. Use when running a complete MeRIP analysis from raw reads, when chaining the constituent epitranscriptomics skills (merip-preprocessing -> m6a-peak-calling -> m6a-differential -> modification-visualization), or when wrapping the pipeline in Snakemake / Nextflow.
tool_type: mixed
primary_tool: exomePeak2
---

## Version Compatibility

Reference examples tested with: STAR 2.7.11+, samtools 1.19+, fastp 0.23+, deepTools 3.5+, PreSeq 3.2+, exomePeak2 1.14+ (Bioconductor 3.18+), MACS3 3.0+, ChIPseeker 1.38+, Guitar 2.18+, BSgenome.Hsapiens.UCSC.hg38 1.4+, TxDb.Hsapiens.UCSC.hg38.knownGene 3.18+, HOMER 4.11+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

exomePeak2 has NO `mode=` or `experiment_design=` argument; differential is triggered by populating `bam_treated_ip` + `bam_treated_input`. MeTPeak defaults are `WINDOW_WIDTH=50, SLIDING_STEP=50, FRAGMENT_LENGTH=100`. MACS3 default `--keep-dup` is 1 and MUST be overridden to `all` for non-UMI MeRIP. Guitar `txTxdb=` is the modern argument name (older releases used `txdb=`).

# MeRIP-seq End-to-End Pipeline

**"Analyze my MeRIP-seq data from FASTQ to differential m6A peaks"** -> Orchestrate read alignment (STAR splice-aware to GENOME), IP-enrichment QC (deepTools plotFingerprint, replicate Spearman, PreSeq saturation), m6A peak calling (exomePeak2 transcript-aware default, MACS3 broad as cross-check), DRACH motif sanity check (HOMER), exomePeak2 differential via the four-BAM-vector interface, ChIPseeker feature annotation, and Guitar transcript-feature metagene confirming canonical stop-codon enrichment. Defer per-skill deep treatment to `epitranscriptomics/merip-preprocessing`, `epitranscriptomics/m6a-peak-calling`, `epitranscriptomics/m6a-differential`, and `epitranscriptomics/modification-visualization`.

## Pipeline Overview

```
FASTQ -> fastp trim -> STAR genome align -> samtools sort/index -> deepTools QC + PreSeq saturation
       -> exomePeak2 peak calling (+ MeTPeak / MACS3 cross-check)
       -> HOMER DRACH sanity check
       -> exomePeak2 differential (bam_ip + bam_treated_ip)
       -> ChIPseeker feature annotation
       -> Guitar metagene (stop-codon QC anchor) + pyGenomeTracks browser figures
```

## Step 1: Adapter Trimming

```bash
fastp \
    --in1 raw/IP_R1.fastq.gz --in2 raw/IP_R2.fastq.gz \
    --out1 trimmed/IP_R1.fq.gz --out2 trimmed/IP_R2.fq.gz \
    --json qc/IP_fastp.json --html qc/IP_fastp.html \
    --length_required 25 --detect_adapter_for_pe --thread 8

fastp \
    --in1 raw/Input_R1.fastq.gz --in2 raw/Input_R2.fastq.gz \
    --out1 trimmed/Input_R1.fq.gz --out2 trimmed/Input_R2.fq.gz \
    --json qc/Input_fastp.json --html qc/Input_fastp.html \
    --length_required 25 --detect_adapter_for_pe --thread 8
```

Standard non-UMI MeRIP: do NOT pass `--umi`. See `epitranscriptomics/merip-preprocessing` for the do-NOT-dedup rationale.

## Step 2: STAR Splice-Aware Genome Alignment

```bash
STAR --runMode alignReads \
    --genomeDir refs/star_index \
    --readFilesIn trimmed/IP_R1.fq.gz trimmed/IP_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 20 \
    --outSAMattributes NH HI AS nM NM MD \
    --outFileNamePrefix aligned/IP_ \
    --runThreadN 12

samtools index aligned/IP_Aligned.sortedByCoord.out.bam
```

Repeat for Input. Align to GENOME (not transcriptome) for downstream MeRIP peak calling. Do NOT deduplicate (no UMI in standard MeRIP).

## Step 3: IP-Enrichment + Replicate-Concordance QC

```bash
multiBamSummary bins \
    --bamfiles aligned/IP_rep*.bam aligned/Input_rep*.bam \
    --binSize 10000 --numberOfProcessors 8 \
    -o qc/cov.npz

plotCorrelation --corData qc/cov.npz --corMethod spearman --skipZeros \
    --whatToPlot heatmap --colorMap RdYlBu_r --plotNumbers \
    -o qc/replicate_correlation.pdf

plotFingerprint \
    --bamfiles aligned/IP_rep*.bam aligned/Input_rep*.bam \
    --skipZeros --numberOfProcessors 8 \
    --outQualityMetrics qc/fingerprint_metrics.tab \
    -o qc/fingerprint.pdf

preseq lc_extrap -B -o qc/IP_rep1_lc_extrap.txt aligned/IP_rep1_Aligned.sortedByCoord.out.bam
```

For peak-count comparison across conditions, rarefy BAMs to a common unique-read depth informed by the saturation curve before calling peaks.

## Step 4: exomePeak2 Peak Calling (Per-Condition)

**Goal:** Produce a transcript-aware set of m6A peaks with FDR and IP/input fold-change from paired IP/Input genome BAM files, suitable as input to differential analysis, motif scanning, or downstream visualisation.

**Approach:** Build a TxDb from the matched GTF; pass paired IP/Input BAM vectors to `exomePeak2()` with `txdb` and `genome` (BSgenome) for GC correction; export BED12 + RDS to `save_dir/experiment_name/`.

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

peaks <- result
length(peaks)
```

`exomePeak2()` writes BED12 + RDS + per-peak fold-change / FDR to `save_dir/experiment_name/`.

## Step 5: MACS3 Broad-Peak Cross-Check (Optional)

```bash
macs3 callpeak \
    --treatment aligned/IP_rep*.bam \
    --control aligned/Input_rep*.bam \
    --format BAMPE --gsize hs \
    --nomodel --extsize 150 \
    --keep-dup all \
    --broad --broad-cutoff 0.1 --qvalue 0.05 \
    --outdir macs3_output --name m6a_run1
```

`--keep-dup all` is non-negotiable for non-UMI MeRIP (default `--keep-dup 1` destroys signal at high-coverage transcripts).

## Step 6: DRACH Motif Sanity Check

```bash
findMotifsGenome.pl \
    exomepeak2_output/m6a_run1/peaks.bed \
    hg38 motif_output \
    -rna -size 100 -len 5,6 -p 8
```

Report DRACH enrichment on the peak set as a sanity check (E-value < 1e-50 expected). NEVER post-hoc filter individual peaks by DRACH.

## Step 7: exomePeak2 Differential (Control vs Treatment)

**Goal:** Identify m6A peaks that change between control and treatment conditions, with per-peak log2FC + FDR, using exomePeak2's integrated peak-calling + differential interface.

**Approach:** Populate `bam_ip` + `bam_input` with the control arm and `bam_treated_ip` + `bam_treated_input` with the treatment arm; populating the treated arms triggers differential mode (there is NO `mode=` argument). Apply effect-size + FDR filters downstream.

```r
library(exomePeak2)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb <- makeTxDbFromGFF('refs/annotation.gtf', format='gtf')

ctrl_ip     <- c('aligned/ctrl_IP1.bam', 'aligned/ctrl_IP2.bam', 'aligned/ctrl_IP3.bam')
ctrl_input  <- c('aligned/ctrl_Input1.bam', 'aligned/ctrl_Input2.bam', 'aligned/ctrl_Input3.bam')
treat_ip    <- c('aligned/treat_IP1.bam', 'aligned/treat_IP2.bam', 'aligned/treat_IP3.bam')
treat_input <- c('aligned/treat_Input1.bam', 'aligned/treat_Input2.bam', 'aligned/treat_Input3.bam')

diff_result <- exomePeak2(
    bam_ip            = ctrl_ip,
    bam_input         = ctrl_input,
    bam_treated_ip    = treat_ip,
    bam_treated_input = treat_input,
    txdb              = txdb,
    genome            = BSgenome.Hsapiens.UCSC.hg38,
    paired_end        = TRUE,
    library_type      = 'unstranded',
    peak_calling_mode = 'exon',
    save_dir          = 'exomepeak2_diff_output',
    experiment_name   = 'ctrl_vs_treat'
)

diff_table <- as.data.frame(diff_result)
sig <- diff_table[diff_table$padj < 0.05 & abs(diff_table$log2FC) > 0.5, ]
nrow(sig)
```

exomePeak2 has NO `mode=` or `experiment_design=` argument. Populating `bam_treated_ip` + `bam_treated_input` triggers differential output. For batch / antibody-lot covariate adjustment, fall through to featureCounts-on-peaks -> DESeq2 (see `epitranscriptomics/m6a-differential`).

## Step 8: Peak Annotation to Transcript Features

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)

peaks <- import('exomepeak2_output/m6a_run1/peaks.bed')
anno <- annotatePeak(peaks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level='transcript')
plotAnnoBar(anno)
plotDistToTSS(anno)
```

Flag peaks within ~50 nt of TSS as m6A-or-m6Am ambiguous (antibody cross-reactivity with PCIF1-deposited cap m6Am).

## Step 9: Guitar Metagene (Biological QC Anchor)

```r
library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

GuitarPlot(
    txTxdb          = TxDb.Hsapiens.UCSC.hg38.knownGene,
    stBedFiles      = list('exomepeak2_output/m6a_run1/peaks.bed'),
    saveToPDFprefix = 'figures/m6a_metagene'
)
```

Expected pattern: peak density rises toward and peaks near the stop codon (3'UTR-proximal end of CDS). If absent, suspect IP failure or wrong antibody; do NOT proceed to downstream interpretation.

## Complete Bash Driver

```bash
#!/usr/bin/env bash
set -euo pipefail

STAR_INDEX=$1
GTF=$2
IP_R1=$3
IP_R2=$4
INPUT_R1=$5
INPUT_R2=$6
OUTPUT_DIR=$7

mkdir -p "${OUTPUT_DIR}"/{qc,trimmed,aligned,peaks,figures}

fastp --in1 "${IP_R1}" --in2 "${IP_R2}" \
    --out1 "${OUTPUT_DIR}/trimmed/IP_R1.fq.gz" --out2 "${OUTPUT_DIR}/trimmed/IP_R2.fq.gz" \
    --json "${OUTPUT_DIR}/qc/IP_fastp.json" --length_required 25 --detect_adapter_for_pe --thread 8

fastp --in1 "${INPUT_R1}" --in2 "${INPUT_R2}" \
    --out1 "${OUTPUT_DIR}/trimmed/Input_R1.fq.gz" --out2 "${OUTPUT_DIR}/trimmed/Input_R2.fq.gz" \
    --json "${OUTPUT_DIR}/qc/Input_fastp.json" --length_required 25 --detect_adapter_for_pe --thread 8

for sample in IP Input; do
    STAR --runMode alignReads --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${OUTPUT_DIR}/trimmed/${sample}_R1.fq.gz" "${OUTPUT_DIR}/trimmed/${sample}_R2.fq.gz" \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 20 \
        --outFileNamePrefix "${OUTPUT_DIR}/aligned/${sample}_" --runThreadN 12
    samtools index "${OUTPUT_DIR}/aligned/${sample}_Aligned.sortedByCoord.out.bam"
done

macs3 callpeak \
    --treatment "${OUTPUT_DIR}/aligned/IP_Aligned.sortedByCoord.out.bam" \
    --control "${OUTPUT_DIR}/aligned/Input_Aligned.sortedByCoord.out.bam" \
    --format BAMPE --gsize hs --nomodel --extsize 150 --keep-dup all \
    --broad --broad-cutoff 0.1 --qvalue 0.05 \
    --outdir "${OUTPUT_DIR}/peaks" --name m6a
```

The full pipeline (incl. exomePeak2 peak calling, DRACH check, ChIPseeker annotation, Guitar metagene) is best orchestrated in Snakemake or Nextflow with the per-skill recipes from the four `epitranscriptomics/` skills.

## QC Checkpoints

| Checkpoint | Expected | Action if Failed |
|------------|----------|------------------|
| Properly-paired rate (samtools flagstat) | >=85% | Check trimming and adapter contamination |
| Replicate Spearman within condition (10 kb bins) | >=0.85 IP-IP | Inspect divergent replicate; consider exclusion |
| plotFingerprint IP-vs-input JS distance | >=0.5 | Suspect failed IP if lower |
| Saturation plateau depth | ~30-60M unique reads | Sequence deeper if not plateaued |
| DRACH motif enrichment (HOMER, peak set) | E-value < 1e-50 | Suspect IP failure or wrong antibody |
| Stop-codon enrichment in Guitar metagene | Clear 3'UTR-proximal peak | Suspect IP failure, wrong antibody, or non-m6A modification |
| 5'UTR peaks fraction | Note ambiguity zone (~50 nt of TSS) | Flag as m6A-or-m6Am ambiguous; PCIF1 cross-reactivity |

## Output Files

| File | Description |
|------|-------------|
| `exomepeak2_output/m6a_run1/peaks.bed` | exomePeak2 peak BED12 |
| `exomepeak2_diff_output/ctrl_vs_treat/...` | Differential peaks with log2FC + FDR |
| `motif_output/` | HOMER DRACH motif enrichment report |
| `figures/m6a_metagene.pdf` | Guitar transcript-feature metagene (stop-codon QC anchor) |
| `qc/replicate_correlation.pdf` | deepTools Spearman heatmap |
| `qc/fingerprint.pdf` | deepTools Lorenz IP-enrichment plot |
| `qc/IP_rep*_lc_extrap.txt` | PreSeq saturation curves |

## Related Skills

- epitranscriptomics/merip-preprocessing - Per-step preprocessing (trim, align, QC, saturation, IP-over-Input bigWig)
- epitranscriptomics/m6a-peak-calling - exomePeak2 / MeTPeak / MACS3 deep treatment, DRACH sanity check, m6A-vs-m6Am 5'UTR flag
- epitranscriptomics/m6a-differential - Differential methods (exomePeak2, QNB, RADAR), batch / lot covariate handling, stoichiometry-vs-expression confound
- epitranscriptomics/modification-visualization - Guitar metagene, peak-centred heatmaps, pyGenomeTracks browser figures
- epitranscriptomics/m6anet-analysis - ONT direct-RNA alternative for orthogonal stoichiometry validation
- chip-seq/peak-calling - Sibling IP-vs-input peak-calling framework
- chip-seq/chipseq-qc - IP enrichment QC concepts that transfer to MeRIP
- read-alignment/star-alignment - General STAR splice-aware alignment
- workflow-management/snakemake-workflows - Snakemake orchestration patterns
- workflow-management/nextflow-pipelines - Nextflow orchestration patterns
- workflows/rnaseq-to-de - General RNA-seq -> DE pipeline patterns
