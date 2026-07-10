---
name: bio-workflows-smrna-pipeline
description: Orchestrates the end-to-end small RNA-seq pipeline from FASTQ to differential miRNAs and expression-filtered targets, chaining kit-aware cutadapt trimming (adapter on every read, UMI/4N handling), miRge3 known+isomiR quantification or miRDeep2 novel discovery, compositionally-aware DESeq2, and miRanda target prediction. Use when committing the library-kit adapter/UMI handling once, choosing the NORMALIZER (which drives which miRNAs are called DE more than the DE model does), deciding known quantification vs novel discovery, handling biofluid/plasma libraries that lack a trustworthy endogenous normalizer, routing tRF/piRNA reads to their own profiling, or feeding RAW (not RPM) counts with size-factor inspection into DE. Hands mechanism to the small-rna-seq component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: miRge3.0
workflow: true
depends_on:
  - small-rna-seq/smrna-preprocessing
  - small-rna-seq/mirge3-analysis
  - small-rna-seq/mirdeep2-analysis
  - small-rna-seq/differential-mirna
  - small-rna-seq/target-prediction
  - small-rna-seq/trf-pirna-profiling
qc_checkpoints:
  - after_trim: "Read-length peak 21-23 nt (30+ smear = degradation/tRNA/rRNA; 26-32 peak = piRNA)"
  - after_quant: "RNA-class composition checked (miRNA vs tRF/rRF/piRNA); abundant non-miRNA = different story"
  - before_de: "RAW counts confirmed (not RPM); size factors inspected for compositional distortion"
  - after_de: "baseMean reported with every call (significant FC on a ~5-count miRNA is noise)"
---

## Version Compatibility

Reference examples tested with: cutadapt 4.4+, miRge3.0 0.1.4+, miRDeep2 2.0.1.3+, DESeq2 1.42+, apeglm 1.24+, miRanda 3.3a, umi_tools 1.1+, miRTrace 1.0+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: the correct 3' adapter sequence and the UMI/4N scheme are KIT-specific (TruSeq vs NEXTflex vs QIAseq) — confirm against the kit before trimming. The normalizer that best fits compositional miRNA data drifts with the method literature; treat the table below as a decision aid to validate, not a fixed rule.

# Small RNA-seq Pipeline

**"Analyze my small RNA-seq data from FASTQ to differential miRNAs"** -> Chain kit-aware trimming, known-miRNA quantification (or novel discovery), a compositionally-aware DE test, and expression-filtered target prediction.
- CLI + R: cutadapt -> (miRge3.0 | miRDeep2) -> DESeq2 (inspect size factors) -> miRanda (filter by anti-correlated mRNA)

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A small RNA-seq result turns on two things the mRNA pipeline never faces — the adapter is on every read, and the counts are compositional — so the trustworthy result is decided at these seams.

1. **The library kit fixes adapter + UMI/4N handling, committed once at trimming and un-reversible.** The insert (~22 nt) is far shorter than the read, so the 3' adapter is on EVERY real read and `--discard-untrimmed` correctly drops no-insert junk (the inverse of genomic DNA). NEXTflex 4N random bases are stripped AFTER adapter removal; QIAseq UMIs are extracted and deduped. NEVER position-dedup a non-UMI small-RNA library — abundant miRNAs legitimately share start positions, so position-dedup destroys real signal.
2. **The normalizer is the dominant analytic commitment — more than the DE model.** miRNA counts are few-featured and compositional; a handful of miRNAs can dominate the library. Global scaling (TMM's 30%/5% trimming, median-of-ratios) is built for ~20k mRNAs and is harmful among hundreds of miRNAs (Garmire 2012). The normalizer choice CHANGES which miRNAs are called DE. Commit and report it; do not default silently.
3. **RAW counts (not RPM) flow into DESeq2, and size factors must be inspected.** A large rise in one abundant miRNA mechanically deflates all others (the closed-composition artifact), manufacturing spurious "down" calls. Confirm the input is raw integer counts and inspect `sizeFactors(dds)` before trusting any call.
4. **Targets are hypotheses, not findings.** miRanda output must be intersected with anti-correlated mRNA DE and validated databases before interpretation.

## Pipeline map

```
FASTQ (single-end small RNA)
  | [1] Trim (kit-aware) --------> cutadapt (+ umi_tools / 4N)   (small-rna-seq/smrna-preprocessing)
  v     ^-- commitment: kit adapter + UMI/4N; --discard-untrimmed; NO position-dedup w/o UMI
  | [2] Quantify OR discover ----> miRge3.0 (known+isomiR) | miRDeep2 (novel)  (small-rna-seq/mirge3-analysis, mirdeep2-analysis)
  v
  | [3] Differential expression -> DESeq2 on RAW counts; INSPECT size factors   (small-rna-seq/differential-mirna)
  v     ^-- commitment: the NORMALIZER (drives which miRNAs are DE)
  | [4] Target prediction -------> miRanda, filtered by anti-correlated mRNA DE  (small-rna-seq/target-prediction)
  v
Differential miRNAs + expression-supported targets
   (tRF/piRNA reads -> small-rna-seq/trf-pirna-profiling)
```

## Made-once commitments

| Commitment | Choice | Consequence inherited downstream |
|------------|--------|----------------------------------|
| Library kit adapter + UMI/4N | TruSeq / NEXTflex-4N / QIAseq-UMI handling, set at trimming | Wrong adapter or missed 4N/UMI corrupts every count; non-UMI position-dedup destroys real miRNAs |
| Quantification target | Known-miRNA + isomiR (miRge3) vs novel discovery (miRDeep2) | Discovery is high-FP and needs a genome + bowtie1; quantification is the common curated case |
| Normalizer | median-of-ratios/TMM vs quantile/loess vs RUVg vs CoDA/ALDEx2 | Which miRNAs are DE (more than the DE model choice) |
| Biofluid vs tissue | plasma/serum has NO trustworthy endogenous normalizer | DE may be uninterpretable without spike-ins + hemolysis modeling |

## The canonical order and why

1. **Trim (kit-aware)** — adapter removal with `--discard-untrimmed`; 4N/UMI handling before or during collapse. Order-trap: position-deduping a non-UMI library removes real high-abundance miRNAs.
2. **miRTrace QC** — length/complexity, RNA-class composition, cross-clade contamination BEFORE quantification (a good mapping rate hides sample swaps and reagent contamination).
3. **Quantify (miRge3) or discover (miRDeep2)** — to raw integer counts. For miRDeep2, choose the score cutoff from `survey.pl` signal-to-noise, not a fixed rule.
4. **Low-count filter** (`rowSums >= 10`, lower than mRNA) BEFORE normalization.
5. **DE with the committed normalizer** — raw counts into DESeq2; inspect `sizeFactors`. Order-trap: feeding RPM/CPM bakes compositional distortion into an invalid model.
6. **Target prediction** — intersect miRanda predictions with anti-correlated mRNA DE (order-trap: enrichment on raw predictions is false-positive-laden).

## Fork/decision points

Pipeline-level selection only; mechanism lives in the component skills.

| Fork | Lean toward | Hand off to |
|------|-------------|-------------|
| Quantify vs discover | miRge3.0 (known + isomiR, curated, common) vs miRDeep2 (novel only; high FP, genome + bowtie1, survey.pl cutoff) | small-rna-seq/mirge3-analysis, small-rna-seq/mirdeep2-analysis |
| Normalizer | median-of-ratios/TMM (risky under composition shift) vs quantile/loess (better on skewed miRNA) vs RUVg (control/empirical miRNAs) vs CoDA/ALDEx2 (compositional sensitivity check) | small-rna-seq/differential-mirna |
| Biofluid/plasma | no trustworthy endogenous normalizer (miR-16 is hemolysis-sensitive, U6 degrades); cel-miR-39 spike controls EXTRACTION not biological scale; model batch/hemolysis explicitly | small-rna-seq/differential-mirna |
| RNA class | miRNA vs tRF/piRNA (26-32 nt peak) -> route to trf-pirna-profiling (MINTmap/unitas) | small-rna-seq/trf-pirna-profiling |

## Primary path: cutadapt -> miRge3.0 -> DESeq2 -> miRanda

**Goal:** turn kit-specific FASTQ into differential miRNAs with expression-supported targets.

**Approach:** trim to the kit, quantify known miRNAs/isomiRs, test RAW counts with size-factor inspection, then filter targets by anti-correlation. (The runnable script in examples/ shows the miRDeep2 novel-discovery route below; the miRge3.0 commands are shown here.)

```bash
# Kit-aware trim: adapter on EVERY read, so --discard-untrimmed drops no-insert junk (inverse of gDNA).
# NEXTflex 4N: cutadapt -u 4 -u -4 AFTER adapter. QIAseq UMIs: umi_tools. NEVER position-dedup non-UMI.
cutadapt -a TGGAATTCTCGGGTGCCAAGG --minimum-length 18 --maximum-length 30 --discard-untrimmed \
    -o trimmed.fastq.gz reads.fastq.gz

# Known-miRNA + isomiR quantification (the common case). NOTE: v3 has NO 'annotate' subcommand
# (that was miRge2); it is `miRge3.0 -s ...`. Reads are already cutadapt-trimmed, so -a is OMITTED
# (v3 has no 'none' keyword; passing one makes cutadapt treat it as an adapter sequence and fail).
# isomiRs are annotated by default; -gff emits the isomiR GFF.
# (-ai would instead compute A-to-I editing, a separate analysis, not isomiR reporting.)
miRge3.0 -s trimmed.fastq.gz -lib /path/to/miRge3_Lib -on human -db mirbase \
    -gff -cpu 8 -o mirge_out
```

```r
library(DESeq2)
# RAW counts (miR.Counts.csv), NOT RPM. A few miRNAs can dominate, so inspect sizeFactors first.
# miRge3.0 writes into a timestamped subfolder (mirge_out/miRge.YYYY-M-D_h-m-s/), not directly into -o
counts <- read.csv(Sys.glob('mirge_out/miRge.*/miR.Counts.csv')[1], row.names = 1)
dds <- DESeqDataSetFromMatrix(round(counts), colData, ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]      # lower prefilter than mRNA
dds <- DESeq(dds)
print(sizeFactors(dds))                       # a dominant-miRNA shift distorts these -> spurious calls
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
```

```bash
# Targets are a hypothesis: intersect with anti-correlated mRNA DE before trusting them.
# -sc 140: miRanda default alignment-score floor; -en -20: keep duplexes with MFE <= -20 kcal/mol (stable binding)
miranda mature_mirnas.fa target_3utrs.fa -sc 140 -en -20 -strict -out targets.txt
```

## Novel discovery alternative: miRDeep2

**Goal:** discover previously unannotated miRNAs (only when that is the question).

**Approach:** collapse reads, map to the genome with bowtie1, run miRDeep2, and pick the score cutoff from the signal-to-noise survey — not a fixed threshold. Full runnable script: `examples/smrna_full_pipeline.sh`.

```bash
gunzip -kc trimmed.fastq.gz > trimmed.fastq            # mapper.pl reads plain-text FASTQ, not gzip
mapper.pl trimmed.fastq -e -h -i -j -l 18 -m -p genome_index \
    -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf
miRDeep2.pl reads_collapsed.fa genome.fa reads_collapsed_vs_genome.arf \
    mature_ref.fa none hairpin_ref.fa -t Human      # score cutoff from survey.pl signal-to-noise
```

## QC checkpoints between steps

| After | Gate | Interpretation |
|-------|------|----------------|
| Trim | Read-length peak 21-23 nt | 30+ nt smear = degradation/tRNA/rRNA; 26-32 nt peak = piRNA (route to trf-pirna-profiling) |
| miRTrace/align | RNA-class composition (miRNA vs tRF/rRF/piRNA); cross-clade contamination | Abundant non-miRNA classes mean this is a tRF/piRNA story, not a miRNA one |
| Before DE | RAW counts (not RPM); size factors inspected | A dominant-miRNA shift distorts size factors and manufactures spurious "down" calls |
| After DE | baseMean reported with every call | A significant fold-change on a ~5-count miRNA is noise |
| Targets | filtered by anti-correlated mRNA DE + validated DBs | Raw miRanda predictions are hypotheses, not findings |

Ligation bias means absolute cross-miRNA abundance WITHIN a sample is untrustworthy; compare the same miRNA across samples, never across kits. For a fully reproducible run, nf-core/smrnaseq chains these steps (workflow-management/nf-core-pipelines).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Real high-abundance miRNAs vanish | Position-deduped a non-UMI library | Dedup ONLY with UMIs (umi_tools); non-UMI libraries are not position-deduped |
| Invalid model / distorted calls | Fed RPM/CPM to DESeq2 | Use RAW integer counts; RPM bakes in composition |
| Many spurious "down" miRNAs | One abundant miRNA rose and deflated the rest (closed composition) | Inspect size factors; consider quantile/RUVg/ALDEx2 as a sensitivity check |
| DE looks strong but is noise | Fold-change on a ~5-count miRNA | Report and gate on baseMean |
| Plasma miRNA DE uninterpretable | No trustworthy endogenous normalizer; hemolysis confound | Spike-ins for extraction control + model hemolysis/batch explicitly |
| "miRNA" library is mostly tRF/piRNA | 26-32 nt peak / broad tRF distribution | Route to trf-pirna-profiling (MINTmap/unitas) |

## Pipeline map (hand-offs)

- small-rna-seq/smrna-preprocessing - kit-specific adapter, UMI, and 4N handling
- small-rna-seq/mirge3-analysis - known-miRNA and isomiR quantification
- small-rna-seq/mirdeep2-analysis - novel miRNA discovery and survey.pl cutoff
- small-rna-seq/differential-mirna - compositionally-aware DE and normalizer selection
- small-rna-seq/target-prediction - seed prediction filtered by expression
- small-rna-seq/trf-pirna-profiling - tRF and piRNA profiling

The runnable miRDeep2 discovery-route script is in this skill's examples/ (`smrna_full_pipeline.sh`).

## Related Skills

- small-rna-seq/smrna-preprocessing - Kit-specific adapter, UMI, and 4N handling
- small-rna-seq/mirdeep2-analysis - Novel miRNA discovery
- small-rna-seq/mirge3-analysis - Known-miRNA and isomiR quantification
- small-rna-seq/differential-mirna - Compositionally-aware DE
- small-rna-seq/target-prediction - Seed prediction filtered by expression
- small-rna-seq/trf-pirna-profiling - tRF and piRNA profiling
- differential-expression/deseq2-basics - DESeq2 model, contrasts, shrinkage
- workflow-management/nf-core-pipelines - Run nf-core/smrnaseq as a curated, reproducible pipeline

## References

- Garmire LX, Subramaniam S (2012) Evaluation of normalization methods in mammalian microRNA-Seq data. *RNA* 18:1279-1288. DOI 10.1261/rna.030916.111. (normalization, not the DE model, drives miRNA DE results.)
- Risso D, Ngai J, Speed TP, Dudoit S (2014) Normalization of RNA-seq data using factor analysis of control genes or samples. *Nature Biotechnology* 32:896-902. DOI 10.1038/nbt.2931. (RUVg for miRNA/unwanted variation.)
- Friedländer MR, Mackowiak SD, Li N, Chen W, Rajewsky N (2012) miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. *Nucleic Acids Research* 40:37-52. DOI 10.1093/nar/gkr688.
