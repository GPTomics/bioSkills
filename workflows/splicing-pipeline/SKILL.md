---
name: bio-workflows-splicing-pipeline
description: Orchestrates the end-to-end bulk short-read alternative-splicing pipeline from FASTQ to differential splicing, chaining fastp QC, cohort-consistent STAR 2-pass alignment (one shared junction DB), junction QC, event-level differential splicing (rMATS-turbo + leafcutter, optional MAJIQ V3), parallel isoform-level DTU (Salmon -> tximport dtuScaledTPM -> DRIMSeq/DEXSeq -> stageR), and sashimi visualization. Use when committing the annotation GTF and a shared 2-pass junction database for the whole cohort, keeping the analysis at splice-aware resolution (never collapsing to gene), choosing event-level vs isoform-level DTU and reconciling them, applying the stageR two-stage gene->transcript FDR, or off-ramping to splice-variant / outlier / long-read / single-cell splicing. Hands mechanism to the alternative-splicing component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: rMATS-turbo
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - read-alignment/star-alignment
  - alternative-splicing/splicing-qc
  - alternative-splicing/splicing-quantification
  - alternative-splicing/differential-splicing
  - rna-quantification/alignment-free-quant
  - rna-quantification/tximport-workflow
  - alternative-splicing/isoform-switching
  - alternative-splicing/sashimi-plots
qc_checkpoints:
  - after_qc: "Q30 >80%, adapter <5%; reads NOT over-trimmed (short reads lose junction-spanning power)"
  - after_align: "Uniquely-mapped >80%; junction-saturation curves plateau (else sequence deeper)"
  - after_diff: "|deltaPSI| >0.1 (lenient) / >0.2 (stringent), FDR <0.05, >=10 junction reads supporting the event"
  - after_dtu: "stageR gene-level screen passed BEFORE trusting any transcript-level q-value"
---

## Version Compatibility

Reference examples tested with: STAR 2.7.11+, fastp 0.23+, rMATS-turbo 4.3+, leafcutter 0.2.9+, Salmon 1.10+, tximport 1.30+, DRIMSeq 1.30+, DEXSeq 1.48+, stageR 1.24+, IsoformSwitchAnalyzeR 2.2+, RSeQC 5.0+, ggsashimi 1.1+ (numpy 1.26+, pandas 2.2+)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: rMATS-turbo `--readLength` must match the trimmed read length (add `--variable-read-length` if trimming produced a range); IsoformSwitchAnalyzeR's `importRdata` argument is genuinely spelled `isoformExonAnnoation` (a package typo, not an error here). Confirm in-tool before quoting.

# Alternative Splicing Analysis Pipeline

**"Find differential alternative splicing between my two conditions"** -> Chain QC/trim, cohort-consistent 2-pass alignment, junction QC, event-level and (parallel) isoform-level differential testing, and sashimi visualization.
- CLI + R: fastp -> STAR 2-pass (shared SJ DB) -> junction QC -> rMATS-turbo + leafcutter -> [Salmon -> tximport dtuScaledTPM -> DRIMSeq/DEXSeq -> stageR] -> ggsashimi

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

Splicing quantification is comparative by construction, so the trustworthy result is decided at four seams, not inside any one caller.

1. **Splice-aware alignment must be cohort-consistent — one shared junction database for the whole experiment.** STAR 2-pass where pass 2 uses the COMBINED `SJ.out.tab` from all samples is a made-once commitment, the splicing analogue of a shared reference. Junctions discovered per-sample and applied unevenly make PSI values non-comparable across samples, so a "differential" event is really a coverage artifact. Never run 2-pass independently per sample for a cohort.
2. **The analysis lives at splice-aware resolution and does NOT collapse to gene.** Feeding the standard gene-level tximport import (`txOut=FALSE`) averages the isoform-switch signal away. The DTU branch uses the OPPOSITE import from workflows/rnaseq-to-de: `countsFromAbundance="dtuScaledTPM"` + `txOut=TRUE`.
3. **Event-level and isoform-level tests answer different questions and are reconciled, not merged.** rMATS/leafcutter test exon-inclusion (ΔPSI); DTU tests within-gene isoform-proportion shifts. A gene can be significant for one and not the other; treating their p-values as interchangeable is a category error.
4. **Transcript-level DTU FDR is only honest through the stageR two-stage test.** Report per-transcript q-values only after a gene-level screen has passed; raw transcript-level FDR is inflated. stageR is the seam that converts gene-screened-then-transcript-confirmed tests into calibrated FDR.

## Pipeline map

```
FASTQ (paired)
  | [1] QC & trim ---------------> fastp                (read-qc/fastp-workflow)
  v
  | [2] STAR 2-pass -------------> pass1 (all) -> COMBINED SJ.out.tab -> pass2 (all)  (read-alignment/star-alignment)
  v     ^-- commitment: annotation GTF + ONE shared junction DB + readLength
  | [3] Junction QC -------------> saturation plateau, entropy   (alternative-splicing/splicing-qc)
  v
  +---------------------------+-----------------------------------+
  | EVENT level               | ISOFORM level (parallel, optional)|
  v                           v
[4a] rMATS-turbo + leafcutter [4b] Salmon -> tximport dtuScaledTPM(txOut=TRUE)
     (alternative-splicing/       -> DRIMSeq/DEXSeq -> stageR
      differential-splicing)      (alternative-splicing/isoform-switching)
  |   ^-- reconcile, don't merge      ^-- gene screen BEFORE transcript q
  v
[5] Sashimi on top events -----> ggsashimi   (alternative-splicing/sashimi-plots)
```

## Made-once commitments

Decided before the first alignment; every downstream PSI/DTU value inherits them.

| Commitment | Choice | Consequence inherited downstream |
|------------|--------|----------------------------------|
| Annotation GTF | ONE GTF used by rMATS, leafcutter, and IsoformSwitchAnalyzeR | Different GTFs make their events irreconcilable; fixes the event universe |
| 2-pass junction DB | Combined `SJ.out.tab` from ALL samples fed into pass 2 | Per-sample junctions -> non-comparable PSI, false "differential" events from coverage |
| Measurement level | Splice-aware (event PSI and/or transcript DTU); NEVER gene collapse | Gene-level tximport (`txOut=FALSE`) averages away the switch signal |
| `--readLength` | Matches the trimmed read length (or `--variable-read-length`) | A wrong value miscomputes inclusion-junction lengths and biases PSI |

## The canonical order and why

1. **QC/trim** — but do NOT over-trim; short reads lose junction-spanning power.
2. **STAR pass 1 (all samples)** -> collect every `SJ.out.tab`.
3. **Build ONE combined junction DB** from the concatenated `SJ.out.tab` (order-trap: skipping this / per-sample DBs makes PSI incomparable).
4. **STAR pass 2 (all samples, same combined DB)** -> coordinate-sorted BAMs.
5. **Junction QC** — saturation curves must plateau (else deeper sequencing), entropy sane.
6. **Event-level differential splicing** — rMATS-turbo plus leafcutter for concordance.
7. **(Parallel) isoform-level DTU** — Salmon transcript quant -> tximport `dtuScaledTPM`+`txOut=TRUE` -> DRIMSeq/DEXSeq -> **stageR** (order-trap: reporting transcript q-values without the stageR gene screen inflates FDR).
8. **Sashimi** on the top reconciled events.

Order-traps that silently produce wrong results: per-sample (not cohort) 2-pass junctions; collapsing to gene before splicing analysis; transcript FDR without stageR; mixing ΔPSI significance with DTU proportion significance as if interchangeable.

## Choosing event-level vs isoform-level (and the caller)

Pipeline-level selection only; mechanism lives in the component skills.

| Fork | Lean toward | Hand off to |
|------|-------------|-------------|
| Event-level vs isoform-level | rMATS/leafcutter (which exon?) vs IsoformSwitchAnalyzeR/DRIMSeq+DEXSeq (which isoform, + NMD/ORF/domain consequences) | alternative-splicing/differential-splicing, alternative-splicing/isoform-switching |
| rMATS vs leafcutter vs MAJIQ | rMATS (known event types, replicate-based) + leafcutter (annotation-free intron clusters) for concordance; MAJIQ V3 HET for complex/heterogeneous cohorts | alternative-splicing/differential-splicing |
| DTU import scale | Salmon -> tximport `dtuScaledTPM` + `txOut=TRUE` (the ONLY correct DTU import) | rna-quantification/tximport-workflow |
| Transcript FDR | stageR two-stage: gene-level screen then transcript confirmation | alternative-splicing/isoform-switching |

## Primary path: STAR 2-pass + rMATS-turbo (+ leafcutter)

**Goal:** produce cohort-comparable PSI and a filtered differential-event table.

**Approach:** align all samples through one shared junction DB, QC junction saturation, then run rMATS-turbo (and leafcutter for concordance). Full runnable script: `examples/splicing_pipeline.sh`.

```bash
# Pass 1 (all samples) -> collect junctions
STAR --runThreadN 8 --genomeDir star_index/ --readFilesIn ${s}_R1.fq.gz ${s}_R2.fq.gz \
    --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${s}_pass1_
cat *_pass1_SJ.out.tab > combined_SJ.out.tab          # ONE shared DB for the whole cohort
# Pass 2 (all samples, same combined DB) -> comparable coordinates
STAR --runThreadN 8 --genomeDir star_index/ --readFilesIn ${s}_R1.fq.gz ${s}_R2.fq.gz \
    --readFilesCommand zcat --sjdbFileChrStartEnd combined_SJ.out.tab \
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${s}_

# Differential splicing. --readLength MUST match the trimmed reads, and --variable-read-length is
# required because the fastp step above trims to a RANGE, not a single length; without it rMATS
# miscomputes inclusion-junction lengths and biases PSI.
rmats.py --b1 cond1_bams.txt --b2 cond2_bams.txt --gtf annotation.gtf \
    -t paired --readLength 150 --variable-read-length --nthread 8 --od rmats_output --tmp rmats_tmp
```

Filter events on `|IncLevelDifference| > 0.1`, `FDR < 0.05`, and >=10 supporting junction reads averaged per replicate (sum `IJC_SAMPLE_1`, `SJC_SAMPLE_1`, `IJC_SAMPLE_2`, `SJC_SAMPLE_2` — each a comma-separated per-replicate list — and divide by the replicate count); rank by `-log10(FDR) * |IncLevelDifference|` (clamp FDR with `max(FDR, 1e-300)`). Resolve every column by header name: MXE carries two extra coordinate columns, so fixed positions silently read the wrong field. The read floor is not cosmetic — PSI is a ratio, so a 2-read event can reach `|dPSI| = 0.9` and pass FDR on noise alone.

## Parallel path: isoform-level DTU (Salmon -> stageR)

**Goal:** detect within-gene isoform-proportion switches with honest transcript-level FDR.

**Approach:** quantify transcripts with Salmon, import at DTU scale, test proportions with DRIMSeq/DEXSeq, and gate transcript q-values through stageR.

```r
library(tximport)
# DTU import is the OPPOSITE of gene-level DGE: keep transcripts, dtuScaledTPM.
# dtuScaledTPM scales by median tx length AMONG a gene's isoforms, so tx2gene is required even with txOut.
txi <- tximport(files, type = 'salmon', txOut = TRUE, countsFromAbundance = 'dtuScaledTPM', tx2gene = tx2gene)

# Canonical two-stage FDR route (Love, Soneson & Patro 2018): DRIMSeq/DEXSeq proportion test, then
# stageR gene-level SCREEN -> transcript-level CONFIRM. This is the "stageR seam" the principle names;
# mechanism lives in alternative-splicing/isoform-switching.

# Alternative route -- IsoformSwitchAnalyzeR adds NMD/ORF/protein-domain consequences. Its gene-level
# q-value comes from DEXSeq's perGeneQValue (min-p aggregation), NOT the stageR package:
library(IsoformSwitchAnalyzeR)
sList <- importRdata(isoformCountMatrix = txi$counts, isoformRepExpression = txi$abundance,
                     designMatrix = design, isoformExonAnnoation = 'annotation.gtf',
                     isoformNtFasta = 'transcripts.fa')      # 'isoformExonAnnoation' is the real (typo'd) arg
sList <- isoformSwitchTestDEXSeq(sList, reduceToSwitchingGenes = TRUE)   # gene q via DEXSeq perGeneQValue
```

## When NOT to use this pipeline (regime off-ramps)

This pipeline targets **bulk short-read differential splicing between two groups**. For other regimes, use the dedicated skill.

| Question | Use instead |
|----------|-------------|
| "Does this DNA variant alter splicing?" | alternative-splicing/splice-variant-prediction (SpliceAI, Pangolin, MMSplice) |
| "What is aberrant in this single rare-disease patient?" | alternative-splicing/outlier-splicing-detection (FRASER 2.0, OUTRIDER, DROP) |
| "Full-isoform analysis from PacBio Iso-Seq / ONT" | alternative-splicing/long-read-splicing (FLAIR, IsoQuant, Bambu, SQANTI3) |
| "Single-cell splicing analysis" | alternative-splicing/single-cell-splicing (chemistry-first; MARVEL, BRIE2) |
| "Heterogeneous cohort, n>=10 vs n>=10" | This pipeline + MAJIQ V3 HET (alternative-splicing/differential-splicing) |
| "Microexon-focused (3-27 nt)" | This pipeline with VAST-TOOLS or MicroExonator (alternative-splicing/splicing-quantification) |

## QC checkpoints between steps

| After | Gate | Interpretation |
|-------|------|----------------|
| QC/trim | Q30 >80%, adapter <5%, reads not over-trimmed | Aggressive trimming below ~75 nt weakens junction-spanning evidence |
| Alignment | Uniquely-mapped >80%; junction-saturation curves plateau | Still-rising curves = under-sequenced for splicing; deeper reads needed |
| Differential | \|ΔPSI\| >0.1 / >0.2, FDR <0.05, >=10 junction reads | Low read support = PSI is noise; require the read floor per event |
| DTU | stageR gene-level screen passed | Transcript q-values are only valid after the gene screen (alternative-splicing/isoform-switching) |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "Differential" events that are really coverage differences | Per-sample 2-pass junctions, not a shared DB | Concatenate all `SJ.out.tab`, feed the combined DB to pass 2 for every sample |
| Isoform-switch signal disappears | Gene-level tximport (`txOut=FALSE`) before splicing | Keep transcripts: `txOut=TRUE` + `dtuScaledTPM` on the DTU branch |
| Too many "significant" transcripts | Transcript q-values reported without stageR | Run the stageR two-stage gene->transcript test |
| Biased PSI across samples | `--readLength` != trimmed length | Set `--readLength` to the real length or add `--variable-read-length` |
| Event- and isoform-level calls "disagree" | Treating ΔPSI and DTU proportion tests as the same question | Reconcile as complementary; they test different things |

## Pipeline map (hand-offs)

- read-qc/fastp-workflow - QC/trim without over-trimming junction-spanning reads
- read-alignment/star-alignment - STAR 2-pass cohort-style, the shared junction DB
- alternative-splicing/splicing-qc - junction saturation/entropy, depth thresholds
- alternative-splicing/splicing-quantification - PSI computation, event taxonomy, sign conventions
- alternative-splicing/differential-splicing - rMATS/leafcutter/MAJIQ selection and reconciliation
- rna-quantification/alignment-free-quant - Salmon transcript quant for the DTU branch
- rna-quantification/tximport-workflow - dtuScaledTPM + txOut import
- alternative-splicing/isoform-switching - DTU + stageR + NMD/ORF/domain consequences
- alternative-splicing/sashimi-plots - ggsashimi/leafviz visualization

The complete runnable script is in this skill's examples/ (`splicing_pipeline.sh`).

## Related Skills

- read-qc/fastp-workflow - QC/trim options
- read-alignment/star-alignment - STAR 2-pass cohort-style configuration
- alternative-splicing/splicing-quantification - PSI computation, event taxonomy, sign conventions
- alternative-splicing/differential-splicing - Tool selection, MAJIQ V3, leafcutter, reconciliation
- alternative-splicing/isoform-switching - DTU + NMD/ORF/domain consequences (IsoformSwitchAnalyzeR v2, stageR)
- alternative-splicing/sashimi-plots - ggsashimi, MAJIQ-VOILA, leafviz visualization
- alternative-splicing/splicing-qc - STAR 2-pass cohort-style, library prep, depth thresholds
- alternative-splicing/single-cell-splicing - 10X chemistry decision; plate-based and long-read SC
- alternative-splicing/splice-variant-prediction - SpliceAI / Pangolin / MMSplice variant interpretation
- alternative-splicing/outlier-splicing-detection - FRASER 2.0 / DROP rare-disease workflow
- alternative-splicing/long-read-splicing - PacBio HiFi / ONT full-isoform analysis
- rna-quantification/alignment-free-quant - Salmon TPM for SUPPA2 and DTU pipelines
- rna-quantification/tximport-workflow - dtuScaledTPM + txOut DTU import

## References

- Shen S, Park JW, Lu ZX, et al (2014) rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *PNAS* 111:E5593-E5601. DOI 10.1073/pnas.1419161111.
- Li YI, Knowles DA, Humphrey J, et al (2018) Annotation-free quantification of RNA splicing using LeafCutter. *Nature Genetics* 50:151-158. DOI 10.1038/s41588-017-0004-9.
- Vitting-Seerup K, Sandelin A (2019) IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. *Bioinformatics* 35:4469-4471. DOI 10.1093/bioinformatics/btz247.
- Van den Berge K, Soneson C, Robinson MD, Clement L (2017) stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. *Genome Biology* 18:151. DOI 10.1186/s13059-017-1277-0.
- Love MI, Soneson C, Patro R (2018) Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification. *F1000Research* 7:952. DOI 10.12688/f1000research.15398.3. (the dtuScaledTPM -> DRIMSeq/DEXSeq -> stageR two-stage DTU workflow.)
