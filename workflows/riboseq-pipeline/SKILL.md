---
name: bio-workflows-riboseq-pipeline
description: End-to-end Ribo-seq analysis from FASTQ through periodicity QC, P-site calibration, ORF detection, translation efficiency, and stalling. Use when orchestrating a full ribosome profiling pipeline and deciding harvest/dedup/alignment options and which downstream analyses the library can support.
tool_type: mixed
primary_tool: STAR
workflow: true
depends_on:
  - ribo-seq/riboseq-preprocessing
  - ribo-seq/ribosome-periodicity
  - ribo-seq/orf-detection
  - ribo-seq/translation-efficiency
  - ribo-seq/ribosome-stalling
  - ribo-seq/initiation-site-mapping
qc_checkpoints:
  - after_align: "EndToEnd footprint alignment; deduplicate ONLY with UMIs; transcriptome BAM emitted"
  - periodicity_gate: "Strong 3-nt frame-0 periodicity REQUIRED before any ORF/stalling analysis; else gene-level counts only"
  - after_psite: "Per-read-length P-site offsets calibrated (not a single hardcoded 28 / read 5' end)"
  - after_te: "TE via count-based GLM (riborex/Xtail/anota2seq), never a ratio of ratios"
---

## Version Compatibility

Reference examples tested with: cutadapt 4.4+, umi_tools 1.1+, STAR 2.7.11+, bowtie2 2.5.3+, plastid 0.6+, riboWaltz 2.0+, RiboCode 1.2+, riborex 2.4+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Ribo-seq Pipeline

**"Analyze my ribosome profiling data from FASTQ to translation efficiency"** -> Orchestrate UMI handling, trimming, rRNA depletion, footprint-aware alignment, periodicity QC, P-site calibration, ORF detection, and differential translation, gating each downstream analysis on library quality.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

A Ribo-seq analysis is decided at four seams, two of them set at the bench before any sequencing.

1. **The harvest method is the deepest commitment and it fixes which analyses are even valid.** Cycloheximide (CHX) pre-treatment freezes elongating ribosomes but distorts codon-level dwell times (it can flip the codon-occupancy/tRNA-abundance correlation, Hussmann 2015); dwell-time / stalling / pausing analyses are only valid on flash-frozen, no-drug libraries. Gene-level footprint counts and TE are robust to CHX; codon-resolution analyses are not.
2. **Deduplicate ONLY when a UMI is present.** Ribosome footprints are ~28-30 nt and massively over-represented at abundant transcripts, so identical reads are mostly real signal; position-deduplicating a non-UMI library destroys it. With UMIs, dedup BOTH the genome and the transcriptome BAM.
3. **3-nt periodicity is a HARD gate, not a QC nicety.** Frame-based ORF calling and P-site analyses require strong frame-0 periodicity; a library that lacks it supports only gene-level counts. Certify it (riboWaltz frame-0 fraction) BEFORE any ORF/stalling step.
4. **P-site offsets are per-read-length, and TE is a count model.** Never hardcode a single 28-nt offset or the read 5' end; calibrate per length. Differential translation efficiency is a count-based GLM (riborex/Xtail/anota2seq) on CDS counts from BOTH assays, never a ratio of ribo/RNA ratios.

## Pipeline overview

```
FASTQ -> UMI extract -> trim -> rRNA remove -> STAR (EndToEnd) -> dedup (UMI only)
      -> periodicity QC (HARD GATE) + per-length P-site offsets -> [ORF detection | translation efficiency | stalling]
```

## Step 1: Preprocess

**Goal:** Produce a clean, footprint-aware alignment.

**Approach:** Extract UMIs first (if present), trim with a permissive floor, deplete rRNA before alignment, align end-to-end, and deduplicate only with UMIs. See riboseq-preprocessing for the decision tables.

```bash
# UMI-extract (if present) -> trim -> rRNA remove -> STAR EndToEnd -> dedup (UMI only)
cutadapt -a CTGTAGGCACCATCAAT --discard-untrimmed -m 15 -M 40 -o trimmed.fq.gz reads.fq.gz
bowtie2 -x contaminant_index -U trimmed.fq.gz --un-gz noncontam.fq.gz -S /dev/null -p 8
STAR --genomeDir STAR_index --readFilesIn noncontam.fq.gz --readFilesCommand zcat \
    --alignEndsType EndToEnd --seedSearchStartLmax 15 --outFilterMismatchNmax 2 \
    --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ribo_
samtools index ribo_Aligned.sortedByCoord.out.bam
```

`--quantMode TranscriptomeSAM` writes a SEPARATE `ribo_Aligned.toTranscriptome.out.bam` alongside the sorted genome BAM. RiboCode and riboWaltz transcriptome paths consume the TRANSCRIPTOME BAM; the sorted genome BAM is for plastid/genome-coordinate steps. With UMIs, deduplicate the transcriptome BAM too (coordinate-sort it first, then plain `umi_tools dedup --method directional`; see riboseq-preprocessing), or its ORF/periodicity inputs stay PCR-inflated. Do NOT add `--per-contig`/`--per-gene` here: they treat all reads on a transcript as one position, collapsing the per-codon footprints periodicity depends on.

## Step 2: Periodicity QC and P-site offsets

**Goal:** Certify the library and obtain per-length P-site offsets.

**Approach:** Run riboWaltz to filter periodic read lengths and calibrate offsets; the frame-0 fraction is the pass/fail metric. See ribosome-periodicity.

```r
library(riboWaltz)
annotation <- create_annotation("annotation.gtf")
reads <- bamtolist("bams", annotation = annotation)
reads <- length_filter(reads, length_filter_mode = "periodicity", periodicity_threshold = 50)
offsets <- psite(reads, extremity = "auto")   # per-length P-site offsets
```

Either riboWaltz (above) or the plastid `metagene generate` + `psite` CLI (used in the example script) is acceptable for offsets; pick one per project.

## Step 3: Detect ORFs

**Goal:** Call translated ORFs once offsets are known.

**Approach:** Run RiboCode; read lengths come from the metaplots config, and `-l` is the longest-ORF toggle. See orf-detection.

```bash
prepare_transcripts -g annotation.gtf -f genome.fa -o annot
metaplots -a annot -r ribo_Aligned.toTranscriptome.out.bam -o metaplots
RiboCode -a annot -c metaplots_pre_config.txt -A CTG,GTG -l no -p 0.05 -o ribocode_result
```

## Step 4: Translation efficiency

**Goal:** Test differential translation with matched RNA-seq.

**Approach:** Count both assays over the CDS and use a count-based GLM; use anota2seq when buffering vs control matters. See translation-efficiency.

```r
library(riborex)
res <- riborex(rnaCntTable = rna_cds_counts, riboCntTable = ribo_cds_counts,
               rnaCond = cond, riboCond = cond, engine = "DESeq2")
sig <- res[which(res$padj < 0.05), ]
```

## Step 5: Optional analyses

Stalling/pausing (only on flash-frozen no-drug data; see ribosome-stalling) and initiation-site mapping (needs a harringtonine/LTM library; see initiation-site-mapping) run off the same aligned BAM and calibrated offsets.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Downstream analyses all noisy | Periodicity QC skipped | Gate ORF/stalling on the frame-0 fraction first |
| P-site offsets look wrong | Single hardcoded offset across lengths | Calibrate per length with riboWaltz |
| RiboCode uses wrong read lengths | `-l` passed read lengths | Read lengths come from metaplots; `-l` is a toggle |
| TE hits dominated by low-count genes | Ratio testing | Use riborex/Xtail/anota2seq count GLMs |

## Related Skills

- ribo-seq/riboseq-preprocessing - UMI handling, trimming, rRNA removal, alignment
- ribo-seq/ribosome-periodicity - Periodicity QC and P-site calibration
- ribo-seq/orf-detection - Translated ORF calling
- ribo-seq/translation-efficiency - Differential TE and buffering
- ribo-seq/initiation-site-mapping - Start-codon mapping from TI-seq
- differential-expression/deseq2-basics - Count-based differential testing

## References

- McGlincy NJ, Ingolia NT. 2017. Transcriptome-wide measurement of translation by ribosome profiling. Methods 126:112-129. doi:10.1016/j.ymeth.2017.05.028
- Lauria F, Tebaldi T, Bernabò P, Groen EJN, Gillingwater TH, Viero G. 2018. riboWaltz: Optimization of ribosome P-site positioning in ribosome profiling data. PLoS Comput Biol 14(8):e1006169. doi:10.1371/journal.pcbi.1006169
- Xiao Z, Huang R, Xing X, Chen Y, Deng H, Yang X. 2018. De novo annotation and characterization of the translatome with ribosome profiling data. Nucleic Acids Res 46(10):e61. doi:10.1093/nar/gky179
- Hussmann JA, Patchett S, Johnson A, Sawyer S, Press WH. 2015. Understanding biases in ribosome profiling experiments reveals signatures of translation dynamics in yeast. PLoS Genet 11(12):e1005732. doi:10.1371/journal.pgen.1005732 (cycloheximide distorts codon-level dwell times)
