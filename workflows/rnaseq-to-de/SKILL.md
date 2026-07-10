---
name: bio-workflows-rnaseq-to-de
description: Orchestrates the end-to-end bulk RNA-seq differential-expression pipeline from FASTQ to an annotated DE gene table, chaining fastp QC/trim, Salmon (decoy-aware) or STAR+featureCounts quantification, tximport gene-level collapse, DESeq2/edgeR/limma-voom testing, apeglm shrinkage, and VST-based visualization. Use when committing the reference release and gene-ID namespace once for the whole run, sequencing steps in the defensible order (tximport before DE, raw counts into the model, VST only for viz/clustering), choosing alignment-free vs align-then-count and the DE engine, setting strandedness correctly, keeping batch in the design instead of correcting-then-testing, or handing the signed ranking statistic to downstream enrichment. Hands mechanism to the component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: DESeq2
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - rna-quantification/alignment-free-quant
  - read-alignment/star-alignment
  - read-qc/rnaseq-qc
  - rna-quantification/tximport-workflow
  - rna-quantification/count-matrix-qc
  - differential-expression/deseq2-basics
  - differential-expression/edger-basics
  - differential-expression/de-results
  - differential-expression/de-visualization
qc_checkpoints:
  - after_qc: "Q30 >80%, adapter content <5% (RNA has a lower quality floor than DNA)"
  - after_quant: "Mapping rate >70%, >10M reads mapped, flat gene-body coverage, low rRNA/intronic"
  - after_import: "tx2gene release matches the Salmon index; ID conversion loses few transcripts"
  - after_de: "Dispersion trend sane, PCA separates condition not batch, no Cook's outliers"
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, tximport 1.30+, apeglm 1.24+, STAR 2.7.11+, Salmon 1.10+, Subread/featureCounts 2.0.2+ (--countReadPairs added in 2.0.2), fastp 0.23+, ggplot2 3.5+ (kallisto 0.50+ as a Salmon alternative)

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: Salmon selective alignment is default since 1.0 (the historical `--validateMappings` is now a no-op); `DESeqDataSetFromTximport` carries the average-transcript-length offset automatically; `lfcShrink(type='apeglm')` requires `coef` to name a `resultsNames(dds)` coefficient and DROPS the `stat` column. Confirm these in-tool before quoting.

# RNA-seq to Differential Expression Workflow

**"Find differentially expressed genes from my RNA-seq FASTQ files"** -> Chain QC/trim, decoy-aware quantification, tximport gene-level collapse, a count-based DE test, shrinkage, and visualization into one annotated DE table.
- CLI + R: fastp -> (salmon | STAR + featureCounts) -> tximport -> DESeq2/edgeR/limma-voom -> lfcShrink -> VST/volcano

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A bulk RNA-seq result is decided at four seams between steps, not inside any one tool.

1. **The reference RELEASE + transcriptome/GTF pair is a pipeline-wide commitment made once at quantification and inherited by everything downstream.** The transcriptome FASTA that builds the Salmon index and the GTF that builds the tx2gene map (and drives featureCounts) must be the SAME Ensembl/GENCODE release. Mixing an index built on release 104 with a tx2gene from release 110 silently drops renamed/removed transcripts — no error, just missing genes. This choice also fixes the gene-ID namespace (ENSG is the safe backbone; convert to symbol/Entrez only at the reporting/enrichment seam). Changing the release later forces re-quantification.
2. **Raw counts flow forward; normalized/transformed values are terminal.** Integer counts (or tximport count-scale output) are the ONLY valid input to DESeq2/edgeR/limma-voom. TPM/CPM are for within-sample ranking only; a VST/rlog matrix is for PCA, clustering, heatmaps, and ML — never for re-running a count-based test. Feeding the wrong scale across a join is the single most common silent corruption.
3. **Collapse transcript->gene through tximport, not by summing counts.** tximport carries the average-transcript-length offset that corrects for isoform-usage shifts; naively summing Salmon `NumReads` biases gene counts whenever isoform usage changes across conditions.
4. **Batch belongs in the design, not "corrected" then tested.** Put known batch in the formula (`~ batch + condition`). Running `removeBatchEffect`/ComBat and then testing on the corrected matrix exaggerates confidence (Nygaard 2016); the corrected matrix is for visualization/clustering/ML input only.

## Pipeline map

```
FASTQ (paired)
  | [1] QC & trim ----------------> fastp              (read-qc/fastp-workflow)
  v
  | [2] Quantify -----------------> salmon (decoy-aware)   (rna-quantification/alignment-free-quant)
  v     |   OR  STAR + featureCounts (need a BAM?)         (read-alignment/star-alignment)
  v     ^-- commitment: reference RELEASE + tx/GTF pair, gene-ID namespace
  | [3] Import & collapse tx->gene -> tximport         (rna-quantification/tximport-workflow)
  v     ^-- carries the length offset; NEVER sum NumReads
  | [4] Pre-DE QC ----------------> PCA / dispersion / outliers  (rna-quantification/count-matrix-qc)
  v
  | [5] DE test ------------------> DESeq2 | edgeR-QL | limma-voom  (differential-expression/deseq2-basics)
  v     ^-- RAW counts in; batch in the design, not corrected-then-tested
  | [6] Shrink & extract ---------> lfcShrink(apeglm); pull Wald `stat` for ranking  (differential-expression/de-results)
  v
  | [7] Visualize ----------------> VST heatmap/PCA, volcano  (differential-expression/de-visualization)
  v
Annotated DE table (ENSG + symbol + biotype, log2FC, stat, pvalue, padj, baseMean)
```

## Reference, IDs, and quantification target: the made-once commitments

Decided before the first `salmon quant`; everything downstream inherits them. Mechanism lives in the component skills; the reasoning below is what a reviewer expects justified.

| Commitment | Options | Consequence inherited downstream |
|------------|---------|----------------------------------|
| Reference release | One Ensembl/GENCODE release for BOTH the transcriptome FASTA (index) and the GTF (tx2gene / featureCounts) | Any mismatch silently drops renamed transcripts; fixes DE row names and the pathway-DB key space |
| Gene-ID namespace | ENSG backbone (convert to symbol/Entrez only at reporting) | Symbol space is lossy (aliases, many-ENSG-one-symbol merges genes); stripping the ENSG `.version` with `\..*` also destroys the GENCODE `_PAR_Y` tag, collapsing chrY-PAR onto chrX (rna-quantification/tximport-workflow) |
| Quantification target | Gene-level DGE (`countsFromAbundance="no"`, offset carried) vs transcript-level DTU | DTU needs a DIFFERENT import (`txOut=TRUE` + `dtuScaledTPM`); switching later is a re-import, not a filter — see workflows/splicing-pipeline |
| 3'-tagged vs full-length | 3'-tagged (QuantSeq/bulk-10x): `countsFromAbundance="no"`, no length offset | Length-bias correction does not apply to 3'-tagged libraries |

## The canonical order and why

Each step assumes the previous; two reorderings silently produce wrong results.

1. **QC/trim before quantification** — adapter/quality tails corrupt pseudo-mapping and duplicate structure.
2. **Quantify to the committed reference** — Salmon decoy-aware (genome as decoy) so intron/pseudogene reads are not misassigned to transcripts; STAR only when a genome BAM is also needed downstream.
3. **Import via tximport** — the tx->gene collapse happens HERE, carrying the length offset (order-trap: summing `NumReads` biases genes under isoform shift).
4. **Pre-filter low-count genes** (`rowSums(counts) >= 10`) — this is a speed/memory step, NOT the FDR filter. Order-trap: it does not replace the baseMean independent filtering that `results()` applies at the FDR step; `filterByExpr(y, design)` (edgeR) is the design-aware version and must run once BEFORE dispersion, never after.
5. **DESeq()** on raw counts with batch in the design — size factors, dispersion, Wald/LRT.
6. **results() then lfcShrink()** — independent filtering happens inside `results()` on baseMean; shrink LFC for effect sizes/ranking, but p-values stay from the unshrunken test. Order-trap: apeglm/ashr objects DROP the `stat` column — pull the Wald `stat` from unshrunk `results()` if a signed ranking metric is needed for GSEA.
7. **Visualize on VST** (heatmaps/PCA); volcano uses shrunken LFC + unshrunken p.

## Choosing the quantifier and the DE engine

Pipeline-level selection only; mechanism lives in the component skills.

| Fork | Lean toward | Hand off to |
|------|-------------|-------------|
| Alignment-free (Salmon/kallisto) vs align-then-count (STAR+featureCounts) | Salmon for gene-level DGE (decoy-aware, GC/seq-bias correction, no BAM); STAR when a genome BAM is also needed (splicing, coverage, novel junctions, variants) | rna-quantification/alignment-free-quant, read-alignment/star-alignment |
| DESeq2 vs edgeR-QL vs limma-voom | limma-voom when library sizes vary >3x or outliers dominate; edgeR-QL for tight finite-sample type-I control; DESeq2 for the apeglm/downstream ecosystem (70-90% top-gene overlap on well-designed data) | differential-expression/deseq2-basics, differential-expression/edger-basics |
| `countsFromAbundance` | `no` (gene DGE via DESeqDataSetFromTximport) / `lengthScaledTPM` (DGE when the tool can't take offsets) / `dtuScaledTPM`+`txOut` (DTU) | rna-quantification/tximport-workflow |
| Strandedness `-s` | Confirm, never assume: STAR `ReadsPerGene.out.tab` cols 3 vs 4, or RSeQC `infer_experiment.py`; dUTP/TruSeq is reverse (`-s 2`) | read-qc/rnaseq-qc |

## Primary path: Salmon + tximport + DESeq2

**Goal:** turn trimmed FASTQ into a shrunken, annotated gene-level DE table.

**Approach:** build a decoy-aware index once, quantify each sample, collapse to genes via tximport (release-matched tx2gene), test raw counts with batch in the design, shrink for ranking. Full runnable script: `examples/salmon_deseq2_workflow.R`.

```bash
# Index once: decoy-aware (genome as decoy) so intron/pseudogene reads are not misassigned
grep "^>" genome.fa | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
cat transcriptome.fa genome.fa > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -i salmon_index -k 31 -p 8

# Quantify (selective alignment is default since 1.0; --gcBias/--seqBias correct known biases)
salmon quant -i salmon_index -l A -1 trimmed/${s}_R1.fq.gz -2 trimmed/${s}_R2.fq.gz \
    -o quants/${s} --gcBias --seqBias -p 8
```

```r
library(tximport); library(DESeq2)
# tx2gene MUST come from the same release as the index (else renamed transcripts drop silently)
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ batch + condition)  # batch in design
dds <- dds[rowSums(counts(dds)) >= 10, ]              # speed filter, NOT the FDR filter
dds$condition <- relevel(dds$condition, ref = 'control')
dds <- DESeq(dds)                                     # RAW counts in
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')  # ranking/effect size
# For GSEA ranking, pull the Wald stat from the UNSHRUNK results (apeglm drops `stat`):
res_stat <- results(dds, name = 'condition_treated_vs_control')$stat
```

## Alternative path: STAR + featureCounts + DESeq2

**Goal:** produce a genome BAM (reused by splicing/coverage/variant steps) alongside gene counts.

**Approach:** align with `--sjdbOverhang` = readlen-1, count with the verified strandedness, then `DESeqDataSetFromMatrix`. Full script: `examples/star_deseq2_workflow.sh`.

```bash
STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles genome.fa \
    --sjdbGTFfile genes.gtf --sjdbOverhang 149 --runThreadN 8   # 149 for 2x150, not a blanket 100
STAR --genomeDir star_index --readFilesIn trimmed/${s}_R1.fq.gz trimmed/${s}_R2.fq.gz \
    --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
    --outFileNamePrefix aligned/${s}_ --runThreadN 8
# -s from ReadsPerGene.out.tab cols 3 vs 4 (or infer_experiment.py); -s 2 = dUTP/TruSeq reverse
featureCounts -T 8 -p --countReadPairs -s 2 -a genes.gtf -o counts.txt aligned/*_Aligned.sortedByCoord.out.bam
```

```r
counts <- read.table('counts.txt', header = TRUE, row.names = 1, skip = 1)[, -(1:5)]
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ batch + condition)
```

## QC checkpoints between steps

| After | Gate | Interpretation |
|-------|------|----------------|
| QC/trim | Q30 >80%, adapter <5% | RNA has a lower quality floor than DNA; sharp Q30 drop = degraded input |
| Quant/align | Mapping >70%, >10M reads mapped; flat gene-body coverage; low rRNA%/intronic% | 3' bias = degradation/oligo-dT; high intronic = pre-mRNA/gDNA; high intergenic = gDNA/annotation gap — all compromise DE BEFORE it runs (read-qc/rnaseq-qc) |
| Import | tx2gene release == index release; few transcripts dropped; report the ID-conversion rate (<0.85 => wrong ID type or organism) | Mismatched release silently loses renamed transcripts |
| Pre-DE | Dispersion trend sane; PCA separates condition not batch; no Cook's outliers | PCA clustering by batch => batch dominates; add it to the design (rna-quantification/count-matrix-qc) |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Many genes missing / low tx conversion | Salmon index and tx2gene from different releases | Rebuild both from ONE release; pin it for the whole cohort |
| Gene counts biased where isoforms switch | Summed Salmon `NumReads` instead of importing | Collapse via tximport (carries the length offset) |
| "invalid class DESeqDataSet" or nonsense LFCs | TPM/VST fed into DESeq2 | Raw counts (or `lengthScaledTPM`) only; VST is for viz/ML |
| Counts collapsed / library looks failed | Wrong featureCounts `-s` strandedness | Infer with `infer_experiment.py` or STAR cols 3 vs 4 before counting |
| Suspiciously many DE genes, tiny p-values | Batch-corrected matrix fed to the test | Keep batch in the design; correct only for visualization |
| `lfcShrink` error / wrong contrast | `coef` not in `resultsNames(dds)`, or ranking off the shrunk object | Use a `resultsNames` coefficient; pull `stat` from unshrunk `results()` |

## Pipeline map (hand-offs)

- read-qc/fastp-workflow - adapter/quality trimming and report interpretation
- rna-quantification/alignment-free-quant - Salmon/kallisto decoy-aware quantification
- read-alignment/star-alignment - STAR index/sjdbOverhang, 2-pass, GeneCounts strandedness
- rna-quantification/tximport-workflow - tx->gene collapse, countsFromAbundance, tx2gene, ID-version traps
- rna-quantification/count-matrix-qc - pre-DE PCA, dispersion, Cook's outliers, batch checks
- differential-expression/deseq2-basics - the DESeq2 model, design, contrasts
- differential-expression/de-results - extracting/annotating results and the signed ranking statistic
- differential-expression/de-visualization - volcano/MA/heatmap on the right scale

The complete runnable scripts for both paths are in this skill's examples/ (`salmon_deseq2_workflow.R`, `star_deseq2_workflow.sh`).

## Related Skills

- database-access/geo-data - Find a GSE on GEO, detect SuperSeries, link to SRA
- database-access/sra-data - Download paired-end FASTQ from SRA / ENA / STRIDES cloud
- sequence-io/fastq-quality - Confirm the FASTQ quality encoding before trimming public or pre-2011 data
- sequence-io/paired-end-fastq - Keep R1/R2 mates synchronized; independent per-mate filtering desyncs pairs
- read-qc/fastp-workflow - Detailed QC options and parameters
- read-qc/rnaseq-qc - Post-alignment RNA QC: strandedness, gene-body coverage, rRNA/intronic
- read-alignment/star-alignment - The align path (BAM for splicing/coverage/variants)
- rna-quantification/alignment-free-quant - Salmon and kallisto details
- rna-quantification/tximport-workflow - tximport options, countsFromAbundance, tx2gene creation
- rna-quantification/count-matrix-qc - Pre-DE QC and diagnostics
- differential-expression/deseq2-basics - Complete DESeq2 reference
- differential-expression/de-results - Results extraction, annotation, ranking statistic
- differential-expression/de-visualization - Advanced visualization options
- alternative-splicing/isoform-switching - Transcript-level DTU when gene-level is not enough (splicing fork)
- pathway-analysis/go-enrichment - Next step: functional enrichment (workflows/expression-to-pathways)

## References

- Soneson C, Love MI, Robinson MD (2015) Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research* 4:1521. DOI 10.12688/f1000research.7563.1. (tximport; the tx->gene length-offset seam.)
- Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550. DOI 10.1186/s13059-014-0550-8.
- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C (2017) Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods* 14:417-419. DOI 10.1038/nmeth.4197.
- Nygaard V, Rødland EA, Hovig E (2016) Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17:29-39. DOI 10.1093/biostatistics/kxv027. (batch belongs in the design.)
- Ewels PA, Peltzer A, Fillinger S, et al (2020) The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology* 38:276-278. DOI 10.1038/s41587-020-0439-x. (nf-core/rnaseq: the reproducible reference orchestration.)
