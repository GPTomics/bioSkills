---
name: bio-rna-quantification-tximport-workflow
description: Import transcript-level quantifications from Salmon/kallisto/RSEM into R for gene-level analysis with DESeq2/edgeR using tximport or tximeta. Use when summarizing transcript abundances to gene counts with the correct length offset, choosing a countsFromAbundance mode (full-length vs 3'-tag vs DTU), resolving transcript-ID version mismatches, or handing off to DESeq2/edgeR without double-applying the offset.
tool_type: r
primary_tool: tximport
---

## Version Compatibility

Reference examples tested with: tximport 1.30+, tximeta 1.20+, DESeq2 1.42+, edgeR 4.0+, txdbmaker 1.0+, Salmon 1.10+, kallisto 0.50+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# tximport Workflow

**"Import Salmon/kallisto results into DESeq2"** -> Summarize transcript-level abundance estimates to gene-level counts AND compute a per-gene, per-sample length offset that the DE model consumes.
- R: `tximport::tximport(files, type='salmon', tx2gene=tx2gene)`

## Why this is not just a sum

tximport does not merely add transcript counts to a gene total. The fragment count a gene produces depends on the average length of the isoforms expressed in that sample, because longer molecules yield more fragments (more start positions). When isoform usage shifts between conditions (differential transcript usage), the gene's average effective length changes, so a naive summed count is length-biased in a condition-correlated way and masquerades as differential expression. tximport corrects this by returning a per-gene, per-sample average-length matrix (`txi$length`) and passing it as a normalization offset to DESeq2/edgeR. A single per-gene length cannot capture this because the bias is sample-specific (Soneson, Love, Robinson 2015).

## Basic tximport

**Goal:** Import transcript-level quantifications into R as gene-level counts plus the length offset for DESeq2 or edgeR.

**Approach:** Build a transcript-to-gene map, then run tximport over the quant files; the returned `txi$counts`/`txi$length` carry both the gene counts and the offset.

```r
library(tximport)

files <- c(sample1 = 'sample1_quant/quant.sf',
           sample2 = 'sample2_quant/quant.sf',
           sample3 = 'sample3_quant/quant.sf')

tx2gene <- read.csv('tx2gene.csv')                 # column order: TXNAME, then GENEID
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
```

`txi` is a list: `$abundance` (TPM), `$counts` (estimated counts), `$length` (the average-length offset source), `$countsFromAbundance`.

## Decision: countsFromAbundance mode

This argument silently determines correctness; nothing errors when it is wrong.

| Mode | What it returns | Use when |
|------|-----------------|----------|
| `'no'` (default) | Estimated counts + separate length offset | Full-length library -> DESeq2/edgeR (they consume the offset). The cleanest path. |
| `'lengthScaledTPM'` | Counts with the length correction baked in, no separate offset | A tool that cannot take an offset (e.g. limma-voom) |
| `'scaledTPM'` | TPM scaled to library size, no length scaling | Transcript-level DTU with `txOut=TRUE` (DRIMSeq/DEXSeq); the established Love et al. workflow input |
| `'dtuScaledTPM'` | Scaled by median isoform length | DTU alternative (tximport >= 1.10), needs `tx2gene`; helps when isoform lengths within a gene differ widely |

For DTU, `scaledTPM` is the established default; `dtuScaledTPM` is the newer purpose-built mode, preferable when a gene's isoforms span very different lengths.

The 3'-tag exception: for 3'-end protocols (10x, QuantSeq, Lexogen) a read count does not scale with transcript length, so there is no length bias to correct, and length-correcting injects one. Do not use the length-scaled modes (`lengthScaledTPM`/`dtuScaledTPM`) for tag-seq. Import with the default, but build the DESeqDataSet from the plain counts so the length offset is NOT auto-applied:

```r
# 3'-tag: bypass the length offset that DESeqDataSetFromTximport would otherwise apply
dds <- DESeqDataSetFromMatrix(round(txi$counts), colData = coldata, design = ~ condition)
```

```r
# Full-length, DESeq2/edgeR (default): keep the offset path
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)

# Transcript-level for DTU (hand off to alternative-splicing/isoform-switching)
txi_tx <- tximport(files, type = 'salmon', txOut = TRUE,
                   countsFromAbundance = 'scaledTPM')
```

## Creating tx2gene

The map is a two-column data frame; column ORDER is load-bearing (TXNAME first, GENEID second), names do not matter.

**Goal:** Map every quantified transcript ID to its gene, with IDs that exactly match the quant files.

**Approach:** Derive from the annotation that built the index (GTF, ensembldb, biomaRt, or the index t2g); strip version suffixes to match.

```r
# From a GTF: makeTxDbFromGFF moved to txdbmaker in Bioconductor >= 3.19
# (defunct in GenomicFeatures >= 1.61.1; on older Bioconductor use GenomicFeatures::makeTxDbFromGFF)
library(txdbmaker)
txdb <- makeTxDbFromGFF('annotation.gtf')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = 'TXNAME',
                                 columns = c('TXNAME', 'GENEID'))

# From biomaRt (useEnsembl; useMart is deprecated)
library(biomaRt)
mart <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl')
tx2gene <- getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id'), mart = mart)
```

## The #1 silent failure: transcript-ID version mismatch

If `quant.sf` IDs carry version suffixes (`ENST00000456328.4`) but `tx2gene` does not (or vice versa), the IDs do not match. Total non-overlap raises an error; partial mismatch silently drops the non-matching transcripts and prints a summary, deflating affected genes toward zero. Fix by stripping versions consistently or with `ignoreTxVersion`:

```r
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene,
                ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
```

## Handoff to DESeq2 (offset applied automatically)

**Goal:** Build a DESeqDataSet that uses the tximport length offset without any manual step.

**Approach:** `DESeqDataSetFromTximport` stores `txi$length` as the `avgTxLength` assay and converts it to per-gene normalization factors inside `DESeq()`.

```r
library(DESeq2)
coldata <- data.frame(condition = factor(c('control', 'control', 'treated', 'treated')),
                      row.names = names(files))
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # light pre-filter (speed); results() does the inferential filter
dds <- DESeq(dds)
res <- results(dds)
```

Passing a `countsFromAbundance='no'` txi prints "using counts and average transcript lengths from tximport"; a length-scaled txi prints "using just counts" and applies no offset. Both are handled correctly by the function.

## Handoff to edgeR (manual offset)

**Goal:** Carry the length offset into an edgeR DGEList.

**Approach:** Geometric-mean-center the length matrix, fold in composition-corrected library sizes, log it, attach via `scaleOffset`.

```r
library(edgeR)
cts <- txi$counts
normMat <- txi$length / exp(rowMeans(log(txi$length)))   # center each gene on its geometric mean
normCts <- cts / normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, '*')
y <- scaleOffset(DGEList(cts), log(normMat))
y <- y[filterByExpr(y, group = coldata$condition), , keep.lib.sizes = FALSE]   # group-aware filter
```

Do not double-apply the offset: if `countsFromAbundance='lengthScaledTPM'` already baked the correction into the counts, do not also attach a length offset. Use `'no'` for the offset path, the scaled modes for the no-offset path, never both.

## Transcript-level uncertainty (DTE/DTU)

Gene-level estimates are robust because per-isoform assignment uncertainty cancels on summation. Transcript-level testing must propagate it: edgeR `catchSalmon` deflates counts by per-transcript overdispersion (differential-expression/edger-basics), swish/fishpond tests across Salmon Gibbs samples (alternative-splicing/isoform-switching), and sleuth uses kallisto bootstraps (expression-matrix/counts-ingest). Generate the replicates at quantification time (rna-quantification/alignment-free-quant).

## tximeta: provenance by checksum

tximeta hashes the index's reference sequences and looks the digest up against known GENCODE/Ensembl/RefSeq releases, attaching transcript ranges and release metadata automatically, so the exact reference becomes a verified property of the object rather than lab lore.

```r
library(tximeta)
makeLinkedTxome(indexDir = 'salmon_index', source = 'Ensembl', organism = 'Homo sapiens',
                release = '110', genome = 'GRCh38', fasta = 'transcripts.fa', gtf = 'annotation.gtf')
coldata <- data.frame(names = names(files), files = files,
                      condition = c('control', 'control', 'treated', 'treated'))
se <- tximeta(coldata)
gse <- summarizeToGene(se)
dds <- DESeqDataSet(gse, design = ~ condition)
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Many genes import as zero or deflated | Transcript-ID version mismatch (partial drop) | `ignoreTxVersion = TRUE`; or strip `\.\d+$` from both sides |
| Error: none of the transcripts present in tx2gene | Total ID mismatch (versions or wrong annotation) | Rebuild tx2gene from the annotation that built the index |
| Summarized at the wrong level, no error | tx2gene columns reversed (GENEID first) | Order as TXNAME, then GENEID |
| Length bias appears in 3'-tag data | DESeqDataSetFromTximport auto-applied the length offset | Build via `DESeqDataSetFromMatrix(round(txi$counts), ...)` so no offset is applied |
| Fold changes inflated near isoform switches with manual edgeR | Offset double-applied or omitted | One path only: `'no'`+offset, or scaled mode without offset |

## Related Skills

- rna-quantification/alignment-free-quant - Upstream Salmon/kallisto and inferential replicates
- differential-expression/deseq2-basics - Gene-level DE from a DESeqDataSet
- differential-expression/edger-basics - edgeR DE and catchSalmon transcript DTE
- alternative-splicing/isoform-switching - DTU and swish from transcript-level import
- expression-matrix/counts-ingest - sleuth and other quantifier ingestion paths
- genome-intervals/gtf-gff-handling - Building tx2gene from a GTF

## References

- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4:1521. doi:10.12688/f1000research.7563
- Love MI, Soneson C, Hickey PF, et al. 2020. Tximeta: Reference sequence checksums for provenance identification in RNA-seq. PLoS Comput Biol 16(2):e1007664. doi:10.1371/journal.pcbi.1007664
