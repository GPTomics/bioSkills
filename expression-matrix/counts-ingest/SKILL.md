---
name: bio-expression-matrix-counts-ingest
description: Imports gene expression count matrices from featureCounts, HTSeq, STAR ReadsPerGene, Salmon/kallisto via tximport or tximeta, RSEM, 10X Genomics MTX/H5, AnnData H5AD, and RDS. Handles silent-miscounting traps (featureCounts -p v2.0.2 API break, STAR strandedness column choice, salmon NumReads-sum without tximport, RSEM non-integer expected_count, GENCODE _PAR_Y suffix, zero-length-transcript TPM divide-by-zero), and encodes the tximport countsFromAbundance decision tree with the "lengthScaledTPM is not TPM" warning. Use when assembling a gene-by-sample count matrix from aligner or quantifier output, importing salmon/kallisto for DESeq2 vs limma-voom, choosing strandedness column for STAR, debugging zero-count panics, or building tx2gene mapping.
tool_type: mixed
primary_tool: tximport
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, numpy 1.26+, scanpy 1.10+, anndata 0.10+, tximport 1.30+, tximeta 1.20+, GenomicFeatures 1.54+, Subread/featureCounts 2.0.6+ (post-v2.0.2 API), STAR 2.7.10+, Salmon 1.10+, kallisto 0.48+, RSEM 1.3.3+, HTSeq 2.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Count Matrix Ingestion

**"Load my featureCounts / Salmon / STAR output into a count matrix"** -> Parse the per-sample quantification, strip metadata, choose the correct strandedness or NumReads column, optionally apply length-bias correction via tximport, and return a gene-by-sample matrix appropriate for the downstream DE tool.

## The Single Most Important Modern Insight -- `lengthScaledTPM` is NOT TPM; the naming has fooled many

`tximport(..., countsFromAbundance='lengthScaledTPM')` returns a **count-scale matrix** (sum to library size) with the gene-length bias removed, intended as input to DE tools that cannot accept offsets (limma-voom). The word "TPM" in the option name has misled users into reporting these values as normalized abundance -- they are not.

The decision tree for `countsFromAbundance` (Soneson, Love, Robinson 2015 *F1000Res* 4:1521):

| Option | Returns | Use with |
|--------|---------|----------|
| `'no'` (default) | Raw NumReads; length matrix passed as DESeq2/edgeR offset | DESeq2 via `DESeqDataSetFromTximport`, edgeR via `DGEList` -- the correct path for DGE |
| `'scaledTPM'` | TPM scaled to library size | DGE when downstream tool cannot accept offsets (rare) |
| `'lengthScaledTPM'` | TPM scaled by avg transcript length AND library size | limma-voom for DGE; recommended modern default when offsets unavailable |
| `'dtuScaledTPM'` | Per-transcript scaling by median isoform length | Differential Transcript Usage (DRIMSeq, DEXSeq) ONLY; requires `txOut=TRUE` |

Two adjacent traps with the same flavor of silent miscount:

1. **featureCounts `-p` since v2.0.2 (March 2021) needs `--countReadPairs`**. Pre-v2.0.2, `-p` flagged paired-end AND counted pairs as fragments. Post-v2.0.2, `-p` only flags paired-end -- pairs are NOT counted as one fragment unless `--countReadPairs` is added. Old scripts run on new installs produce ~2x the expected counts (one count per mate). No warning. Always pass `-p --countReadPairs` together for paired-end.

2. **STAR `ReadsPerGene.out.tab` column choice depends on library protocol**. Column 2 is unstranded; column 3 is forward-stranded; column 4 is reverse-stranded. Illumina TruSeq Stranded (the dominant kit, dUTP-based) is reverse-stranded -- column 4. Reading column 3 for TruSeq throws away ~95% of reads. Reading column 2 conflates antisense expression. Verify with RSeQC `infer_experiment.py` on a few BAMs before assembling the matrix.

## Algorithmic Taxonomy

| Source | Output structure | Read into | Caveats |
|--------|------------------|-----------|---------|
| featureCounts (Liao 2014) | TSV with 6 metadata cols + N sample cols | `read.delim` (R), `pd.read_csv(comment='#')` (Python) | -p/--countReadPairs v2.0.2 break; -O double-counts; -s strandedness |
| HTSeq-count (Anders 2015) | Per-sample 2-col TSV with `__` summary lines | Per-sample read, drop `__` rows, concat | `--mode` matters (union vs intersection-strict vs intersection-nonempty); `-s` strandedness |
| STAR `--quantMode GeneCounts` | Per-sample 4-col `ReadsPerGene.out.tab` | Skip first 4 summary rows; pick column by strandedness | Column 4 = TruSeq Stranded |
| Salmon `quant.sf` | Per-sample 5-col TSV (Name, Length, EffectiveLength, TPM, NumReads) | tximport (recommended) or manual NumReads sum (length-biased) | Selective alignment is the default since Salmon 1.0.0 (Srivastava 2020); the `--validateMappings` flag is now a no-op |
| kallisto `abundance.tsv` / `.h5` | Per-sample 5-col TSV; bootstraps in `.h5` | tximport (gene-level) or sleuth (transcript-level with bootstrap variance) | `kallisto quant -b 100` for sleuth |
| RSEM `*.genes.results` / `*.isoforms.results` | Per-sample TSV; expected_count is non-integer | tximport (`type='rsem'`) with `round()` via DESeq2; or manual | Zero-length transcripts cause `lengths > 0 is not TRUE` |
| 10X Genomics CellRanger | filtered_feature_bc_matrix/ MTX dir OR `.h5` | scanpy `read_10x_mtx` / `read_10x_h5`; Seurat `Read10X` | Single-cell convention: cells in rows |
| AnnData `.h5ad` | scverse single-file binary | scanpy `read_h5ad`; in R via zellkonverter | `.X` vs `.layers['counts']` vs `.raw.X` semantics |
| RDS (from R) | R-serialized object | `pyreadr` (Python); base R `readRDS` | For Seurat objects, convert first |

## Decision Tree by Scenario

| Scenario | Recommended approach |
|----------|---------------------|
| Bulk RNA-seq with featureCounts paired-end | `featureCounts -p --countReadPairs -s 2 -t exon -g gene_id` (TruSeq Stranded) |
| Bulk RNA-seq with STAR | `--quantMode GeneCounts`; read column 4 for TruSeq Stranded |
| Salmon/kallisto -> gene-level DGE with DESeq2 | `tximport(type='salmon', tx2gene)` + `DESeqDataSetFromTximport()` -- offsets handled |
| Salmon -> gene-level DGE with limma-voom | `tximport(..., countsFromAbundance='lengthScaledTPM')` -- no offset |
| Salmon/kallisto -> DTE (differential transcript expression) | `tximport(..., txOut=TRUE, countsFromAbundance='dtuScaledTPM')` for DRIMSeq/DEXSeq; OR `edgeR::catchSalmon()` for the Baldoni 2024 framework |
| kallisto with bootstraps -> sleuth (uncertainty-aware DE) | `sleuth_prep()` directly; do not go through tximport |
| RSEM -> DESeq2 | `tximport(files, type='rsem', txIn=FALSE)` + `DESeqDataSetFromTximport()` |
| Want automatic annotation provenance | `tximeta` (Love 2020 *PLoS Comp Biol* 16:e1007664) |
| 3'-tagged library (10x bulk, QuantSeq) | `countsFromAbundance='no'` WITHOUT length offset -- length bias negligible |
| 10X single-cell | `sc.read_10x_h5()` or `Read10X_h5()` |
| Strandedness unknown | RSeQC `infer_experiment.py` on 1-2 BAMs BEFORE re-running quantification |

## featureCounts -- The `-p` and `-O` Traps

**Goal:** Run featureCounts correctly on paired-end stranded RNA-seq, then read the output without inheriting the metadata columns.

**Approach:** CLI invocation with `-p --countReadPairs -s <strandedness>`, plus optional `-O --fraction` for overlapping-gene handling; parse with pandas/read.delim stripping the 6 metadata columns.

```bash
featureCounts -T 8 -p --countReadPairs -s 2 \
    -t exon -g gene_id \
    -a annotation.gtf \
    -o featurecounts.txt \
    sample1.bam sample2.bam sample3.bam
```

```python
import pandas as pd

fc = pd.read_csv('featurecounts.txt', sep='\t', comment='#')
counts = fc.set_index('Geneid').iloc[:, 5:]
counts.columns = [c.replace('.bam', '').split('/')[-1] for c in counts.columns]
```

```r
fc <- read.delim('featurecounts.txt', comment.char = '#', row.names = 1)
counts <- fc[, 6:ncol(fc)]
colnames(counts) <- gsub('.*/|\\.bam$', '', colnames(counts))
```

| Flag | Meaning | Default | When to flip |
|------|---------|---------|--------------|
| `-p` | Input is paired-end | off | Always for paired-end -- AND add `--countReadPairs` |
| `--countReadPairs` | Count fragment (pair) as one | off in v2.0.2+ | Always for paired-end (post-v2.0.2 API change) |
| `-s 0|1|2` | Strandedness | 0 (unstranded) | `-s 2` for TruSeq Stranded; `-s 1` for forward kits |
| `-O` | Allow multi-overlap | off | Off for typical DGE; on creates double-counts in overlapping genes |
| `-M` | Count multi-mappers | off | Off for DGE; on with `--fraction` for fractional counting |
| `-t` | Feature type | `exon` | Almost always `exon` |
| `-g` | Group attribute | `gene_id` | `gene_id` for gene-level; `transcript_id` is wrong (use Salmon for transcript-level) |

## STAR `--quantMode GeneCounts`

**Goal:** Build a count matrix from STAR's per-sample 4-column output, choosing the correct strandedness column.

**Approach:** Skip the first 4 summary rows; index column 0 (gene ID); pick column for strandedness; concat across samples.

```python
import pandas as pd
from pathlib import Path

def load_star_genecounts(filepaths, strandedness='reverse'):
    '''Load STAR ReadsPerGene.out.tab files.
    File columns (1-indexed): 1=gene_id, 2=unstranded, 3=forward, 4=reverse.
    After read_csv with index_col=0, the three remaining columns are 0=unstranded, 1=forward, 2=reverse.
    Illumina TruSeq Stranded is 'reverse'.
    '''
    col_map = {'unstranded': 0, 'forward': 1, 'reverse': 2}
    col_idx = col_map[strandedness]
    dfs = {}
    for fp in filepaths:
        sample = Path(fp).name.replace('_ReadsPerGene.out.tab', '')
        df = pd.read_csv(fp, sep='\t', header=None, index_col=0)
        dfs[sample] = df.iloc[4:, col_idx]
    return pd.DataFrame(dfs)
```

Strandedness verification: for a stranded library, the "wrong" strand column total should be <5% of the "right" strand column total. If comparable, the library was unstranded or there was a kit / config mix-up.

```bash
infer_experiment.py -r annotation.bed -i sample.bam
```

## Salmon / kallisto via tximport

**Goal:** Import transcript-level Salmon or kallisto quantifications to gene-level counts with length-bias correction.

**Approach:** Build a `tx2gene` mapping; call `tximport()` with the right `type` and `countsFromAbundance`; hand the result to DESeq2 or limma-voom.

```r
library(tximport)
library(DESeq2)

tx2gene <- read.csv('tx2gene.csv')

files <- file.path('salmon_out', samples$id, 'quant.sf')
names(files) <- samples$id

txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
```

The naive alternative -- summing NumReads to gene level -- is wrong when isoform usage varies across samples. NumReads is normalized against EffectiveLength per transcript; sum-then-DE introduces a length-by-condition bias indistinguishable from differential expression. tximport handles this by carrying the per-sample average transcript length matrix as a DESeq2/edgeR offset.

For limma-voom (no offset mechanism):

```r
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene,
                countsFromAbundance = 'lengthScaledTPM')
y <- DGEList(counts = txi$counts)
v <- voom(y, design, plot = TRUE)
```

For DTU (DRIMSeq, DEXSeq), use `dtuScaledTPM` with `txOut=TRUE`.

For 3'-tagged libraries (10x Chromium 3', QuantSeq), length bias is negligible; use `countsFromAbundance='no'` WITHOUT the length offset (`tximeta`'s argument or manual disable).

## tximeta -- Automatic Provenance

**Goal:** Import Salmon/kallisto with automatic linkage to the exact annotation release used in the index.

**Approach:** `tximeta` (Love, Soneson, Hickey et al. 2020 *PLoS Comp Biol* 16:e1007664) inspects the Salmon index hash and pulls matching annotation and metadata; the resulting `SummarizedExperiment` carries the provenance.

```r
library(tximeta)

coldata <- data.frame(
    names = samples$id,
    files = file.path('salmon_out', samples$id, 'quant.sf'),
    condition = samples$condition
)

se <- tximeta(coldata)
gse <- summarizeToGene(se)

library(DESeq2)
dds <- DESeqDataSet(gse, design = ~ condition)
```

For new projects, `tximeta` is the modern preference over hand-managed `tx2gene` -- it eliminates a class of "wrong annotation" bugs.

## HTSeq-count

```python
import pandas as pd
from pathlib import Path

def load_htseq_counts(filepaths):
    '''Load HTSeq count files, dropping the __no_feature etc. summary rows.'''
    dfs = {}
    for fp in filepaths:
        sample = Path(fp).stem.replace('_counts', '')
        df = pd.read_csv(fp, sep='\t', header=None, index_col=0,
                         names=['gene', 'count'])
        df = df[~df.index.str.startswith('__')]
        dfs[sample] = df['count']
    return pd.DataFrame(dfs)
```

HTSeq overlap modes (Anders, Pyl, Huber 2015 *Bioinformatics* 31:166):

| Mode | Behavior |
|------|----------|
| `union` (default) | Read assigned to union of overlapping features; >1 gene -> `__ambiguous` |
| `intersection-strict` | Every base of read must overlap the same single feature |
| `intersection-nonempty` | Intersection across positions; nonempty -> that gene wins |

`-s yes|no|reverse`: TruSeq Stranded is `reverse`. The `--mode` and `-s` flags must match the library; mismatches are silent miscounts.

The `__no_feature` and `__ambiguous` totals are QC signals. If `__no_feature` > 30%: annotation incomplete, wrong reference, or strandedness misspecified. If `__ambiguous` > 10%: many overlapping gene annotations (common with comprehensive GTFs); consider `intersection-nonempty`.

## RSEM expected_count

```r
library(tximport)
library(DESeq2)

files <- file.path('rsem_out', paste0(samples$id, '.genes.results'))
names(files) <- samples$id
txi <- tximport(files, type = 'rsem', txIn = FALSE, txOut = FALSE)

txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
```

RSEM `expected_count` is non-integer (EM-derived). DESeq2 requires integers, but `DESeqDataSetFromTximport` rounds appropriately. The `txi$length == 0` substitution handles the "Error: all(lengths > 0) is not TRUE" panic from rRNA-filtered or zero-length transcripts.

Avoid `round()` + `DESeqDataSetFromMatrix` -- that path loses the length correction.

## kallisto + sleuth (Uncertainty-Aware DE)

```bash
kallisto quant -i index -o sample1_out -b 100 -t 8 read1_1.fq.gz read1_2.fq.gz
```

```r
library(sleuth)

s2c <- data.frame(sample = samples$id,
                  condition = samples$condition,
                  path = file.path('kallisto_out', samples$id))

so <- sleuth_prep(s2c, ~ condition,
                  target_mapping = tx2gene,
                  aggregation_column = 'gene_id',
                  gene_mode = TRUE)
so <- sleuth_fit(so, ~ condition, 'full')
so <- sleuth_fit(so, ~ 1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results <- sleuth_results(so, 'reduced:full')
```

Sleuth (Pimentel et al. 2017 *Nat Methods* 14:687) uses the 100 bootstrap replicates from kallisto to decouple BIOLOGICAL variance (between-replicate) from INFERENTIAL variance (within-replicate, from EM uncertainty). For transcript-level analyses where quantification uncertainty matters (similar isoforms, low coverage), sleuth's response error linear model is more conservative than DESeq2 on plain counts. For well-quantified gene-level DGE, sleuth and DESeq2-via-tximport converge.

## Salmon Selective Alignment

Selective alignment (Srivastava 2020 *Genome Biol* 21:239) -- which combines fast pseudo-mapping with traditional alignment of seed extensions to improve quantification accuracy for transcripts with sequence similarity -- has been the **default** since Salmon 1.0.0. The historical `--validateMappings` flag that explicitly requested it is now a deprecated no-op; passing it does nothing.

```bash
salmon quant -i index -l A \
    -1 read_1.fq.gz -2 read_2.fq.gz \
    --gcBias --seqBias \
    -p 8 -o sample_out
```

`--gcBias` corrects fragment GC bias; `--seqBias` corrects random hexamer priming bias. Both are recommended for any Salmon run from RNA-seq with biological condition variation. They affect EffectiveLength estimates and propagate through tximport.

## 10X Genomics

```python
import scanpy as sc

adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')
adata = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
```

```r
library(Seurat)
mat <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
mat <- Read10X_h5('filtered_feature_bc_matrix.h5')
```

10X convention: cells in rows (AnnData) or cells in columns (Seurat after Read10X). The R/Python convention differs by transpose.

## Annotation Pre-Filtering -- rRNA, Mt-rRNA, Pseudogenes

**Goal:** Drop biotypes that distort downstream normalization or are irrelevant to the question.

**Approach:** Filter on `gene_biotype` from the GTF (Ensembl) or `gene_type` (GENCODE).

```r
library(rtracklayer)
gtf <- import('Homo_sapiens.GRCh38.110.gtf.gz')
gene_info <- as.data.frame(gtf[gtf$type == 'gene',
                                c('gene_id', 'gene_name', 'gene_biotype')])

keep_biotypes <- c('protein_coding', 'lncRNA',
                   'IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene',
                   'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene')
keep_genes <- gene_info$gene_id[gene_info$gene_biotype %in% keep_biotypes]
counts_filt <- counts[rownames(counts) %in% keep_genes, ]
```

Mt-encoded protein-coding genes (MT-CO1, MT-ND1, ...) are a judgment call: include for tissue-specific work (heart, muscle); drop when mitochondrial fraction varies with cell stress and could confound normalization.

In single-cell, mitochondrial percentage is a STANDARD QC metric (cells with high %mito are stressed/dying); see `single-cell/preprocessing`.

## GENCODE vs Ensembl -- What Differs at the File Level

| Difference | Ensembl | GENCODE |
|------------|---------|---------|
| chromosome naming | `1`, `2`, ..., `MT` | `chr1`, `chr2`, ..., `chrM` |
| PAR gene encoding (releases 25-43) | Once on chrX | Both chrX and chrY with `_PAR_Y` suffix on chrY copies |
| PAR gene encoding (releases 44+ / Ensembl 110+) | Once on chrX | chrY copies get their own ENSG accessions; `_PAR_Y` retired |
| Subset releases | Single release | `basic` (high-confidence subset) and `comprehensive` (everything) |

CRITICAL: code that strips Ensembl version suffixes via `sub('\\..*', '', x)` ALSO strips the `_PAR_Y` tag in GENCODE 25-43, collapsing the chrY duplicate onto the chrX gene -- silently introducing duplicate row indices. Use `sub('\\.[0-9]+(_PAR_Y)?$', '\\1', x)` to preserve.

BAM files aligned to GENCODE (`chr1`) cannot be quantified against Ensembl GTF (`1`) without rename. Mismatched naming causes `__no_feature` to dominate.

## Filter Low-Count Genes

```r
library(edgeR)

y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y, design = model.matrix(~ condition, coldata))
y <- y[keep, , keep.lib.sizes = FALSE]
```

```python
min_counts, min_samples = 10, 3
expressed = (counts >= min_counts).sum(axis=1) >= min_samples
counts_filt = counts.loc[expressed]
```

See `differential-expression/edger-basics` for `filterByExpr` semantics. DESeq2 has automatic independent filtering at `results()` time; manual pre-filter is speed-only.

## Per-Method Failure Modes

### featureCounts paired-end double-counted

**Trigger:** Pipeline written against featureCounts pre-v2.0.2 used `-p` only; rerun on Subread 2.0.6 produces counts ~2x expected.

**Mechanism:** Post-v2.0.2 `-p` flags paired-end but does NOT count pairs as fragments. Each mate is counted separately.

**Symptom:** Library sizes ~2x what RNA-seq QC reports; downstream CPM compressed; "more reads than expected" panic.

**Fix:** Add `--countReadPairs`. Re-run featureCounts.

### STAR wrong strandedness column

**Trigger:** TruSeq Stranded library quantified via STAR; user reads column 3 (forward); ~95% of reads dropped.

**Mechanism:** TruSeq Stranded is reverse-stranded (column 4). Column 3 is forward-stranded; for a reverse library, almost no reads align in the forward orientation.

**Symptom:** Library sizes ~5% of expected; almost no DE detectable; very low gene detection rate.

**Fix:** Read column 4. Verify with `infer_experiment.py`.

### Salmon NumReads summed naively -> length-by-condition bias

**Trigger:** Salmon output read without tximport; per-transcript NumReads summed to gene level via groupby.

**Mechanism:** NumReads is normalized against per-transcript EffectiveLength. Summing ignores that. If treatment shifts isoform usage from short to long, the gene appears upregulated even with constant total mRNA.

**Symptom:** DE genes overlap with known isoform-switching genes (e.g., during development); fold changes don't replicate at the protein level.

**Fix:** Use `tximport` (or `catchSalmon` in edgeR) to carry the length matrix as a DE offset.

### RSEM zero-length transcript breaks tximport

**Trigger:** `tximport(files, type='rsem')` errors out with "all(lengths > 0) is not TRUE".

**Mechanism:** RSEM reports `effective_length = 0` for certain very-short or filtered transcripts.

**Symptom:** Pipeline halts at import.

**Fix:** `txi$length[txi$length == 0] <- 1` after `tximport`; or filter `tx2gene` to exclude affected transcripts.

### `_PAR_Y` stripped, chrY duplicates lost

**Trigger:** GENCODE v40 quantification; `rownames(counts) <- sub('\\..*', '', rownames(counts))`; duplicate row indices.

**Mechanism:** Default regex strips `_PAR_Y` along with the version suffix; chrY PAR copies collapse onto chrX gene IDs.

**Symptom:** Duplicate row warnings; sample-specific counts double for affected genes; downstream `aggregate(... ~ rownames)` adds them.

**Fix:** Use the regex that preserves `_PAR_Y`: `sub('\\.[0-9]+(_PAR_Y)?$', '\\1', x)`. Or upgrade to GENCODE 44+ where the issue is gone.

### HTSeq mode mismatched to library

**Trigger:** Mouse RNA-seq with `--mode intersection-strict`; many overlapping-gene reads dropped; library size 30% lower than featureCounts on the same data.

**Mechanism:** `intersection-strict` requires every base to overlap the same feature; reads crossing exon-intron boundaries or overlapping gene boundaries get dropped.

**Symptom:** Lower count totals than other quantifiers; `__no_feature` and `__ambiguous` high.

**Fix:** `--mode union` (default) or `--mode intersection-nonempty` for richer recovery.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `Error: all(lengths > 0) is not TRUE` | RSEM zero-length transcripts | `txi$length[txi$length == 0] <- 1` |
| Duplicate row indices after version strip | `_PAR_Y` collapsed | Use the preserving regex |
| `__no_feature` >30% | Wrong reference; wrong strandedness; chrN naming mismatch | Verify reference + strandedness; check chr vs N naming |
| `counts matrix should be integers` | Salmon/RSEM raw counts to DESeq2 | Use `DESeqDataSetFromTximport`; or round (loses length correction) |
| Library sizes ~2x expected for paired-end | featureCounts v2.0.2 `-p` without `--countReadPairs` | Add `--countReadPairs` |
| TPM column has Inf values | Zero-length features | Filter `Length > 0` before TPM computation |
| Salmon TPM doesn't sum to 1e6 | Pre-summed across samples or filtered subset | Recompute TPM after filter; or use original abundance column |

## References

- Liao Y, Smyth GK, Shi W. 2014. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30(7):923-930. doi:10.1093/bioinformatics/btt656
- Anders S, Pyl PT, Huber W. 2015. HTSeq -- a Python framework to work with high-throughput sequencing data. *Bioinformatics* 31(2):166-169. doi:10.1093/bioinformatics/btu638
- Dobin A et al. 2013. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29(1):15-21. doi:10.1093/bioinformatics/bts635
- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. 2017. Salmon provides fast and bias-aware quantification of transcript expression. *Nat Methods* 14(4):417-419. doi:10.1038/nmeth.4197
- Srivastava A, Malik L, Sarkar H, Zakeri M, Almodaresi F, Soneson C, Love MI, Kingsford C, Patro R. 2020. Alignment and mapping methodology influence transcript abundance estimation. *Genome Biol* 21:239. doi:10.1186/s13059-020-02151-8
- Bray NL, Pimentel H, Melsted P, Pachter L. 2016. Near-optimal probabilistic RNA-seq quantification. *Nat Biotechnol* 34(5):525-527. doi:10.1038/nbt.3519
- Pimentel H, Bray NL, Puente S, Melsted P, Pachter L. 2017. Differential analysis of RNA-seq incorporating quantification uncertainty. *Nat Methods* 14(7):687-690. doi:10.1038/nmeth.4324
- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Res* 4:1521. doi:10.12688/f1000research.7563.2
- Love MI, Soneson C, Hickey PF, Johnson LK, Pierce NT, Shepherd L, Morgan M, Patro R. 2020. Tximeta: Reference sequence checksums for provenance identification in RNA-seq. *PLoS Comput Biol* 16(2):e1007664. doi:10.1371/journal.pcbi.1007664
- Frankish A et al. 2021. GENCODE 2021. *Nucleic Acids Res* 49(D1):D916-D923. doi:10.1093/nar/gkaa1087

## Related Skills

- gene-id-mapping - Building tx2gene; Ensembl version stripping with PAR_Y preservation
- normalization - Composition bias, TMM/RLE, the "biotype filter before normalize" rationale
- metadata-joins - Sample name reconciliation across counts and metadata
- sparse-handling - Single-cell sparse matrix formats
- differential-expression/deseq2-basics - tximport -> DESeqDataSetFromTximport workflow
- differential-expression/edger-basics - catchSalmon for transcript-level DTE
- differential-expression/batch-correction - Batch covariate vs subtraction
- rna-quantification/featurecounts-counting - Detailed featureCounts run patterns
- rna-quantification/alignment-free-quant - Salmon/kallisto invocation details
- rna-quantification/tximport-workflow - Detailed tximport workflow
- read-alignment/star-alignment - STAR upstream of GeneCounts
- read-qc/quality-reports - RIN, DV200 inputs to consider
