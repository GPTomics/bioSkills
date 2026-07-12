---
name: bio-single-cell-splicing
description: Analyzes alternative splicing at single-cell resolution. The first decision is library chemistry — 10X 3' is fundamentally limited (RT primes from poly-A, R2 falls in 3' UTR, <0.1 junction read per cell per AS event). Plate-based full-length methods (Smart-seq3, FLASH-seq, VASA-seq, STORM-seq) and single-cell long-read (MAS-Iso-seq, scISOr-Seq2) are the chemistries that give per-cell isoform structure. Tools include MARVEL (R, Smart-seq integrated), BRIE2 (Bayesian PSI with regulatory features and ELBO_gain test), scQuint (junction-cluster, plate-based; not for 10X), SpliZ (annotation-free Z-score), Psix (graph-smoothness regulated AS), and Sierra (alternative polyadenylation, often confused with AS). Use when analyzing isoform usage in scRNA-seq, identifying cell-type-specific splicing, or determining whether scRNA-seq chemistry supports splicing analysis at all.
tool_type: python
primary_tool: MARVEL
---

## Version Compatibility

Reference examples tested with: MARVEL 2.0+, BRIE2 0.2.4+, scQuint 0.1+, SpliZ 0.0.1+, Sierra 1.0+, Psix 0.1+, anndata 0.10+, scanpy 1.10+, pandas 2.2+, scipy 1.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Single-Cell Splicing Analysis

The fundamental decision is **chemistry**, not tool. Most droplet 3' scRNA-seq cannot support transcriptome-wide splicing inference because reverse transcription primes from the poly(A) tail and most reads land in the 3' UTR — far from CDS-region splicing events. Plate-based full-length methods and single-cell long-read sequencing are the chemistries that give per-cell isoform structure across the gene body.

## The 10X 3' Problem (Quantified)

Three compounding mechanisms make 10X Chromium 3' (v3.1, GEM-X, v4) hostile to splicing:

1. **3' enrichment**: median fragment <1 kb from poly(A); >70% of unique reads fall within 3' UTR.
2. **Short R2 (~91 nt)**: each read straddles at most one junction; usually none, because R2 lands in 3' UTR.
3. **PCR concatemers and TSO artifacts**: pollute junction detection; UMI collapse is gene-level, not isoform-level.

**Quantitative estimate:** Only a small fraction of cassette exons sit close enough to the polyA site to be sampled by 3' chemistry (empirical estimates from APA/3'-end atlases — see Tian & Manley 2017 *Nat Rev Mol Cell Biol* for the 3' UTR isoform landscape). Effective junction read yield from 10X 3' is **<0.1 per cell per AS event** — vs the 5-10 needed for stable per-cell PSI. Most splicing analyses on 10X 3' data report artifacts.

**The 5' kit (10X 5' GEX) does not solve this** — it shifts capture from 3' UTR to 5' UTR / TSS-proximal regions. Marginal improvement; not a transcriptome-wide solution. Note that V(D)J recovery requires the **10X Chromium Single Cell Immune Profiling kit** (with TCR/BCR-specific enrichment), not 5' GEX alone — postdocs designing immune-repertoire experiments must use the dedicated V(D)J kit.

## Decision: Does the Chemistry Support Splicing Analysis?

| Chemistry | Splicing analysis viable? | Best alternative if no |
|-----------|----------------------------|--------------------------|
| 10X 3' (Chromium v3, GEM-X, v4, Flex) | No (transcriptome-wide); maybe near-3'-end events | Sierra for APA |
| 10X 5' GEX | Limited; near-5'-end events only | Sierra for alternative TSS; switch to MAS-Iso-seq |
| Smart-seq2 | Yes (full transcript) | MARVEL or BRIE2 |
| Smart-seq3 / Smart-seq3xpress | Yes + UMI molecule counting | MARVEL or BRIE2 |
| FLASH-seq | Yes (faster, cheaper Smart-seq3) | MARVEL or BRIE2 |
| VASA-seq | Yes + total RNA (incl. nascent, IR) | MARVEL with IR analysis |
| STORM-seq | Yes + total RNA + ribodepletion | MARVEL with IR analysis |
| MAS-Iso-seq + 10X 5' (PacBio Kinnex) | Yes — full isoforms per cell | FLAMES, scNanoGPS, IsoQuant, see long-read-splicing |
| scISOr-Seq2 (PacBio + 10X) | Yes — full isoforms with cell-typing | FLAMES, IsoQuant |
| ONT direct cDNA scRNA | Yes | FLAMES |
| ONT direct RNA scRNA | Yes + native modifications | FLAMES |

## Tool Selection Matrix

| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| MARVEL | Smart-seq plate-based and (v2+) 10X droplet unified workflow | Plate or droplet BAMs + Seurat | SE/A5SS/A3SS/MXE/RI/AFE/ALE; modality classification; native Seurat integration; v2 droplet support | R-only |
| BRIE2 | Plate-based with regulatory feature prior | Plate BAM + GFF3 events | Bayesian variational PSI + ELBO_gain test; principled uncertainty; CLI-driven (`brie-count`, `brie-quant`) | TensorFlow dependency; slow at scale |
| scQuint | Plate-based annotation-free junction-cluster quantification (validated on Smart-seq2) | STAR junctions across cells | Cluster-level junction usage; latent Dirichlet | Authors recommend AGAINST use on 10X 3'/5' data (3'-bias confounds); plate-based only |
| SpliZ | Annotation-free discovery of cell-state-associated splicing | STAR-aligned BAMs | Per-gene Z-score; no event database needed | Annotation-free = power tradeoff |
| Psix | Regulated AS along trajectories | PSI matrix + kNN graph | Tests graph smoothness; robust to dropout | Needs cell-state graph upstream |
| Sierra | APA in 10X 3' (NOT splicing) | 10X BAM + GTF | Peak-calling 3' ends; DEXSeq DTU on UTR isoforms | APA only; not for cassette exons |
| pseudobulk leafcutter / rMATS | Between-cell-type differential splicing | Aggregated BAMs | Bulk-level statistical power | Loses within-cluster heterogeneity |
| MAS-Iso-seq + FLAMES | Full-length single-cell isoforms | 10X 5' + PacBio Kinnex | Full isoforms per cell at scale | Cost; complex pipeline |

## Decision Tree by Goal

| Goal | Recommended approach |
|------|----------------------|
| "Will my 10X 3' data support splicing?" | No transcriptome-wide; consider Sierra for APA. Note: scQuint authors recommend against use on 10X data |
| Cassette exon analysis in cell types from Smart-seq2 | MARVEL with `ComputePSI` + `AssignModality` + `CompareValues` |
| Discover cell-state-associated splicing without an event database | SpliZ |
| Test regulated AS along developmental pseudotime | Psix |
| Per-cell PSI with uncertainty in low-coverage cells | BRIE2 |
| Differential splicing between two well-defined cell types | Pseudobulk leafcutter or rMATS on aggregated BAMs |
| APA (alternative polyadenylation, often confused with AS) | Sierra |
| Full-length single-cell isoforms at scale | MAS-Iso-seq + FLAMES (long-read) |
| Microexons (3-27 nt) | Long-read or aligner with low overhang (uLTRA, deSALT) |
| snRNA-seq (nuclei) — IR question | Library captures nuclear RNA enriched for incomplete splicing — interpret IR cautiously |

## MARVEL Plate-Based Workflow

**Goal:** Run a unified workflow from STAR junctions to cell-type-specific splicing calls.

**Approach:** Build a wide splice-junction count matrix (rows = junctions keyed by `coord.intron`, columns = cells), assemble per-event feature tables, then construct MARVEL object with named slots (`SpliceJunction`, `SplicePheno`, `SpliceFeature`, `IntronCounts`, `GeneFeature`, `Exp`, `GTF`). Quantify PSI per event class, classify modality, test differential splicing.

```r
library(MARVEL); library(Seurat); library(data.table)

seurat_obj <- readRDS('cells.rds')

# Build wide SJ matrix: first column 'coord.intron' (e.g. 'chr1:100007082:100022621'),
# subsequent columns are per-cell sample IDs with junction counts as values.
# This is constructed from STAR SJ.out.tab files (one per cell) merged on intron coord.
sj_files <- list.files('star_pass2/', pattern='SJ.out.tab$', full.names=TRUE)
sj_long <- rbindlist(lapply(sj_files, function(f) {
    d <- fread(f, sep='\t', header=FALSE,
               col.names=c('chr','start','end','strand','motif','annot','unique','multi','overhang'))
    d$coord.intron <- paste(d$chr, d$start, d$end, sep=':')
    d$sample <- gsub('_SJ.out.tab$', '', basename(f))
    d[, .(coord.intron, sample, unique)]
}))
sj <- dcast(sj_long, coord.intron ~ sample, value.var='unique', fill=0)

# SpliceFeature is a NAMED LIST keyed by event class
df.feature.list <- list(
    SE   = read.table('events_SE.txt',   header=TRUE, sep='\t'),
    A5SS = read.table('events_A5SS.txt', header=TRUE, sep='\t'),
    A3SS = read.table('events_A3SS.txt', header=TRUE, sep='\t'),
    MXE  = read.table('events_MXE.txt',  header=TRUE, sep='\t'),
    RI   = read.table('events_RI.txt',   header=TRUE, sep='\t')
)

# SplicePheno: per-cell metadata; sample.id column maps to SpliceJunction column names
df.pheno <- seurat_obj@meta.data
df.pheno$sample.id <- rownames(df.pheno)

marvel <- CreateMarvelObject(
    SpliceJunction = sj,
    SplicePheno    = df.pheno,
    SpliceFeature  = df.feature.list,
    GeneFeature    = read.table('gene_features.tsv', header=TRUE, sep='\t'),
    Exp            = read.table('tpm.tsv', header=TRUE, sep='\t', row.names=1),
    GTF            = rtracklayer::import('annotation.gtf')
)

marvel <- ComputePSI(marvel, CoverageThreshold=10, EventType='SE')
marvel <- AssignModality(marvel, EventType='SE')
marvel <- CompareValues(
    marvel,
    cell.group.g1 = neurons, cell.group.g2 = glia,
    method = 'wilcox', n.cells = 25, psi.delta = 0.1
)
```

For 10X droplet data, MARVEL v2+ provides `CreateMarvelObject.10x()` and `AnnotateSJ.10x()` constructors. Verify the exact API via `?CreateMarvelObject.10x` in installed MARVEL.

MARVEL classifies events into modalities (Song 2017 *Mol Cell*): included (PSI~1), excluded (PSI~0), bimodal (mixture at 0/1), middle (peaked ~0.5), multimodal. Bimodality usually reflects mixed cell states or stochastic monoallelic-like bursting. Mid-modality (peaked at 0.5) can be technical (mixed cells in a droplet) — confirm with full-length data.

## BRIE2 Bayesian PSI

**Goal:** Estimate per-cell PSI with informative regulatory-feature prior; test cell-state association via likelihood-ratio testing on covariate effects.

**Approach:** BRIE2 is a CLI-driven workflow (`brie-count` for read counting, `brie-quant` for variational inference + LRT). Prepare a GFF3 of splicing events, count cell-barcoded junction reads, then fit the model with covariate testing.

```bash
# 1. Count splicing events per cell
brie-count \
    -a splicing_events.gff3 \
    -S sample_list.tsv \
    -o brie_counts/ \
    -p 16

# 2. Fit BRIE2 with LRT against the cell-type covariate
brie-quant \
    -i brie_counts/brie_count.h5ad \
    -c cell_metadata.tsv \
    -o brie_quant.h5ad \
    --interceptMode gene \
    --LRTindex All \
    --testBase null \
    --MCsize 3 \
    --batchSize 1000000 \
    -p 16
```

`--interceptMode gene` fits a gene-specific intercept (recommended); `--LRTindex All` tests all covariates; `--testBase null` uses the null model as the LRT reference. Verify exact flag set via `brie-quant -h` in installed BRIE2.

```python
import scanpy as sc

adata_splice = sc.read_h5ad('brie_quant.h5ad')
# Per-event covariate effects, ELBO values, and LRT statistics live in
# adata_splice.varm and adata_splice.var; column names depend on BRIE2 version.
# Inspect with: print(adata_splice); print(adata_splice.varm.keys())
# Per-event significance is typically derived from LRT delta-ELBO.
```

BRIE2 (Huang & Sanguinetti 2021 *Genome Biol*) uses a sequence-derived feature prior (exon length, GC content, splice site strength, motif counts) to regularize PSI estimates in low-coverage cells. The LRT-based covariate test answers "is this event associated with cell state?" without requiring per-cell PSI accuracy. Threshold the delta-ELBO at ~3 (analogous to log-Bayes-factor); confirm against version-specific output keys via the brie-tutorials repo.

## SpliZ for Annotation-Free Discovery

**Goal:** Identify splicing-defined cell populations without an event database.

**Approach:** Compute per-gene splicing Z-score across cells; test for cell-state association via permutation.

```bash
# SpliZ is a Nextflow pipeline (not a standalone CLI). Configure inputs in a .config
# file (dataname, input_file, libraryType, grouping_level_1/2) - either SICILIAN
# output (SICILIAN=true) or BAMs via a samplesheet CSV + metadata + GTF (SICILIAN=false).
nextflow run salzmanlab/spliz -r main -latest -c spliz.config
```

SpliZ (Olivieri 2022 *Nat Methods*) is robust to dropout because it pools junction information across the gene; particularly useful for discovering splicing diversity in heterogeneous tumor samples.

## Psix for Regulated AS Along Trajectories

**Goal:** Detect AS that varies coherently with cell state along a developmental trajectory, robust to dropout.

**Approach:** Score whether observed PSI is smooth on the cell-cell kNN graph from expression-space embedding.

```python
import psix
import scanpy as sc

adata = sc.read_h5ad('cells.h5ad')
sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca')

psix_obj = psix.Psix(adata, psi_matrix_path='psi_matrix.tsv')
psix_obj.run_psix()

regulated = psix_obj.psix_results.query('psix_score > 1.5 and pvalue < 0.05')
```

Psix (Buen Abad Najar 2022 *Genome Res* 32:1385) is the principled alternative to imputing PSI: do not impute (it obliterates heterogeneity); test for graph smoothness instead.

## Sierra for APA (Not Splicing)

**Goal:** Detect alternative polyadenylation in 10X 3' data — frequently confounded with AS.

**Approach:** Peak-call read pile-ups at 3' ends, then DEXSeq-style DTU on 3' UTR isoforms.

```r
library(Sierra)

peak_file <- FindPeaks(
    output.file = 'peaks.txt',
    gtf.file = 'annotation.gtf',
    bam.file = 'possorted_genome_bam.bam'
)

counts <- CountPeaks(
    peak.sites.file = 'peaks.txt',
    gtf.file = 'annotation.gtf',
    bamfile = 'possorted_genome_bam.bam',
    whitelist.file = 'barcodes.tsv'
)

# CountPeaks returns a peak x cell matrix; annotate it and build a peak Seurat
# object before differential-usage testing.
peak.annotations <- AnnotatePeaksFromGTF(
    peak.sites.file = 'peaks.txt',
    gtf.file = 'annotation.gtf',
    output.file = 'peak_annotations.txt'
)

peaks.seurat <- NewPeakSeurat(
    peak.data = counts,
    annot.info = peak.annotations,
    cell.idents = cell_identities
)

apa_results <- DUTest(peaks.seurat, population.1 = ctrl_cells, population.2 = trt_cells)
```

If only 10X 3' data is available, this is often what is actually wanted. Distinct UTRs change miRNA targeting, RBP binding, and stability — biologically meaningful but not splicing.

## Pseudobulk for Statistical Power

**Goal:** Recover bulk-level statistical power for differential splicing between cell types.

**Approach:** Sum junction counts across cells of the same cluster, then run leafcutter / rMATS on aggregated counts.

```python
import pandas as pd
import numpy as np

def pseudobulk_junctions(junction_counts, cell_metadata, groupby='cell_type'):
    out = {}
    for group, cells in cell_metadata.groupby(groupby).groups.items():
        mask = junction_counts.columns.isin(cells)
        out[group] = junction_counts.loc[:, mask].sum(axis=1)
    return pd.DataFrame(out)
```

Use pseudobulk for differential splicing **between** well-defined cell types; use per-cell methods for **within-population heterogeneity** (graded splicing along pseudotime, bimodal cell-state mixtures).

## Single-Cell Long-Read = Future of Single-Cell Splicing

In 2024-2026, full-length single-cell long-read sequencing has become practical and is the recommended chemistry for splicing-focused single-cell experiments:

- **MAS-Iso-seq / PacBio Kinnex**: concatenated full-length cDNA arrays, ~16x throughput vs plain Iso-Seq, compatible with 10X 5' libraries (Al'Khafaji 2024 *Nat Biotech*)
- **scISOr-Seq2**: hybrid 10X + PacBio for cell typing + isoform structure (Joglekar et al 2024 *Nat Neurosci* 27:1051-1063, single-cell long-read brain isoform mapping)
- **ONT direct cDNA + 10X**: lower cost, similar information content
- **FLAMES**: barcode demultiplexing + isoform quantification + SNV calling for ONT scRNA (Tian 2021 *Genome Biol* 22:310)

For splicing-specific full-length single-cell analysis, see `long-read-splicing` skill.

## Per-Tool Failure Modes

### MARVEL: SpliceJunction Matrix Format

**Trigger:** Building the SpliceJunction matrix from STAR SJ.out.tab incorrectly (e.g. long-format instead of wide).

**Mechanism:** MARVEL plate-based `CreateMarvelObject(SpliceJunction = ...)` expects a **wide matrix** with first column `coord.intron` (formatted `chr:start:end`) and subsequent columns being per-cell sample IDs with integer junction counts. Long-format data.frames or missing `coord.intron` column cause runtime errors.

**Symptom:** "no `coord.intron` column found" errors; or empty PSI tables despite junction reads being present.

**Fix:** Verify wide-matrix structure; ensure SJ.out.tabs are merged on the `chr:start:end` key with cells as columns. Use `data.table::dcast` for the long->wide reshape.

### BRIE2: TensorFlow Memory

**Trigger:** Large cohort (>10k cells) with deep coverage.

**Mechanism:** Variational inference loads full count matrix; TensorFlow allocates GPU memory aggressively.

**Symptom:** OOM kills; training stalls.

**Fix:** Reduce `--batchSize` from default (500000) to 100000 or 50000; train per-chromosome batch; use CPU mode for very small cohorts. Note flag is camelCase `--batchSize`, not `--batch_size`.

### scQuint: 3' Data Sparsity

**Trigger:** Running scQuint on 10X 3' v3 data hoping for splicing signal.

**Mechanism:** scQuint's latent Dirichlet model needs junction counts; 10X 3' yields too few junction reads to fit the model robustly.

**Symptom:** All cells assign to one cluster; no informative splicing signal.

**Fix:** Pivot to APA analysis with Sierra; or upgrade chemistry to MAS-Iso-seq.

### Psix: Missing kNN Graph

**Trigger:** Running Psix without precomputed cell-cell graph.

**Mechanism:** Psix tests PSI smoothness on a pre-existing cell-cell graph; without one, no smoothness statistic.

**Symptom:** Empty results or error about missing `connectivities`.

**Fix:** Run `sc.pp.neighbors(adata)` before Psix; ensure `connectivities` is in `adata.obsp`.

### Sierra: Annotation Gaps

**Trigger:** GTF missing 3'UTR annotations.

**Mechanism:** Sierra peak-calls within annotated 3'UTRs; missing annotations mean missed peaks.

**Symptom:** Few peaks detected; gene-level coverage but no APA calls.

**Fix:** Use comprehensive GENCODE annotation; or run de-novo peak calling first.

## Reconciliation: When Single-Cell Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MARVEL sig, BRIE2 not | Per-cell PSI noise (BRIE2 conservative); MARVEL pseudobulk-like | Trust MARVEL for cell-type comparisons; BRIE2 for within-cluster |
| BRIE2 sig, MARVEL not | Cell-state effect smoother than cell-type boundary | Test along trajectory with Psix |
| SpliZ sig, MARVEL not | Annotation-free SpliZ catches novel events | Investigate junction structure manually |
| Sierra sig, MARVEL not | Sierra is APA, MARVEL is splicing — different biology | Distinguish in interpretation |
| Pseudobulk sig, per-cell not | Power issue; effect averaged out per-cell | Report at cluster level, not per-cell |

## Quantitative Concepts Unique to Single-Cell

**Per-cell PSI vs pseudobulk PSI:**
- Per-cell PSI: meaningful only when junction coverage exceeds ~10-20 reads per cell per event (plate-based or long-read).
- Pseudobulk PSI: aggregate, recovers bulk-level statistical power, discards within-cluster heterogeneity.

**Modality detection in PSI distributions** (Song 2017 *Mol Cell*):
| Modality | PSI distribution | Biology |
|----------|------------------|---------|
| Included | Peaked at 1 | Constitutive inclusion |
| Excluded | Peaked at 0 | Constitutive skipping |
| Bimodal | Mixture at 0 and 1 | Mixed cell states or monoallelic-like bursting |
| Middle | Peaked ~0.5 | Often technical (well-contamination, doublets, or low-coverage shrinkage to prior); confirm with full-length |
| Multimodal | Multiple peaks | Complex regulation; deserves follow-up |

**Beta-binomial vs binomial models:** with sparse counts, binomial PSI is overdispersed. Beta-binomial models (BRIE2; leafcutter2 as Dirichlet-multinomial cluster-level) handle this. For very sparse droplet data, even beta-binomial fits poorly per cell — collapse to pseudobulk.

**Imputation pitfalls:** naive imputation (MAGIC, scImpute, ALRA) of expression matrices is **not** appropriate for PSI: imputing missing junction counts averages over neighboring cells and obliterates the very heterogeneity under study. Psix's approach — testing smoothness of observed PSI on the kNN graph — is the principled alternative.

## Cell-Type-Specific Splicing Biology

| System | Event | Regulator |
|--------|-------|-----------|
| Neural microexons | 3-27 nt exons enriched in brain | SRRM4/nSR100 (Irimia 2014 *Cell*); SRRM3 in retina/photoreceptors (Ciampi 2022 *PNAS*) |
| Neural differentiation | PTBP1 -> PTBP2 switch | miR-124 represses PTBP1; derepresses neural exons (Boutz 2007 *Genes Dev*) |
| T-cell activation | CD45 RA -> RO | hnRNP-L, ESRP-mediated |
| Erythropoiesis | EPB41 exon 16 | Splicing factor switching during maturation |
| Cardiac development | TTN N2BA -> N2B | MBNL1/CELF1 antagonism |
| EMT | FGFR2 IIIb -> IIIc, ENAH exon 11a | ESRP1/2 loss in mesenchymal state (Warzecha 2009 *Mol Cell*) |
| Activated T cell | CD45 isoform shift | Multiple SR/hnRNP regulators |

## Quality Thresholds

| Metric | Recommendation |
|--------|----------------|
| Cells per event with reads | >=50 (per-cell PSI); >=200 cells per cluster (pseudobulk) |
| Junction reads per event per cell | >=5 with coverage; <=1 = unreliable |
| PSI variance for cell-type call | <0.1 within cluster, >0.2 between clusters |
| Library | full-length plate or long-read for transcriptome-wide; 3' for APA only |
| Doublet filtering | Required before splicing analysis (DoubletFinder, Scrublet) |
| Cells per cluster (pseudobulk) | >=100 ideal; >=50 minimum |
| nuclear vs whole-cell | snRNA-seq enriches IR; treat with caution |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `MARVEL: ComputePSI returns empty` | STAR SJ.out.tab missing strand info | Re-run STAR with `--outSJtype Standard` |
| `brie.tl.fit: NaN loss` | Insufficient junction reads per cell | Filter cells with `min_reads=20`; raise threshold |
| `scQuint: convergence not reached` | LDA model fit on too-few junctions | Aggregate by chromosome; or switch chemistry |
| `Psix: missing connectivities` | Neighbors graph not computed | Run `sc.pp.neighbors(adata)` first |
| `Sierra: no peaks called` | GTF missing 3'UTR annotations | Use comprehensive GENCODE; or de-novo peak-call |
| `MARVEL: ggplot error` | Seurat version mismatch | Match MARVEL and Seurat versions |
| `FLAMES: barcode rescue failed` | Short-read 10X output not in expected directory | Verify cellranger output structure |

## Common Pitfalls

- **Treating 10X 3' splicing analysis as legitimate** — the chemistry doesn't support it. Use Sierra for APA or upgrade to MAS-Iso-seq.
- **Imputing PSI matrices** — destroys the heterogeneity to be detected. Use Psix or BRIE2 instead.
- **Per-cell PSI on droplet data** — typically too sparse for stable estimates. Use pseudobulk first, then drill down to per-cell.
- **Confusing APA with splicing** — Sierra results look like AS but are 3' UTR isoforms. Different machinery, different biology.
- **snRNA-seq IR signal misinterpreted as splicing dysregulation** — nuclear RNA is enriched for incompletely spliced transcripts; baseline IR is high.
- **Trusting per-cell PSI from BRIE2 without ELBO_gain test** — BRIE2's per-cell point estimates are noisy; the principled output is the ELBO_gain cell-state-association statistic.
- **Microexon analysis with default short-read aligners** — anchors >=20 nt miss most microexons; use VAST-TOOLS, MicroExonator, or long-read.
- **Skipping doublet filtering before splicing** — doublets create artificial PSI mid-modality.

## Related Skills

- single-cell/preprocessing - QC and normalization (must run before splicing)
- single-cell/clustering - Cell type annotation prerequisite
- single-cell/doublet-detection - Doublet filtering critical for splicing
- single-cell/data-io - h5ad / Seurat I/O
- splicing-quantification - Bulk RNA-seq comparison context
- long-read-splicing - Full-isoform analysis from MAS-Iso-seq, scISOr-Seq2; future of single-cell splicing

## References

- Huang & Sanguinetti 2021 *Genome Biol* - BRIE2
- Wen et al 2023 *Nucleic Acids Research* 51:e29 - MARVEL
- Benegas, Fischer & Song 2022 *eLife* - scQuint (annotation-free single-cell splicing analysis, validated on Smart-seq2)
- Olivieri et al 2022 *Nat Methods* - SpliZ
- Buen Abad Najar et al 2022 *Genome Research* 32:1385 - Psix
- Patrick et al 2020 *Genome Biol* - Sierra
- Song et al 2017 *Mol Cell* - splicing modality classification
- Picelli et al 2014 *Nat Protoc* - Smart-seq2
- Hagemann-Jensen et al 2020 *Nat Biotech* - Smart-seq3
- Hagemann-Jensen et al 2022 *Nat Biotech* - Smart-seq3xpress
- Hahaut et al 2022 *Nat Biotech* - FLASH-seq
- Salmen et al 2022 *Nat Biotech* - VASA-seq
- Johnson et al 2022 *bioRxiv* 10.1101/2022.03.14.484332 - STORM-seq (preprint)
- Al'Khafaji et al 2024 *Nat Biotech* - MAS-Iso-seq / Kinnex
- Tian et al 2021 *Genome Biology* 22:310 - FLAMES
- Joglekar et al 2024 *Nat Neurosci* 27:1051-1063 - scISOr-Seq2 single-cell brain isoform mapping
- Irimia et al 2014 *Cell* - neural microexons / SRRM4
- Ciampi et al 2022 *PNAS* 119:e2117090119 - SRRM3-dependent photoreceptor microexons
- Boutz et al 2007 *Genes Dev* - PTBP1/PTBP2 neural switch
- Tian & Manley 2017 *Nat Rev Mol Cell Biol* - alternative polyadenylation and 3' UTR isoforms
