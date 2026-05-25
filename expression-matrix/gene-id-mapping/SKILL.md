---
name: bio-expression-matrix-gene-id-mapping
description: Maps between gene identifier systems (Ensembl, Entrez, HGNC symbol, UniProt, RefSeq, MANE) using AnnotationDbi, biomaRt, mygene, pyensembl, and Ensembl REST. Encodes Ensembl version stripping with GENCODE _PAR_Y preservation, the Ziemann 2016 Excel autocorrect debacle and Bruford 2020 HGNC renames (SEPT*->SEPTIN*, MARCH*->MARCHF*, MARC*->MTARC*, DEC1->DELEC1), OCT4/POU5F1 alias resolution, biomaRt archive endpoints for release pinning, the `filters` (plural) gotcha, MANE Select for clinical reporting, cross-species orthology via Ensembl Compara / OMA / OrthoDB, and tx2gene construction for tximport. Use when converting gene IDs across systems, handling renamed symbols, building tx2gene, pinning to a specific Ensembl release for reproducibility, or mapping cross-species orthologs.
tool_type: mixed
primary_tool: biomaRt
---

## Version Compatibility

Reference examples tested with: biomaRt 2.58+, AnnotationDbi 1.66+, org.Hs.eg.db 3.18+, org.Mm.eg.db 3.18+, GenomicFeatures 1.54+, mygene 1.38+ (Python), pyensembl 2.3+, pandas 2.2+, rtracklayer 1.62+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Gene ID Mapping

**"Convert gene IDs from X to Y"** -> Query the appropriate annotation source (local org.db for speed, biomaRt for Ensembl-specific attributes, mygene for cross-database aliases, Ensembl REST for low-level access), with version pinning for reproducibility and explicit handling of one-to-many mappings, withdrawn symbols, and species-specific naming.

## The Single Most Important Modern Insight -- Excel autocorrect renamed real genes

Ziemann, Eren, El-Osta 2016 *Genome Biol* 17:177 scanned 18 leading genomics journals and found ~20% of papers with Excel-attached supplementary gene lists had silently mangled symbols (`SEPT2` -> `2-Sep`, `MARCH1` -> `1-Mar`, ...). Five years later the problem persisted. HGNC's response (Bruford, Braschi, Denny, Jones, Seal, Tweedie 2020 *Nat Genet* 52:754) was to rename the affected genes:

| Old | New | Affected |
|-----|-----|----------|
| `SEPT#` | `SEPTIN#` | SEPT1 - SEPT14 -> SEPTIN1 - SEPTIN14 |
| `MARCH#` | `MARCHF#` | MARCH1 - MARCH11 -> MARCHF1 - MARCHF11 |
| `MARC#` | `MTARC#` | MARC1, MARC2 -> MTARC1, MTARC2 |
| `DEC1` | `DELEC1` | DEC1 -> DELEC1 |

Code that hard-codes old symbols silently drops these genes when joined against post-2020 annotations. Detection on import: if a gene column contains `^\d{1,2}-(Jan|Feb|Mar|...|Dec)$` patterns, the file was Excel-corrupted. Always `read.csv(colClasses=c(gene='character'))` (R) or `pd.read_csv(dtype={'gene': str})` (Python) -- but the damage is at Excel-save time, not import time.

Two related insights that determine half the practical work:

1. **Ensembl version suffixes matter sometimes and not others.** `ENSG00000123456.7` is release-specific; the unversioned `ENSG00000123456` is the stable cross-release ID. STRIP for cross-release joins, MSigDB lookups, gene-set databases. KEEP for intra-release reproducibility and clinical reports. CRITICAL: the naive `sub('\\..*', '', x)` regex ALSO strips the GENCODE `_PAR_Y` suffix in releases 25-43, collapsing chrY PAR duplicates onto their chrX counterparts. Use `sub('\\.[0-9]+(_PAR_Y)?$', '\\1', x)`.

2. **Never use HGNC symbols as the primary computational key.** Symbols change. Use Ensembl or Entrez as keys; carry symbols only as display labels in the final results table.

## Algorithmic Taxonomy

| Tool | Source | Speed | Strength | Use for |
|------|--------|-------|----------|---------|
| AnnotationDbi + `org.Hs.eg.db` / `org.Mm.eg.db` | NCBI Gene snapshot, pinned at Bioc install | Fast, local | Stable, version-pinned | Default for Ensembl <-> Entrez <-> Symbol within Bioconductor |
| biomaRt | Ensembl BioMart over HTTP | Slow for >5k queries; timeouts | Ensembl-specific attributes (biotype, transcript versions, paralogs, orthologs) | Need Ensembl-specific fields; archive endpoints for release pinning |
| mygene.info / mygene (Python) | REST API to a curated meta-database | Server-side batching of 1000 IDs | Best for symbol/alias/prev_symbol resolution | Cross-database; HGNC withdrawn symbol resolution; non-R environments |
| Ensembl REST | Direct REST API to Ensembl | Rate-limited (15 req/sec) | Low-level access to variant consequence, sequence, etc. | Specialized queries not covered by biomaRt |
| pyensembl | Local Ensembl database (Python) | Fast, local, version-pinned | Reproducible offline; gene objects with transcript and exon access | Python pipelines needing rich annotation |
| HGNC API direct | https://rest.genenames.org | REST | Authoritative source for HGNC | Symbol provenance, prev/alias detection |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| R Bioconductor pipeline, Ensembl <-> Entrez <-> Symbol | `AnnotationDbi::mapIds(org.Hs.eg.db, ...)` | Fastest, version-pinned, stable |
| Need Ensembl-only attributes (biotype, paralog, ortholog) | `biomaRt::useEnsembl(version=N)` | Only biomaRt exposes these |
| Cross-database with alias and withdrawn-symbol fallback | mygene `querymany(scopes='symbol,alias,prev_symbol')` | Designed for this case |
| Python pipeline, reproducible | `pyensembl` with pinned release | Offline, version-locked |
| Clinical report needing canonical transcript per gene | MANE Select (Morales 2022 *Nature* 604:310) | Cross-database consensus (RefSeq + Ensembl) |
| Cross-species mouse <-> human | Ensembl Compara `getLDS` filtered to `one2one` | Compara has best coverage; one2one most defensible |
| Building tx2gene for tximport | `GenomicFeatures::makeTxDbFromGFF` on the SAME GTF used in quantification | Annotation pinning matters |
| Need to reproduce a 2023 analysis exactly | `useEnsembl(version=109)` (or whichever release was used) | Without `version=`, biomaRt floats to current release |
| GRCh37 (legacy clinical) | `useEnsembl(GRCh=37)` dedicated permanent endpoint | GRCh37 -> GRCh38 mappings are not 1:1 |

## AnnotationDbi + org.db

**Goal:** Map Ensembl gene IDs to symbols, Entrez IDs, or descriptions using a local Bioconductor annotation package.

**Approach:** `mapIds()` with the source keytype and target column; handle one-to-many via `multiVals`.

```r
library(org.Hs.eg.db)
library(AnnotationDbi)

ensembl_ids <- sub('\\.[0-9]+(_PAR_Y)?$', '\\1', rownames(counts))

symbols  <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                    keytype = 'ENSEMBL', column = 'SYMBOL',
                    multiVals = 'first')
entrez   <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                    keytype = 'ENSEMBL', column = 'ENTREZID',
                    multiVals = 'first')
descrips <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                    keytype = 'ENSEMBL', column = 'GENENAME',
                    multiVals = 'first')

keytypes(org.Hs.eg.db)
```

`multiVals` options: `'first'` (silent), `'asNA'` (NA for ambiguous), `'list'` (preserve all). For DE results tables, `'first'` is typical but the mapping rate should be reported.

For mouse: `org.Mm.eg.db`. For other organisms: check Bioconductor `AnnotationData -> OrgDb` list.

## biomaRt with Version Pinning

**Goal:** Query Ensembl BioMart with the EXACT release version, for reproducibility.

**Approach:** `useEnsembl(version=N)` pins; `listEnsemblArchives()` lists available archives.

```r
library(biomaRt)

ensembl <- useEnsembl(biomart = 'genes',
                       dataset = 'hsapiens_gene_ensembl',
                       version = 110)

ensembl_grch37 <- useEnsembl(biomart = 'genes',
                              dataset = 'hsapiens_gene_ensembl',
                              GRCh = 37)

mapping <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id',
                   'gene_biotype', 'description'),
    filters    = 'ensembl_gene_id',
    values     = ensembl_ids,
    mart       = ensembl
)
```

The `filters=` argument is PLURAL. The singular `filter=` may work via R's partial matching but breaks unpredictably if another argument starts with `f`. Always spell `filters=` and `values=` fully.

Multiple filters:

```r
genes_in_region <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters    = c('chromosome_name', 'start', 'end'),
    values     = list('16', 1100000, 1250000),
    mart       = ensembl
)
```

Without `version=`, biomaRt floats to the current release -- a script written in 2023 against Ensembl 109 produces different mappings in 2026 against Ensembl 113. ALWAYS pin for any published analysis. Cache the mapping table alongside the analysis for reproducibility.

`listEnsemblArchives()` shows the available historical releases.

## mygene (Python)

**Goal:** Map between any identifier systems using the curated MyGene.info meta-database with alias fallback.

**Approach:** `MyGeneInfo().querymany(ids, scopes, fields, species)`; auto-batches at 1000 IDs server-side.

```python
import mygene

mg = mygene.MyGeneInfo()

results = mg.querymany(['ENSG00000141510', 'ENSG00000012048', 'ENSG00000141736'],
                        scopes='ensembl.gene', fields='symbol,entrezgene,uniprot',
                        species='human')

mapping = {r['query']: r.get('symbol', None) for r in results}

results = mg.querymany(['SEPT1', 'MARCH1', 'OCT4'],
                        scopes='symbol,alias,prev_symbol',
                        fields='symbol,entrezgene,ensembl.gene',
                        species='human')
```

For paper-derived gene lists where symbols may be old or aliases (OCT4 vs POU5F1, MARCH1 vs MARCHF1, SEPT2 vs SEPTIN2), `scopes='symbol,alias,prev_symbol'` handles the resolution. The MyGene database aggregates HGNC's prev/alias columns.

OCT4 is the common usage; POU5F1 is the official HGNC symbol; in MSigDB the gene is POU5F1; in a Western blot legend it's "Oct4". For mapping a stem-cell paper to an Ensembl-quantified matrix, scope to aliases.

## pyensembl

```python
from pyensembl import EnsemblRelease

ensembl = EnsemblRelease(110, species='human')

gene = ensembl.gene_by_id('ENSG00000141510')
gene.gene_name

gene = ensembl.genes_by_name('TP53')[0]
gene.gene_id

mapping = {}
for eid in ensembl_ids:
    try:
        gene = ensembl.gene_by_id(eid.split('.')[0])
        mapping[eid] = gene.gene_name
    except ValueError:
        mapping[eid] = None
```

pyensembl downloads and caches the release database on first use; thereafter offline and version-locked.

## Apply Mapping to Count Matrix -- Handling One-to-Many

**Goal:** Convert the gene index of a count matrix to a different ID type, summing reads from multiple source IDs that map to the same target.

**Approach:** Look up mapping, replace index, aggregate duplicates by SUM (not mean -- counts add).

```python
import pandas as pd
import mygene

def map_count_matrix_ids(counts, from_type='ensembl.gene', to_type='symbol',
                         species='human'):
    '''Map gene IDs in count matrix index, summing reads when multiple source map to one target.'''
    mg = mygene.MyGeneInfo()
    clean = [g.split('.')[0] for g in counts.index]
    results = mg.querymany(clean, scopes=from_type, fields=to_type, species=species)
    mapping = {r['query']: r[to_type] for r in results if to_type in r}
    new_index = [mapping.get(g.split('.')[0], g) for g in counts.index]
    counts_mapped = counts.copy()
    counts_mapped.index = new_index
    counts_mapped = counts_mapped.groupby(counts_mapped.index).sum()
    return counts_mapped

mapped = map_count_matrix_ids(counts, 'ensembl.gene', 'symbol')
```

Counts ADD when collapsing multiple source genes to one target. Means or medians would be wrong (they understate library size for the merged target).

```r
library(biomaRt)

ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 110)

clean <- sub('\\.[0-9]+(_PAR_Y)?$', '\\1', rownames(counts))

mapping <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters    = 'ensembl_gene_id',
    values     = clean,
    mart       = ensembl
)

counts_df <- as.data.frame(counts)
counts_df$ensembl <- clean
merged <- merge(counts_df, mapping, by.x = 'ensembl', by.y = 'ensembl_gene_id')
counts_by_symbol <- aggregate(. ~ hgnc_symbol,
                              data = merged[, setdiff(colnames(merged), 'ensembl')],
                              FUN = sum)
rownames(counts_by_symbol) <- counts_by_symbol$hgnc_symbol
counts_by_symbol$hgnc_symbol <- NULL
```

## Handle Unmapped IDs

```python
def robust_id_mapping(gene_ids, from_type, to_type, species='human'):
    import mygene
    mg = mygene.MyGeneInfo()
    clean = [g.split('.')[0] for g in gene_ids]
    results = mg.querymany(clean, scopes=from_type, fields=to_type, species=species)
    mapping, unmapped = {}, []
    for r in results:
        original = gene_ids[clean.index(r['query'])]
        if to_type in r:
            mapping[original] = r[to_type]
        else:
            mapping[original] = original
            unmapped.append(original)
    print(f'Mapped: {len(gene_ids) - len(unmapped)}/{len(gene_ids)}')
    return mapping, unmapped
```

Unmapped fraction is a QC signal:
- <5% unmapped: normal (rare genes, recent deprecations)
- 5-20% unmapped: check Ensembl release alignment; check for HGNC renames
- >20% unmapped: wrong annotation release, wrong species, or wrong source ID type

## MANE Select for Clinical Reporting

**Goal:** Use the single representative transcript per gene with identical exon/CDS in RefSeq AND Ensembl for clinical variant reporting.

**Approach:** Download the MANE TSV; join on `Ensembl_Gene` -> `Ensembl_nuc` (transcript) and `RefSeq_nuc`.

Morales J, Pujar S, Loveland JE et al. 2022 *Nature* 604:310-315 established MANE Select. ~19,000+ protein-coding genes have a single agreed transcript with matched coordinates across RefSeq (NM_xxxxxx) and Ensembl/GENCODE (ENST00000xxxxxxx). MANE Plus Clinical adds extra transcripts at loci where Select misses clinical variants.

For clinical reports with HGVS notation like `NM_000546.6:c.215C>G`, use the MANE Select RefSeq accession. The MANE TSV (downloadable from NCBI) provides the Ensembl crosswalk.

## Cross-Species Orthologs

**Goal:** Map mouse <-> human (or any pair) for cross-species integration or pathway transfer.

**Approach:** Ensembl Compara via biomaRt `getLDS`; filter to orthology type appropriate to use.

```r
library(biomaRt)

human <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 110)
mouse <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl', version = 110)

orthologs <- getLDS(
    attributes  = c('hgnc_symbol', 'ensembl_gene_id'),
    filters     = 'ensembl_gene_id',
    values      = human_gene_ids,
    mart        = human,
    attributesL = c('mgi_symbol', 'ensembl_gene_id', 'mmusculus_homolog_orthology_type'),
    martL       = mouse
)
```

| Strategy | When | Trade-off |
|----------|------|-----------|
| `one2one` orthologs only | Cross-species scRNA-seq integration; conservative DE comparison | Loses genes with paralog expansions; lower coverage |
| Include `one2many` | Broader gene coverage needed | Must select within group (highest confidence; highest expression) |
| Include `many2many` | Maximum inclusivity | Introduces ambiguity; use with caution |

The "homology threshold" problem: no automatic threshold reliably separates true orthologs from paralogs across all gene families. For pathway transfer (mouse signature -> human), filter to one2one and accept the coverage loss.

Alternative sources: OMA (Hierarchical Orthologous Groups, cleaner one2one when present, smaller coverage); OrthoDB (hierarchical at multiple taxonomic levels). OrthoFinder for custom genomes.

## PAR Gene Complications

Pseudo-autosomal region (PAR) genes exist on both X and Y with identical sequences. In GENCODE 25-43, the chrY copy has a `_PAR_Y` suffix. In GENCODE 44+ (Ensembl 110+), chrY PAR genes get their own ENSG accessions.

```python
par_genes_human = ['SHOX', 'IL3RA', 'SLC25A6', 'P2RY8', 'AKAP17A', 'ASMT', 'DHRSX']
dup_ids = counts.index[counts.index.duplicated()].unique()
if len(dup_ids) > 0:
    print(f'Duplicate gene entries: {len(dup_ids)}')
    counts = counts.groupby(counts.index).sum()
```

Reads from PAR regions cannot be unambiguously assigned to X or Y. Some references mask the Y-chromosome PAR to avoid double-counting; verify what the alignment reference does before building the matrix.

## Build tx2gene for tximport

**Goal:** Create the transcript-to-gene mapping needed by tximport for gene-level summarization.

**Approach:** Build from the SAME GTF used to construct the Salmon/kallisto index, OR pull from biomaRt with version pinning.

```r
library(GenomicFeatures)

txdb <- makeTxDbFromGFF('annotation.gtf.gz')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- AnnotationDbi::select(txdb, k, 'GENEID', 'TXNAME')
```

```r
library(biomaRt)

mart <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 110)
tx2gene <- getBM(
    attributes = c('ensembl_transcript_id_version', 'ensembl_gene_id_version'),
    mart       = mart
)
colnames(tx2gene) <- c('TXNAME', 'GENEID')
```

```python
import pandas as pd

def tx2gene_from_gtf(gtf_path):
    records = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#') or '\ttranscript\t' not in line:
                continue
            attrs = line.strip().split('\t')[8]
            gene_id = [a.split('"')[1] for a in attrs.split(';') if 'gene_id' in a][0]
            tx_id   = [a.split('"')[1] for a in attrs.split(';') if 'transcript_id' in a][0]
            records.append({'TXNAME': tx_id, 'GENEID': gene_id})
    return pd.DataFrame(records).drop_duplicates()
```

CRITICAL: the tx2gene MUST use the same versioning convention as the Salmon/kallisto index. If the index used `ENST00000269305.9` and tx2gene has `ENST00000269305` (unversioned), tximport drops the transcripts. Mismatched versions silently lose data.

## ID Type Reference

| Type | Example | Stability | Use case |
|------|---------|-----------|----------|
| Ensembl Gene | ENSG00000141510 | Stable across releases; versioned | RNA-seq, GTFs, primary computational key |
| Ensembl Transcript | ENST00000269305 | Stable; versioned | Transcript-level analysis |
| Entrez Gene | 7157 | Stable; never reused | NCBI databases, KEGG pathways |
| HGNC Symbol | TP53 | Changes (see SEPT/MARCH renames) | Display labels only |
| UniProt | P04637 | Stable; versioned releases | Protein databases |
| RefSeq mRNA | NM_000546 | Stable; versioned | Clinical reports, HGVS notation |
| MANE Select | NM_000546.6 / ENST00000269305.9 | Stable consensus | Clinical variant reporting |

## Per-Method Failure Modes

### `_PAR_Y` stripped, chrY duplicates collapsed

**Trigger:** GENCODE v40 count matrix; `rownames(counts) <- sub('\\..*', '', rownames(counts))`; duplicate row indices and inflated chrY PAR gene counts.

**Mechanism:** Default regex strips `_PAR_Y` along with the version suffix. Two distinct rows (chrX and chrY copies) become the same ENSG ID; `aggregate` sums them.

**Symptom:** Counts for PAR genes double; sex check shows females expressing chrY genes; downstream `rowGroupBy` returns warnings.

**Fix:** Use the preserving regex: `sub('\\.[0-9]+(_PAR_Y)?$', '\\1', x)`. Or upgrade quantification to GENCODE 44+ where `_PAR_Y` is retired.

### `biomaRt` returned 0 rows without warning

**Trigger:** `getBM(attributes=..., filter='ensembl_gene_id', values=ids, mart=mart)` -- note singular `filter`.

**Mechanism:** R's partial matching usually resolves `filter` -> `filters`, but in some package versions or with conflicting argument names, the call silently passes nothing.

**Symptom:** Empty result data frame; no error.

**Fix:** Always spell `filters=` and `values=` fully.

### HGNC SEPT/MARCH symbols silently dropped

**Trigger:** Code copies a pre-2020 list of septin genes (`SEPT1`, `SEPT2`, ...); current org.db / biomaRt returns no matches.

**Mechanism:** HGNC renamed all `SEPT#` to `SEPTIN#` in 2020.

**Symptom:** 0% mapping rate for septin genes; functional analyses missing septin pathways.

**Fix:** Use mygene `scopes='symbol,alias,prev_symbol'`; or update the input list to current symbols.

### biomaRt drift between runs

**Trigger:** A 2023 analysis used `useEnsembl()` without `version=`; rerun in 2026 produces 200 fewer significant genes.

**Mechanism:** Without `version=`, biomaRt floats to the current release. Symbols, biotypes, and gene boundaries change between releases.

**Symptom:** Non-reproducible results across runs of the same script.

**Fix:** Pin `useEnsembl(version=N)` where N is the release used in the original analysis. Cache the mapping table.

### tx2gene version mismatch with Salmon index

**Trigger:** `tximport(files, type='salmon', tx2gene)` runs but the gene-level counts have far fewer genes than expected.

**Mechanism:** Salmon index built with versioned transcript IDs (`ENST00000269305.9`) but tx2gene has unversioned IDs (`ENST00000269305`). Transcripts silently drop during the mapping step.

**Symptom:** Lower-than-expected gene count; warning from tximport about missing transcript IDs.

**Fix:** Match versioning convention: rebuild tx2gene with the same versioning as the index. `GenomicFeatures::makeTxDbFromGFF` on the same GTF as the index is the safest path.

### Cross-species mapping reports many2many, user picks one arbitrarily

**Trigger:** Mouse-to-human mapping returns 1.3 mouse genes per human gene on average; user takes the first row of each duplicate.

**Mechanism:** Many2many orthology is genuinely ambiguous; "first row" is unprincipled and irreproducible across biomaRt API versions.

**Symptom:** Different mappings on rerun; conflicting downstream gene sets.

**Fix:** Either filter to `mmusculus_homolog_orthology_type == 'ortholog_one2one'` (conservative) or aggregate via highest homology confidence score (`mmusculus_homolog_perc_id_r1`).

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `filters` returns empty | Singular `filter=` partial-matched against another argument | Spell `filters=` fully |
| `1-Mar` in gene column | Excel autocorrected `MARCH1` | Re-import with explicit string type; map back to MARCHF1 |
| pyensembl `ValueError: gene not found` | ID not in pinned release; or unversioned ID against versioned database | Strip version before lookup; verify release |
| Duplicate rownames after aggregate | Collapsed multiple source IDs to one target; OR `_PAR_Y` stripped | Sum-collapse expected; for PAR_Y use preserving regex |
| biomaRt timeout for >5k IDs | Query too large | Chunk into batches of 1000 |
| Wrong species mapping | Default `species='human'` in mygene; mouse query returns nothing | Pass `species='mouse'` explicitly |
| `ENSEMBL` keytype not available | Older org.db package or non-human/mouse | `keytypes(orgdb)` to verify |

## References

- Ziemann M, Eren Y, El-Osta A. 2016. Gene name errors are widespread in the scientific literature. *Genome Biol* 17:177. doi:10.1186/s13059-016-1044-7
- Bruford EA, Braschi B, Denny P, Jones TEM, Seal RL, Tweedie S. 2020. Guidelines for human gene nomenclature. *Nat Genet* 52:754-758. doi:10.1038/s41588-020-0669-3
- Morales J, Pujar S, Loveland JE, et al. 2022. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. *Nature* 604:310-315. doi:10.1038/s41586-022-04558-8
- Durinck S, Spellman PT, Birney E, Huber W. 2009. Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. *Nat Protoc* 4(8):1184-1191. doi:10.1038/nprot.2009.97
- Carlson M, Falcon S, Pages H, Li N. 2019. org.Hs.eg.db: Genome wide annotation for Human. R package. Bioconductor.
- Wu C et al. 2013. BioGPS and MyGene.info: organizing online, gene-centric information. *Nucleic Acids Res* 41(D1):D561-D565. doi:10.1093/nar/gks1114
- Frankish A et al. 2021. GENCODE 2021. *Nucleic Acids Res* 49(D1):D916-D923. doi:10.1093/nar/gkaa1087
- Howe KL et al. 2021. Ensembl 2021. *Nucleic Acids Res* 49(D1):D884-D891. doi:10.1093/nar/gkaa942
- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Res* 4:1521. doi:10.12688/f1000research.7563.2

## Related Skills

- counts-ingest - Building count matrices and tx2gene
- metadata-joins - Joining annotation with sample tables
- normalization - Biotype filtering before normalization
- sparse-handling - Single-cell row metadata in AnnData
- differential-expression/de-results - Annotating DE results; gene-symbol display
- rna-quantification/tximport-workflow - Detailed tximport + tx2gene workflow
- pathway-analysis/go-enrichment - Entrez IDs required
- pathway-analysis/kegg-pathways - Entrez IDs; strain-specific organism codes
- database-access/biomart-queries - General biomaRt patterns
- database-access/uniprot-access - UniProt mapping details
- database-access/ortholog-inference - De novo ortholog inference for custom genomes
