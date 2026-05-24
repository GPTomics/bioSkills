---
name: bio-copy-number-cnv-annotation
description: Annotate copy number variant segments with overlapping genes, dosage-sensitivity scores, cancer driver databases, population frequencies, and clinical-variant content. Covers bedtools/pybedtools interval intersection, AnnotSV comprehensive annotation and ranking, ClinGen haploinsufficiency/triplosensitivity scoring, gnomAD-SV/DGV frequency filtering, COSMIC Cancer Gene Census, and ClinVar overlap. Use when interpreting which genes a CNV affects, distinguishing the driver gene of a focal event from passengers, filtering against population CNVs, separating whole-gene from partial-gene overlap, or preparing CNVs for clinical classification.
tool_type: mixed
primary_tool: bedtools
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, AnnotSV 3.4+, Python 3.10+ with
pybedtools 0.9+, pandas 2.2+, pysam 0.22+; R 4.3+ with clusterProfiler 4.10+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `bedtools --version`, `AnnotSV --version`
- Python: `pip show pybedtools pandas pysam`
- R: `packageVersion('clusterProfiler')`

If code throws an error, introspect the installed package and adapt the example. AnnotSV
output column names change between major versions — verify against the installed version.

# CNV Annotation

**"Annotate my CNV calls with the genes they affect"** -> Overlap CNV segments with gene models, dosage-sensitivity maps, and clinical databases. The hard part is not the intersection — it is deciding *which* genes matter. A focal amplification overlapping 30 genes usually has one driver (the peak gene); a deletion's consequence depends on whether each gene is dosage-sensitive and whether the whole gene or only part is removed.

- CLI: `bedtools intersect -a cnvs.bed -b genes.bed -wa -wb`; `AnnotSV` for full annotation
- Python: `pybedtools` for interval logic; `pysam` for VCF database queries

## Annotation Strategy — Pick the Database for the Question

| Question | Resource | What it answers |
|----------|----------|-----------------|
| Which genes does this CNV span? | RefSeq/GENCODE gene BED | Raw overlap (not yet consequence) |
| Is loss of this gene damaging? | ClinGen haploinsufficiency (HI) score | Dosage sensitivity to deletion |
| Is gain of this gene damaging? | ClinGen triplosensitivity (TS) score | Dosage sensitivity to duplication |
| Is this CNV common in the population? | gnomAD-SV, DGV, 1000G CNV | Benign-frequency filtering |
| Is this a known recurrent disorder locus? | ClinGen dosage regions, DECIPHER | Genomic-disorder context |
| Is this a cancer driver? | COSMIC Cancer Gene Census, OncoKB | Oncogene vs tumor-suppressor role |
| Is there pathogenic small-variant content? | ClinVar | Coincident SNV/indel pathogenicity |
| One-shot comprehensive annotation + ranking | AnnotSV | Aggregates most of the above |

For constitutional CNV *classification* (assigning pathogenic/VUS/benign), the annotated
output feeds the ACMG/ClinGen points framework — see germline-cnv-interpretation. For
cohort-level recurrence and driver-peak identification, see recurrent-cnv.

## The Core Distinction: Overlap Is Not Consequence

A CNV overlapping a gene does not necessarily change that gene's dosage in a way that
matters. Three refinements separate annotation from interpretation:

1. **Whole-gene vs partial overlap.** A deletion spanning an entire gene removes one copy
   (clean haploinsufficiency test). A deletion removing only the last two exons creates a
   truncated allele — a different, often more damaging, consequence. Always record the
   fraction of each gene covered and whether coding exons or only introns/UTRs are hit.
2. **Dosage sensitivity.** Most genes tolerate single-copy loss. ClinGen HI/TS scores
   (3 = sufficient evidence for dosage sensitivity, 0 = no evidence, 30 = gene associated
   with an autosomal-recessive phenotype, 40 = dosage sensitivity unlikely) indicate which
   genes' loss/gain is actually consequential.
3. **Driver vs passenger in focal events.** A focal amplification carries many genes; the
   driver is the one under selection, typically at the recurrence peak across a cohort
   (GISTIC) and a known oncogene. Annotating all 30 genes as "amplified" overstates.

## Gene Overlap with bedtools

**Goal:** Find genes overlapping each CNV segment, recording overlap extent.

**Approach:** Convert segments to BED, intersect with a gene model, keep both feature
sets (`-wo` reports the overlap length) so partial vs whole-gene overlap is recoverable.

```bash
# Segments to BED (CNVkit .cns example; columns chrom/start/end/log2)
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' sample.cns > sample.cnv.bed

# Intersect; -wo appends the number of overlapping bases
bedtools intersect -a sample.cnv.bed -b gencode.genes.bed -wo > cnv_gene_overlap.txt
```

## Comprehensive Annotation with AnnotSV

**Goal:** Annotate CNVs against genes, dosage maps, population frequency, and clinical
databases in one pass, with a built-in pathogenicity ranking.

**Approach:** Export CNVs to VCF or BED and run AnnotSV; it returns a "full" line per
gene plus a "split" summary, with an ACMG-aligned rank (1-5).

```bash
AnnotSV \
    -SVinputFile sample.cnv.vcf \
    -genomeBuild GRCh38 \
    -annotationMode both \
    -outputFile sample_annotated.tsv

# Output includes: overlapped genes, ClinGen HI/TS, gnomAD-SV/DGV frequency, OMIM,
# ClinVar, DECIPHER, and an ACMG-class rank per SV.
```

AnnotSV's rank is a useful triage signal, not a final classification — confirm against
the ClinGen points framework for clinical reporting.

## Dosage-Sensitivity and Driver Annotation

**Goal:** Tag each affected gene with its dosage sensitivity and, for tumors, its driver
role, so passengers can be separated from drivers.

**Approach:** Join the gene-overlap table to the ClinGen dosage map (HI/TS scores) and to
the COSMIC Cancer Gene Census; flag CNVs whose direction matches a known mechanism
(oncogene amplified, tumor suppressor deleted).

```python
import pandas as pd

def annotate_dosage_and_drivers(overlap_tsv, clingen_dosage, cgc_file):
    '''Tag overlapped genes with ClinGen HI/TS and COSMIC driver role.'''
    cols = ['cnv_chrom', 'cnv_start', 'cnv_end', 'log2',
            'gene_chrom', 'gene_start', 'gene_end', 'gene', 'overlap_bp']
    df = pd.read_csv(overlap_tsv, sep='\t', names=cols)
    df['gene_len'] = df['gene_end'] - df['gene_start']
    df['gene_frac_covered'] = (df['overlap_bp'] / df['gene_len']).clip(upper=1.0)
    df['whole_gene'] = df['gene_frac_covered'] >= 0.99

    dosage = pd.read_csv(clingen_dosage, sep='\t')  # gene, HI_score, TS_score
    df = df.merge(dosage, on='gene', how='left')

    cgc = pd.read_csv(cgc_file, sep='\t')
    role = dict(zip(cgc['Gene Symbol'], cgc['Role in Cancer']))
    df['driver_role'] = df['gene'].map(role)

    # Direction-consistent driver hits: oncogene amplified or TSG deleted.
    df['driver_hit'] = (
        ((df['log2'] > 0.3) & df['driver_role'].fillna('').str.contains('oncogene')) |
        ((df['log2'] < -0.3) & df['driver_role'].fillna('').str.contains('TSG')))
    return df
```

## Population-Frequency Filtering

**Goal:** Remove common, presumed-benign CNVs before clinical interpretation.

**Approach:** Reciprocal-overlap match each CNV against a population SV catalog
(gnomAD-SV, DGV); a CNV with high reciprocal overlap to a common population CNV of the
same type is likely benign.

```bash
# 50% reciprocal overlap (-f 0.5 -r): same-type, similar-extent population match.
# Reciprocal overlap, not one-sided, prevents a tiny CNV inside a huge population CNV
# (or vice versa) from being wrongly matched.
bedtools intersect -a sample.cnv.bed -b gnomad_sv.bed -f 0.5 -r -wa -wb \
    > cnv_population_match.txt
```

## Pathway Enrichment of Affected Genes

**Goal:** Test whether genes in amplified (or deleted) regions are enriched for pathways.

**Approach:** Extract genes by CNV direction, map to Entrez IDs, run GO/KEGG enrichment.
Caveat: CNVs are large and gene-dense, so enrichment is biased toward whatever pathways
cluster in CNV-prone genomic regions — interpret as hypothesis-generating.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

amp_genes <- unique(cnv_annot$gene[cnv_annot$log2 > 0.3])
entrez <- na.omit(mapIds(org.Hs.eg.db, keys = amp_genes,
                         keytype = 'SYMBOL', column = 'ENTREZID'))
go_bp <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = 'BP',
                  pAdjustMethod = 'BH', qvalueCutoff = 0.05)
```

## Failure Modes

### Genome-build mismatch between CNVs and annotation

**Trigger:** CNV coordinates on GRCh37 intersected with a GRCh38 gene model (or vice versa).

**Mechanism:** Coordinates silently shift; the intersection succeeds and returns wrong genes.

**Symptom:** Implausible gene assignments; a known driver locus annotated with the wrong gene; systematic offset.

**Fix:** Confirm both inputs are the same build. If not, liftOver the CNVs (note that liftOver can split or drop segments across assembly gaps) and verify a known landmark.

### Annotating all overlapped genes as the "affected" genes

**Trigger:** Reporting every gene a focal amplification spans as amplified/driver.

**Mechanism:** Focal events are megabases wide and gene-dense; only the selected gene is the driver.

**Symptom:** A 2 Mb amplicon "amplifies" 40 genes; the report cannot distinguish ERBB2 from its passengers.

**Fix:** For focal events, prioritize the gene at the cohort recurrence peak (GISTIC, see recurrent-cnv) and known drivers (CGC/OncoKB). Report passengers separately or not at all.

### ClinVar CLNSIG parsing errors

**Trigger:** Naive string matching on the ClinVar `CLNSIG` INFO field.

**Mechanism:** `CLNSIG` is multi-valued, mixes terms ("Conflicting_classifications", "Pathogenic/Likely_pathogenic", "Benign/Likely_benign"), and is per-small-variant — not per-CNV. A substring match for "pathogenic" silently captures "Likely_pathogenic" (intended) but a careless match also fires on records that are conflicting or benign once underscores and slashes are involved.

**Symptom:** Benign or conflicting variants reported as pathogenic; CNV flagged on incidental nearby SNVs.

**Fix:** Parse `CLNSIG` against the controlled vocabulary; exclude "Conflicting" and benign terms explicitly. Remember ClinVar SNV/indel pathogenicity does not transfer to a CNV — use it as context, and use ClinVar's own CNV records or ClinGen dosage regions for the CNV itself.

### Equating overlap with consequence

**Trigger:** Treating any gene-overlapping CNV as functionally significant.

**Mechanism:** Most single-copy losses are tolerated; partial overlaps may hit only introns/UTRs.

**Symptom:** Long lists of "affected" dosage-insensitive genes; benign CNVs over-called as significant.

**Fix:** Require dosage evidence (ClinGen HI/TS) and record coding-exon overlap and whole-gene-vs-partial status before calling a gene affected.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Population-CNV reciprocal overlap | >= 50% (`-f 0.5 -r`) | Standard reciprocal-overlap match for benign filtering (convention traces to gnomAD-SV / DGV workflows; Collins SR et al 2020 *Nature* 581:444 uses comparable reciprocal-overlap thresholds for benign-population matching) |
| Common-CNV benign frequency | > 1% population frequency | ACMG/ClinGen: high frequency supports benign |
| ClinGen HI/TS dosage-sensitive | score = 3 | ClinGen: sufficient evidence for dosage sensitivity |
| Whole-gene overlap | >= 99% gene length covered | Distinguishes clean haploinsufficiency from partial/truncating |
| AnnotSV pathogenic ranks | rank 4-5 | AnnotSV ACMG-aligned ranking (1 benign - 5 pathogenic) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Wrong genes assigned to a CNV | hg19/hg38 build mismatch | Match builds; liftOver and verify a landmark |
| 40 genes called "amplified" for one amplicon | All overlapped genes reported as drivers | Prioritize recurrence-peak + known-driver genes |
| Benign variants flagged pathogenic | Substring match on ClinVar CLNSIG | Parse the controlled vocabulary explicitly |
| Tiny CNV matched to a huge population CNV | One-sided overlap used for frequency filter | Use reciprocal overlap (`-f 0.5 -r`) |
| AnnotSV columns not found | Column names differ across AnnotSV versions | Check `head` of the output for the installed version |
| Enrichment dominated by gene-dense loci | CNVs span gene clusters | Treat CNV-gene enrichment as hypothesis-generating |

## References

- Geoffroy V et al 2018. AnnotSV: an integrated tool for structural variations annotation. Bioinformatics 34:3572
- Riggs ER et al 2020. Technical standards for the interpretation and reporting of constitutional copy-number variants: ACMG and ClinGen. Genet Med 22:245
- Collins RL et al 2020. A structural variation reference for medical and population genetics (gnomAD-SV). Nature 581:444
- Sondka Z et al 2018. The COSMIC Cancer Gene Census. Nat Rev Cancer 18:696

## Related Skills

- copy-number/germline-cnv-interpretation - ACMG/ClinGen points-based CNV classification
- copy-number/recurrent-cnv - GISTIC2 recurrence peaks and driver-gene identification
- copy-number/cnvkit-analysis - Generates the CNV segments to annotate
- copy-number/cnv-visualization - Visualizing annotated CNVs
- pathway-analysis/go-enrichment - GO/KEGG enrichment methodology and caveats
- genome-intervals/bed-file-basics - BED interval operations
- clinical-databases/clinvar-lookup - Querying ClinVar for variant pathogenicity
