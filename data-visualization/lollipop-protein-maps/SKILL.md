---
name: bio-data-visualization-lollipop-protein-maps
description: Plot per-gene mutation distributions on a protein-domain map (lollipop / needle plots) showing mutation position, recurrence count, and variant classification with maftools, g3-lollipop, trackViewer, and ProteinPaint. Use when visualizing recurrent mutation hotspots on a single gene's protein, marking domain boundaries from UniProt/Pfam, comparing missense vs truncating distributions, or contrasting two cohorts on the same lollipop.
tool_type: mixed
primary_tool: maftools
---

## Version Compatibility

Reference examples tested with: maftools 2.18+, trackViewer 1.38+, g3-lollipop (JavaScript via R `g3viz` 1.2+), Bio.PDB 1.83+ (for domain coordinates). ProteinPaint is a hosted service.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Lollipop / Needle Protein Maps

**"Plot mutations on a gene's protein"** -> Render a horizontal protein backbone with colored domain rectangles (from UniProt/Pfam/InterPro), then stack vertical lines ("stems") at mutated amino-acid positions, capped with circles ("lollipops") whose size reflects mutation count and whose color encodes variant class. The biological story is hotspot identification — a tall stack of recurrences at a single residue (e.g., KRAS G12, PIK3CA E545/H1047) is the visual signature of a driver mutation.

- R: `maftools::lollipopPlot`, `trackViewer::lolliplot`, `g3viz::g3Lollipop`
- Python: `pyLollipop` (limited maintenance); ProteinPaint via API
- Web: cBioPortal, ProteinPaint, MutationMapper

## The Single Most Important Modern Insight -- Hotspot Recurrence Drives the Plot

A lollipop plot exists to identify hotspots — residues with disproportionate recurrence. The MutSig hotspot test (Lawrence 2014 *Nature* 505:495) and statisticalhotspot methods (Chang 2016 *Nat Biotechnol* 34:155) formalize this: a residue's mutation count should exceed the gene-wide background rate × residue count. Visualizing this on a domain map IS the diagnostic.

Key practical consequences:
- **Stack height ≠ frequency**: a tall lollipop at residue 600 means recurrence, not population frequency. Annotate the count.
- **Domain colors should encode functional class** (kinase, SH2, binding), not random hue.
- **Mark known activating/inactivating residues** (G12 for KRAS, R175 for TP53) with bold labels.

## Decision Tree by Question

| Question | Approach |
|----------|----------|
| Where are the hotspots? | Lollipop with size = count; label top 5 recurrent residues |
| Missense vs truncating distribution? | Color stems by class; tumor suppressors show truncating spread; oncogenes show missense hotspots |
| Compare two cohorts | Stacked lollipops (one cohort up, one down) on shared domain map |
| 3D-cluster hotspot detection? | Use HotMAPS / 3D Hotspots — beyond linear lollipop |
| Druggable position? | Add ClinVar / OncoKB level annotation at the residue |

## maftools::lollipopPlot

**Goal:** Render per-gene mutation distribution on Pfam domain map with count-sized lollipops and class-colored stems.

**Approach:** Pass MAF and gene to `lollipopPlot`; maftools queries Pfam for domain coordinates automatically; outputs ggplot2 object.

```r
library(maftools)
maf <- read.maf(maf = 'cohort.maf')

# Default lollipop
lollipopPlot(maf = maf, gene = 'TP53',
             AACol = 'HGVSp_Short',
             labelPos = c(175, 248, 273),                   # mark canonical hotspots
             labPosSize = 1.0,
             showMutationRate = TRUE,
             domainLabelSize = 1,
             printCount = TRUE,
             colors = c(Missense_Mutation = '#D55E00',
                        Nonsense_Mutation = '#000000',
                        Frame_Shift_Del   = '#0072B2',
                        Frame_Shift_Ins   = '#56B4E9',
                        Splice_Site       = '#CC79A7',
                        In_Frame_Del      = '#009E73'))
```

```r
# Compare two cohorts -- one up, one down
lollipopPlot2(m1 = cohort_a, m2 = cohort_b,
              gene = 'TP53',
              m1_name = 'Cohort A',
              m2_name = 'Cohort B',
              AACol1 = 'HGVSp_Short', AACol2 = 'HGVSp_Short',
              colors = my_palette)
```

## trackViewer::lolliplot -- Fine Control over Track Layout

```r
library(trackViewer)
library(GenomicRanges)

# Build SNP (lollipop) and feature (domain) GRanges
snps <- GRanges('chr17', IRanges(c(175, 248, 273), width = 1, names = c('R175H', 'R248Q', 'R273H')),
                color = c('#D55E00', '#D55E00', '#D55E00'),
                score = c(45, 38, 29))                       # mutation count
features <- GRanges('chr17',
                    IRanges(c(102, 323, 363), width = c(190, 30, 30),
                            names = c('DNA-binding', 'Tetramerization', 'Regulatory')),
                    fill = c('#0072B2', '#009E73', '#CC79A7'),
                    height = 0.04)

lolliplot(snps, features, ylab = 'Mutation count',
          xaxis = TRUE, yaxis = TRUE)
```

trackViewer is more flexible than maftools for non-standard layouts (custom domain sources, multi-protein stacking, integration with genome coordinates).

## g3viz / g3-lollipop -- Interactive HTML

```r
library(g3viz)
mutation_data <- hgvspChange2protein(maf, gene = 'TP53')
g3Lollipop(mutation_data,
           gene.symbol = 'TP53',
           protein.change.col = 'AA_Change',
           plot.options = g3Lollipop.theme(theme.name = 'nature'),
           output.filename = 'TP53_lollipop.html')
```

g3-lollipop produces an interactive HTML — hover tooltips, click-to-filter, exportable. Suitable for supplementary HTML supplement; not for journal figure submission directly.

## Domain Annotation Sources

| Source | Format | Stability | Caveat |
|--------|--------|-----------|--------|
| Pfam (via maftools) | Pfam-A domain coordinates | Updated occasionally | maftools caches local; may lag Pfam release |
| UniProt | Domain + Region features (varied types) | Daily updates | API-driven; rate limits |
| InterPro | Integrated multi-database | More inclusive than Pfam | Different sub-classifications |
| Custom | Hand-curated for specific paper | Reproducible | Cite source |

For canonical isoform: maftools uses the canonical UniProt isoform by default. For specific isoform: pass `refSeqID` or `proteinID` explicitly. Mutations annotated against a different isoform will be off-by-residue.

## Per-Method Failure Modes

### Mutations not labeled with AA position

**Trigger:** MAF column `HGVSp_Short` missing or malformed.

**Mechanism:** maftools expects `HGVSp_Short` (e.g., 'p.R175H'); falls back to other columns inconsistently.

**Symptom:** "No mutations to plot" or wrong positions.

**Fix:** Verify `HGVSp_Short` column exists; reformat from HGVSp if needed. Use `AACol` argument to specify which column.

### Isoform mismatch

**Trigger:** Mutations called against ENST00000269305 but plotted against canonical ENST00000288602 (TP53).

**Mechanism:** Residue numbering differs across isoforms.

**Symptom:** Known R175H plotted at R177H or in a different domain.

**Fix:** Annotate the isoform in the figure caption; pass `proteinID` to `lollipopPlot` to force a specific isoform.

### Domain map outdated

**Trigger:** maftools' cached Pfam annotation is older than the protein's current Pfam release.

**Mechanism:** Domain coordinates can shift across Pfam versions.

**Symptom:** Domain boundaries off by a few residues; published-figure mismatch.

**Fix:** Pull domain coordinates from UniProt directly (current); pass via `trackViewer::lolliplot` features.

### Recurrence at low-coverage region overinterpreted

**Trigger:** "Hotspot" identified at a residue with high coverage variance — looks recurrent but is a sequencing artifact.

**Mechanism:** Capture-bait coverage variability; some residues sequenced more deeply.

**Symptom:** "Hotspot" in untargeted region; not validated in WGS.

**Fix:** Verify recurrence in independent cohort (TCGA Pan-Cancer + ICGC); use MutSig hotspot test (Lawrence 2014) for formal hotspot calling.

### Counts encoded only as size; no actual numbers shown

**Trigger:** Default `printCount = FALSE`.

**Mechanism:** Size-encoded counts beyond ~10 saturate visually.

**Symptom:** Reader cannot tell whether the top lollipop is 30 vs 300 mutations.

**Fix:** `printCount = TRUE` annotates each lollipop with its count.

### Domain colors random; no functional grouping

**Trigger:** Default rainbow domain colors.

**Mechanism:** Domains colored by accident, not by function class.

**Symptom:** Reader cannot quickly identify which domain is the kinase.

**Fix:** Manually map domain colors by functional class (kinase = blue, binding = green, regulatory = purple).

## Reconciliation: When Hotspots Disagree

| Pattern | Cause | Action |
|---------|-------|--------|
| Hotspot in cohort A absent in B | Cohort A enriched for a subtype OR small N | Stratify by subtype; cite both N |
| 3D hotspot test calls residues not on lollipop | Linear adjacency misses 3D proximity | Use HotMAPS / 3D Hotspots for spatial clusters |
| Recurrent residue lacks OncoKB evidence | Novel hotspot OR sequencing artifact | Confirm via independent cohort + WGS |
| Frame-shift indels not aligned to expected codon | Different annotation tool (VEP vs SnpEff) | Standardize annotation; verify HGVSp |

**Operational rule:** annotate the isoform; show absolute counts on lollipops; verify hotspots against TCGA Pan-Cancer + ICGC before novel-hotspot claims.

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Hotspot recurrence cutoff | depends on gene length + cohort size | Lawrence 2014 — formal MutSig test |
| Display all mutations vs filter | Recurrent (count ≥ 2) for clarity; show all in supplement | Visualization practical |
| Domain source default | Pfam (maftools default); UniProt for current | Tool-specific |
| Cohort N for credible hotspot | ≥200 for a single gene; pan-cancer for novel | Standard practice |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No mutations on plot | HGVSp_Short column missing | Verify / reformat from HGVSp |
| Mutations at wrong position | Isoform mismatch | Specify proteinID; document isoform |
| Domain boundaries slightly off | maftools Pfam cache outdated | Pull from UniProt; use trackViewer |
| Hotspot size saturates | Counts >10 indistinguishable by size | `printCount = TRUE` to annotate numbers |
| Random domain colors | Default rainbow | Manual mapping by functional class |
| Novel hotspot from one cohort | Insufficient N | Verify in TCGA + ICGC |

## References

- Chang MT, Asthana S, Gao SP, et al. 2016. Identifying recurrent mutations in cancer reveals widespread lineage diversity and mutational specificity. *Nat Biotechnol* 34(2):155-163.
- Gao J, Aksoy BA, Dogrusoz U, et al. 2013. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. *Sci Signal* 6(269):pl1.
- Lawrence MS, Stojanov P, Mermel CH, et al. 2014. Discovery and saturation analysis of cancer genes across 21 tumour types. *Nature* 505:495-501.
- Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Res* 28(11):1747-1756.
- Ou J, Zhu LJ. 2019. trackViewer: a Bioconductor package for interactive and integrative visualization of multi-omics data. *Nat Methods* 16:453-454.
- Zhou X, Edmonson MN, Wilkinson MR, et al. 2016. Exploring genomic alteration in pediatric cancer using ProteinPaint. *Nat Genet* 48(1):4-6.

## Related Skills

- data-visualization/oncoprint-mutation-matrices - Cohort-wide mutation matrix
- variant-calling/variant-annotation - Annotate HGVSp upstream
- clinical-databases/variant-prioritization - Filter variants before lollipop
- data-visualization/color-palettes - CVD-safe class palettes
- structural-biology/structure-navigation - 3D protein structure for hotspot interpretation
