---
name: bio-metagenomics-abundance
description: Turns shotgun classifier output into a defensible abundance table with Bracken Bayesian re-estimation, then compositional treatment (CLR, zero handling), library-size normalization, reference-frame differential abundance, and optional absolute quantification. Covers why a relative-abundance change is not a change, why Bracken read fractions and MetaPhlAn percentages are different physical quantities, the silent -r read-length bias, the genome-size confound no library-size method fixes, and the rarefaction debate. Use when estimating species abundance from a Kraken2 report, normalizing a community count table, choosing a compositional transform, or converting relative to absolute load. For classification see kraken-classification; for diversity/ordination/DA mechanics see metagenome-visualization.
tool_type: mixed
primary_tool: Bracken
---

## Version Compatibility

Reference examples tested with: Bracken 2.9+, Kraken2 2.1.3+, pandas 2.2+, scikit-bio 0.6+, R zCompositions 1.5+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `bracken -h`, `bracken-build -h` to confirm flags and defaults
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('zCompositions')` then `?cmultRepl` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The Bracken `databaseRLENmers.kmer_distrib` is built per database at a fixed read length and MUST match both the Kraken2 database and the actual (post-trim) read length; `-r` is not auto-detected and a mismatch silently biases every species fraction. Bracken needs the default Kraken2 report format, not mpa-style.

# Abundance Estimation

**"How much of each taxon is in my sample?"** -> Re-estimate species reads with Bracken, then treat the table as a composition - because the sequencer fixed the total, so the numbers are relative and a change in one forces apparent changes in the rest.
- CLI: `bracken -d DB -i kraken.kreport -o out.bracken -w out.bracken.kreport -r 150 -l S -t 10`

Scope: Bracken mechanics plus what happens to the numbers afterward - estimand choice, compositional transforms, normalization, reference-frame DA logic, absolute conversion. Read classification -> kraken-classification, metaphlan-profiling. Diversity/ordination/DA tool mechanics -> metagenome-visualization. Generic 16S diversity -> the microbiome category.

## The Single Most Important Modern Insight -- A Relative-Abundance Change Is Not a Change

A shotgun abundance table is a composition: the sequencer fixes the total number of reads, not the sample's microbial load. So every number is relative to every other, and "taxon X went up 2-fold" is undefined without a reference frame or an external load anchor. Without one, a single blooming taxon makes every other taxon look depleted (the blooming-taxon illusion), and a Pearson correlation of two taxa's proportions is biased negative by arithmetic, not biology. Bracken is step one of about six: classify -> re-estimate -> compositional transform -> normalize -> reference-frame test -> (optional) absolute conversion. Two corollaries:

1. **Bracken percent and MetaPhlAn percent are different physical quantities.** Bracken is a fraction of READS (large-genome biased, approximately fraction of DNA); MetaPhlAn is a genome-size-normalized marker abundance (approximately fraction of cells). They legitimately disagree 2-3x for the same sample. Picking one chooses an estimand (DNA vs cells), not the "more accurate tool."
2. **No library-size normalizer fixes the genome-size confound.** TSS/CSS/TMM/GMPR all operate on the count table and inherit the large-genome bias; only a coverage-based estimand (CoverM) or a genome-normalized profiler (MetaPhlAn) removes it.

## What Quantity Is Being Estimated?

| Estimand | Definition | Bias / use |
|----------|------------|-----------|
| Read count | raw reads to a taxon | library-size and genome-size dependent; never compare raw across samples |
| Relative abundance | read count / total | compositional (sums to 1); still genome-size biased |
| Coverage abundance | reads x readlen / genome size | removes large-genome bias; ~ fraction of genomes/cells (CoverM) |
| Cell fraction | fraction of organisms | needs coverage + an assumption of one genome per cell; polyploidy/growth bias it |
| Absolute load | composition x external total | the only basis for "increased/decreased" (flow/qPCR/spike-in) |

## Bracken: Step One, Not the Finish Line

Bracken redistributes reads stranded at the genus/family node (where shared k-mers stopped Kraken) down to species, using a database-derived expectation of where length-L reads classify. It does not classify and it does not add precision.

```bash
bracken-build -d "$KRAKEN_DB" -t 8 -k 35 -l 150   # one-time; -k MUST equal the Kraken2 build k (35)
bracken -d "$KRAKEN_DB" -i kraken.kreport -o out.bracken -w out.bracken.kreport \
    -r 150 \   # MUST equal the bracken-build -l AND the actual post-trim read length (not auto-detected)
    -l S -t 10 # species level; -t drops taxa with fewer than 10 clade-level reads (strict <) before redistribution
```

Output columns: `name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads`. `fraction_total_reads` is a fraction of classified-and-retained reads, not of all input - so two samples with different unclassified (e.g. host) fractions have non-comparable denominators. `combine_bracken_outputs.py --files *.bracken -o matrix.tsv` builds a taxa-by-sample matrix.

### Bracken's defining failure mode
Bracken can only redistribute among species already in the database. A true organism absent from the database has its reads parked by Kraken at the shared genus node, and Bracken hands them to the database-present congeners - confidently fabricating or inflating those species. Distrust any species with high `added_reads` but tiny `kraken_assigned_reads`; gate presence on upstream unique-minimizer evidence (kraken-classification) and a coverage breadth check.

## Treat the Table as a Composition

**Goal:** Make the abundance table valid for multivariate stats and correlation by removing closure with a centered log-ratio, after handling zeros (log of zero is undefined).

**Approach:** Impute count zeros with Bayesian-multiplicative replacement (preserves ratios), then CLR-transform; use Aitchison distance (Euclidean on CLR) downstream, never raw-proportion Pearson or Bray-Curtis for correlation.

```python
import pandas as pd
import numpy as np
from skbio.stats.composition import clr, multi_replace   # multiplicative_replacement was renamed multi_replace in skbio 0.6

counts = pd.read_csv('matrix.tsv', sep='\t', index_col=0)   # taxa x samples (new_est_reads)
mat = counts.T.values.astype(float)                          # samples x taxa for transform
mat_nozero = multi_replace(mat / mat.sum(axis=1, keepdims=True))
clr_mat = clr(mat_nozero)                                    # closure removed; rows are CLR coordinates
clr_df = pd.DataFrame(clr_mat, index=counts.columns, columns=counts.index)
```

For sparse tables prefer R `zCompositions::cmultRepl()` (Bayesian-multiplicative, posterior-imputed) over a fixed +1 pseudocount, which is arbitrary, distorts ratios, and re-opens closure. Structural zeros (taxon genuinely absent) and sampling zeros (present below detection) are usually indistinguishable from the table - disclose the assumption rather than pretend otherwise.

## Library-Size Normalization

| Method | What it does | Shotgun applicability |
|--------|--------------|-----------------------|
| TSS (proportions) | divide by library size | IS the compositional closure; fine for viz, biased for DA/correlation |
| Rarefaction | subsample to common depth | discards data; defensible for diversity, contested for DA (see below) |
| CSS (metagenomeSeq) | scale by a cumulative-sum quantile | robust to dominant taxa; designed for marker surveys, usable on counts |
| TMM (edgeR) | trimmed mean of M-values | "most features unchanged" often violated in microbiome; use with caution |
| RLE (DESeq median-of-ratios) | geometric-mean reference | breaks on zeros (geometric mean -> 0); needs poscounts workaround |
| GMPR | pairwise median ratios then geometric-mean | purpose-built for zero-inflated counts; good default size factor |

None of these fixes the genome-size confound - that needs a coverage estimand or a genome-normalized profiler.

### The rarefaction debate (do not take a side; decide per analysis)
McMurdie & Holmes 2014 (*PLoS Comput Biol* 10:e1003531) showed rarefying for DIFFERENTIAL ABUNDANCE is statistically wasteful versus modeling library size. Schloss 2024 (*mSphere* 9:e00354-23 and the companion e00355-23) argues rarefaction is currently the best control of uneven effort for RICHNESS and COMMUNITY-DISTANCE analyses, where scaling alone does not remove the depth effect. They concern different downstream analyses: for DA do not rarefy (use a CoDA/reference-frame method or a modeled size factor); for alpha/beta diversity rarefaction is defensible. Treating "rarefy: yes/no" as one global switch is the mistake.

## Reference-Frame Differential Abundance (logic here; mechanics in metagenome-visualization)

"Taxon X increased" needs a reference frame because the total is fixed. Methods differ chiefly by their implicit frame: total-sum (naive, wrong), the geometric mean of all taxa (CLR / ALDEx2), or an estimated per-sample sampling fraction (ANCOM-BC). ALDEx2 (Fernandes 2014 *Microbiome* 2:15) Monte-Carlo samples a Dirichlet posterior, CLR-transforms, and tests each draw; ANCOM-BC (Lin & Peddada 2020 *Nat Commun* 11:3514) estimates and corrects each sample's sampling fraction. Do not run naive t-tests or Wilcoxon on TSS proportions. Run the tools and read their output in metagenome-visualization.

## Absolute Quantification: the Only Basis for "Increased/Decreased"

Relative methods recover the composition; absolute load = composition x an external total. Anchors: flow-cytometry cell counts (QMP, Vandeputte 2017 *Nature* 551:507, which showed apparent Crohn's "increases" were a microbial-load artifact - the absolute trajectory was opposite); per-taxon flow density (Props 2017 *ISME J* 11:584); a known spike-in organism added before extraction (SCML, Stammler 2016 *Microbiome* 4:28; a cellular spike also captures extraction bias a DNA spike misses); or total 16S qPCR copies (cheap, copy-number biased). Any "X increased" claim without an anchor or a reference-frame method is unsupported.

## Per-Method Failure Modes

### Bracken `-r` read-length mismatch
**Trigger:** `-r 150` on trimmed ~120 bp reads, or a database built only for a different length. **Mechanism:** the redistribution prior is fragment-length specific and not auto-detected. **Symptom:** biased species fractions with no error if the `.kmer_distrib` exists, hard crash if not. **Fix:** `-r` = actual post-trim read length and a matching built distribution.

### Pearson/Bray-Curtis on proportions
**Trigger:** correlating two taxa's relative abundances, or Bray-Curtis distance for a "who-co-occurs" claim. **Mechanism:** closure biases proportion correlations negative and makes Euclidean/Bray-Curtis sub-compositionally incoherent. **Symptom:** spurious negative interactions; unstable clustering. **Fix:** CLR + Aitchison distance, or SparCC/proportionality for co-occurrence.

### Relative fold-change reported as absolute
**Trigger:** "taxon X doubled" from a relative table. **Mechanism:** one bloom deflates every other proportion. **Symptom:** whole-community "depletion" that is really one taxon rising. **Fix:** anchor to load (flow/qPCR/spike-in) or use ANCOM-BC; state which frame.

### RLE/DESeq on a sparse table
**Trigger:** DESeq2/edgeR median-of-ratios on a species count table. **Mechanism:** the geometric-mean reference collapses to 0 when any feature has a zero. **Symptom:** degenerate size factors, errors, or nonsense fold-changes. **Fix:** GMPR or CSS; or DESeq2 with a poscounts estimator.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Bracken `-k` = 35 | Lu 2017 *PeerJ Comput Sci* 3:e104 | must equal the Kraken2 database k-mer length |
| Bracken `-r` = post-trim read length | Bracken docs | redistribution prior is fragment-length specific; silent bias otherwise |
| Bracken `-t` 10 default | Bracken docs | redistribution floor; raise to suppress noise, too high deletes real rare taxa |
| Multiplicative/Bayesian zero replacement over +1 | Martin-Fernandez 2015 *Stat Modelling* 15:134 | preserves ratios; fixed pseudocount distorts them and re-opens closure |
| Rarefy for diversity, not for DA | McMurdie 2014; Schloss 2024 *mSphere* 9:e00354-23 | decision is per-analysis, not global |
| Coverage breadth (`covered_fraction`) presence gate | Aroney 2025 *Bioinformatics* 41:btaf147 | high mean coverage over a few % of a genome is a conserved-region artifact |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "kmer_distrib not found" | `-r` has no matching built distribution | `bracken-build -l <readlen>` or pick a DB shipping it |
| Bracken rejects input | mpa-style or per-read file passed to `-i` | give the default Kraken2 `--report` |
| Species with huge added_reads, tiny assigned | redistribution into a DB-present relative of an absent taxon | gate presence on unique minimizers + coverage breadth |
| CLR returns -inf / NaN | zeros not replaced before log-ratio | multiplicative or `cmultRepl` replacement first |
| Cross-sample fractions not comparable | differing unclassified (host) fractions in the denominator | track classified fraction as a covariate; host-deplete upstream |
| MetaPhlAn and Bracken disagree | different estimands (cells vs reads) | do not merge; pick one estimand and state it |

## References

- Lu J, Breitwieser FP, Thielen P, Salzberg SL. 2017. Bracken: estimating species abundance in metagenomics data. *PeerJ Comput Sci* 3:e104.
- Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. 2017. Microbiome datasets are compositional: and this is not optional. *Front Microbiol* 8:2224.
- Quinn TP, Erb I, Richardson MF, Crowley TM. 2018. Understanding sequencing data as compositions: an outlook and review. *Bioinformatics* 34:2870-2878.
- Fernandes AD, Reid JN, Macklaim JM, et al. 2014. Unifying the analysis of high-throughput sequencing datasets. *Microbiome* 2:15.
- Lin H, Peddada SD. 2020. Analysis of compositions of microbiomes with bias correction. *Nat Commun* 11:3514.
- Martin-Fernandez JA, Hron K, Templ M, Filzmoser P, Palarea-Albaladejo J. 2015. Bayesian-multiplicative treatment of count zeros in compositional data sets. *Stat Modelling* 15:134-158.
- McMurdie PJ, Holmes S. 2014. Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Comput Biol* 10:e1003531.
- Schloss PD. 2024. Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere* 9:e00354-23.
- Weiss S, Xu ZZ, Peddada S, et al. 2017. Normalization and microbial differential abundance strategies depend upon data characteristics. *Microbiome* 5:27.
- Chen L, Reeve J, Zhang L, et al. 2018. GMPR: a robust normalization method for zero-inflated count data. *PeerJ* 6:e4600.
- Vandeputte D, Kathagen G, D'hoe K, et al. 2017. Quantitative microbiome profiling links gut community variation to microbial load. *Nature* 551:507-511.
- Stammler F, Glasner J, Hiergeist A, et al. 2016. Adjusting microbiome profiles for differences in microbial load by spike-in bacteria. *Microbiome* 4:28.
- Aroney STN, Newell RJP, Nissen JN, Camargo AP, Tyson GW, Woodcroft BJ. 2025. CoverM: read alignment statistics for metagenomics. *Bioinformatics* 41:btaf147.

## Related Skills

- kraken-classification - Generates the Kraken2 report and owns the read-count-is-not-abundance reframe
- metaphlan-profiling - Genome-size-normalized cell-fraction abundance; a different estimand
- metagenome-visualization - Diversity, ordination, and differential-abundance tool mechanics
- contamination-controls - Host-fraction and blank handling that affect the denominator
- genome-assembly/metagenome-assembly - Coverage-based MAG abundance via read mapping
- workflows/metagenomics-pipeline - End-to-end profiling and abundance
