---
name: bio-microbiome-amplicon-processing
description: Infers exact amplicon sequence variants (ASVs) from demultiplexed 16S rRNA or ITS amplicon FASTQ with DADA2 - removing primers with cutadapt (--discard-untrimmed), learning a per-run error model (filterAndTrim -> learnErrors -> dada -> mergePairs), merging run-level tables with mergeSequenceTables, then one removeBimeraDenovo. Covers why primers come OFF before truncation, why the error model is per-run, truncLen as a merge-overlap detection budget (V4 vs V3-V4), DADA2 vs Deblur and q2-dada2 (denoise-paired/single/pyro/ccs), ASV vs OTU, NovaSeq binned-quality error-fit breakage, ITSxpress for variable-length ITS, and decontam removal of reagent/kit contaminants. Use when turning demultiplexed amplicon reads into an ASV/feature table, choosing truncation lengths, handling multi-run studies, or ITS. For shotgun reads see metagenomics/kraken-classification; for QIIME2 CLI mechanics see qiime2-workflow; for primer trimming theory see read-qc/adapter-trimming.
tool_type: mixed
primary_tool: DADA2
---

## Version Compatibility

Reference examples tested with: DADA2 1.30+, cutadapt 4.6+, ITSxpress 2.0+, QIIME2 2024.2+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The error model is a PER-RUN artifact, not a version: `learnErrors` is fit to one sequencing run (flowcell/chemistry/instrument). Multi-run studies run the per-run inference separately, then `mergeSequenceTables`, then a single chimera removal - never pool FASTQs across runs before `learnErrors`. DADA2 `dada()` defaults (`OMEGA_A` 1e-40, `mergePairs` `minOverlap` 12) and QIIME2 plugin flag spellings drift between releases; confirm with `?dada` and `qiime dada2 --help`.

# Amplicon Processing with DADA2

**"Process my 16S amplicon data to get ASVs"** -> Strip primers, learn a per-run error model, denoise into exact amplicon sequence variants, merge pairs, and remove chimeras - because an ASV is a model-inferred sequence conditioned on one run, not a clustered consensus or an organism.
- R: `dada(filtFs, err=learnErrors(filtFs, multithread=TRUE), multithread=TRUE)`
- CLI: `cutadapt -g FWD -G REV --discard-untrimmed ...` then DADA2, or `qiime dada2 denoise-paired`

Scope: demultiplexed amplicon reads -> chimera-free ASV/feature table + representative sequences. Shotgun reads -> metagenomics/kraken-classification. Taxonomy of the ASVs -> taxonomy-assignment. Diversity/DA of the table -> diversity-analysis, differential-abundance. Compositional/normalization theory (shared) -> metagenomics/abundance-estimation. QIIME2 artifact/provenance/demux mechanics -> qiime2-workflow. Primer-trimming theory -> read-qc/adapter-trimming.

## The Single Most Important Modern Insight -- An ASV Is a Denoiser's Output on One Run, Not a Ground-Truth Organism

The feature table is not an observation of the community; it is the residue of modeling decisions made BEFORE any result exists - which primers were stripped, where reads were truncated, what error model the run's quality scores supported, what was called a chimera. Turn the knobs differently and the table changes. Three corollaries each common misuse violates:

1. **Primers and truncLen decide what is detectable, silently.** Leftover primers corrupt the error model (mismatches read as sequencing error) and masquerade as chimeras; truncating reads below the merge-overlap budget erases taxa by arithmetic, not biology. These knobs are set before the answer exists - declare them.
2. **The error model is fit PER RUN.** Illumina error rates are run-specific. Concatenating runs before `learnErrors` fits one model to a mixture of error structures and denoises wrong. Infer each run separately, then `mergeSequenceTables` (the exact-sequence string is the join key), then one chimera removal.
3. **An ASV is an exact sequence, not a cell, genome, or species.** One genome carries multiple, often divergent 16S copies, so one organism becomes several ASVs and inflates richness (Schloss 2021). Reads are not cells (16S copy number varies); a species-level 16S call is usually overconfident.

Organize the work around declaring and defending these knobs - not around running `dada()` and calling the columns "species."

## ASV vs OTU -- the Methodological Fork

An ASV (DADA2/Deblur) is an exact inferred sequence at single-nucleotide resolution; a 97% OTU is a centroid of a 3%-identity cluster. Both sides are live (present both, do not declare a winner):

- **ASVs replace OTUs** (Callahan 2017 *ISME J* 11:2639): the sequence IS the identity, so ASVs are portable across studies without re-clustering, higher-resolution, and reproducible (no clustering-order/abundance dependence). This is the field default.
- **ASVs over-split genomes** (Schloss 2021 *mSphere* 6:e00191-21; Pan 2023 *Appl Environ Microbiol* 89:e02108-22): intragenomic 16S copy heterogeneity (in ~60% of prokaryotes; *E. coli* K-12 has 7 copies in ~5 sequence types) makes one organism appear as several ASVs, inflating richness; 97% OTUs lump those copies back. Defensible practice: use ASVs, but treat ASV count as an upper bound on richness and collapse to a taxonomic rank (taxonomy-assignment) before richness claims.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| DADA2 | Callahan 2016 *Nat Methods* 13:581 | per-run parametric error model, abundance-partition denoising, merge, chimera | the default; variable length, ITS, singleton sensitivity via pseudo-pooling |
| q2-dada2 | (DADA2 engine; Bolyen 2019 *Nat Biotechnol* 37:852) | QIIME2 wrapper: `denoise-paired`/`single`/`pyro`/`ccs` | DADA2 inside a QIIME2 artifact/provenance workflow -> qiime2-workflow |
| Deblur | Amir 2017 *mSystems* 2:e00191-16 | static upper-bound Illumina error profile (positive filter), one fixed length | fast, per-sample-independent, trivially combinable runs; 16S only |
| cutadapt | Martin 2011 *EMBnet J* 17:10 | primer/adapter trimming (`-g`/`-G`, linked adapters) | MUST run before filterAndTrim; primer removal -> read-qc/adapter-trimming |
| ITSxpress | Rivers 2018 *F1000Research* 7:1418 | HMM-trims the variable-length ITS spacer, keeping quality scores | ITS only; ITS has no valid fixed truncLen |
| VSEARCH | Rognes 2016 *PeerJ* 4:e2584 | open-source 97% OTU clustering, dereplication, chimera | the OTU path, if a 97% clustering is required (legacy) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| 16S V4 (~253 bp), 2x250 | DADA2 paired, truncLen with comfortable overlap | huge merge slack; truncate to quality freely |
| 16S V3-V4 (~460 bp), 2x250 | DADA2 paired, protect overlap; loosen `maxEE` R | only ~28 bp slack - the merge budget dominates quality |
| ITS (variable length) | cutadapt + ITSxpress + DADA2 `truncLen=0` | fixed truncation slices real biology and breaks merging |
| Full-length 16S (PacBio HiFi/CCS) | DADA2 / `qiime dada2 denoise-ccs` | resolves to species/strain; single-end CCS, not paired |
| Multiple sequencing runs | per-run inference -> `mergeSequenceTables` -> one chimera removal | error model is per-run; never pool FASTQs first |
| Want speed, fixed length, many runs, 16S only | Deblur (`denoise-16S`) | static positive filter; per-sample independent |
| Need singleton/rare-ASV sensitivity | DADA2 `dada(..., pool='pseudo')` | pseudo-pooling approximates full pooling in linear time |
| NovaSeq/NextSeq/iSeq (binned Q) | inspect `plotErrors`; enforce monotonic error fit | ~4 quality bins starve the loess fit -> wrong denoising |
| Shotgun (random WGS) reads, not amplicon | -> metagenomics/kraken-classification | no primers/per-run denoising; different category |

## Remove Primers First (cutadapt)

**Goal:** Strip synthetic, often-degenerate primer sequence before any quality/error step.

**Approach:** Match the forward primer as a 5' adapter on R1 and the reverse primer on R2, discarding pairs where the primer is absent. The order primers -> filter -> learn-errors is non-negotiable: leftover primers corrupt the error model, shift the truncLen frame, and inflate chimeras.

```bash
# -g = 515F forward primer (5' adapter on R1); -G = 806R reverse primer (5' adapter on R2);
# --discard-untrimmed drops pairs lacking the primer (a primerless read is suspect).
cutadapt \
    -g GTGYCAGCMGCCGCGGTAA \
    -G GGACTACNVGGGTWTCTAAT \
    --discard-untrimmed \
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
    sample_R1.fastq.gz sample_R2.fastq.gz
```

The QIIME2 equivalent is `qiime cutadapt trim-paired --p-front-f FWD --p-front-r REV --p-discard-untrimmed`.

## The Per-Run DADA2 Pipeline

**Goal:** Turn one run's primer-trimmed FASTQs into a denoised, merged sequence table.

**Approach:** Filter on expected errors and truncate within the merge budget, learn the run's error model, denoise each read set against it, merge pairs, then tabulate. Run this block once PER sequencing run.

```r
library(dada2)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(240, 160),     # region/read-length specific; subject to the merge budget below
                     maxEE=c(2, 2), truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
errF <- learnErrors(filtFs, multithread=TRUE)  # fit THIS run only
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)                # observed points must track the fitted line and fall with Q
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)   # pool='pseudo' for rare-ASV sensitivity
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab_run <- makeSequenceTable(mergers)
```

### truncLen Is a Detection Budget, Not a Quality Setting

Paired-end merging needs `truncLen_F + truncLen_R >= amplicon_length + ~12` (DADA2 `minOverlap` default is 12). truncLen is jointly constrained by quality (cut where median Q drops below ~Q30 on `plotQualityProfile`) AND this overlap budget; the two fight, and for long amplicons the budget wins.

- **V4 (515F/806R, ~253 bp), 2x250:** 250+250 vs 253+12 leaves huge slack - truncate to quality freely (e.g. `c(240, 200)`).
- **V3-V4 (341F/805R, ~460 bp), 2x250:** 250+250 vs 460+12 leaves only ~28 bp slack. Aggressive truncation of both reads kills the overlap and the merge rate collapses to near zero. Preserve length: barely truncate the reverse and loosen `maxEE` to `c(2, 5)` to keep low-Q reverse reads.

A merge cliff in the read-tracking table is a budget problem, not bad data - the taxa were erased by arithmetic.

## Combine Runs, Then Remove Chimeras

**Goal:** Merge per-run sequence tables into one study table and remove PCR chimeras once.

**Approach:** Join run-level tables by exact sequence string, then detect bimeras (an ASV reconstructable from two more-abundant parents) across the combined table.

```r
st_all <- mergeSequenceTables(seqtab_run1, seqtab_run2)   # exact-sequence string is the join key
seqtab_nochim <- removeBimeraDenovo(st_all, method='consensus', multithread=TRUE, verbose=TRUE)
sum(seqtab_nochim) / sum(st_all)   # chimeras = many ASVs but few READS (~0.8-0.99 retained)
```

Carry "run" forward as a batch covariate into differential abundance. A large READ fraction removed as chimeric is a leftover-primer smell (degenerate bases look chimeric), not a real chimera storm.

## Decontamination and Controls (low-biomass)

**Goal:** Identify and remove reagent/kit ("kitome") contaminant ASVs before any downstream analysis - decisive for low-biomass samples, where contaminants can outnumber real signal.

**Approach:** Sequence negative controls (extraction blanks, no-template PCR) and a positive mock community alongside the samples, then classify contaminant ASVs with decontam (Davis 2018): the prevalence method when only controls are available, the frequency method when per-sample DNA concentration was measured, combined when both.

```r
library(decontam)
# seqtab_nochim is samples (rows) x ASVs (cols) - decontam's expected orientation.
# is_control: logical, TRUE for negative-control samples; dna_conc: per-sample DNA concentration (qPCR/Qubit).
# prevalence-only threshold 0.1 default; 0.5 = aggressive (ASV more prevalent in controls than samples = contaminant).
contam <- isContaminant(seqtab_nochim, neg = meta$is_control, conc = meta$dna_conc, method = 'combined', threshold = 0.1)
seqtab_clean <- seqtab_nochim[, !contam$contaminant]
```

Low-biomass samples (skin, biopsy, BAL, sterile-site swabs) can be dominated by the kitome, so a "community" there may be mostly contamination - never interpret a low-biomass result without controls. The shotgun analogue is metagenomics/contamination-controls.

## ITS: Never Fixed-Truncate

**Goal:** Isolate the biologically variable-length ITS spacer without slicing real sequence.

**Approach:** Strip primers with cutadapt, then HMM-trim the conserved SSU/5.8S/LSU flanks with ITSxpress (preserving quality scores), then denoise with `truncLen=0`, filtering on `maxEE`/`minLen` only.

```bash
itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz \
    --region ITS2 --taxa Fungi \   # ITS1/ITS2/ALL; --taxa selects the HMM model
    --outfile trimmed.fastq.gz --threads 4
```

```r
out_its <- filterAndTrim(trimmed, filtered, truncLen=0,   # NEVER fix-truncate ITS (variable length)
                         maxEE=2, minLen=50, maxN=0, rm.phix=TRUE, multithread=TRUE)
```

## QIIME2 and Deblur Equivalents

DADA2 inside QIIME2: `qiime dada2 denoise-paired --p-trunc-len-f --p-trunc-len-r` (also `denoise-single`, `denoise-pyro` for 454/Ion Torrent, `denoise-ccs` with `--p-front`/`--p-adapter`/`--p-min-len`/`--p-max-len` for PacBio CCS). Deblur (static positive filter, one fixed length, 16S only):

```bash
qiime deblur denoise-16S --i-demultiplexed-seqs qc.qza \
    --p-trim-length 250 --p-sample-stats \   # ONE fixed length; Deblur cannot handle variable length
    --o-representative-sequences rep-seqs.qza --o-table table.qza --o-stats stats.qza
```

Do not merge a DADA2 ASV table with a Deblur sOTU table - different feature definitions.

## Per-Method Failure Modes

### Primers left on before truncation
**Trigger:** running filterAndTrim/learnErrors on reads that still carry primers. **Mechanism:** synthetic, often-degenerate primer bases are read as sequencing error and create spurious split points. **Symptom:** wrong error fit, a huge READ fraction removed as chimeric, inflated ASV count. **Fix:** cutadapt `--discard-untrimmed` first; order is primers -> filter -> learnErrors.

### Pooling runs before learnErrors
**Trigger:** concatenating multiple runs' FASTQs into one pipeline. **Mechanism:** one error model is fit to a mixture of run-specific error structures. **Symptom:** distorted denoising; ASVs that vanish or appear when runs are split. **Fix:** per-run inference, then `mergeSequenceTables`, then one chimera removal; carry run as a batch covariate.

### Merge cliff from over-truncation
**Trigger:** truncLen_F + truncLen_R below amplicon length + 12. **Mechanism:** denoised pairs no longer overlap enough to merge. **Symptom:** near-zero `merged` column in read tracking; misread as "low diversity"/"bad data". **Fix:** compute the budget from amplicon and read length first; for long amplicons keep length and loosen `maxEE` R.

### Fixed-truncating ITS
**Trigger:** any `truncLen` on ITS. **Mechanism:** ITS length is biological (ITS1 ~200-600 bp), so a fixed cut slices real sequence off long variants and merge-fails short ones. **Symptom:** lost long fungal taxa, poor merging. **Fix:** cutadapt + ITSxpress, then `truncLen=0`, filter on `maxEE`/`minLen`.

### NovaSeq/NextSeq binned-quality error fit
**Trigger:** default `learnErrors` on ~4-bin quality data. **Mechanism:** the loess error-vs-Q fit is starved and can become non-monotonic (error rising at high Q). **Symptom:** in `plotErrors` the fitted line diverges from observed points. **Fix:** enforce monotonicity in the error matrix (nf-core/ampliseq `--illumina_novaseq`, or set sub-max-Q entries to the max-Q error); never trust the default fit on binned Q.

### ASV count read as species richness
**Trigger:** reporting ASV count as richness or each ASV as one organism. **Mechanism:** intragenomic 16S copy divergence splits one genome into several ASVs (Schloss 2021); reads are not cells (copy number 1-15+). **Symptom:** inflated richness, "species" that are copies of one organism. **Fix:** collapse to genus/species (taxonomy-assignment) before richness claims; treat ASV count as an upper bound.

### Low-biomass contamination ignored (no controls / no decontam)
**Trigger:** analysing low-biomass samples (skin, biopsy, BAL, sterile site) without sequencing controls or running decontam. **Mechanism:** reagent/kit DNA (the kitome) is amplified alongside scarce template and can dominate the reads. **Symptom:** a plausible "community" in a near-sterile sample; reagent-associated genera prominent; results track DNA yield. **Fix:** sequence extraction-blank + no-template-PCR negatives (and a positive mock), run decontam (prevalence or combined), report what was removed (Davis 2018; metagenomics/contamination-controls).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `maxEE` c(2,2) (loosen R to 5 for long amplicons) | Callahan 2016 *Nat Methods* 13:581 | expected-errors filter beats a hard Q cutoff; computed on the TRUNCATED read, so it interacts with truncLen |
| truncLen budget: truncLen_F + truncLen_R >= amplicon_len + 12 | DADA2 `mergePairs` `minOverlap` default | below this, denoised pairs cannot merge; the merge cliff is arithmetic, not data |
| truncLen cut where median Q < ~25-30 | DADA2 docs | quality target, secondary to the merge budget for long amplicons |
| `maxN` = 0 | DADA2 docs | DADA2 cannot model ambiguous bases; mandatory |
| chimera retained-read fraction ~0.8-0.99 | DADA2 docs | chimeras are many ASVs but few reads; a large read loss flags leftover primers |
| `pool='pseudo'` for rare ASVs | DADA2 docs | approximates full pooling (quadratic) in linear time; default FALSE misses cross-sample singletons |
| Deblur `--p-trim-length` one fixed value | Amir 2017 *mSystems* 2:e00191-16 | the positive filter requires a single read length |
| 16S copy-number correction: report, do not assume | Louca 2018 *Microbiome* 6:41 | predictable only near reference genomes; correction can ADD error ("unsolved problem") |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Near-zero merge rate | truncLen below the overlap budget | recompute budget; keep length, loosen `maxEE` R |
| Large read fraction "chimeric" | primers not trimmed (degenerate bases) | cutadapt `--discard-untrimmed` before filtering |
| `plotErrors` fitted line diverges from points | binned quality (NovaSeq/NextSeq) | enforce monotonic error matrix; nf-core/ampliseq `--illumina_novaseq` |
| ASVs vanish/appear when runs split | one error model fit across runs | per-run `learnErrors`, then `mergeSequenceTables` |
| Few reads pass filter | `maxEE` too strict or truncLen too long (low-Q tail) | loosen `maxEE`, shorten truncLen within the budget |
| ITS taxa lost / poor merging | fixed `truncLen` on ITS | cutadapt + ITSxpress, then `truncLen=0` |

## References

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: high-resolution sample inference from Illumina amplicon data. *Nat Methods* 13:581-583.
- Callahan BJ, McMurdie PJ, Holmes SP. 2017. Exact sequence variants should replace operational taxonomic units in marker-gene data analysis. *ISME J* 11:2639-2643.
- Callahan BJ, Wong J, Heiner C, Oh S, Theriot CM, Gulati AS, McGill SK, Dougherty MK. 2019. High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution. *Nucleic Acids Res* 47:e103.
- Amir A, McDonald D, Navas-Molina JA, Kopylova E, Morton JT, Zech Xu Z, Kightley EP, Thompson LR, Hyde ER, Gonzalez A, Knight R. 2017. Deblur rapidly resolves single-nucleotide community sequence patterns. *mSystems* 2:e00191-16.
- Martin M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet J* 17:10-12.
- Rivers AR, Weber KC, Gardner TG, Liu S, Armstrong SD. 2018. ITSxpress: software to rapidly trim internally transcribed spacer sequences with quality scores for marker gene analysis. *F1000Research* 7:1418.
- Rognes T, Flouri T, Nichols B, Quince C, Mahe F. 2016. VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4:e2584.
- Bolyen E, Rideout JR, Dillon MR, et al. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nat Biotechnol* 37:852-857.
- Schloss PD. 2021. Amplicon sequence variants artificially split bacterial genomes into separate clusters. *mSphere* 6:e00191-21.
- Pan P, et al. 2023. Microbial diversity biased estimation caused by intragenomic heterogeneity and interspecific conservation of 16S rRNA genes. *Appl Environ Microbiol* 89:e02108-22.
- Louca S, Doebeli M, Parfrey LW. 2018. Correcting for 16S rRNA gene copy numbers in microbiome surveys remains an unsolved problem. *Microbiome* 6:41.
- Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. 2018. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. *Microbiome* 6:226.

## Related Skills

- taxonomy-assignment - Assign taxonomy to the ASVs produced here
- diversity-analysis - Alpha/beta diversity of the resulting community table
- differential-abundance - Compositional DA on the ASV/feature table
- qiime2-workflow - The QIIME2 CLI equivalent of this R workflow
- read-qc/adapter-trimming - cutadapt primer removal before DADA2
- metagenomics/kraken-classification - Shotgun (not amplicon) read classification
- metagenomics/abundance-estimation - Shared compositional/normalization theory
- metagenomics/contamination-controls - Negative/positive controls and decontam for low-biomass (shotgun analogue)
- phylogenetics/tree-io - Phylogenetic tree for UniFrac / Faith PD
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
