---
name: bio-ecological-genomics-edna-metabarcoding
description: Processes eDNA metabarcoding from raw paired-end reads to species tables, navigating ASV (DADA2, UNOISE3) vs OTU (swarm v2) decision (Callahan 2017 vs Schloss multi-copy-16S critique), marker/primer choice (Leray COI, MiFish 12S, 515F/806R 16S, ITS2) with primer-specific bias, OBITools3 v3 command-name break (obi stats plural; .tar.gz taxonomy), tag-jumping with dual-indexing (Schnell 2015; NovaSeq 10x MiSeq), decontam as screening-not-classifier (Davis 2018), read-counts-not-abundance critique (Lamb 2019), site-occupancy modeling (Ficetola 2015), Naive-Bayes calibration limits (Bokulich 2018), and eDNA decay (Strickler 2015). Use when going from raw eDNA FASTQ to species tables, picking marker + denoising pipeline, deciding whether read counts represent abundance, applying occupancy modeling, configuring OBITools3 v3, or interpreting decontam output. Not for clinical 16S microbiome (see microbiome/amplicon-processing).
tool_type: mixed
primary_tool: dada2
---

## Version Compatibility

Reference examples tested with: DADA2 1.30+, cutadapt 4.7+, OBITools3 (Python 3), decontam 1.20+, microDecon 1.0+, occumb 1.0+, vsearch 2.27+, swarm 3.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# eDNA Metabarcoding

**"Process eDNA samples to identify species present"** -> Trim primers, denoise to ASVs (or cluster to OTUs), detect chimeras, assign taxonomy, filter contamination with negative controls AND DNA concentration, decompose tag-jumping artifacts, and quantify detection uncertainty via site-occupancy modeling. For the foundational eDNA-for-wildlife review, see Bohmann et al. 2014 *Trends Ecol Evol* 29:358-367.
- CLI: `cutadapt` for primer removal (linked-adapter mode)
- R: `dada2::filterAndTrim()` -> `dada()` -> `assignTaxonomy()` for ASV pipeline
- CLI: `obi stats` / `obi clean` / `obi ecotag` for OBITools3 (NOTE: v3 plural commands)
- R: `decontam::isContaminant()` for contamination screening
- R: `occumb::occumb()` for detection-corrected occurrence

## The Single Most Important Modern Insight -- Read Counts Are NOT Abundance

Elbrecht & Leese 2015 *PLoS One* 10:e0130324 and Lamb et al. 2019 *Mol Ecol* 28:420-430 (meta-analysis) established that metabarcoding read counts have weak-to-moderate, taxon-specific, NONLINEAR correlation with biomass or DNA input. Primer-binding bias dominates; PCR replicates introduce stochasticity. **Reporting read counts as abundance without mock-community calibration is malpractice.** Modern practice: report PRESENCE/ABSENCE or relative abundance with explicit calibration; use multiple PCR replicates; apply site-occupancy models for detection correction.

A second cornerstone: the ASV-vs-OTU debate is taxon-specific, not universal. Callahan, McMurdie, Holmes 2017 *ISME J* 11:2639-2643 argued ASVs replace OTUs because modern denoising resolves single-nucleotide differences. Schloss 2021 *mSphere* 6:e00191-21 showed that for bacterial 16S with 1-15 intra-genomic rRNA copies, a single E. coli strain produces ~7 distinct ASVs, splitting bacterial genomes across artificial clusters. **For COI metazoan metabarcoding, ASVs (DADA2/UNOISE3) are recommended; for bacterial 16S, ASVs inflate alpha-diversity and OTUs may be appropriate.**

A third: **decontam (Davis 2018) is a SCREENING tool, not a deterministic classifier.** It flags candidates; biological plausibility check is required before deletion. The default `threshold=0.1` over-flags in low-biomass data.

## Algorithmic Taxonomy

| Method | Output | Strength | Fails when |
|--------|--------|----------|------------|
| DADA2 | Single-nucleotide ASVs | High resolution; learned error model; standard for COI/12S/18S/fungal-ITS | Small datasets (< 100 samples) for error learning; multi-copy bacterial rRNA |
| UNOISE3 (USEARCH/VSEARCH; Edgar 2016) | zOTUs (essentially ASVs) | Fast; algorithmic simplicity | Limited Linux/Mac binary distribution under license |
| Swarm v2 `-d 1 --fastidious` (Mahé 2015) | Abundance-weighted single-linkage OTUs | Modern OTU pipeline; better than legacy 97% UCLUST | OTUs by design (not single-nt resolution) |
| 97% UCLUST | Classical OTUs | Legacy familiarity | Biologically arbitrary threshold; supersedes by DADA2/swarm |
| VSEARCH global pairwise | Taxonomic assignment via best-hit | Fast, transparent, no training | Conservative; mis-assigns sister species when ref incomplete |
| Naive Bayes (q2-feature-classifier, RDP) | Probabilistic taxonomic assignment | Probabilistic confidence; standard for 16S | Confidence values are scikit-learn calibrated, not true probabilities (Bokulich 2018) |
| SINTAX (Edgar) | Bootstrap-supported taxonomy | Fast; no training | Less accurate than Naive Bayes for divergent sequences |
| LCA (BASTA, MEGAN-LCA) | Lowest common ancestor of multiple hits | Conservative; never over-confident | Can over-merge to high taxonomic ranks |
| Phylogenetic placement (EPA-ng + gappa) | Position on reference tree | Most rigorous; phylogenetically explicit | 10-100x slower; emerging not yet standard |
| decontam | Flagged contaminant candidates | Statistical screening of negative controls and DNA concentration patterns | Output is screening, not classification; needs biological-plausibility check |
| UCHIME3 (in DADA2/VSEARCH) | Chimera detection | Standard for de novo chimera removal | Some divergent chimeras escape |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Metazoan COI metabarcoding (water, gut content) | mlCOIintF/jgHCO2198 (Leray 2013) primers; DADA2 ASVs | Standard primer set; ASVs preserve single-nt resolution |
| Fish eDNA from water | MiFish-U/E (Miya 2015) 12S primers; DADA2 ASVs | Dominant eDNA fish marker globally |
| Freshwater macroinvertebrate bioassessment | BF1/BR1 freshwater-optimized COI primers | Higher primer-binding inclusivity for aquatic insects |
| Bacterial community 16S | 515F/806R (V4) Parada modified; ASVs OR Swarm v2 | Schloss 2021 caveat applies; ASVs may oversplit multi-copy rRNA |
| Fungal community ITS | ITS2 primers; DADA2 or UNITE pipeline | UNITE is curated for fungal ITS |
| Plant community DNA | trnL P6 loop (Taberlet 2007) for degraded DNA | Robust to degradation |
| Deciding ASV vs OTU | ASVs for COI/12S/18S/fungi; OTU consideration for 16S with multi-copy concern | Taxon-specific |
| NovaSeq library (patterned flow cell) | Heavier tag-jumping correction; expect 10x higher rates than MiSeq | Patterned-cell index hopping |
| Low-biomass eDNA (deep ocean, ancient) | decontam frequency + prevalence methods; explicit reagent-contamination check | Reagent contamination dominates |
| Quantitative comparison across samples | Mock-community calibration BEFORE reporting read counts | Without mock, read counts are biased estimators of biomass |
| Detection probability with replication | Site-occupancy models (occumb, eDNAoccupancy; Ficetola 2015) | Read counts alone underestimate occurrence; replicates correct |
| Taxonomic assignment for marker > 80% covered | Naive Bayes (q2-feature-classifier) | Probabilistic; well-supported |
| Taxonomic assignment for sparse reference | Phylogenetic placement (EPA-ng) | Robust to incomplete references |
| OBITools3 pipeline | `obi stats` (NOTE: plural), DMS-based, `.tar.gz` taxonomy | v3 syntax differs from v1 |

## Primer Trimming with cutadapt

**Goal:** Remove primer sequences while discarding reads that lack primers, before quality filtering.

**Approach:** Use cutadapt linked-adapter mode with marker-specific 5' and 3' primer pairs. `--discard-untrimmed` removes reads lacking expected primers; `min_overlap` prevents false primer detection in random sequence regions.

```bash
# COI metazoan (Leray mlCOIintF / jgHCO2198 -> 313 bp)
cutadapt -g 'GGWACWGGWTGAACWGTWTAYCCYCC;min_overlap=20' \
         -G 'TAIACYTCIGGRTGICCRAARAAYCA;min_overlap=20' \
         --discard-untrimmed --pair-filter=any \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz

# Fish 12S (MiFish-U -> 163-185 bp)
cutadapt -g 'GTCGGTAAAACTCGTGCCAGC;min_overlap=18' \
         -G 'CATAGTGGGGTATCTAATCCCAGTTTG;min_overlap=18' \
         --discard-untrimmed --pair-filter=any \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz

# Fungal ITS2
cutadapt -g 'GTGAATCATCGAATCTTTGAAC;min_overlap=18' \
         -G 'TCCTCCGCTTATTGATATGC;min_overlap=18' \
         --discard-untrimmed --pair-filter=any \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz
```

## DADA2 ASV Pipeline

**Goal:** Denoise paired-end amplicon reads into exact amplicon sequence variants (ASVs) with chimera removal and reference-based taxonomy assignment, per Callahan et al. 2016 *Nat Methods* 13:581-583.

**Approach:** Filter to length/quality thresholds, learn error rates per dataset, run dada() to denoise, merge pairs, build sequence table, remove chimeras with UCHIME3-equivalent in DADA2, then assign taxonomy against the marker-appropriate reference DB. CRITICAL: primers must be removed (cutadapt) BEFORE filterAndTrim, OR the error model is corrupted.

```r
library(dada2)

# CRITICAL: primer removal MUST precede filterAndTrim
# DADA2's error model assumes primer-free reads
fwd_reads <- sort(list.files('primer_trimmed/', pattern = '_R1', full.names = TRUE))
rev_reads <- sort(list.files('primer_trimmed/', pattern = '_R2', full.names = TRUE))
filt_fwd <- file.path('filtered', basename(fwd_reads))
filt_rev <- file.path('filtered', basename(rev_reads))

# Filter and trim
# maxEE=c(2,2): expected errors per read; tradeoff sensitivity/specificity
# truncLen: set from quality profile inspection; do not guess
out <- filterAndTrim(fwd_reads, filt_fwd, rev_reads, filt_rev,
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     truncLen = c(220, 180),     # data-dependent; inspect plotQualityProfile()
                     minLen = 100, rm.phix = TRUE, multithread = TRUE)

# Learn error rates
# For small datasets (< 100 samples), pool aggressively or use pre-learned model
err_fwd <- learnErrors(filt_fwd, multithread = TRUE)
err_rev <- learnErrors(filt_rev, multithread = TRUE)

# Denoise
dada_fwd <- dada(filt_fwd, err = err_fwd, multithread = TRUE)
dada_rev <- dada(filt_rev, err = err_rev, multithread = TRUE)

# Merge pairs with minimum overlap
merged <- mergePairs(dada_fwd, filt_fwd, dada_rev, filt_rev, minOverlap = 12)

# Build sequence table
seqtab <- makeSequenceTable(merged)

# Remove chimeras
# method='consensus': per-sample then consensus; conservative (default)
# method='pooled': pooled across samples; aggressive; can over-merge real diversity
# Chimera rate >30% typically indicates library prep problems
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus',
                                     multithread = TRUE)
cat('Chimera rate:', round(1 - sum(seqtab_nochim) / sum(seqtab), 3), '\n')

# Taxonomy assignment
# minBoot=80: standard genus-level confidence; 50 for family-level
# IMPORTANT: pair the marker with the appropriate reference DB
# COI -> MIDORI2 LONGEST_NUC_GB259_CO1 (or BOLD with curation)
# 12S -> MitoFish (Miya lab)
# 16S V4 -> SILVA 138.1+
# 18S V4/V9 -> SILVA 138.1+ or PR2
# Fungal ITS -> UNITE 9.0+
taxa <- assignTaxonomy(seqtab_nochim,
                       'MIDORI2_LONGEST_NUC_GB259_CO1_DADA2.fasta.gz',
                       minBoot = 80, multithread = TRUE)
```

## OBITools3 Pipeline — The v1 -> v3 Command Break

**Goal:** Process eDNA reads through the Unix-style OBITools v3 pipeline (Boyer et al. 2016 *Mol Ecol Resour* 16:176-182 introduced OBITools v1; v3 is the post-2018 Python 3 rewrite) with DMS-based sequence management.

**Approach:** v3 introduces a Database Management System (DMS) abstraction; sequences are imported into a DMS rather than read directly from FASTQ. Commands use spaces (e.g., `obi stats` plural, not `obistat`). Taxonomy import expects `.tar.gz` archive, not a directory.

```bash
# v1 -> v3 command-name changes (critical):
# v1: obistat       -> v3: obi stats
# v1: obigrep       -> v3: obi grep
# v1: obiuniq       -> v3: obi uniq
# v1: obitab        -> v3: obi annotate / obi export --tab-output (different semantics)
# v1: ngsfilter     -> v3: obi ngsfilter
# v1: taxdump dir   -> v3: .tar.gz archive

# Import paired FASTQ into DMS
obi import --fastq-input raw_R1.fastq.gz EDNA/reads1
obi import --fastq-input raw_R2.fastq.gz EDNA/reads2

# Paired-end alignment
obi alignpairedend -R EDNA/reads2 EDNA/reads1 EDNA/aligned

# Filter by alignment score and length
obi grep -p 'sequence["score"] >= 50' EDNA/aligned EDNA/filtered
obi grep -p 'len(sequence) >= 100 and len(sequence) <= 500' \
    EDNA/filtered EDNA/length_filtered

# Demultiplex (NGS filter file maps barcodes -> samples)
obi ngsfilter -t ngsfilter.txt -u EDNA/unassigned \
    EDNA/length_filtered EDNA/demux

# Dereplicate (obi uniq creates merged_sample attribute automatically)
obi uniq EDNA/demux EDNA/derep

# Remove suspected error singletons
obi grep -p 'sequence["count"] >= 2' EDNA/derep EDNA/no_singletons

# Denoise via obi clean
obi clean -s merged_sample -r 0.05 -H EDNA/no_singletons EDNA/denoised

# Taxonomy assignment against reference database
obi ecotag -R EDNA/refdb --taxonomy EDNA/taxonomy EDNA/denoised EDNA/assigned

# Export tab-separated species table
obi export --tab-output EDNA/assigned > species_table.tsv
```

Across most metabarcoding studies, 50-85% of ASVs cannot be assigned to species level due to incomplete references (Wangensteen et al. 2018 *PeerJ* 6:e4705 documented this for marine COI + 18S). Report this gap honestly; do not infer ecology from "unassigned" reads.

## Tag-Jumping Mitigation — Schnell 2015 + NovaSeq Caveat

**Goal:** Detect and remove sequence-to-sample misassignments arising from chimeric library molecules with mismatched indices.

**Approach:** Use dual-indexing (different indices at both ends; cross-jumped pairs are discarded). Quantify residual tag-jumping rate from per-ASV cross-sample appearance and apply per-ASV abundance threshold filtering with `metabaR::tagjumpslayer`. For NovaSeq libraries, expect ~10x higher tag-jumping than MiSeq due to patterned flow cells.

```r
library(metabaR)

# metabaR expects an metabarlist object (asv table + sample info + ngsfilter)
# tagjumpslayer applies per-ASV abundance-threshold filter
# threshold: 0.01 (1% of ASV total) is conservative; 0.001 for aggressive removal
# Adjust threshold higher for NovaSeq (~0.005-0.01) than MiSeq (~0.001-0.005)

# Quantify residual tag-jumping rate before filtering:
# Count reads in sample x ASV combinations that should be 0 by experimental design
# (e.g., samples explicitly excluded from a particular condition)
# That rate / total reads = empirical tag-jumping rate
# Report this rate in methods section
```

## Contamination Screening with decontam

**Goal:** Identify candidate contaminant ASVs from negative controls and DNA-concentration patterns.

**Approach:** Use `decontam::isContaminant` with `method='combined'` when both DNA concentration AND negative controls are available. Treat flagged ASVs as SCREENING CANDIDATES; verify biological plausibility before deletion. The default `threshold=0.1` is over-aggressive in low-biomass data.

```r
library(decontam)

# Frequency method: contaminants more frequent at LOW DNA concentration
# Prevalence method: contaminants more frequent in negative controls
# Combined: uses both signals (most robust)

contam <- isContaminant(seqtab_nochim,
                        conc = dna_concentration,        # qPCR or Qubit per sample
                        neg = is_negative_control,       # logical: which samples are controls
                        method = 'combined',
                        threshold = 0.1)                 # default; lower for high-confidence calls

# CRITICAL: decontam output is SCREENING, not classification
# Manually inspect each flagged ASV: is the taxonomic assignment plausibly a reagent contaminant?
# Common reagent contaminants: Delftia, Sphingomonas, Burkholderia, Propionibacterium
flagged <- which(contam$contaminant)
cat('Decontam flagged', length(flagged), 'ASVs as candidates\n')

# After manual review, remove confirmed contaminants
confirmed_contam <- intersect(flagged, biological_plausibility_check_result)
seqtab_clean <- seqtab_nochim[, !(colnames(seqtab_nochim) %in% confirmed_contam)]
```

## Site-Occupancy Modeling — Correcting for Imperfect Detection

**Goal:** Estimate true species occurrence probabilities from replicated eDNA samples, accounting for false negatives in any single PCR replicate.

**Approach:** Fit a multi-species occupancy model via MCMC (Ficetola 2015 *Mol Ecol Resour* 15:543-556) on a 3D array of replicated read counts. Output: per-site, per-species occupancy probabilities corrected for detection.

```r
library(occumb)

# y: 3D array [species, sites, replicates] of read counts
# spec_cov: species covariates (traits)
# site_cov: site covariates (env)
data_obj <- occumbData(y = count_array, spec_cov = species_covariates,
                       site_cov = site_covariates)

# Fit hierarchical occupancy model
# Requires JAGS installation
# n.iter >= 10000, n.burn >= 2500 for publication-quality posteriors
fit <- occumb(data = data_obj, n.chains = 4, n.iter = 10000,
              n.thin = 5, n.burn = 2500)

# Extract detection-corrected occupancy
summary(fit)
```

## Per-Method Failure Modes

### Reporting read counts as biomass without mock-community calibration

**Trigger:** Comparing read counts of two ASVs and reporting the ratio as a biomass / abundance estimate.

**Mechanism:** Primer-template binding affinity varies systematically across taxa; PCR amplification is non-linear (saturates); read counts have weak-to-moderate, NONLINEAR correlation with biomass (Elbrecht 2015; Lamb 2019).

**Symptom:** Reviewer asks "how is it known that reads = biomass?"; cross-study quantitative comparisons fail to replicate.

**Fix:** Either (a) restrict reporting to presence/absence; (b) report read counts as relative abundances with explicit caveat; or (c) include mock-community of known composition for primer-specific calibration. Do not silently equate reads with biomass.

### NovaSeq tag-jumping with MiSeq-tuned filtering

**Trigger:** Applying tag-jumping filters calibrated on MiSeq libraries to NovaSeq data.

**Mechanism:** NovaSeq patterned flow cells have ~10x higher index hopping than MiSeq. MiSeq-calibrated thresholds (often ~0.001 fraction) are too permissive on NovaSeq data.

**Symptom:** Apparent rare-species detections in NovaSeq libraries do not replicate; per-ASV cross-sample appearance is unusually broad.

**Fix:** Use NovaSeq-appropriate tag-jumping thresholds (~0.005-0.01) and report the empirical tag-jumping rate from explicit-zero combinations.

### decontam threshold over-aggressive in low-biomass data

**Trigger:** Applying default `threshold = 0.1` to ASVs from open-ocean water, ancient sediments, or other dilute samples.

**Mechanism:** In low-biomass samples, the contaminant signal/background ratio approaches 1; decontam over-flags real but dilute biology as "contaminant" because the statistical pattern looks similar.

**Symptom:** Many ASVs flagged from low-biomass samples; taxonomic profile of "flagged contaminants" looks biologically realistic.

**Fix:** Lower threshold (0.05 or 0.01); always manually review flagged ASVs for biological plausibility; cite Salter 2014 *BMC Biol* 12:87 for the low-biomass reagent-contamination caveat.

### Skipping primer removal before DADA2 filterAndTrim

**Trigger:** Running DADA2's `filterAndTrim()` on FASTQ files that still contain primer sequences.

**Mechanism:** DADA2 learns sequencing error from the empirical data; if primer sequences are present, they look like "perfect agreement" and corrupt the error model. ASVs are inferred with primer artifacts attached.

**Symptom:** DADA2 reports "phix-like contamination" (false; it's primers); ASVs start with the primer sequence; chimera rate elevated.

**Fix:** Always run cutadapt (or similar) BEFORE filterAndTrim. Verify with `head` of trimmed FASTQ that primer sequences are gone.

### OBITools v3 commands with v1 syntax

**Trigger:** Running `obistat` or `obigrep` on a v3 install.

**Mechanism:** v1 used concatenated command names (`obistat`); v3 uses subcommand syntax with a space (`obi stats` — note plural).

**Symptom:** Bash error `obistat: command not found`; tutorial documentation does not match installed version.

**Fix:** Use `obi <subcommand>` syntax; consult `obi --help` for current command list. Taxonomy import requires `.tar.gz` archive, not unpacked directory.

## Quantitative Thresholds

| Threshold | Value | Source / rationale |
|-----------|-------|-------------------|
| DADA2 maxEE per read | 2 | Standard sensitivity/specificity balance |
| DADA2 chimera rate alarm | > 30% suggests library issues | Empirical convention |
| DADA2 minBoot for taxonomy | 80 for genus; 50 for family | Standard confidence cutoffs |
| Tag-jumping filter MiSeq | 0.001-0.005 fraction of ASV total | Schnell 2015 |
| Tag-jumping filter NovaSeq | 0.005-0.01 fraction of ASV total | Patterned-cell index hopping ~10x higher |
| decontam threshold | 0.1 default; 0.05 for low-biomass | Davis 2018; reduce for dilute samples |
| Per-sample minimum reads | 1000 (after filtering) | Below this rare-species detection unreliable |
| Singleton removal | count >= 2 | Singletons often error-driven |
| Bootstrap nperm for tests | 999 | Standard permutation count |
| Occupancy model iterations | n.iter >= 10000, n.burn >= 2500 | occumb default for stable posteriors |
| eDNA decay (20 deg C surface water) | half-life ~4-15 hours | Strickler 2015 Biol Conserv 183:85-92 |

## Common errors

| Error | Cause | Solution |
|-------|-------|----------|
| `obistat: command not found` | OBITools v3 uses `obi stats` (plural) | Use v3 syntax |
| DADA2 error rate plot looks pathological | Primer sequences still in reads | Re-run cutadapt before filterAndTrim |
| Chimera rate > 30% | Library-prep issue or primer dimers | Inspect raw FASTQ; check PCR conditions |
| decontam flags many real species | Default threshold too aggressive for low-biomass | Lower threshold; manual review |
| Naive Bayes confidence 0.95 but species is wrong | scikit-learn-calibrated "confidence" not true probability | Use phylogenetic placement for borderline assignments |
| occumb JAGS not found error | JAGS not installed system-wide | Install JAGS (CRAN page has platform instructions) |
| eDNA detections do not replicate | Read counts treated as abundance | Switch to presence/absence; use mock-community calibration |
| MIDORI2 download path expired | Database updated; old URL gone | Check current MIDORI2 / MitoFish download page |

## References

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016) DADA2. *Nat Methods* 13(7):581-583. doi:10.1038/nmeth.3869
- Callahan BJ, McMurdie PJ, Holmes SP (2017) Exact sequence variants should replace OTUs. *ISME J* 11(12):2639-2643. doi:10.1038/ismej.2017.119
- Edgar RC (2016) UNOISE2 / UNOISE3. *bioRxiv* preprint. doi:10.1101/081257
- Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015) Swarm v2. *PeerJ* 3:e1420. doi:10.7717/peerj.1420
- Leray M, Yang JY, Meyer CP et al. (2013) mlCOIintF/jgHCO2198 metazoan COI primer. *Front Zool* 10:34. doi:10.1186/1742-9994-10-34
- Miya M, Sato Y, Fukunaga T et al. (2015) MiFish 12S fish eDNA primer. *R Soc Open Sci* 2(7):150088. doi:10.1098/rsos.150088
- Schnell IB, Bohmann K, Gilbert MTP (2015) Tag jumps illuminated. *Mol Ecol Resour* 15(6):1289-1303. doi:10.1111/1755-0998.12402
- Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ (2018) decontam. *Microbiome* 6:226. doi:10.1186/s40168-018-0605-2
- Boyer F, Mercier C, Bonin A, Le Bras Y, Taberlet P, Coissac E (2016) OBITools. *Mol Ecol Resour* 16(1):176-182. doi:10.1111/1755-0998.12428
- Elbrecht V, Leese F (2015) DNA-based ecosystem quantification critique. *PLoS One* 10(7):e0130324. doi:10.1371/journal.pone.0130324
- Lamb PD, Hunter E, Pinnegar JK, Creer S, Davies RG, Taylor MI (2019) How quantitative is metabarcoding: meta-analysis. *Mol Ecol* 28(2):420-430. doi:10.1111/mec.14920
- Ficetola GF, Pansu J, Bonin A et al. (2015) Replication levels and false presences in eDNA. *Mol Ecol Resour* 15(3):543-556. doi:10.1111/1755-0998.12338
- Wangensteen OS, Palacin C, Guardiola M, Turon X (2018) COI + 18S marine metabarcoding. *PeerJ* 6:e4705. doi:10.7717/peerj.4705
- Strickler KM, Fremier AK, Goldberg CS (2015) eDNA degradation kinetics. *Biol Conserv* 183:85-92. doi:10.1016/j.biocon.2014.11.038
- Bokulich NA, Kaehler BD, Rideout JR et al. (2018) Optimizing taxonomic classification with q2-feature-classifier. *Microbiome* 6:90. doi:10.1186/s40168-018-0470-z
- Bohmann K, Evans A, Gilbert MTP et al. (2014) eDNA for wildlife and biodiversity. *Trends Ecol Evol* 29(6):358-367. doi:10.1016/j.tree.2014.04.003
- Schloss PD (2021) Amplicon sequence variants artificially split bacterial genomes into separate clusters. *mSphere* 6(4):e00191-21. doi:10.1128/mSphere.00191-21
- Salter SJ, Cox MJ, Turek EM et al. (2014) Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biol* 12:87. doi:10.1186/s12915-014-0087-z
- Taberlet P, Coissac E, Pompanon F et al. (2007) Power and limitations of the chloroplast trnL (UAA) intron for plant DNA barcoding. *Nucleic Acids Res* 35(3):e14. doi:10.1093/nar/gkl938

## Related Skills

- ecological-genomics/biodiversity-metrics - Diversity analysis from species occurrence tables (Hill numbers, beta partition)
- ecological-genomics/community-ecology - Environmental gradient analysis of community composition (PERMANOVA + PERMDISP, ordination)
- microbiome/amplicon-processing - 16S clinical microbiome alternative pipeline
- read-qc/quality-reports - Upstream read-quality assessment before primer trimming
- database-access/entrez-fetch - Retrieve reference sequences for custom taxonomy databases
