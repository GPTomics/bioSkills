---
name: bio-microbiome-taxonomy-assignment
description: Assigns taxonomy to amplicon ASVs/OTUs (16S, ITS, 18S) with a classifier conditioned on a reference database and primer region - DADA2 assignTaxonomy + addSpecies (RDP naive Bayes), DECIPHER IDTAXA, and QIIME2 q2-feature-classifier (classify-sklearn naive Bayes, classify-consensus-vsearch alignment-consensus). Covers region-specific training (extract-reads, fit-classifier-naive-bayes), why a full-length classifier fabricates calls on a V4 read, the scikit-learn version-pinning trap on pre-trained .qza classifiers, confidence thresholds (classify-sklearn 0.7, assignTaxonomy minBoot 50), and choosing SILVA/GTDB/Greengenes2/UNITE/PR2/RDP. Use when classifying ASVs after DADA2, picking a reference database, training a region-matched classifier, setting a confidence threshold, or deciding whether a 16S species call is defensible (usually not - genus at best). For shotgun read classification see metagenomics/kraken-classification and metagenomics/metaphlan-profiling.
tool_type: mixed
primary_tool: DADA2
---

## Version Compatibility

Reference examples tested with: DADA2 1.30+, DECIPHER 2.30+, QIIME2 2024.10+ (q2-feature-classifier), vsearch 2.22+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The REFERENCE DATABASE and the CLASSIFIER TRAINING REGION are the versions that matter most. Results track the database release and format - SILVA 138.1/138.2, GTDB r220, Greengenes2 2024.09, UNITE (its own release scheme for ITS), PR2 5.x. A DADA2-formatted training FASTA, a DECIPHER trainingSet `.RData`, and a QIIME2 pre-trained classifier `.qza` are NOT interchangeable: each must match the marker and the primer region. A pre-trained naive-Bayes `.qza` is a pickled scikit-learn model tied to the QIIME2 release that built it (see the version-pinning failure mode). Record the database release AND the training region (full-length vs 515-806/V4) for every label.

# Taxonomy Assignment

**"What is this ASV?"** -> Return the most probable lineage in a reference database at a stated confidence - because a label is a classifier + database + primer-region-conditioned hypothesis, not an identification, and a short 16S read licenses genus at best.
- R: `assignTaxonomy(seqs, refFasta, minBoot=50, multithread=TRUE)` then `addSpecies(taxa, speciesFasta)` (exact match only)
- CLI: `qiime feature-classifier classify-sklearn --i-classifier region-matched.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza`

Scope: classification of per-feature amplicon sequences (ASVs/OTUs) from 16S/ITS/18S. ASV inference is upstream -> amplicon-processing. Diversity/DA on the classified table -> diversity-analysis, differential-abundance. Shotgun read classification (raw reads, not ASVs) -> metagenomics/kraken-classification (k-mer LCA), metagenomics/metaphlan-profiling (clade markers). Primer removal before assignment -> read-qc/adapter-trimming.

## The Single Most Important Modern Insight -- A Taxonomic Label Is a Classifier + Database + Region-Conditioned Hypothesis, Not an Identification

A classifier never identifies an organism; it returns the most probable lineage in THIS database, under a model that assumes the true source taxon is represented. The label is a function of three choices made before seeing the answer - the classifier, the reference database (and its naming authority), and the primer region the reference was trimmed to. Change any one and the genus or species call can change. Three corollaries each common misuse violates:

1. **Confidence is certainty WITHIN the database, not proof the answer is in it.** A bootstrap of 95 means the model is sure GIVEN the taxon is in the reference. If the true organism is absent, the classifier confidently returns the nearest WRONG relative. The score cannot detect absence.
2. **Assigning a label and the label being correct are different events.** A ~250 bp single hypervariable region (V4) is identical across many species and often across genera - the discriminating substitutions are elsewhere in the gene. The tool will still EMIT a species name (over-classification); that name is not licensed by the data. Report the rank the data supports - genus at best for many bacteria, family for poorly resolved clades.
3. **The reference must match the primer region.** A full-length-trained classifier applied to a V4 read mismatches k-mer composition and both fabricates and erases calls (Werner 2012; Bokulich 2018). Use a region-matched (515-806) classifier or extract-reads -> fit-classifier-naive-bayes.

ITS (the fungal barcode) is the exception: it resolves to species far more reliably than 16S. 18S (eukaryotes) resolves coarsely, like 16S.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| classify-sklearn (q2-feature-classifier) | Bokulich 2018 *Microbiome* 6:90 | multinomial naive Bayes over 7-mers; pre-trained pickled model | QIIME2 default; train once, classify many; fast |
| classify-consensus-vsearch | Bokulich 2018 *Microbiome* 6:90 | global alignment, consensus over top hits; no trained model | immune to sklearn pinning; transparent hits; slower |
| assignTaxonomy + addSpecies (DADA2) | Wang 2007 *Appl Environ Microbiol* 73:5261 | RDP 8-mer naive Bayes + bootstrap; addSpecies = exact 100% match | R workflow; genus via NB, species via exact match only |
| IDTAXA (DECIPHER) | Murali 2018 *Microbiome* 6:140 | tree-descent with learned per-node confidence; refuses to over-descend | conservative; minimizes over-classification; novelty-aware |
| weighted/clawback classifier | Kaehler 2019 *Nat Commun* 10:4643 | naive Bayes with habitat-specific abundance prior | known habitat (gut, soil); raises species accuracy |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| 16S V4 (515F/806R), standard | SILVA region-matched (515-806) pre-trained NB classifier OR DADA2 assignTaxonomy on SILVA | region-matched, fast, well-benchmarked |
| Have primers but no pre-trained artifact for them | extract-reads -> fit-classifier-naive-bayes (RESCRIPt) | builds the matched reference; full-length underperforms |
| Known habitat, best species accuracy | weighted/clawback classifier | habitat prior cuts species error (Kaehler 2019) |
| scikit-learn version error / no retraining wanted | classify-consensus-vsearch | stores no pickled model; immune to version pinning |
| Conservative, novelty-aware, minimize over-calls | IDTAXA (DECIPHER) | learned per-node refusal to over-descend |
| Unify 16S with shotgun on one tree | Greengenes2 (q2-greengenes2) | single genome-backbone tree, GTDB-harmonized; closed-reference |
| Environmental / under-named bacteria | GTDB SSU reference (via RESCRIPt) | genome-based, rank-normalized, polyphyly-pruned |
| Fungal ITS | UNITE (species hypotheses) + NB or vsearch | formal fungal barcode; species-resolved; never position-trim ITS |
| Protist / eukaryote 18S | PR2 (or SILVA 18S) | curated microeukaryote SSU |
| Species name from one 16S region | DO NOT (or addSpecies exact-match only) | region lacks the information; report genus |
| Raw shotgun reads, not ASVs | -> metagenomics/kraken-classification, metaphlan-profiling | different input artifact and output semantics |

Method choice is contested (see References). Bokulich 2018 found naive Bayes and vsearch-consensus broadly comparable and both top-tier; verify current best practice against the latest q2-feature-classifier docs rather than hard-coding one classifier.

## DADA2 assignTaxonomy + addSpecies

**Goal:** Assign each ASV to the deepest rank the data and reference support, reporting NA below the confidence floor rather than forcing a label.

**Approach:** Run the RDP naive Bayes (8-mer, 100 bootstraps) against a DADA2-formatted reference to a genus-level call, leaving ranks NA below `minBoot`; then attempt species ONLY by exact match against a species reference - never inferring species from the noisy read.

```r
library(dada2)
seqtab_nochim <- readRDS('seqtab_nochim.rds')

# minBoot 50 = the DADA2 default and the RDP recommendation for reads <=250 nt; tutorials
# often use 80 (a stricter CHOICE, not the default). Raising it truncates to shallower but
# more reliable ranks; ranks below the threshold are returned as NA, not guessed.
taxa <- assignTaxonomy(seqtab_nochim, 'silva_nr99_v138.1_train_set.fa.gz', minBoot = 50, tryRC = TRUE, multithread = TRUE)

# addSpecies assigns species by EXACT (100%) match against a species reference. It does NOT
# license a species name from a noisy read - it reports a species only when the ASV is
# identical to a reference over the amplicon, else leaves it NA. This is the honest 16S path.
taxa <- addSpecies(taxa, 'silva_species_assignment_v138.1.fa.gz')
```

The reference FASTA must be DADA2-formatted (rank-labelled headers) AND ideally trimmed to the amplicon region; a full-length SILVA training set on V4 reads is the Trap-1 failure mode below. For ITS, use a UNITE DADA2 reference and do NOT trim ITS to a fixed length (it is variable-length).

## DECIPHER IDTAXA

**Goal:** Get a conservative, novelty-aware classification that refuses to descend into a clade the query likely does not belong to.

**Approach:** Convert ASV sequences to a DNAStringSet, classify with a pre-trained DECIPHER trainingSet, then flatten the per-rank output to a matrix, mapping IDTAXA's "unclassified_" placeholders to NA.

```r
library(DECIPHER)
load('SILVA_SSU_r138_2019.RData')  # provides the trainingSet object
dna <- DNAStringSet(getSequences(seqtab_nochim))

# threshold 60 = DECIPHER default confidence cutoff; raise for stricter calls. IDTAXA's
# tree-descent stops (leaves the rank unclassified) when the query likely belongs to a taxon
# absent from the reference - this is the intended anti-over-classification behaviour.
ids <- IdTaxa(dna, trainingSet, strand = 'both', threshold = 60, processors = NULL)

ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
taxa_idtaxa <- t(sapply(ids, function(x) {
    out <- x$taxon[match(ranks, x$rank)]
    out[startsWith(replace(out, is.na(out), ''), 'unclassified_')] <- NA
    out
}))
colnames(taxa_idtaxa) <- ranks
```

## QIIME2 classify-sklearn + Region-Specific Training

**Goal:** Classify ASVs with a naive-Bayes classifier whose reference is trimmed to the EXACT primer region, avoiding the full-length-on-V4 accuracy loss.

**Approach:** In-silico PCR the reference to the primer-bounded region with extract-reads, train a naive-Bayes classifier on the extracted reads, then classify at the default confidence (0.7), which truncates each lineage to the deepest rank clearing the threshold.

```bash
# 1. Extract the V4 (515F/806R) region from a full-length reference (matched k-mer composition)
qiime feature-classifier extract-reads \
    --i-sequences silva-138-99-seqs.qza \
    --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT \
    --p-min-length 50 --p-max-length 0 \
    --o-reads ref-seqs-515-806.qza

# 2. Train the naive-Bayes classifier on the EXTRACTED region (or download the region-matched
#    pre-trained .qza built for THIS QIIME2 release - never a different release, see below)
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ref-seqs-515-806.qza \
    --i-reference-taxonomy silva-138-99-tax.qza \
    --o-classifier silva-138-99-515-806-nb-classifier.qza

# 3. Classify. --p-confidence default 0.7: below it the lineage is truncated to a shallower,
#    more confident rank. 0 = compute but never truncate (deepest always); 'disable' = skip.
qiime feature-classifier classify-sklearn \
    --i-classifier silva-138-99-515-806-nb-classifier.qza \
    --i-reads rep-seqs.qza \
    --p-confidence 0.7 --p-read-orientation auto --p-n-jobs 1 \
    --o-classification taxonomy.qza
```

`--p-n-jobs >1` multiplies memory (each job holds a copy of the classifier); a full-length SILVA classifier is multi-GB, so reduce `--p-n-jobs` / `--p-reads-per-batch` if OOM-killed, or use the smaller region-extracted classifier.

## QIIME2 classify-consensus-vsearch (no pickled model)

**Goal:** Classify without a trained, sklearn-pinned model - by aligning each query to the reference and taking a consensus taxonomy across top hits.

**Approach:** VSEARCH global alignment keeps the top `maxaccepts` hits above `perc-identity`; a rank is reported only if at least `min-consensus` of those hits agree on it, else it is dropped.

```bash
qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seqs.qza \
    --i-reference-reads silva-138-99-seqs.qza \
    --i-reference-taxonomy silva-138-99-tax.qza \
    --p-maxaccepts 10 --p-perc-identity 0.8 --p-min-consensus 0.51 \
    --p-threads 8 \
    --o-classification taxonomy-vsearch.qza --o-search-results vsearch-hits.qza
```

## Filtering Host Organelle and Off-Target Features

**Goal:** Remove host mitochondrial 16S, chloroplast/plastid 16S, and domain-unassigned features BEFORE diversity and differential abundance - universal 16S primers amplify host organelle rRNA, and leaving it in inflates the feature table and deflates every real taxon by compositional closure.

**Approach:** Use the taxonomy just assigned to exclude the Mitochondria and Chloroplast lineages (the labels enable the filter), then carry the filtered table and sequences forward.

```bash
# QIIME2: exclude mitochondria + chloroplast (case-insensitive substring match on the lineage)
qiime taxa filter-table --i-table table.qza --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table table-no-organelle.qza
qiime taxa filter-seqs --i-data rep-seqs.qza --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-data rep-seqs-no-organelle.qza
```

```r
# phyloseq (SILVA ranks: Order 'Chloroplast', Family 'Mitochondria'); the is.na guard keeps unranked taxa
ps <- subset_taxa(ps, is.na(Order)  | Order  != 'Chloroplast')
ps <- subset_taxa(ps, is.na(Family) | Family != 'Mitochondria')
```

Organelle contamination is heaviest in plant, rhizosphere, and host-tissue/biopsy samples (often the majority of reads); inspect per-sample read retention after filtering. Exception: in phototroph-focused, aquatic, or microbial-mat communities, Chloroplast-binned 16S can be the signal of interest (cyanobacterial vs algal-plastid 16S are hard to separate) - inspect what falls in the Chloroplast bin before excluding it. The phyloseq rank-equality form is SILVA-138-specific (Chloroplast at Order, Mitochondria at Family); for GTDB/Greengenes2/RDP verify the rank (`get_taxa_unique(ps, 'Order')`) or use the QIIME2 substring exclude, which is reference-robust.

## Per-Method Failure Modes

### Full-length classifier on a sub-region
**Trigger:** classifying V4 (~250 bp) ASVs with a classifier trained on full-length (~1500 bp) 16S "because it was the one downloaded." **Mechanism:** conserved-region k-mers dominate the model; the discriminating signal is outside the amplified window, so k-mer composition mismatches the query. **Symptom:** systematically shallower or wrong genus/species calls vs a region-matched run. **Fix:** use the 515-806 (or primer-matched) classifier, or extract-reads -> fit-classifier-naive-bayes (Werner 2012; Bokulich 2018).

### Over-reading species from 16S
**Trigger:** reporting a species name (e.g. from naive Bayes descending past the region's resolution) as an identification. **Mechanism:** one hypervariable region is identical across many species/genera; the read lacks the information regardless of classifier quality. **Symptom:** confident species labels that an exact-match (addSpecies) call would leave NA. **Fix:** report genus; use addSpecies (exact match) or IDTAXA's conservative descent for any species claim; treat 16S species as a hypothesis. ITS may legitimately reach species.

### Mixed or mismatched database
**Trigger:** comparing a SILVA genus to a GTDB genus, merging two tables built on different references, or using a database whose scope excludes the marker (GTDB has no Eukarya -> cannot classify 18S/ITS). **Mechanism:** NCBI/SILVA/GTDB assign different names, and GTDB normalizes ranks; labels are not comparable across authorities. **Symptom:** spurious genus mismatches across cohorts; empty/garbled ITS calls against a bacteria-only DB. **Fix:** pick one database matched to the marker (16S->SILVA/GTDB/GG2; ITS->UNITE; 18S->PR2/SILVA), state its release, never merge labels across authorities without a crosswalk.

### scikit-learn version break on a pre-trained classifier
**Trigger:** loading a pre-trained NB `.qza` built for a different QIIME2 release. **Mechanism:** the classifier is a pickled scikit-learn model; QIIME2 records the sklearn version and refuses to run under a different one. **Symptom:** the error "The scikit-learn version (X) used to generate this artifact does not match the current version of scikit-learn installed (Y). Please retrain..." **Fix:** download the classifier built for the exact installed QIIME2 release, OR retrain locally with fit-classifier-naive-bayes, OR use classify-consensus-vsearch (no pickled model). Pinning to a hard-coded classifier URL is exactly this brittleness.

### Confidence default left unexamined
**Trigger:** running classify-sklearn at 0.7 or assignTaxonomy at the default, then reporting whatever rank comes out without stating the threshold. **Mechanism:** the threshold trades sensitivity for specificity - lowering it over-classifies (deeper but wronger), raising it truncates to shallower-but-reliable ranks. **Symptom:** either an over-deep label list or silently dropped/force-filled Unassigned features. **Fix:** state the threshold, tune to region/DB, and keep the truncated ("unassigned at rank X") output honestly - do not drop or force-fill it.

### Host organelle reads not filtered
**Trigger:** running diversity/DA on a host-associated or plant sample without removing Mitochondria/Chloroplast features. **Mechanism:** universal 16S primers amplify host mitochondrial and plastid 16S; classifiers label them `f__Mitochondria`/`o__Chloroplast`, and compositional closure then deflates every real taxon. **Symptom:** a large read fraction labelled Mitochondria/Chloroplast; diversity/DA tracks host content. **Fix:** filter them (the Filtering section above) after assignment, before diversity/DA.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| classify-sklearn `--p-confidence` 0.7 | Bokulich 2018 *Microbiome* 6:90 | benchmarked default for 16S/ITS; below it the lineage is truncated to a more confident rank |
| assignTaxonomy `minBoot` 50 (default); 80 common | Wang 2007 *Appl Environ Microbiol* 73:5261 | RDP recommended 50 for reads <=250 nt, 80 generally; below the floor the rank is NA |
| IDTAXA `threshold` 60 (default) | Murali 2018 *Microbiome* 6:140 | per-node confidence cutoff; raise for stricter, novelty-conservative calls |
| classify-consensus-vsearch `--p-min-consensus` 0.51 | QIIME2 docs | a rank is reported only if a majority of accepted hits agree on it |
| classify-consensus-vsearch `--p-perc-identity` 0.8 | QIIME2 docs | minimum query-reference identity for an accepted hit |
| classify-consensus-vsearch `--p-maxaccepts` 10 | QIIME2 docs | top hits kept per query for the consensus vote |
| addSpecies match 100% (exact) | DADA2 docs | species licensed only by an exact amplicon match; high precision, low recall |
| extract-reads `--p-min-length` 50 | QIIME2 docs | drops too-short in-silico amplicons that would mistrain the classifier |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "scikit-learn version ... does not match" | pre-trained `.qza` built under a different QIIME2/sklearn release | retrain locally, download the release-matched classifier, or use classify-consensus-vsearch |
| Mostly genus/family, few species | correct behaviour for short 16S | report the supported rank; use addSpecies/IDTAXA for species; do not lower confidence to force species |
| OOM kill during classify-sklearn | large classifier x `--p-n-jobs` | lower `--p-n-jobs`/`--p-reads-per-batch`; use a region-extracted (smaller) classifier |
| Empty or garbage ITS taxonomy | bacteria-only DB (GTDB) or position-trimmed ITS | use UNITE; remove ITS primers but never fixed-length-trim ITS |
| All Unassigned at domain level | off-target ASVs (host, chimera, primer artifact) or wrong-orientation reads | filter off-target; leave read-orientation on `auto`; document, do not force-fill |
| Large read fraction labelled Mitochondria/Chloroplast | host organelle 16S amplified by universal primers | `qiime taxa filter-table --p-exclude mitochondria,chloroplast` (or phyloseq subset_taxa) before diversity/DA |
| Genus mismatch across cohorts | labels from different databases (SILVA vs GTDB) | use one database+release for all samples |

## References

- Bokulich NA, Kaehler BD, Rideout JR, Dillon M, Bolyen E, Knight R, Huttley GA, Caporaso JG. 2018. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. *Microbiome* 6:90.
- Kaehler BD, Bokulich NA, McDonald D, Knight R, Caporaso JG, Huttley GA. 2019. Species abundance information improves sequence taxonomy classification accuracy. *Nat Commun* 10:4643.
- Werner JJ, Koren O, Hugenholtz P, DeSantis TZ, Walters WA, Caporaso JG, Angenent LT, Knight R, Ley RE. 2012. Impact of training sets on classification of high-throughput bacterial 16S rRNA gene surveys. *ISME J* 6:94-103.
- Wang Q, Garrity GM, Tiedje JM, Cole JR. 2007. Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Appl Environ Microbiol* 73:5261-5267.
- Murali A, Bhargava A, Wright ES. 2018. IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. *Microbiome* 6:140.
- Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glockner FO. 2013. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Res* 41:D590-D596.
- Parks DH, Chuvochina M, Waite DW, Rinke C, Skarshewski A, Chaumeil PA, Hugenholtz P. 2018. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nat Biotechnol* 36:996-1004.
- McDonald D, Jiang Y, Balaban M, et al. 2024. Greengenes2 unifies microbial data in a single reference tree. *Nat Biotechnol* 42:715-718.
- Nilsson RH, Larsson KH, Taylor AFS, et al. 2019. The UNITE database for molecular identification of fungi: handling dark taxa and parallel taxonomic classifications. *Nucleic Acids Res* 47:D259-D264.
- Guillou L, Bachar D, Audic S, et al. 2013. The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. *Nucleic Acids Res* 41:D597-D604.
- Robeson MS 2nd, O'Rourke DR, Kaehler BD, Ziemski M, Dillon MR, Foster JT, Bokulich NA. 2021. RESCRIPt: Reproducible sequence taxonomy reference database management. *PLoS Comput Biol* 17:e1009581.

## Related Skills

- amplicon-processing - Generate the ASV table that is classified here
- diversity-analysis - Alpha/beta diversity of the classified community table
- differential-abundance - Compositional DA on the classified feature table
- qiime2-workflow - The QIIME2 CLI workflow this classification step plugs into
- read-qc/adapter-trimming - cutadapt primer removal before ASV inference and assignment
- metagenomics/kraken-classification - Shotgun (raw-read, not ASV) k-mer classification
- metagenomics/metaphlan-profiling - Shotgun marker-gene profiling; a different input artifact
- phylogenetics/tree-io - Phylogenetic tree for UniFrac / Faith PD on the classified table
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
