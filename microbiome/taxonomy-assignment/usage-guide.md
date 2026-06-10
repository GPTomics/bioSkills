# Taxonomy Assignment - Usage Guide

## Overview

Taxonomy assignment classifies amplicon ASVs/OTUs (16S, ITS, 18S) against a reference database. The output label is a classifier + database + primer-region-conditioned hypothesis at a stated confidence, NOT an identification. The single most important honest message: a short 16S read (one hypervariable region, ~250 bp for V4) licenses GENUS at best - frequently only family for poorly resolved clades - and rarely species. The confidence/bootstrap score measures the model's certainty assuming the true taxon is in the database; it cannot detect when the true taxon is absent, in which case the classifier returns the nearest wrong relative. ITS (the fungal barcode) is the exception and does resolve to species; 18S resolves coarsely like 16S.

Three classifier families are covered: DADA2 `assignTaxonomy` + `addSpecies` (RDP naive Bayes, R), DECIPHER IDTAXA (conservative tree-descent, R), and QIIME2 q2-feature-classifier (`classify-sklearn` naive Bayes and `classify-consensus-vsearch` alignment-consensus, CLI). Shotgun read classification is a different category - see metagenomics/kraken-classification and metagenomics/metaphlan-profiling.

## Prerequisites

```bash
# R classifiers (DADA2 + DECIPHER IDTAXA)
# BiocManager::install(c('dada2', 'DECIPHER'))

# QIIME2 installs as its own conda environment; the release tag defines the plugin API and the
# .qza classifier format (e.g. qiime2-amplicon-2024.10). vsearch ships with the QIIME2 env.
# conda env create -n qiime2-amplicon-2024.10 --file <release env file>
```

Conceptual prerequisites:
- Inputs are per-feature SEQUENCES (a DADA2 `seqtab_nochim`, or a QIIME2 `FeatureData[Sequence]` of representative sequences) - ASV inference happens upstream in amplicon-processing.
- The reference database is matched to the marker: 16S -> SILVA / GTDB / Greengenes2; ITS -> UNITE; 18S -> PR2 or SILVA. GTDB has no Eukarya and cannot classify ITS/18S.
- The reference should be trimmed to the SAME primer region as the reads. A full-length classifier on V4 reads fabricates and erases calls.
- Reference database downloads are large (a full-length SILVA NB classifier is multi-GB).
- A pre-trained naive-Bayes `.qza` is a pickled scikit-learn model tied to the QIIME2 release that built it; a mismatched release errors out.
- Pick ONE database+release for all samples in a study; never merge labels across SILVA/GTDB/Greengenes2.

## Quick Start

Tell your AI agent what you want to do:
- "Assign taxonomy to my 16S V4 ASVs against SILVA and report the rank the data supports"
- "Train a region-matched naive-Bayes classifier for my 515F/806R primers, then classify"
- "Classify my fungal ITS sequences against UNITE"
- "I am getting a scikit-learn version error on a pre-trained classifier - what are my options?"
- "Use IDTAXA for a conservative, novelty-aware classification"
- "Filter host mitochondria and chloroplast features out of my table before diversity and DA"

## Example Prompts

### Region-matched classification
> "I have 16S V4 (515F/806R) ASVs from DADA2. Assign taxonomy with SILVA, but make sure the classifier is matched to the V4 region rather than full-length, and report genus where species is not supported."

### Database selection
> "These are environmental samples with many under-named bacteria. Should I use SILVA or GTDB, and what changes about the names if I switch?"

### Method selection
> "Classify my ASVs two ways - naive Bayes and vsearch alignment-consensus - and tell me where they disagree and why one is immune to the scikit-learn version error."

### Confidence and honest reporting
> "Run assignTaxonomy and tell me what minBoot you used, what it trades, and which ASVs you left unassigned at genus rather than force-filling."

### Species-level question
> "Can I report species for this 16S ASV? If not, what would it take to defend a species call?"

### Fungal / eukaryote markers
> "My amplicon is fungal ITS. Classify against UNITE and explain why ITS can reach species when my 16S could not."

### Filtering host organelle reads
> "My samples are plant-associated and a big fraction of ASVs are labelled Chloroplast or Mitochondria. Filter the host organelle features out of the feature table after assignment, before diversity and differential abundance, and tell me what fraction of reads that removed."

## What the Agent Will Do

1. Confirm the marker (16S/ITS/18S) and primer region, and load the per-feature sequences.
2. Select a reference database matched to the marker, and state its release.
3. Match the classifier reference to the primer region (use a region-matched pre-trained classifier or run extract-reads -> fit-classifier-naive-bayes).
4. Choose a classifier (naive Bayes, vsearch-consensus, or IDTAXA) appropriate to the goal, and state the confidence threshold.
5. Classify, leaving ranks unassigned (NA) below the confidence floor rather than forcing a label.
6. Attempt species ONLY by exact match (addSpecies) for 16S, or accept species for ITS.
7. Filter host organelle (Mitochondria, Chloroplast) and domain-unassigned/off-target features out of the table before any diversity or DA step (unless chloroplast is the study target in a phototroph community).
8. Report the rank the data supports and the three conditioning choices (classifier, database+release, region); hand the labelled, filtered table to diversity-analysis / differential-abundance.

## Tips

- Report the rank the data supports. For 16S, that is genus at best for many taxa, family for poorly resolved clades - not species. The tool will emit a species name; that does not make it licensed.
- A confidence of 0.95 means the model is sure GIVEN the taxon is in the database. It is not evidence the taxon is present. If the true organism is absent, the call is the nearest wrong relative.
- Match the reference to the primer region. Pointing a full-length classifier at V4 reads is the single most common silent accuracy loss.
- For a scikit-learn version error: retrain locally with fit-classifier-naive-bayes, download the classifier built for the exact QIIME2 release, or switch to classify-consensus-vsearch (no pickled model).
- Keep `Unassigned`/truncated features and label them honestly (e.g. `g__; s__`); do not silently drop or force-fill them. `Unassigned` at domain level usually means off-target (host, chimera, primer artifact) - filter, but document it.
- Never position-trim ITS to a fixed length - ITS is variable-length. Remove ITS primers, then classify.
- Use one database+release for an entire study; SILVA, GTDB, and Greengenes2 disagree on names and ranks.
- Filter host mitochondria and chloroplast features after assignment and before diversity/DA: universal 16S primers amplify host organelle rRNA, which otherwise inflates the table and deflates every real taxon by closure. Use `qiime taxa filter-table --p-exclude mitochondria,chloroplast` or the phyloseq subset_taxa equivalent (the rank-equality form is SILVA-138-specific). Exception: in phototroph/aquatic/mat communities, chloroplast 16S can be the signal - inspect before excluding.

## Classification Methods

| Method | Strength | Limitation |
|--------|----------|------------|
| Naive Bayes (classify-sklearn / assignTaxonomy) | fast (train once), well-benchmarked default | over-classifies if confidence left low; pre-trained model is sklearn-pinned |
| Alignment-consensus (classify-consensus-vsearch) | no trained model, immune to sklearn pinning, transparent hits | slower; consensus thresholds are extra knobs |
| IDTAXA (DECIPHER) | conservative, novelty-aware, refuses to over-descend | needs a pre-trained DECIPHER trainingSet |
| addSpecies (DADA2, exact match) | high-precision species calls | exact-match only; low recall (most ASVs left NA) |

## Reference Databases

| Database | Marker / scope | Best for |
|----------|----------------|----------|
| SILVA 138.x | 16S + 18S; Bacteria, Archaea, Eukarya; curated rRNA taxonomy | general-purpose 16S/18S default |
| GTDB r220 | 16S; Bacteria + Archaea only; genome-based, rank-normalized | environmental / under-named bacteria (names differ from SILVA) |
| Greengenes2 | 16S (V4-focused) on a genome-backbone tree, GTDB-harmonized | unifying 16S with shotgun on one tree (closed-reference) |
| UNITE | fungal ITS; species hypotheses | fungi; resolves to species |
| PR2 5.x | protist/eukaryote 18S; curated | microeukaryote 18S |
| RDP | 16S; legacy | historical reproducibility only |

## Related Skills

- amplicon-processing - Generate the ASV table that is classified here
- diversity-analysis - Alpha/beta diversity of the classified community table
- differential-abundance - Compositional differential abundance on the classified feature table
- qiime2-workflow - The QIIME2 CLI workflow this classification step plugs into
- read-qc/adapter-trimming - cutadapt primer removal before ASV inference and assignment
- metagenomics/kraken-classification - Shotgun (raw-read, not ASV) k-mer classification
- metagenomics/metaphlan-profiling - Shotgun marker-gene profiling; a different input artifact
- phylogenetics/tree-io - Phylogenetic tree for UniFrac / Faith PD on the classified table
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
