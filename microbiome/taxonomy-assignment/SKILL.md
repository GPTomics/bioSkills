---
name: bio-microbiome-taxonomy-assignment
description: Taxonomic classification of ASVs using reference databases like SILVA, GTDB, or UNITE. Covers naive Bayes classifiers (DADA2, IDTAXA) and exact matching approaches. Use when assigning taxonomy to ASVs after DADA2 amplicon processing.
tool_type: mixed
primary_tool: dada2
---

# Taxonomy Assignment

## DADA2 Naive Bayes Classifier

```r
library(dada2)

seqtab_nochim <- readRDS('seqtab_nochim.rds')

# SILVA for 16S (download from https://zenodo.org/record/4587955)
taxa <- assignTaxonomy(seqtab_nochim, 'silva_nr99_v138.1_train_set.fa.gz',
                       multithread = TRUE)

# Add species-level (exact matching)
taxa <- addSpecies(taxa, 'silva_species_assignment_v138.1.fa.gz')

# Check results
head(taxa)
```

## GTDB for 16S

```r
# GTDB-formatted database (better for environmental samples)
taxa_gtdb <- assignTaxonomy(seqtab_nochim, 'GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz',
                            multithread = TRUE)
```

## UNITE for ITS (Fungi)

```r
# UNITE database for fungal ITS
taxa_its <- assignTaxonomy(seqtab_nochim, 'sh_general_release_dynamic_25.07.2023.fasta',
                           multithread = TRUE)
```

## QIIME2 Feature Classifier

```bash
# Train classifier (one-time)
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138-99-seqs.qza \
    --i-reference-taxonomy silva-138-99-tax.qza \
    --o-classifier silva-138-99-nb-classifier.qza

# Classify ASVs
qiime feature-classifier classify-sklearn \
    --i-classifier silva-138-99-nb-classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza
```

## VSEARCH Exact Matching

```bash
# Faster but requires exact or near-exact matches
vsearch --usearch_global asv_seqs.fasta \
    --db silva_138_SSURef_NR99.fasta \
    --id 0.97 \
    --blast6out taxonomy_vsearch.tsv \
    --top_hits_only
```

## RDP Classifier

```r
library(dada2)

# RDP training set (less detailed than SILVA)
taxa_rdp <- assignTaxonomy(seqtab_nochim, 'rdp_train_set_18.fa.gz',
                           multithread = TRUE)
```

## IDTAXA (DECIPHER) - Often More Accurate

```r
library(DECIPHER)

# Load IDTAXA training set (download from http://www2.decipher.codes/Downloads.html)
load('SILVA_SSU_r138_2019.RData')  # Creates 'trainingSet' object

# Convert ASV sequences to DNAStringSet
dna <- DNAStringSet(getSequences(seqtab_nochim))

# Classify with IDTAXA
ids <- IdTaxa(dna, trainingSet, strand = 'top', processors = NULL, verbose = TRUE)

# Convert to matrix format like assignTaxonomy
ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
taxa_idtaxa <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, 'unclassified_')] <- NA
    taxa
}))
colnames(taxa_idtaxa) <- ranks
```

## Confidence Filtering

```r
# assignTaxonomy returns bootstrap confidence
# Filter low-confidence assignments
taxa_filtered <- taxa
taxa_filtered[taxa_filtered < 80] <- NA  # If using minBoot output

# Or use confidence threshold during assignment
taxa <- assignTaxonomy(seqtab_nochim, 'silva_nr99_v138.1_train_set.fa.gz',
                       minBoot = 80, multithread = TRUE)
```

## Combine into phyloseq

```r
library(phyloseq)

# Create phyloseq object
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE),
               tax_table(taxa))

# Add sample metadata
sample_data(ps) <- read.csv('sample_metadata.csv', row.names = 1)

# Rename ASVs for readability
taxa_names(ps) <- paste0('ASV', seq(ntaxa(ps)))
```

## Database Comparison

| Database | Organisms | Taxonomy | Updated |
|----------|-----------|----------|---------|
| SILVA 138.1 | Bacteria, Archaea, Eukaryotes | 7 ranks | 2024 |
| GTDB R220 | Bacteria, Archaea | 7 ranks (genome-based) | 2024 |
| RDP 18 | Bacteria, Archaea | 6 ranks | 2016 |
| UNITE 10.0 | Fungi | 7 ranks | 2024 |
| PR2 5.0 | Protists | 8 ranks | 2024 |

## Related Skills

- amplicon-processing - Generate ASV table for classification
- diversity-analysis - Analyze classified communities
- metagenomics/kraken-classification - Read-level taxonomic classification
