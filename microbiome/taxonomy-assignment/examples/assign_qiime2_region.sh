#!/usr/bin/env bash
# Reference: QIIME2 2024.10+ (q2-feature-classifier) | Verify API if version differs
# Region-specific naive-Bayes training + classification for 16S V4 (515F/806R) ASVs.
# A full-length classifier on V4 reads fabricates and erases calls (Werner 2012; Bokulich 2018):
# match the reference to the primer region by extracting the amplicon in silico, then training.
set -euo pipefail

REF_SEQS='silva-138-99-seqs.qza'   # full-length SILVA reference sequences (FeatureData[Sequence])
REF_TAX='silva-138-99-tax.qza'     # matching taxonomy (FeatureData[Taxonomy])
REP_SEQS='rep-seqs.qza'            # the ASV representative sequences to classify

# 515F / 806R primers (Earth Microbiome Project V4); r-primer is given 5'->3', NOT reverse-complemented.
F_PRIMER='GTGYCAGCMGCCGCGGTAA'
R_PRIMER='GGACTACNVGGGTWTCTAAT'

# 1. In-silico PCR: extract the V4 window so reference k-mer composition matches the reads.
#    --p-min-length 50 drops too-short in-silico amplicons that would mistrain the classifier.
qiime feature-classifier extract-reads \
    --i-sequences "${REF_SEQS}" \
    --p-f-primer "${F_PRIMER}" --p-r-primer "${R_PRIMER}" \
    --p-min-length 50 --p-max-length 0 \
    --o-reads ref-seqs-515-806.qza

# 2. Train naive Bayes on the EXTRACTED region. The output .qza is a pickled scikit-learn model
#    tied to THIS QIIME2 release; a classifier from another release errors on load. Retraining
#    here (vs downloading a pre-trained .qza) guarantees the sklearn versions match.
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ref-seqs-515-806.qza \
    --i-reference-taxonomy "${REF_TAX}" \
    --o-classifier silva-138-99-515-806-nb-classifier.qza

# 3. Classify. --p-confidence 0.7 is the benchmarked default (Bokulich 2018): below it the lineage
#    is truncated to a shallower, more confident rank (honest unassigned), rather than over-called.
#    --p-n-jobs 1 keeps memory low; raising it multiplies memory (one classifier copy per job).
qiime feature-classifier classify-sklearn \
    --i-classifier silva-138-99-515-806-nb-classifier.qza \
    --i-reads "${REP_SEQS}" \
    --p-confidence 0.7 --p-read-orientation auto --p-n-jobs 1 \
    --o-classification taxonomy.qza

# Alternative immune to the scikit-learn version trap (no pickled model): alignment-consensus.
# A rank is reported only if --p-min-consensus (0.51) of the top --p-maxaccepts (10) hits agree.
qiime feature-classifier classify-consensus-vsearch \
    --i-query "${REP_SEQS}" \
    --i-reference-reads ref-seqs-515-806.qza \
    --i-reference-taxonomy "${REF_TAX}" \
    --p-maxaccepts 10 --p-perc-identity 0.8 --p-min-consensus 0.51 \
    --p-threads 8 \
    --o-classification taxonomy-vsearch.qza --o-search-results vsearch-hits.qza
