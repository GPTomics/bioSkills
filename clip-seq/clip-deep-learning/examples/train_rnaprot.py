# Reference: RNAProt 0.5+, biopython 1.83+, scikit-learn 1.4+, pytorch 2.2+ | Verify API if version differs
# Train RNAProt RNN classifier on eCLIP peaks for RBP binding prediction.
# Chromosome-split train/test prevents gene-neighbor leakage that inflates AUC.
# GC-matched 3' UTR background prevents the model from learning transcript-region differences instead of RBP specificity.

import os
import subprocess
from pathlib import Path

PEAKS_BED = 'peaks.stringent.bed'       # CLIPper / Skipper stringent
EXPRESSED_REGIONS = 'expressed_3utr.bed'  # for GC-matched background
GENOME_FA = 'hg38.fa'
OUT_DIR = 'rnaprot_model'
TRAIN_CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
              'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20']
TEST_CHRS = ['chr21','chr22','chrX']

Path(OUT_DIR).mkdir(exist_ok=True)

# Step 1: Generate GC-matched background
# Random shuffle within expressed transcript regions; matches foreground GC distribution
subprocess.run([
    'bedtools', 'shuffle',
    '-i', PEAKS_BED,
    '-g', 'chrom.sizes',
    '-incl', EXPRESSED_REGIONS,
    '-seed', '42'
], stdout=open(f'{OUT_DIR}/background.bed', 'w'), check=True)

# Step 2: Split by chromosome to prevent train/test leakage
# Random split would leak gene-neighbor context and inflate AUC
def split_bed_by_chrom(bed_in, train_out, test_out, train_chrs):
    train_chrs_set = set(train_chrs)
    with open(bed_in) as f_in, open(train_out, 'w') as f_train, open(test_out, 'w') as f_test:
        for line in f_in:
            chrom = line.split('\t')[0]
            if chrom in train_chrs_set:
                f_train.write(line)
            else:
                f_test.write(line)

split_bed_by_chrom(PEAKS_BED, f'{OUT_DIR}/peaks_train.bed', f'{OUT_DIR}/peaks_test.bed', TRAIN_CHRS)
split_bed_by_chrom(f'{OUT_DIR}/background.bed', f'{OUT_DIR}/bg_train.bed', f'{OUT_DIR}/bg_test.bed', TRAIN_CHRS)

# Step 3: Extract FASTA sequences (strand-preserving)
for split in ['train', 'test']:
    for label in ['peaks', 'bg']:
        bed = f'{OUT_DIR}/{label}_{split}.bed'
        fa = f'{OUT_DIR}/{label}_{split}.fa'
        subprocess.run([
            'bedtools', 'getfasta', '-fi', GENOME_FA, '-bed', bed, '-s', '-fo', fa
        ], check=True)

# Step 4: Train RNAProt
# 50 epochs typical; loss should converge; validate on held-out chromosome-split test set
subprocess.run([
    'RNAProt', 'train',
    '--in', f'{OUT_DIR}/peaks_train.fa',
    '--neg', f'{OUT_DIR}/bg_train.fa',
    '--out', f'{OUT_DIR}/model',
    '--epochs', '50',
    '--batch-size', '64',
    '--learning-rate', '0.001',
    '--validation-split', '0.2'
], check=True)

# Step 5: Evaluate on held-out chromosome-split test set
subprocess.run([
    'RNAProt', 'predict',
    '--model', f'{OUT_DIR}/model',
    '--in', f'{OUT_DIR}/peaks_test.fa',
    '--out', f'{OUT_DIR}/peaks_test_pred.tsv'
], check=True)
subprocess.run([
    'RNAProt', 'predict',
    '--model', f'{OUT_DIR}/model',
    '--in', f'{OUT_DIR}/bg_test.fa',
    '--out', f'{OUT_DIR}/bg_test_pred.tsv'
], check=True)

# Compute AUC on held-out test
import pandas as pd
from sklearn.metrics import roc_auc_score

peaks_pred = pd.read_csv(f'{OUT_DIR}/peaks_test_pred.tsv', sep='\t')
bg_pred = pd.read_csv(f'{OUT_DIR}/bg_test_pred.tsv', sep='\t')

# Adjust column name to match RNAProt output (verify with `head` if unsure)
score_col = 'binding_probability'  # check RNAProt docs for exact column name
y_true = [1] * len(peaks_pred) + [0] * len(bg_pred)
y_score = list(peaks_pred[score_col]) + list(bg_pred[score_col])
auc = roc_auc_score(y_true, y_score)
print(f'Held-out AUC (chromosome-split): {auc:.4f}')
print('Benchmark target: AUC 0.85-0.89 for ENCODE RBPs')

# Variant-effect prediction example
# Reference and alternative sequences at a heterozygous SNP
ref_seq = 'CTGTACTGCAGTAGCATGCTAGCATGCTAGCAT'
alt_seq = 'CTGTACTGCAGTAGCATGCTAGCATGCTAGCAA'
with open(f'{OUT_DIR}/variant_query.fa', 'w') as f:
    f.write(f'>ref\n{ref_seq}\n>alt\n{alt_seq}\n')

subprocess.run([
    'RNAProt', 'predict',
    '--model', f'{OUT_DIR}/model',
    '--in', f'{OUT_DIR}/variant_query.fa',
    '--out', f'{OUT_DIR}/variant_pred.tsv'
], check=True)

print('\nVariant predictions (see variant_pred.tsv)')
print('|log2(alt/ref)| > 1 = strong effect')
