# Reference: alakazam 1.3+ shazam 1.2+ scoper 1.3+ dowser 2.x tigger 1.x | Verify API if version differs
# Decision-grade BCR repertoire analysis with Immcantation.
# The threshold is DERIVED from the data (distToNearest -> findThreshold), never hardcoded.
# Order matters: TIGGER genotype -> germline -> threshold -> clones -> SHM -> selection -> trees.

library(alakazam)
library(shazam)
library(scoper)
library(tigger)
library(dplyr)

db <- readChangeoDb('bcr_data_airr.tsv')
cat('Loaded', nrow(db), 'sequences\n')

# Step 0: personalize the germline. An unrecorded personal allele otherwise reads as
# recurrent SHM at a fixed position and adds spurious distance to distToNearest.
ighv <- readIgFasta('IMGT_Human_IGHV.fasta')
novel <- findNovelAlleles(db, germline_db = ighv, v_call = 'v_call', nproc = 1)
genotype <- inferGenotypeBayesian(db, germline_db = ighv, novel = novel, find_unmutated = TRUE)
gt_seqs <- genotypeFasta(genotype, germline_db = ighv, novel = novel)
db <- reassignAlleles(db, genotype_db = gt_seqs)

# Step 1: derive the clonal threshold from the bimodal distance-to-nearest distribution.
# normalize='len' puts the threshold on a 0-1 scale comparable across CDR3 lengths.
db <- distToNearest(db, sequenceColumn = 'junction', vCallColumn = 'v_call',
                    jCallColumn = 'j_call', model = 'ham', normalize = 'len', nproc = 1)

thr_obj <- findThreshold(db$dist_nearest, method = 'density')
threshold <- thr_obj@threshold
cat('Derived clonal threshold:', round(threshold, 3), '\n')
# A unimodal dist_nearest gives NA here -> switch to spectralClones(method='novj').

# Step 2: cluster sequences into clonal families at the derived threshold.
# Nucleotide (not AA) junction distance because SHM is a nucleotide process.
results <- hierarchicalClones(db, threshold = threshold, method = 'nt', linkage = 'single')
db <- as.data.frame(results)
n_clones <- length(unique(db$clone_id))
cat('Identified', n_clones, 'clonal lineages\n')

clone_sizes <- db %>% group_by(clone_id) %>% summarize(size = n()) %>% arrange(desc(size))
cat('Largest clone has', max(clone_sizes$size), 'sequences\n')

# Step 3: quantify SHM against the D-masked germline, restricted to the V region.
# regionDefinition=IMGT_V stops before CDR3 (junctional bases have no germline template).
# frequency=TRUE reports mutations per informative position, comparable across coverage.
# NOTE: this assumes germline_alignment_d_mask exists (dowser::createGermlines upstream).
db <- observedMutations(db, sequenceColumn = 'sequence_alignment',
                        germlineColumn = 'germline_alignment_d_mask',
                        regionDefinition = IMGT_V, frequency = TRUE, nproc = 1)
# Columns added: mu_freq_cdr_r, mu_freq_cdr_s, mu_freq_fwr_r, mu_freq_fwr_s

mutation_summary <- db %>%
    summarize(mean_cdr_r = mean(mu_freq_cdr_r, na.rm = TRUE),
              mean_fwr_r = mean(mu_freq_fwr_r, na.rm = TRUE))
cat('\nMutation summary (replacement frequency):\n')
print(mutation_summary)

# Step 4: test for selection. calcBaseline builds per-sequence posteriors from a
# codon+motif-aware null; groupBaseline convolves them per group. Raw R/S is NOT selection.
baseline <- calcBaseline(db, testStatistic = 'focused', regionDefinition = IMGT_V, nproc = 1)
grouped <- groupBaseline(baseline, groupBy = 'clone_id')
cat('\nSelection analysis complete (sigma>0 in CDR = positive selection)\n')

# Step 5: compare diversity at equal depth. uniform=TRUE (default) resamples every
# group to equal N; comparing raw diversity across unequal depths measures depth, not biology.
# min_n floor of 30 is a small demo value; raise it toward the smallest real sample size.
div <- alphaDiversity(db, group = 'clone_id', clone = 'clone_id',
                      min_q = 0, max_q = 2, step_q = 0.5, min_n = 30, ci = 0.95, nboot = 200)
cat('\nDiversity curve computed\n')

# V gene usage counted once per clone to avoid clonal-expansion bias.
v_usage <- countGenes(db, gene = 'v_call', clone = 'clone_id', mode = 'gene')
cat('\nTop V genes:\n')
print(head(v_usage, 10))

cat('\nAnalysis complete.\n')
