#!/usr/bin/env Rscript
# Reference: TitanCNA 1.40+ | Verify API if version differs
#
# Subclonal copy number with TITAN: jointly infers copy number, LOH, and the cellular
# prevalence of clonal clusters from allele counts AND corrected read depth.
#
# TITAN needs BOTH inputs. The read-depth chain -- correctReadDepth -> getPositionOverlap
# -> log transform -> filterData -- is mandatory: runEMclonalCN requires data$logR and
# fails without it. correctReadDepth takes tumour + normal coverage WIGs (from HMMcopy
# readCounter) plus GC and mappability WIGs (reference data).
#
# The number of clonal clusters is not known a priori -- this script sweeps it and
# selects by the S_Dbw validity index.

suppressMessages(library(TitanCNA))

args <- commandArgs(trailingOnly = TRUE)
allele_counts <- ifelse(length(args) >= 1, args[1], 'tumour.allelicCounts.tsv')
sample_id     <- ifelse(length(args) >= 2, args[2], 'tumour')
max_clusters  <- ifelse(length(args) >= 3, as.integer(args[3]), 5)

# Read-depth WIG files: tumour/normal coverage (HMMcopy readCounter) + GC/mappability.
tum_wig  <- 'tumour.wig'
norm_wig <- 'normal.wig'
gc_wig   <- 'gc.wig'
map_wig  <- 'map.wig'

# Allele counts at het SNPs.
data <- loadAlleleCounts(allele_counts, genomeStyle = 'UCSC')

# Mandatory read-depth correction: GC/mappability bias correction, overlay logR onto the
# het positions, natural-log transform, then filter on depth/mappability.
# genomeStyle MUST match loadAlleleCounts ('UCSC' here) -- correctReadDepth defaults to
# 'NCBI', and a mismatch makes getPositionOverlap match no chromosomes (logR all NA).
cnData <- correctReadDepth(tum_wig, norm_wig, gc_wig, map_wig, genomeStyle = 'UCSC')
data$logR <- log(2 ^ getPositionOverlap(data$chr, data$posn, cnData))
data <- filterData(data, 1:24, minDepth = 10, maxDepth = 200, map = NULL)

# Sweep the number of clonal clusters; TITAN does not infer it -- model selection does.
fits <- list()
for (k in seq_len(max_clusters)) {
    params <- loadDefaultParameters(copyNumber = 8, numberClonalClusters = k,
                                    symmetric = TRUE, data = data)
    conv <- runEMclonalCN(data, params, maxiter = 20, txnExpLen = 1e15,
                          useOutlierState = FALSE)
    optimalPath <- viterbiClonalCN(data, conv)
    # computeSDbwIndex needs the CORRECTED results data frame, not the Viterbi path.
    # outputTitanResults(correctResults = TRUE) returns $results and $corrResults.
    titan <- outputTitanResults(data, conv, optimalPath, filename = NULL,
                                posteriorProbs = FALSE, correctResults = TRUE,
                                proportionThreshold = 0.05,
                                proportionThresholdClonal = 0.05)
    # S_Dbw validity index: lower indicates a better-supported cluster number.
    # computeSDbwIndex returns a list with dens.bw, scat, and S_DbwIndex (their sum).
    sdbw <- computeSDbwIndex(titan$corrResults, data.type = 'LogRatio',
                             centroid.method = 'median', S_Dbw.method = 'Tong')
    fits[[k]] <- list(k = k, conv = conv, optimalPath = optimalPath,
                      sdbw = sdbw$S_DbwIndex, loglik = tail(conv$loglik, 1))
    cat(sprintf('clusters=%d  loglik=%.1f  S_Dbw=%.4f\n', k, fits[[k]]$loglik, fits[[k]]$sdbw))
}

# Select the cluster number minimizing the S_Dbw validity index.
best <- fits[[which.min(sapply(fits, function(f) f$sdbw))]]
cat(sprintf('\nSelected %d clonal cluster(s) by S_Dbw.\n', best$k))

outprefix <- paste0(sample_id, '.titan.cluster', best$k)
outputTitanResults(data, best$conv, best$optimalPath,
                   filename = paste0(outprefix, '.segs.txt'), posteriorProbs = FALSE)
cat('Output:', paste0(outprefix, '.segs.txt'), '\n')
cat('Reminder: confirm whole-genome-doubling status from absolute copy number before\n',
    'interpreting subclonal states; a missed WGD halves every copy number.\n')
