# Model-corrected distance tree in ape. Deliverable: a TN93+gamma distance matrix and a FastME
# balanced-minimum-evolution tree with bootstrap -- the correction is where the biology lives,
# unlike Biopython's identity-only DistanceCalculator. Outputs go to tempdir(), no strays.
# Reference: ape 5.8+ | Verify API if version differs

library(ape)

aln <- as.DNAbin(matrix(c(
  strsplit('ATGCATGCATGCATGCATGC', '')[[1]],
  strsplit('ATGCATGCATGAATGCATGC', '')[[1]],
  strsplit('ATGCATGAATGCATGCATGC', '')[[1]],
  strsplit('ATGAATGCATGCATGCATGC', '')[[1]],
  strsplit('ATGAATGAATGCATGCATGC', '')[[1]]),
  nrow = 5, byrow = TRUE,
  dimnames = list(c('Human', 'Chimp', 'Gorilla', 'Mouse', 'Rat'), NULL)))

# TN93 = two transition rates + unequal base freqs; gamma applies ASRV (alpha < 1 = strong heterogeneity)
d <- dist.dna(aln, model = 'TN93', gamma = 0.5)

# Saturation pre-flight: plot transitions vs a corrected distance; a PLATEAU means signal is erased.
ts <- dist.dna(aln, model = 'TS')
jc <- dist.dna(aln, model = 'JC69')
# plot(jc, ts) would show roughly linear (unsaturated) vs bent-over (saturated). Xia Iss test lives in DAMBE.

# FastME balanced minimum evolution: the modern best distance tree (searches, not a single greedy pass).
tree <- fastme.bal(d, nni = TRUE, spr = TRUE)

# Bootstrap: 100-1000 reps standard; this is sampling PRECISION, not accuracy.
bs <- boot.phylo(tree, as.matrix(aln), function(x) fastme.bal(dist.dna(x, model = 'TN93')), B = 100)

out <- file.path(tempdir(), 'distance_tree.nwk')
write.tree(tree, out)
cat('Wrote model-corrected FastME tree to', out, '\n')
print(tree)
