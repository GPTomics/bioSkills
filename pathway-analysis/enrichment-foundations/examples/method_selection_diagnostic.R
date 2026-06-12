# Reference: clusterProfiler 4.18.4+, org.Hs.eg.db 3.22+ | Verify API if version differs
# Method-selection and trustworthiness diagnostic for functional enrichment.
# Offline-runnable: uses only the local org.Hs.eg.db OrgDb (no KEGG/WikiPathways network calls).
# Decides ORA vs GSEA from the input shape, builds and inspects a defensible background
# universe against the all-genome default, demonstrates the named decreasing GSEA vector,
# and notes the inter-gene-correlation issue behind every preranked GSEA.

library(clusterProfiler)
library(org.Hs.eg.db)

out_dir <- tempdir()

choose_generation <- function(named_ranked_vector=NULL, gene_list=NULL) {
  if (!is.null(named_ranked_vector)) {
    is_named      <- !is.null(names(named_ranked_vector))
    is_decreasing <- !is.unsorted(rev(named_ranked_vector))
    if (is_named && is_decreasing) {
      return('FCS / GSEA: a per-gene statistic exists for (nearly) all genes; no cutoff needed')
    }
    return('FCS intended but vector is not a NAMED DECREASING numeric - fix before gseGO/GSEA')
  }
  if (!is.null(gene_list)) {
    return('ORA: a pre-selected list with no ranking; needs a defensible background universe')
  }
  'No usable input'
}

set.seed(123)   # fix the seed so any downstream permutation p-value is reproducible

# A coherent gene list: ENTREZ genes annotated to one BP term, so enrichGO returns real hits and
# the background contrast is visible rather than empty noise.
immune_genes <- unique(unlist(AnnotationDbi::mget('GO:0006955', org.Hs.egGO2ALLEGS)))   # immune response
hits <- head(immune_genes, 200)

all_entrez   <- keys(org.Hs.eg.db, keytype='ENTREZID')
testable_n   <- 12000   # genes that passed expression/independent filtering and entered the DE test
universe_set <- union(hits, sample(setdiff(all_entrez, hits), testable_n - length(hits)))   # testable genes incl. the hits

generation_for_list   <- choose_generation(gene_list=hits)
generation_for_list

# The background decides ORA significance: contrast the testable-gene universe with the
# silent all-genome default. The difference is the count of genes that were NEVER eligible to be DE.
genome_default_n   <- length(all_entrez)
universe_defensible_n <- length(universe_set)
ineligible_in_default <- genome_default_n - universe_defensible_n
background_comparison <- data.frame(
  universe        = c('all-genome default (silent)', 'defensible testable-gene set'),
  size            = c(genome_default_n, universe_defensible_n),
  extra_ineligible = c(ineligible_in_default, 0)
)
background_comparison

ora_default <- enrichGO(gene=hits, OrgDb=org.Hs.eg.db, keyType='ENTREZID', ont='BP',
                        pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.2,
                        minGSSize=10, maxGSSize=500)   # universe= omitted -> all-annotated-genes default

ora_defensible <- enrichGO(gene=hits, OrgDb=org.Hs.eg.db, keyType='ENTREZID', ont='BP',
                           universe=universe_set, pAdjustMethod='BH', pvalueCutoff=0.05,
                           qvalueCutoff=0.2, minGSSize=10, maxGSSize=500)

n_default    <- if (is.null(ora_default)) 0L else nrow(as.data.frame(ora_default))
n_defensible <- if (is.null(ora_defensible)) 0L else nrow(as.data.frame(ora_defensible))
background_effect <- data.frame(setting=c('all-genome default', 'testable-gene universe'),
                                significant_terms=c(n_default, n_defensible))
background_effect   # the universe choice, not the gene list, shifts the result

# The GSEA branch needs a NAMED numeric vector sorted DECREASING (a filtered list errors or mis-ranks).
ranking_stat <- rnorm(length(universe_set))          # stand-in for the signed -log10 p or t-statistic
names(ranking_stat) <- universe_set
gsea_vector <- sort(ranking_stat, decreasing=TRUE)
generation_for_vector <- choose_generation(named_ranked_vector=gsea_vector)
generation_for_vector

head(gsea_vector)   # named, decreasing - the shape gseGO/GSEA require

# Every clusterProfiler/fgsea GSEA is preranked = gene-permutation = competitive null with the
# WRONG variance under inter-gene correlation. A correlated set has a larger true variance than
# the independence assumption allows, so the screen is anti-conservative.
rho_bar  <- 0.1   # illustrative average pairwise correlation within a co-regulated set
set_n    <- 50
var_inflation <- 1 + (set_n - 1) * rho_bar   # Var(mean of n correlated scores) scales by this factor
camera_note <- sprintf('Variance inflation at rho-bar=%.2f, n=%d is %.1fx; CAMERA corrects this.',
                       rho_bar, set_n, var_inflation)
camera_note

provenance <- data.frame(
  clusterProfiler = as.character(packageVersion('clusterProfiler')),
  org.Hs.eg.db    = as.character(packageVersion('org.Hs.eg.db')),
  pAdjustMethod   = 'BH',
  universe        = 'testable genes (same filter as DE list)',
  seed            = 123,
  access_date     = as.character(Sys.Date())
)

write.csv(background_effect, file.path(out_dir, 'background_universe_effect.csv'), row.names=FALSE)
write.csv(provenance, file.path(out_dir, 'enrichment_provenance.csv'), row.names=FALSE)
cat('Diagnostic complete. Outputs written to', out_dir, '\n')
