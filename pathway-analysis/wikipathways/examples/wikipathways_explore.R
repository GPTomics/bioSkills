# Reference: clusterProfiler 4.18+, rWikiPathways 1.26+ | Verify API if version differs
# NOTE: every rWikiPathways query and download hits the live WikiPathways API/archive over the
# network (needs internet). downloadPathwayArchive pins a dated MONTHLY release = the reproducible
# pattern; enrichWP/gseWP without a pin pull current/, which changes month to month.
library(rWikiPathways)
library(clusterProfiler)
library(tidyr)

listOrganisms()                       # supported species; ~30+ full scientific names
get_wp_organisms()                    # plural accessor; the organism string must match exactly

human_pathways <- listPathways('Homo sapiens')
head(human_pathways, 20)

# searchPathways() is NOT a current function (errors); use findPathwaysByText
cancer_pathways <- findPathwaysByText('cancer')
head(cancer_pathways)

pathway_info <- getPathwayInfo('WP554')   # metadata incl. last-edit; check before trusting a hit
pathway_info

# genes by BridgeDb system code: 'L'=Entrez, 'H'=HGNC symbol, 'En'=Ensembl
pathway_entrez <- getXrefList('WP554', 'L')
pathway_entrez

# Reproducible pattern: pin a dated GMT, split the compound term, run enricher on the pinned sets.
# format defaults to gpml -> pass format='gmt'; organism=NULL would open the index, not download.
gmt <- downloadPathwayArchive(date = '20240310', organism = 'Homo sapiens',
                              format = 'gmt', destpath = tempdir())
wp2gene <- read.gmt(file.path(tempdir(), gmt))
wp2gene <- separate(wp2gene, term, c('name', 'version', 'wpid', 'org'), sep = '%')
t2g <- wp2gene[, c('wpid', 'gene')]
t2n <- wp2gene[, c('wpid', 'name')]

# report date='20240310' in methods; entrez_ids/all_entrez built upstream (see wikipathways_ora.R)
wp_pinned <- enricher(entrez_ids, universe = all_entrez, TERM2GENE = t2g, TERM2NAME = t2n)
as.data.frame(wp_pinned)
