library(rWikiPathways)

available_organisms <- listOrganisms()
available_organisms

human_pathways <- listPathways('Homo sapiens')
head(human_pathways, 20)

cancer_pathways <- searchPathways('cancer', 'Homo sapiens')
cancer_pathways

pathway_info <- getPathwayInfo('WP554')
pathway_info

pathway_genes_symbol <- getXrefList('WP554', 'H')
pathway_genes_entrez <- getXrefList('WP554', 'L')
pathway_genes_symbol

downloadPathwayArchive(organism = 'Homo sapiens', format = 'gmt', destpath = '.')
