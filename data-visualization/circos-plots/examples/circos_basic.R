#!/usr/bin/env Rscript
# Basic circos plot with circlize

library(circlize)

# Initialize with human genome
circos.initializeWithIdeogram(species = 'hg38', plotType = c('axis', 'labels'))

# Create example data
set.seed(42)
bed <- data.frame(
    chr = paste0('chr', sample(1:22, 200, replace = TRUE)),
    start = sample(1:1e8, 200),
    value = runif(200)
)
bed$end <- bed$start + 1e6

# Add scatter track
circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = 'steelblue')
}, track.height = 0.1)

# Add histogram track
circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                       col = 'coral', border = NA)
}, track.height = 0.1)

# Add links between regions
link_data <- data.frame(
    chr1 = c('chr1', 'chr3', 'chr7'),
    start1 = c(1e7, 5e7, 8e7),
    end1 = c(2e7, 6e7, 9e7),
    chr2 = c('chr5', 'chr10', 'chr15'),
    start2 = c(3e7, 4e7, 2e7),
    end2 = c(4e7, 5e7, 3e7)
)

for (i in seq_len(nrow(link_data))) {
    circos.link(link_data$chr1[i], c(link_data$start1[i], link_data$end1[i]),
                link_data$chr2[i], c(link_data$start2[i], link_data$end2[i]),
                col = adjustcolor('purple', alpha.f = 0.5))
}

circos.clear()

# Save to PDF
pdf('circos_output.pdf', width = 10, height = 10)
circos.initializeWithIdeogram(species = 'hg38', plotType = c('axis', 'labels'))
circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = 'steelblue')
}, track.height = 0.1)
circos.clear()
dev.off()
cat('Saved circos_output.pdf\n')
