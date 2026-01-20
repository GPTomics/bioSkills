# Circos Plots Usage Guide

Create circular genome visualizations with multiple data tracks.

## Prerequisites

```bash
# Circos (Perl)
conda install -c bioconda circos

# pyCircos (Python)
pip install pyCircos

# circlize (R)
install.packages('circlize')
```

## Quick Start

### circlize (R) - Fastest

```r
library(circlize)

# Initialize with human genome
circos.initializeWithIdeogram(species = 'hg38')

# Add data track
bed <- data.frame(
    chr = paste0('chr', sample(1:22, 100, replace=TRUE)),
    start = sample(1:1e8, 100),
    end = sample(1:1e8, 100) + 1e6,
    value = runif(100)
)

circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, pch=16, cex=0.5, col='red')
})

circos.clear()
```

### pyCircos (Python)

```python
from pycircos import Gcircle, Garc

circle = Gcircle()

chromosomes = [('chr1', 248956422), ('chr2', 242193529), ('chr3', 198295559)]
for name, length in chromosomes:
    arc = Garc(arc_id=name, size=length, interspace=2, raxis_range=(850, 900))
    circle.add_garc(arc)
circle.set_garcs()

# Add scatter plot
import numpy as np
for name, length in chromosomes:
    positions = np.random.randint(0, length, 50)
    values = np.random.random(50)
    circle.scatterplot(name, data=values, positions=positions,
                       raxis_range=(700, 800), facecolor='red')

circle.figure.savefig('circos.png', dpi=300)
```

## Common Use Cases

### CNV Visualization

Show copy number gains (red) and losses (blue) as bars or heatmap tracks.

### Gene Fusions

Display fusion breakpoints as arcs connecting two genomic positions.

### Hi-C Contacts

Visualize chromatin interactions as links between genomic regions.

### Multi-omics Summary

Layer multiple data types: variants, expression, methylation around the genome.

## Tips

- **Track order**: Outermost = chromosome ideograms, inner tracks = data
- **Spacing**: Use `interspace` parameter to separate chromosomes
- **Colors**: Use consistent color schemes across tracks
- **Export**: SVG for publication, PNG for quick viewing

## See Also

- [Circos documentation](http://circos.ca/documentation/)
- [circlize book](https://jokergoo.github.io/circlize_book/)
