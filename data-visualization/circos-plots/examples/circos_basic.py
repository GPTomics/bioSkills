#!/usr/bin/env python3
'''Basic genome circos with pyCirclize (Shimoyama 2024). Matches SKILL.md pyCirclize API.'''
# Reference: pyCirclize 1.4+, matplotlib 3.8+ | Verify API if version differs

from pycirclize import Circos
import matplotlib.pyplot as plt

# Sector sizes = chromosome lengths (hg38, autosomes only for clarity)
sectors = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983, 'chr22': 50818468,
}

circos = Circos(sectors, space=2)                       # 2 degrees gap between sectors

# Outer ring: ideogram + labels
for sector in circos.sectors:
    sector.text(sector.name, r=110, size=8)
    sector.axis(r_lim=(95, 100), fc='lightgrey')

# Inter-sector chord (e.g., BCR-ABL1 t(9;22) translocation)
circos.link(('chr9', 133600000, 133700000),
            ('chr22', 23200000, 23300000),
            color='#D55E00', alpha=0.5)

fig = circos.plotfig()
fig.savefig('circos_basic.pdf', bbox_inches='tight', dpi=300)
plt.close(fig)
