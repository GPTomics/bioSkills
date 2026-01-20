#!/usr/bin/env python3
'''Basic circos plot with pyCircos'''

from pycircos import Gcircle, Garc
import numpy as np

chromosomes = [
    ('chr1', 248956422), ('chr2', 242193529), ('chr3', 198295559),
    ('chr4', 190214555), ('chr5', 181538259), ('chr6', 170805979),
    ('chr7', 159345973), ('chr8', 145138636), ('chr9', 138394717),
    ('chr10', 133797422), ('chr11', 135086622), ('chr12', 133275309)
]

circle = Gcircle()

for name, length in chromosomes:
    arc = Garc(arc_id=name, size=length, interspace=3, raxis_range=(850, 900),
               labelposition=60, label_visible=True)
    circle.add_garc(arc)

circle.set_garcs()

# Add scatter track
for name, length in chromosomes:
    positions = np.random.randint(0, length, 30)
    values = np.random.random(30)
    circle.scatterplot(name, data=values, positions=positions,
                       raxis_range=(700, 800), facecolor='steelblue', markersize=4)

# Add bar track
for name, length in chromosomes:
    positions = np.linspace(0, length, 50, dtype=int)
    values = np.random.random(50) * 100
    circle.barplot(name, data=values, positions=positions,
                   raxis_range=(550, 680), facecolor='coral')

# Add a link between regions
circle.chord_plot(('chr1', 50000000, 60000000), ('chr5', 100000000, 110000000),
                  raxis_range=(0, 500), facecolor='purple', alpha=0.5)

fig = circle.figure
fig.savefig('circos_output.png', dpi=300, bbox_inches='tight')
print('Saved circos_output.png')
