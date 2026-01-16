'''Reading and writing motif files in various formats'''
from Bio import motifs
from Bio.Seq import Seq
import tempfile
import os

# Create a sample motif
print('=== Create Sample Motif ===')
instances = [Seq('TGACCA'), Seq('TGACTA'), Seq('TGACGA'), Seq('TGACAA')]
m = motifs.create(instances)
m.name = 'MyMotif'
m.matrix_id = 'MA0001.1'
print(f'Created motif: {m.name}')
print(f'Consensus: {m.consensus}')
print(f'Degenerate: {m.degenerate_consensus}')

# Write to JASPAR format
print('\n=== JASPAR Format ===')
jaspar_str = m.format('jaspar')
print(jaspar_str)

# Write to TRANSFAC format
print('\n=== TRANSFAC Format ===')
transfac_str = m.format('transfac')
print(transfac_str[:500])

# Demonstrate reading back from file
print('\n=== Round-trip Test ===')
with tempfile.NamedTemporaryFile(mode='w', suffix='.jaspar', delete=False) as f:
    f.write(jaspar_str)
    temp_path = f.name

with open(temp_path) as f:
    m_read = motifs.read(f, 'jaspar')
print(f'Read back motif: {m_read.name}')
print(f'Consensus matches: {m_read.consensus == m.consensus}')
os.unlink(temp_path)

# Example JASPAR format parsing (inline for demo)
print('\n=== Parse JASPAR String ===')
jaspar_example = '''>MA0001.1 AGL3
A [  0   3  79  40  66  48  65  11  65  ]
C [ 94  75   4   3   1   2   5   2   3  ]
G [  1   0   3   4   1   0   5   3  28  ]
T [  2  19  11  50  29  47  22  81   1  ]
'''
from io import StringIO
m_parsed = motifs.read(StringIO(jaspar_example), 'jaspar')
print(f'Parsed: {m_parsed.name}')
print(f'Consensus: {m_parsed.consensus}')
print(f'Degenerate: {m_parsed.degenerate_consensus}')

# Information about the parsed motif
print('\n=== Motif Analysis ===')
pwm = m_parsed.counts.normalize(pseudocounts=0.5)
pssm = pwm.log_odds()
print(f'Length: {len(m_parsed)} positions')
print(f'Mean IC: {pssm.mean():.3f} bits')
print(f'Score range: [{pssm.min:.2f}, {pssm.max:.2f}]')
