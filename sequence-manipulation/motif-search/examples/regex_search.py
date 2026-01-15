'''Using regex for flexible pattern matching'''
from Bio.Seq import Seq
import re

seq = Seq('ATGCTATAAAAGATCTATAAATGCTATATAGCATGC')
print(f'Sequence: {seq}')

# Basic regex search
print('\n=== Regex Search ===')
def regex_search(seq, pattern):
    return [(m.start(), m.group()) for m in re.finditer(pattern, str(seq))]

# Find TATA box variants
print('TATA box variants (TATA[AT]A[AT]):')
matches = regex_search(seq, 'TATA[AT]A[AT]')
for pos, match in matches:
    print(f'  Position {pos}: {match}')

# IUPAC ambiguity code conversion
print('\n=== IUPAC Pattern Search ===')
IUPAC_DNA = {
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def iupac_to_regex(pattern):
    regex = ''
    for char in pattern.upper():
        regex += IUPAC_DNA.get(char, char)
    return regex

# Search for degenerate pattern
iupac_pattern = 'TATAWAW'  # W = A or T
regex_pattern = iupac_to_regex(iupac_pattern)
print(f'IUPAC pattern: {iupac_pattern}')
print(f'Regex pattern: {regex_pattern}')
matches = regex_search(seq, regex_pattern)
for pos, match in matches:
    print(f'  Position {pos}: {match}')

# Find start codons with Kozak context
print('\n=== Kozak Consensus ===')
seq2 = Seq('ATGCGCCATGGCTACCATGGATGC')
kozak_pattern = '[AG]CC(ATG)G'  # Capture group for ATG
for m in re.finditer(kozak_pattern, str(seq2)):
    print(f'Position {m.start()}: {m.group()} (ATG at {m.start(3)})')
