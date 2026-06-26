'''Build a CAI index from highly expressed reference genes, score a query, and codon-optimize it.

CodonAdaptationIndex replaces the Bio.SeqUtils.CodonUsage module removed in Biopython 1.82.
CAI is only meaningful against an expression-biased reference (highly expressed genes of the
target host), so the index is built here from a small reference set, not a bundled table.
'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex
from Bio.Data.CodonTable import standard_dna_table

# Stand-in for highly expressed E. coli genes (ribosomal proteins, EFs), which strongly
# prefer GCT/CTG/GAC/AAA/GGT. In practice these are parsed from a FASTA with SeqIO.parse.
reference_genes = [
    Seq('ATGGCTGCTGCTCTGCTGCTGGACGACAAAAAAGGTGGTGCTCTGGACAAAGGTGCTTAA'),
    Seq('ATGGCTCTGCTGGACAAAGGTGCTGCTCTGCTGGACGACAAAAAAGGTGGTGCTCTGTAA'),
    Seq('ATGCTGGCTGACGACAAAGGTGGTGCTGCTCTGCTGAAAAAAGACGGTGCTCTGGCTTAA'),
]

cai = CodonAdaptationIndex(reference_genes, table=standard_dna_table)
print(f'Index built from {len(reference_genes)} reference genes ({len(cai)} codon weights)')
print(f'Relative adaptiveness of GCT (Ala): {cai["GCT"]:.3f}')

# Score a query against the reference. The query MUST be in frame: calculate steps codons
# from position 0 without checking, so an out-of-frame sequence is silently mis-scored.
# This query uses the host-dispreferred synonyms, so its CAI is low.
query = Seq('ATGGCAGCATTAGATGATAAAGGATAA')
print(f'\nQuery: {query} ({len(query) // 3} codons)')
print(f'Query CAI: {cai.calculate(query):.3f}')

# Codon-optimize: replace each amino acid with the host most-preferred synonymous codon.
# strict=False avoids a ValueError when two codons are equally preferred in a small reference.
optimized = cai.optimize(query, seq_type='DNA', strict=False)
print(f'\nOptimized: {optimized}')
print(f'Optimized CAI: {cai.calculate(optimized):.3f}')

# The optimization must preserve the protein. This is the only correctness guarantee
# optimize() makes; it does not screen for ramp, 5' structure, or cryptic elements.
assert optimized.translate() == query.translate()
print(f'Protein preserved: {optimized.translate()}')

changes = [(i // 3, str(query)[i:i+3], str(optimized)[i:i+3])
           for i in range(0, len(query) - 2, 3) if str(query)[i:i+3] != str(optimized)[i:i+3]]
print(f'\nCodons changed: {len(changes)}')
for codon_index, old, new in changes:
    print(f'  codon {codon_index}: {old} -> {new}')
