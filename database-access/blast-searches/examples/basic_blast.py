'''Basic remote BLAST searches with Biopython'''
from Bio.Blast import NCBIWWW, NCBIXML

# Example nucleotide sequence (human beta-globin partial)
sequence = '''ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG'''

print("Running BLASTN against nt database...")
print("(This may take 30-60 seconds)\n")

result_handle = NCBIWWW.qblast('blastn', 'nt', sequence, hitlist_size=10)
blast_record = NCBIXML.read(result_handle)
result_handle.close()

print(f"Query: {blast_record.query[:50]}...")
print(f"Query length: {blast_record.query_length}")
print(f"Database: {blast_record.database}")
print(f"Hits found: {len(blast_record.alignments)}\n")

print("Top hits:")
for i, alignment in enumerate(blast_record.alignments[:5], 1):
    hsp = alignment.hsps[0]
    identity_pct = 100 * hsp.identities / hsp.align_length
    print(f"\n{i}. {alignment.accession}")
    print(f"   {alignment.title[:70]}...")
    print(f"   E-value: {hsp.expect:.2e}")
    print(f"   Score: {hsp.bits:.1f} bits")
    print(f"   Identity: {hsp.identities}/{hsp.align_length} ({identity_pct:.1f}%)")
    print(f"   Gaps: {hsp.gaps}")
