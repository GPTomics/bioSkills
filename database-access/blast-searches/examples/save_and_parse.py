'''Save BLAST results and parse later'''
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

def run_and_save_blast(sequence, output_file, program='blastn', database='nt'):
    '''Run BLAST and save XML results'''
    print(f"Running {program} against {database}...")
    result_handle = NCBIWWW.qblast(program, database, sequence, hitlist_size=50)

    with open(output_file, 'w') as out:
        out.write(result_handle.read())
    result_handle.close()
    print(f"Results saved to {output_file}")

def parse_blast_xml(xml_file):
    '''Parse saved BLAST XML file'''
    with open(xml_file) as f:
        blast_record = NCBIXML.read(f)

    hits = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            hits.append({
                'accession': alignment.accession,
                'title': alignment.title,
                'length': alignment.length,
                'evalue': hsp.expect,
                'bits': hsp.bits,
                'identities': hsp.identities,
                'align_length': hsp.align_length,
                'identity_pct': 100 * hsp.identities / hsp.align_length,
                'query_start': hsp.query_start,
                'query_end': hsp.query_end,
                'sbjct_start': hsp.sbjct_start,
                'sbjct_end': hsp.sbjct_end
            })
    return hits

# Example sequence
sequence = '''ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG'''

# Run and save
run_and_save_blast(sequence, 'blast_results.xml')

# Parse saved results
print("\nParsing saved results...")
hits = parse_blast_xml('blast_results.xml')

print(f"\nFound {len(hits)} HSPs")
print("\nTop hits by E-value:")
for hit in sorted(hits, key=lambda x: x['evalue'])[:10]:
    print(f"  {hit['accession']}: {hit['identity_pct']:.1f}% identity, E={hit['evalue']:.2e}")
