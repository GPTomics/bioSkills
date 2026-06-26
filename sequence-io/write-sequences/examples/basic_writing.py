#!/usr/bin/env python3
'''Write sequences to FASTA, FASTQ, and GenBank, showing the field requirements of each format'''
# Reference: biopython 1.83+ | Verify API if version differs
import os
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

outdir = tempfile.mkdtemp(prefix='write_seq_')

# FASTA writes '>' + description, so the description starts with the id to keep the header faithful
fasta_records = [
    SeqRecord(Seq('ATGCGATCGATCGATCGATCG'), id='seq1', description='seq1 first sequence'),
    SeqRecord(Seq('GCTAGCTAGCTAGCTAGCTA'), id='seq2', description='seq2 second sequence'),
    SeqRecord(Seq('TTAATTAATTAATTAATTAA'), id='seq3', description='seq3 third sequence')
]
fasta_path = os.path.join(outdir, 'output.fasta')
count = SeqIO.write(fasta_records, fasta_path, 'fasta')
print(f'Wrote {count} records to {fasta_path}')

# FASTQ requires a phred_quality per base; letter_annotations is length-locked to the sequence
fastq_record = SeqRecord(Seq('ATGCGATCG'), id='read1', description='read1')
fastq_record.letter_annotations['phred_quality'] = [40] * len(fastq_record.seq)
fastq_path = os.path.join(outdir, 'output.fastq')
SeqIO.write(fastq_record, fastq_path, 'fastq')
print(f'Wrote {fastq_path}')

# GenBank/EMBL require annotations['molecule_type'] since the BioPython 1.78 alphabet removal
gb_record = SeqRecord(Seq('ATGCGATCGATCG'), id='SEQ001', name='example')
gb_record.annotations['molecule_type'] = 'DNA'
gb_record.annotations['topology'] = 'linear'
gb_record.annotations['organism'] = 'Example organism'
gb_path = os.path.join(outdir, 'output.gb')
SeqIO.write(gb_record, gb_path, 'genbank')
print(f'Wrote {gb_path}')

# The header trap: a stale id-like token in the description is duplicated into the FASTA header
trap = SeqRecord(Seq('ACGT'), id='geneX', description='oldname kinase domain')
print('Header with stale id in description:', trap.format('fasta').splitlines()[0])

for path in (fasta_path, fastq_path, gb_path):
    os.remove(path)
os.rmdir(outdir)
