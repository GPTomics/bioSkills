'''Reading FASTQ files and working with quality scores'''
from Bio import SeqIO

print('=== Parsing FASTQ with Quality Scores ===')
for record in SeqIO.parse('sample.fastq', 'fastq'):
    qualities = record.letter_annotations['phred_quality']
    avg_qual = sum(qualities) / len(qualities)
    min_qual = min(qualities)
    max_qual = max(qualities)

    print(f'{record.id}:')
    print(f'  Length: {len(record.seq)} bp')
    print(f'  Avg Quality: {avg_qual:.1f}')
    print(f'  Min Quality: {min_qual}')
    print(f'  Max Quality: {max_qual}')
    print()

# Filter reads by quality
print('=== Filtering by Quality (avg >= 30) ===')
high_quality = [r for r in SeqIO.parse('sample.fastq', 'fastq')
                if sum(r.letter_annotations['phred_quality']) / len(r.seq) >= 30]
print(f'High quality reads: {len(high_quality)}')
for r in high_quality:
    print(f'  {r.id}')
