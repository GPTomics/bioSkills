---
name: bio-batch-processing
description: Process multiple sequence files in batch using Biopython. Use when working with many files, merging/splitting sequences, or automating file operations across directories.
tool_type: python
primary_tool: Bio.SeqIO
---

# Batch Processing

Process multiple sequence files efficiently using Biopython.

## Required Imports

```python
from pathlib import Path
from Bio import SeqIO
```

## Process Multiple Files

### Iterate Over Files in Directory
```python
from pathlib import Path

for fasta_file in Path('data/').glob('*.fasta'):
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    print(f'{fasta_file.name}: {len(records)} sequences')
```

### Process All FASTQ Files
```python
for fq_file in Path('.').glob('*.fastq'):
    count = sum(1 for _ in SeqIO.parse(fq_file, 'fastq'))
    print(f'{fq_file.name}: {count} reads')
```

### Recursive File Search
```python
for gb_file in Path('data/').rglob('*.gb'):
    print(f'Found: {gb_file}')
```

## Merge Files

### Merge All FASTA Files
```python
from pathlib import Path

def all_records(directory, pattern, format):
    for filepath in Path(directory).glob(pattern):
        yield from SeqIO.parse(filepath, format)

records = all_records('data/', '*.fasta', 'fasta')
count = SeqIO.write(records, 'merged.fasta', 'fasta')
print(f'Merged {count} records')
```

### Merge with Source Tracking
```python
def records_with_source(directory, pattern, format):
    for filepath in Path(directory).glob(pattern):
        for record in SeqIO.parse(filepath, format):
            record.description = f'{record.description} [source={filepath.name}]'
            yield record

records = records_with_source('data/', '*.fasta', 'fasta')
SeqIO.write(records, 'merged_tracked.fasta', 'fasta')
```

### Merge Specific Files
```python
files = ['sample1.fasta', 'sample2.fasta', 'sample3.fasta']

def merge_files(file_list, format):
    for filepath in file_list:
        yield from SeqIO.parse(filepath, format)

SeqIO.write(merge_files(files, 'fasta'), 'combined.fasta', 'fasta')
```

## Split Files

### Split by Number of Records
```python
from itertools import islice

def split_file(input_file, format, records_per_file, output_prefix):
    records = SeqIO.parse(input_file, format)
    file_num = 1
    while True:
        batch = list(islice(records, records_per_file))
        if not batch:
            break
        output_file = f'{output_prefix}_{file_num}.{format}'
        SeqIO.write(batch, output_file, format)
        print(f'Wrote {len(batch)} records to {output_file}')
        file_num += 1

split_file('large.fasta', 'fasta', 1000, 'split')
```

### Split by Sequence ID Prefix
```python
from collections import defaultdict

records_by_prefix = defaultdict(list)
for record in SeqIO.parse('input.fasta', 'fasta'):
    prefix = record.id.split('_')[0]
    records_by_prefix[prefix].append(record)

for prefix, records in records_by_prefix.items():
    SeqIO.write(records, f'{prefix}.fasta', 'fasta')
```

### One Sequence Per File
```python
for record in SeqIO.parse('multi.fasta', 'fasta'):
    SeqIO.write(record, f'{record.id}.fasta', 'fasta')
```

## Batch Convert

### Convert All Files in Directory
```python
from pathlib import Path

for gb_file in Path('genbank/').glob('*.gb'):
    fasta_file = Path('fasta/') / gb_file.with_suffix('.fasta').name
    count = SeqIO.convert(str(gb_file), 'genbank', str(fasta_file), 'fasta')
    print(f'{gb_file.name} -> {fasta_file.name}: {count} records')
```

### Batch Convert with Summary
```python
from pathlib import Path

results = []
for input_file in Path('input/').glob('*.gb'):
    output_file = Path('output/') / input_file.with_suffix('.fasta').name
    count = SeqIO.convert(str(input_file), 'genbank', str(output_file), 'fasta')
    results.append({'file': input_file.name, 'records': count})

print(f'Converted {len(results)} files, {sum(r["records"] for r in results)} total records')
```

## Parallel Processing

### Using multiprocessing
```python
from multiprocessing import Pool
from pathlib import Path

def process_file(filepath):
    records = list(SeqIO.parse(filepath, 'fasta'))
    return {'file': filepath.name, 'count': len(records), 'total_bp': sum(len(r.seq) for r in records)}

files = list(Path('data/').glob('*.fasta'))
with Pool(4) as pool:
    results = pool.map(process_file, files)

for r in results:
    print(f'{r["file"]}: {r["count"]} seqs, {r["total_bp"]} bp')
```

### Using concurrent.futures
```python
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

def count_records(filepath):
    return filepath.name, sum(1 for _ in SeqIO.parse(filepath, 'fasta'))

files = list(Path('data/').glob('*.fasta'))
with ThreadPoolExecutor(max_workers=4) as executor:
    results = executor.map(count_records, files)

for name, count in results:
    print(f'{name}: {count}')
```

## Summary Statistics

### Aggregate Stats Across Files
```python
from pathlib import Path

total_seqs = 0
total_bp = 0
file_count = 0

for fasta_file in Path('data/').glob('*.fasta'):
    for record in SeqIO.parse(fasta_file, 'fasta'):
        total_seqs += 1
        total_bp += len(record.seq)
    file_count += 1

print(f'Files: {file_count}')
print(f'Sequences: {total_seqs}')
print(f'Total bp: {total_bp}')
print(f'Average length: {total_bp / total_seqs:.0f}')
```

### Per-File Summary Report
```python
from pathlib import Path
import csv

summaries = []
for fasta_file in Path('data/').glob('*.fasta'):
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    lengths = [len(r.seq) for r in records]
    summaries.append({
        'file': fasta_file.name,
        'sequences': len(records),
        'total_bp': sum(lengths),
        'min_len': min(lengths) if lengths else 0,
        'max_len': max(lengths) if lengths else 0,
        'avg_len': sum(lengths) / len(lengths) if lengths else 0
    })

with open('summary.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=summaries[0].keys())
    writer.writeheader()
    writer.writerows(summaries)
```

## File Organization

### Organize by Criteria
```python
from pathlib import Path
from Bio.SeqUtils import gc_fraction

Path('high_gc').mkdir(exist_ok=True)
Path('low_gc').mkdir(exist_ok=True)

for fasta_file in Path('input/').glob('*.fasta'):
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    avg_gc = sum(gc_fraction(r.seq) for r in records) / len(records)

    if avg_gc >= 0.5:
        dest = Path('high_gc') / fasta_file.name
    else:
        dest = Path('low_gc') / fasta_file.name

    SeqIO.write(records, dest, 'fasta')
```

## Common Patterns

| Task | Approach |
|------|----------|
| Merge files | Generator yielding from each file |
| Split file | islice with batch size |
| Convert all | Loop with SeqIO.convert |
| Parallel processing | multiprocessing.Pool or ThreadPoolExecutor |
| Summary stats | Accumulate while iterating |

## Related Skills

- read-sequences - Core parsing functions for each file
- write-sequences - Write processed outputs
- sequence-statistics - Generate per-file statistics
- format-conversion - Batch format conversion
- compressed-files - Handle compressed files in batch
- database-access - Batch download sequences from NCBI
