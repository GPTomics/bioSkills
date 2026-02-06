---
name: bio-genome-annotation-prokaryotic-annotation
description: Annotate bacterial and archaeal genomes with Bakta for comprehensive structural and functional annotation, or Prokka for lightweight annotation. Generates GFF3, GenBank, and FASTA outputs with NCBI-compatible locus tags. Use when annotating a newly assembled prokaryotic genome or preparing annotations for NCBI submission.
tool_type: cli
primary_tool: Bakta
---

# Prokaryotic Genome Annotation

Annotate prokaryotic genomes with Bakta (preferred) or Prokka (legacy). Bakta provides more comprehensive functional annotation through up-to-date databases and NCBI-compatible output formatting.

## Bakta

### Database Setup

```bash
# Download the full database (~30 GB, recommended for comprehensive annotation)
bakta_db download --output /path/to/bakta_db --type full

# Lightweight database (~1.5 GB, faster but less comprehensive)
bakta_db download --output /path/to/bakta_db --type light

# Update existing database
bakta_db update --db /path/to/bakta_db
```

### Basic Annotation

```bash
bakta \
    --db /path/to/bakta_db \
    --output bakta_out \
    --prefix my_genome \
    --locus-tag MYORG \
    --threads 8 \
    assembly.fasta
```

### Key Options

| Option | Description |
|--------|-------------|
| `--db` | Path to Bakta database |
| `--output` | Output directory |
| `--prefix` | Output file prefix |
| `--locus-tag` | NCBI-compatible locus tag prefix |
| `--genus` / `--species` | Organism taxonomy |
| `--strain` | Strain designation |
| `--complete` | Flag for complete genomes (enables oriC/oriV detection) |
| `--gram` | Gram type (+ or -) for signal peptide prediction |
| `--threads` | CPU threads |
| `--min-contig-length` | Minimum contig length to annotate (default: 1) |
| `--translation-table` | Genetic code (default: 11 for bacteria) |

### With Organism Metadata

```bash
bakta \
    --db /path/to/bakta_db \
    --output bakta_out \
    --prefix ecoli_k12 \
    --locus-tag ECK12 \
    --genus Escherichia --species coli --strain K-12 \
    --gram - \
    --complete \
    --threads 16 \
    assembly.fasta
```

### Output Files

```
bakta_out/
├── my_genome.gff3       # GFF3 annotation (primary output)
├── my_genome.gbff       # GenBank format
├── my_genome.ffn        # Nucleotide CDS sequences
├── my_genome.faa        # Protein sequences
├── my_genome.fna        # Annotated genome sequence
├── my_genome.embl       # EMBL format
├── my_genome.tsv        # Tab-separated feature table
├── my_genome.json       # Machine-readable JSON
└── my_genome.txt        # Summary statistics
```

## Prokka (Legacy Alternative)

Prokka is lighter weight and faster but uses older databases. Prefer Bakta for new projects.

```bash
prokka \
    --outdir prokka_out \
    --prefix my_genome \
    --locustag MYORG \
    --genus Escherichia --species coli \
    --cpus 8 \
    --rfam \
    assembly.fasta
```

### Prokka vs Bakta

| Feature | Bakta | Prokka |
|---------|-------|--------|
| Database updates | Active (2024+) | Unmaintained since 2021 |
| Functional annotation | Comprehensive (UniProt, COG, Pfam) | Basic (UniProt) |
| ncRNA detection | Infernal + Rfam 14.x | Infernal + Rfam 12.x |
| NCBI compatibility | Full SQN output | Requires tbl2asn |
| Speed | Moderate | Fast |

## Parsing Annotations with Python

```python
import gffutils

def load_annotation(gff_file):
    '''Load GFF3 into a queryable database.'''
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy='merge')
    return db

def extract_cds_features(db):
    '''Extract all CDS features with product annotations.'''
    features = []
    for cds in db.features_of_type('CDS'):
        features.append({
            'id': cds.id,
            'seqid': cds.seqid,
            'start': cds.start,
            'end': cds.end,
            'strand': cds.strand,
            'product': cds.attributes.get('product', ['unknown'])[0],
            'locus_tag': cds.attributes.get('locus_tag', [''])[0]
        })
    return features

def compute_coding_density(db, genome_length):
    '''Compute fraction of genome encoding proteins.

    Typical prokaryotic coding density: 85-95%.
    Values below 80% may indicate pseudogenes or annotation gaps.
    Values above 95% may indicate overlapping annotations.
    '''
    coding_bp = sum(cds.end - cds.start + 1 for cds in db.features_of_type('CDS'))
    return coding_bp / genome_length

db = load_annotation('bakta_out/my_genome.gff3')
cds_features = extract_cds_features(db)
print(f'Total CDSs: {len(cds_features)}')
```

## Annotation QC

### Expected Metrics by Genome Size

| Genome Size | Expected Genes | Coding Density |
|-------------|---------------|----------------|
| 1-2 Mb | 900-2,000 | 85-92% |
| 2-5 Mb | 1,800-5,000 | 85-90% |
| 5-10 Mb | 4,500-9,000 | 82-88% |

### QC Checks

```bash
# Count annotated features
grep -c $'\tCDS\t' bakta_out/my_genome.gff3
grep -c $'\ttRNA\t' bakta_out/my_genome.gff3
grep -c $'\trRNA\t' bakta_out/my_genome.gff3

# Check for hypothetical proteins (ideally <40% of total CDSs)
grep -c 'hypothetical protein' bakta_out/my_genome.tsv
```

### BUSCO on Predicted Proteins

```bash
busco -i bakta_out/my_genome.faa -m proteins -l bacteria_odb10 -o busco_proteins
```

## Troubleshooting

### Low Gene Count
- Check assembly completeness with BUSCO (genome mode)
- Verify correct translation table (--translation-table 4 for Mycoplasma)
- Inspect minimum contig length filter

### Many Hypothetical Proteins
- Normal for novel organisms (30-50% is common)
- Try running InterProScan on the .faa file for additional annotations
- Consider eggNOG-mapper for orthology-based functional assignment

### NCBI Submission
- Use `--compliant` flag for NCBI-ready output
- Ensure locus tags follow NCBI format (3-12 uppercase alphanumeric)
- Review .tsv output for annotation warnings

## Related Skills

- functional-annotation - Add GO/KEGG/Pfam to predicted proteins
- ncrna-annotation - Detailed ncRNA identification with Infernal
- genome-assembly/assembly-qc - Assess assembly quality before annotation
- genome-intervals/gtf-gff-handling - Parse and manipulate GFF3 output
