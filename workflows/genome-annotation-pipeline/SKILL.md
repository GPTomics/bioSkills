---
name: bio-workflows-genome-annotation-pipeline
description: End-to-end genome annotation pipeline from assembled contigs to functional annotation, covering repeat masking, gene prediction, and functional assignment for both prokaryotic and eukaryotic genomes. Use when annotating a newly assembled genome from scratch.
tool_type: mixed
primary_tool: Bakta
workflow: true
depends_on:
  - genome-annotation/prokaryotic-annotation
  - genome-annotation/eukaryotic-gene-prediction
  - genome-annotation/functional-annotation
  - genome-assembly/assembly-qc
qc_checkpoints:
  - after_repeat_masking: "Repeat content within expected range for taxon"
  - after_gene_prediction: "Gene count plausible, BUSCO completeness >90%"
  - after_functional_annotation: ">60% of genes with functional assignment"
---

# Genome Annotation Pipeline

Complete workflow from assembled contigs to functional annotation for prokaryotic or eukaryotic genomes.

## Pipeline Overview

```
Assembled contigs
    |
    v
[0. Assembly QC] ----------> QUAST, BUSCO (confirm assembly quality)
    |
    +----- Prokaryotic? -----> Path A: Bakta (one-step annotation)
    |                                |
    |                                v
    |                          Annotated genome (GFF3, GenBank, FASTA)
    |
    +----- Eukaryotic? ------> Path B: Multi-step pipeline
                                    |
                                    v
                              [1. Repeat Masking] ----> RepeatModeler + RepeatMasker
                                    |
                                    v
                              [2. Gene Prediction] ---> BRAKER3 (RNA-seq + protein evidence)
                                    |
                                    v
                              [3. Functional Annotation] -> eggNOG-mapper + InterProScan
                                    |
                                    v
                              [4. ncRNA Annotation] ---> Infernal + tRNAscan-SE
                                    |
                                    v
                              Annotated genome (GFF3, proteins, functional tables)
```

## Path A: Prokaryotic Annotation (Bakta)

Bakta provides comprehensive one-step annotation for bacteria and archaea. Preferred over Prokka for new projects.

### Database Setup

```bash
bakta_db download --output /path/to/bakta_db --type full
```

### Run Bakta

```bash
bakta \
    --db /path/to/bakta_db \
    --output bakta_out \
    --prefix my_genome \
    --locus-tag MYORG \
    --genus Escherichia --species "coli" \
    --strain K12 \
    --gram - \
    --complete \
    --threads 8 \
    assembly.fasta
```

### Prokaryotic QC Checkpoint

```python
import subprocess
import json

def validate_prokaryotic_annotation(bakta_dir, prefix, expected_cds_range=(500, 8000)):
    '''
    QC gates for prokaryotic annotation.
    - CDS count in expected range for genome size
    - tRNA count >= 20 (typical minimum for free-living bacteria)
    - rRNA operons detected
    '''
    gff_file = f'{bakta_dir}/{prefix}.gff3'

    feature_counts = {'CDS': 0, 'tRNA': 0, 'rRNA': 0, 'tmRNA': 0, 'ncRNA': 0}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3 and fields[2] in feature_counts:
                feature_counts[fields[2]] += 1

    qc_pass = True
    if not (expected_cds_range[0] <= feature_counts['CDS'] <= expected_cds_range[1]):
        print(f'WARNING: CDS count {feature_counts["CDS"]} outside expected range {expected_cds_range}')
        qc_pass = False
    if feature_counts['tRNA'] < 20:
        print(f'WARNING: Only {feature_counts["tRNA"]} tRNAs detected (expect >= 20)')
        qc_pass = False

    print(f'Feature summary: {feature_counts}')
    return qc_pass, feature_counts
```

## Path B: Eukaryotic Annotation

### Step 1: Repeat Masking

```bash
# Build species-specific repeat library
RepeatModeler -database mygenome -pa 8 -LTRStruct

# Combine with known repeats
cat mygenome-families.fa /path/to/RepeatMasker/Libraries/RepeatMaskerLib.h5 > combined_lib.fa

# Mask the genome
RepeatMasker \
    -lib combined_lib.fa \
    -pa 8 \
    -xsmall \
    -gff \
    -dir repeat_out \
    assembly.fasta
```

#### Repeat Masking QC Checkpoint

```python
def check_repeat_content(repeatmasker_tbl, taxon='vertebrate'):
    '''
    Verify repeat content is within expected range for taxon.
    Typical ranges:
    - Vertebrate: 30-60%
    - Insect: 15-45%
    - Plant: 20-85%
    - Fungus: 3-20%
    '''
    expected_ranges = {
        'vertebrate': (30, 60), 'insect': (15, 45),
        'plant': (20, 85), 'fungus': (3, 20)
    }
    low, high = expected_ranges.get(taxon, (5, 80))

    with open(repeatmasker_tbl) as f:
        for line in f:
            if 'total interspersed' in line.lower():
                pct = float(line.strip().split()[-1].replace('%', ''))
                break

    qc_pass = low <= pct <= high
    if not qc_pass:
        print(f'WARNING: Repeat content {pct:.1f}% outside expected range ({low}-{high}%) for {taxon}')
    return qc_pass, pct
```

### Step 2: Gene Prediction with BRAKER3

```bash
# BRAKER3 combines GeneMark-ETP, AUGUSTUS, and TSEBRA
# Uses both RNA-seq and protein evidence for best results
braker.pl \
    --genome=assembly.fasta.masked \
    --bam=rnaseq_sorted.bam \
    --prot_seq=proteins.fa \
    --softmasking \
    --threads 8 \
    --species=my_species \
    --gff3 \
    --workingdir=braker_out

# If only RNA-seq evidence available
braker.pl \
    --genome=assembly.fasta.masked \
    --bam=rnaseq_sorted.bam \
    --softmasking \
    --threads 8 \
    --species=my_species \
    --gff3

# If only protein evidence available (use OrthoDB proteins)
braker.pl \
    --genome=assembly.fasta.masked \
    --prot_seq=orthodb_proteins.fa \
    --softmasking \
    --threads 8 \
    --species=my_species \
    --gff3
```

#### Gene Prediction QC Checkpoint

```bash
# BUSCO completeness check on predicted proteins
busco \
    -i braker_out/braker.aa \
    -l eukaryota_odb10 \
    -o busco_annotation \
    -m protein \
    --cpu 8
```

```python
def check_gene_prediction(braker_gff, busco_summary, expected_genes_range=(15000, 35000)):
    '''
    QC gates after gene prediction.
    - Gene count within expected range for genome
    - BUSCO completeness > 90%
    - Mean exons per gene > 1 (spliced genes expected in eukaryotes)
    '''
    gene_count = 0
    exon_count = 0
    with open(braker_gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            feature = line.strip().split('\t')[2] if len(line.strip().split('\t')) >= 3 else ''
            if feature == 'gene':
                gene_count += 1
            elif feature == 'exon':
                exon_count += 1

    mean_exons = exon_count / gene_count if gene_count > 0 else 0

    with open(busco_summary) as f:
        for line in f:
            if line.strip().startswith('C:'):
                completeness = float(line.strip().split('C:')[1].split('%')[0])
                break

    issues = []
    if not (expected_genes_range[0] <= gene_count <= expected_genes_range[1]):
        issues.append(f'Gene count {gene_count} outside expected range {expected_genes_range}')
    if completeness < 90:
        issues.append(f'BUSCO completeness {completeness:.1f}% < 90%')
    if mean_exons < 2:
        issues.append(f'Mean exons/gene {mean_exons:.1f} is low for eukaryote')

    print(f'Genes: {gene_count}, Mean exons/gene: {mean_exons:.1f}, BUSCO: {completeness:.1f}%')
    return len(issues) == 0, issues
```

### Step 3: Functional Annotation

```bash
# eggNOG-mapper for comprehensive functional annotation
emapper.py \
    -i braker_out/braker.aa \
    --output eggnog_results \
    --cpu 8 \
    -m diamond \
    --tax_scope auto \
    --go_evidence non-electronic \
    --target_orthologs all \
    --seed_ortholog_evalue 1e-5 \
    --override

# InterProScan for domain annotation (complementary to eggNOG)
interproscan.sh \
    -i braker_out/braker.aa \
    -o interpro_results.tsv \
    -f tsv,gff3 \
    -goterms \
    -pa \
    -cpu 8
```

#### Functional Annotation QC Checkpoint

```python
import pandas as pd

def check_functional_annotation(eggnog_annotations, total_genes):
    '''
    QC gate: > 60% of genes should have functional assignment.
    Below 50% suggests database issues or highly divergent organism.
    '''
    df = pd.read_csv(eggnog_annotations, sep='\t', comment='#')
    annotated = len(df[df['Description'] != '-'])
    pct_annotated = annotated / total_genes * 100

    has_go = len(df[df['GOs'] != '-'])
    has_kegg = len(df[df['KEGG_ko'] != '-'])

    print(f'Annotated: {annotated}/{total_genes} ({pct_annotated:.1f}%)')
    print(f'With GO terms: {has_go}, With KEGG: {has_kegg}')

    if pct_annotated < 60:
        print('WARNING: <60% annotated. Check database version or use broader taxonomy scope.')
    return pct_annotated >= 60
```

### Step 4: ncRNA Annotation

```bash
# tRNAscan-SE for tRNA genes
tRNAscan-SE \
    -E \
    --thread 8 \
    -o trna_results.txt \
    --gff trna.gff \
    assembly.fasta

# Infernal for Rfam-based ncRNA annotation
# Download Rfam covariance models first
cmscan \
    --cpu 8 \
    --tblout rfam_results.tbl \
    --fmt 2 \
    --clanin Rfam.clanin \
    Rfam.cm \
    assembly.fasta
```

## Merging Annotations

```python
def merge_annotations(braker_gff, trna_gff, rfam_tbl, eggnog_tsv, output_gff):
    '''Merge gene predictions, ncRNAs, and functional annotations into final GFF3.'''
    import subprocess

    # Use AGAT for GFF merging and validation
    subprocess.run([
        'agat_sp_merge_annotations.pl',
        '--gff', braker_gff,
        '--gff', trna_gff,
        '-o', output_gff
    ], check=True)

    # Validate final GFF3
    subprocess.run([
        'agat_sp_statistics.pl',
        '--gff', output_gff,
        '-o', output_gff.replace('.gff3', '_stats.txt')
    ], check=True)

    print(f'Merged annotations written to {output_gff}')
```

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| Low gene count | Repeat masking too aggressive | Use `-xsmall` (soft-masking) not `-x` (hard-masking) |
| BUSCO < 80% | Poor assembly or missing evidence | Add RNA-seq data; check assembly contiguity |
| Many partial genes | Fragmented assembly | Scaffold first; use `--min_contig` in BRAKER |
| < 60% annotated | Divergent organism | Use broader `--tax_scope`; try InterProScan |
| Too many genes | Gene prediction artifacts | Increase `--min_intron_len`; filter short ORFs |
| Missing ncRNAs | Wrong Rfam models | Verify Rfam version matches genome build |

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

GENOME="assembly.fasta"
RNASEQ_BAM="rnaseq_sorted.bam"
PROTEINS="orthodb_proteins.fa"
BAKTA_DB="/path/to/bakta_db"
THREADS=8

# Determine organism type
ORGANISM_TYPE="${1:-eukaryotic}"  # prokaryotic or eukaryotic

if [ "$ORGANISM_TYPE" == "prokaryotic" ]; then
    echo "Running prokaryotic annotation with Bakta"
    bakta --db $BAKTA_DB --output bakta_out --prefix genome \
          --locus-tag MYORG --threads $THREADS $GENOME
    echo "Done. Results in bakta_out/"

else
    echo "Running eukaryotic annotation pipeline"

    echo "Step 1: Repeat masking"
    RepeatModeler -database mygenome -pa $THREADS -LTRStruct
    RepeatMasker -lib mygenome-families.fa -pa $THREADS -xsmall -gff -dir repeat_out $GENOME

    echo "Step 2: Gene prediction with BRAKER3"
    braker.pl --genome=repeat_out/$(basename $GENOME).masked \
              --bam=$RNASEQ_BAM --prot_seq=$PROTEINS \
              --softmasking --threads $THREADS --gff3 --workingdir=braker_out

    echo "Step 3: BUSCO QC"
    busco -i braker_out/braker.aa -l eukaryota_odb10 -o busco_check -m protein --cpu $THREADS

    echo "Step 4: Functional annotation"
    emapper.py -i braker_out/braker.aa --output eggnog_out --cpu $THREADS -m diamond

    echo "Step 5: ncRNA annotation"
    tRNAscan-SE -E --thread $THREADS -o trna_out.txt --gff trna.gff $GENOME
    cmscan --cpu $THREADS --tblout rfam.tbl --fmt 2 Rfam.cm $GENOME

    echo "Done. Check braker_out/, eggnog_out*, trna.gff, rfam.tbl"
fi
```

## Related Skills

- genome-annotation/prokaryotic-annotation - Bakta and Prokka details
- genome-annotation/eukaryotic-gene-prediction - BRAKER3 and AUGUSTUS options
- genome-annotation/functional-annotation - eggNOG-mapper and InterProScan
- genome-assembly/assembly-qc - Pre-annotation assembly quality checks
- genome-intervals/gtf-gff-handling - GFF3/GTF parsing and manipulation
