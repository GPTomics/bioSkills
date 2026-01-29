---
name: bio-clinical-databases-hla-typing
description: Call HLA alleles from NGS data using OptiType, HLA-HD, or arcasHLA for immunogenomics applications. Use when determining HLA genotype for transplant matching, neoantigen prediction, or pharmacogenomic screening.
tool_type: cli
primary_tool: OptiType
---

# HLA Typing

## OptiType (HLA Class I)

### From DNA-seq

```bash
# Extract HLA reads from BAM
samtools view -h input.bam chr6:28000000-34000000 | \
    samtools fastq -1 hla_R1.fq -2 hla_R2.fq -

# Run OptiType
OptiTypePipeline.py \
    -i hla_R1.fq hla_R2.fq \
    -d \
    -o optitype_output \
    -c config.ini

# Output: optitype_output/sample_result.tsv
# Contains HLA-A, HLA-B, HLA-C alleles (4-field resolution)
```

### From RNA-seq

```bash
# RNA mode
OptiTypePipeline.py \
    -i rna_R1.fq rna_R2.fq \
    -r \
    -o optitype_rna_output \
    -c config.ini
```

### OptiType Config

```ini
# config.ini
[mapping]
razers3=/path/to/razers3
threads=4

[ilp]
solver=glpk
threads=4

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
```

## HLA-HD (Full Resolution)

```bash
# HLA-HD for high-resolution typing
# Supports Class I and Class II

# Extract HLA reads
samtools view -b input.bam chr6:28000000-34000000 > hla_region.bam
samtools sort -n hla_region.bam -o hla_sorted.bam
samtools fastq -1 hla_R1.fq -2 hla_R2.fq hla_sorted.bam

# Run HLA-HD
hlahd.sh \
    -t 8 \
    -m 100 \
    -f freq_data \
    hla_R1.fq \
    hla_R2.fq \
    gene_split_filt \
    dictionary \
    sample_name \
    output_dir

# Output includes HLA-A, B, C, DRB1, DQB1, DPB1 at 4-field resolution
```

## arcasHLA (RNA-seq)

```bash
# Fast HLA typing from RNA-seq
# Extracts and genotypes in one step

# From STAR-aligned BAM
arcasHLA extract sample.bam -o output_dir
arcasHLA genotype output_dir/sample.extracted.fq.gz -o output_dir

# Output: sample.genotype.json
# {
#   "A": ["A*02:01", "A*24:02"],
#   "B": ["B*35:01", "B*44:03"],
#   "C": ["C*04:01", "C*05:01"]
# }
```

### arcasHLA Merge

```bash
# Merge multiple samples
arcasHLA merge output_dir/*.genotype.json -o merged_hla.tsv
```

## HLA Nomenclature

```
HLA-A*02:01:01:01
     |  |  |  |
     |  |  |  +-- Non-coding variation (optional)
     |  |  +----- Synonymous variation (optional)
     |  +-------- Protein sequence (usually reported)
     +----------- Allele group

Resolution levels:
- 2-field: A*02:01 (protein sequence - clinical standard)
- 4-field: A*02:01:01 (includes synonymous changes)
- Full: A*02:01:01:01 (includes non-coding)
```

## HLA and Pharmacogenomics

```python
# Key HLA-drug associations
HLA_DRUG_ASSOCIATIONS = {
    'B*57:01': {
        'drug': 'Abacavir',
        'reaction': 'Hypersensitivity syndrome',
        'screening': 'Required before prescribing'
    },
    'B*15:02': {
        'drug': 'Carbamazepine',
        'reaction': 'SJS/TEN',
        'populations': 'High risk in Han Chinese, Southeast Asian'
    },
    'B*58:01': {
        'drug': 'Allopurinol',
        'reaction': 'SJS/TEN',
        'populations': 'High risk in Han Chinese, Korean, Thai'
    },
    'A*31:01': {
        'drug': 'Carbamazepine',
        'reaction': 'DRESS',
        'populations': 'European, Japanese'
    }
}

def check_hla_drug_risk(hla_alleles, drug):
    '''Check if patient HLA poses drug reaction risk'''
    risks = []
    for allele in hla_alleles:
        if allele in HLA_DRUG_ASSOCIATIONS:
            assoc = HLA_DRUG_ASSOCIATIONS[allele]
            if assoc['drug'].lower() == drug.lower():
                risks.append({
                    'allele': allele,
                    'drug': drug,
                    'reaction': assoc['reaction']
                })
    return risks
```

## Parse OptiType Results

```python
import pandas as pd

def parse_optitype(result_file):
    '''Parse OptiType TSV output'''
    df = pd.read_csv(result_file, sep='\t')
    # Columns: A1, A2, B1, B2, C1, C2, Reads, Objective
    hla_calls = {
        'HLA-A': [df['A1'].iloc[0], df['A2'].iloc[0]],
        'HLA-B': [df['B1'].iloc[0], df['B2'].iloc[0]],
        'HLA-C': [df['C1'].iloc[0], df['C2'].iloc[0]]
    }
    return hla_calls

def format_hla_report(hla_calls):
    '''Format HLA calls for clinical report'''
    report = []
    for gene, alleles in hla_calls.items():
        allele_str = '/'.join(sorted(set(alleles)))
        report.append(f'{gene}: {allele_str}')
    return '\n'.join(report)
```

## Class I vs Class II

| Class | Genes | Function | Typing Priority |
|-------|-------|----------|-----------------|
| Class I | HLA-A, B, C | Present intracellular peptides | Neoantigen, PGx |
| Class II | HLA-DR, DQ, DP | Present extracellular peptides | Transplant, autoimmune |

## Tool Comparison

| Tool | Input | Classes | Resolution | Speed |
|------|-------|---------|------------|-------|
| OptiType | WGS/WES/RNA | I only | 4-field | Fast |
| HLA-HD | WGS/WES | I and II | 4-field | Moderate |
| arcasHLA | RNA-seq | I and II | 4-field | Fast |
| HLA-LA | WGS | I and II | 4-field | Slow |

## Related Skills

- clinical-databases/pharmacogenomics - HLA-drug interactions
- variant-calling/clinical-interpretation - Clinical reporting
- single-cell/cell-type-annotation - HLA expression
