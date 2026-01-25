# CNV Annotation Usage Guide

## Overview

CNV annotation adds biological context to copy number calls by identifying affected genes, pathways, and clinical significance. This is essential for interpreting CNV findings in research and clinical contexts.

## Prerequisites

```bash
conda install -c bioconda annotsv pybedtools
pip install pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate my CNV segments with overlapping genes and cancer gene flags"
- "Filter out common CNVs using population frequency data from DGV"
- "Run AnnotSV on my structural variant VCF for comprehensive annotation"
- "Add ClinVar pathogenicity annotations to my CNV calls"

## Annotation Resources

| Resource | Content | Format |
|----------|---------|--------|
| RefSeq/Ensembl | Gene coordinates | GTF/GFF |
| COSMIC CGC | Cancer genes | TSV |
| ClinVar | Clinical variants | VCF |
| DGV | Common CNVs | BED |
| gnomAD-SV | Population frequencies | VCF |
| OMIM | Disease associations | TSV |

## Complete Annotation Workflow

```python
import pandas as pd
import pybedtools
from pathlib import Path

class CNVAnnotator:
    def __init__(self, gene_gtf, cgc_file=None, clinvar_vcf=None):
        self.genes = self.load_genes(gene_gtf)
        self.cancer_genes = self.load_cgc(cgc_file) if cgc_file else set()
        self.clinvar = clinvar_vcf

    def load_genes(self, gtf_file):
        '''Load gene annotations from GTF.'''
        genes = []
        with open(gtf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] == 'gene':
                    attrs = dict(x.strip().split(' ', 1) for x in fields[8].split(';') if x.strip())
                    gene_name = attrs.get('gene_name', '').strip('"')
                    genes.append({
                        'chrom': fields[0],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'gene_name': gene_name
                    })
        return pd.DataFrame(genes)

    def load_cgc(self, cgc_file):
        '''Load COSMIC Cancer Gene Census.'''
        cgc = pd.read_csv(cgc_file, sep='\t')
        return set(cgc['Gene Symbol'].dropna().tolist())

    def annotate(self, cns_file):
        '''Annotate CNV segments.'''
        cns = pd.read_csv(cns_file, sep='\t')

        # BedTools intersection
        cns_bed = pybedtools.BedTool.from_dataframe(
            cns[['chromosome', 'start', 'end', 'log2']])
        genes_bed = pybedtools.BedTool.from_dataframe(
            self.genes[['chrom', 'start', 'end', 'gene_name']])

        intersect = cns_bed.intersect(genes_bed, wa=True, wb=True)

        # Parse results
        results = []
        for interval in intersect:
            results.append({
                'chrom': interval[0],
                'start': int(interval[1]),
                'end': int(interval[2]),
                'log2': float(interval[3]),
                'gene_name': interval[7]
            })

        df = pd.DataFrame(results)

        # Aggregate
        annotated = df.groupby(['chrom', 'start', 'end', 'log2']).agg({
            'gene_name': lambda x: ','.join(sorted(set(x)))
        }).reset_index()

        # Add cancer gene flags
        annotated['is_cancer_gene'] = annotated['gene_name'].apply(
            lambda x: any(g in self.cancer_genes for g in x.split(',')))

        annotated['cancer_genes'] = annotated['gene_name'].apply(
            lambda x: ','.join(g for g in x.split(',') if g in self.cancer_genes))

        # Add CNV type
        annotated['cnv_type'] = annotated['log2'].apply(self.classify_cnv)

        return annotated

    @staticmethod
    def classify_cnv(log2):
        if log2 < -1.0:
            return 'deep_deletion'
        elif log2 < -0.3:
            return 'deletion'
        elif log2 > 0.8:
            return 'high_amplification'
        elif log2 > 0.3:
            return 'amplification'
        return 'neutral'

# Usage
annotator = CNVAnnotator(
    gene_gtf='genes.gtf',
    cgc_file='cancer_gene_census.tsv'
)
annotated = annotator.annotate('sample.cns')
annotated.to_csv('sample_annotated.tsv', sep='\t', index=False)
```

## Using AnnotSV

AnnotSV provides comprehensive annotation with a single command:

```bash
# Install
conda install -c bioconda annotsv

# Download annotations (first time)
AnnotSV -SVinputFile dummy.vcf -genomeBuild GRCh38

# Annotate CNVs
AnnotSV \
    -SVinputFile sample.cnv.vcf \
    -genomeBuild GRCh38 \
    -outputFile sample_annotsv \
    -SVminSize 1000 \
    -annotationMode full

# Key output columns:
# - Gene_name: overlapping genes
# - DGV_GAIN/LOSS_Freq: population frequency
# - ClinVar_clndn: disease associations
# - ACMG_class: pathogenicity classification
```

## Cancer-Specific Annotation

```python
def annotate_for_cancer(cnv_annotated, cancer_type='breast'):
    '''Add cancer-specific annotations.'''
    # Known driver genes by cancer type
    cancer_drivers = {
        'breast': ['BRCA1', 'BRCA2', 'TP53', 'ERBB2', 'PIK3CA', 'MYC'],
        'lung': ['EGFR', 'KRAS', 'ALK', 'MET', 'ROS1', 'TP53'],
        'colon': ['APC', 'KRAS', 'TP53', 'SMAD4', 'PIK3CA'],
    }

    drivers = set(cancer_drivers.get(cancer_type, []))

    cnv_annotated['driver_genes'] = cnv_annotated['gene_name'].apply(
        lambda x: ','.join(g for g in x.split(',') if g in drivers))

    cnv_annotated['has_driver'] = cnv_annotated['driver_genes'].str.len() > 0

    # Prioritize
    cnv_annotated['priority'] = cnv_annotated.apply(
        lambda x: 'high' if x['has_driver'] or abs(x['log2']) > 1
            else ('medium' if x['is_cancer_gene'] else 'low'), axis=1)

    return cnv_annotated
```

## Population Frequency Filtering

```python
def filter_by_frequency(cnv_annotated, dgv_file, max_freq=0.01):
    '''Filter common CNVs using DGV data.'''
    dgv = pybedtools.BedTool(dgv_file)

    cnv_bed = pybedtools.BedTool.from_dataframe(
        cnv_annotated[['chrom', 'start', 'end']])

    # Find CNVs overlapping common variants (50% reciprocal overlap)
    common = cnv_bed.intersect(dgv, f=0.5, r=True, u=True)
    common_coords = set([(str(i[0]), int(i[1]), int(i[2])) for i in common])

    cnv_annotated['is_common'] = cnv_annotated.apply(
        lambda x: (x['chrom'], x['start'], x['end']) in common_coords, axis=1)

    # Filter
    rare_cnvs = cnv_annotated[~cnv_annotated['is_common']]

    return rare_cnvs
```

## Generate Summary Report

```python
def generate_report(annotated, output_prefix):
    '''Generate annotation summary report.'''
    # Statistics
    stats = {
        'Total CNVs': len(annotated),
        'Amplifications': (annotated['cnv_type'].str.contains('amplification')).sum(),
        'Deletions': (annotated['cnv_type'].str.contains('deletion')).sum(),
        'Affecting cancer genes': annotated['is_cancer_gene'].sum(),
        'Total genes affected': annotated['gene_name'].str.split(',').explode().nunique(),
    }

    with open(f'{output_prefix}_summary.txt', 'w') as f:
        for k, v in stats.items():
            f.write(f'{k}: {v}\n')

    # High-priority CNVs
    high_priority = annotated[annotated['priority'] == 'high']
    high_priority.to_csv(f'{output_prefix}_high_priority.tsv', sep='\t', index=False)

    # Gene list for pathway analysis
    genes = annotated['gene_name'].str.split(',').explode().unique()
    with open(f'{output_prefix}_genes.txt', 'w') as f:
        f.write('\n'.join(sorted(genes)))

    return stats
```

## Example Prompts

> "Annotate my CNV segments with overlapping genes and cancer gene census flags"

> "Filter out common CNVs from my callset using DGV population data"

> "Add ClinVar annotations to my CNV calls and prioritize pathogenic variants"

> "Run AnnotSV on my structural variant VCF to get comprehensive annotations"

## Integration with Pathway Analysis

```r
library(clusterProfiler)

# Load annotated CNV genes
amp_genes <- read.delim('cnv_annotated.tsv')
amp_genes <- amp_genes[amp_genes$log2 > 0.5, ]
genes <- unique(unlist(strsplit(amp_genes$gene_name, ',')))

# Convert to Entrez
entrez <- bitr(genes, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')

# Pathway enrichment
kegg <- enrichKEGG(gene=entrez$ENTREZID, organism='hsa')
go <- enrichGO(gene=entrez$ENTREZID, OrgDb='org.Hs.eg.db', ont='BP')

# Visualize
dotplot(kegg, showCategory=15)
dotplot(go, showCategory=15)
```
