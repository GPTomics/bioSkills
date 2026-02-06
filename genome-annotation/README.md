# genome-annotation

## Overview

Annotate assembled genomes with gene predictions, functional assignments, repeat elements, non-coding RNAs, and annotation transfer between assemblies.

**Tool type:** cli | **Primary tools:** Bakta, BRAKER3, eggNOG-mapper, InterProScan, RepeatMasker, Liftoff

## Skills

| Skill | Description |
|-------|-------------|
| prokaryotic-annotation | Bacterial/archaeal annotation with Bakta or Prokka |
| eukaryotic-gene-prediction | Gene prediction with BRAKER3, GALBA, Augustus |
| functional-annotation | GO, KEGG, Pfam assignment with eggNOG-mapper and InterProScan |
| ncrna-annotation | Non-coding RNA detection with Infernal and tRNAscan-SE |
| repeat-annotation | Repeat and TE annotation with RepeatModeler and RepeatMasker |
| annotation-transfer | Liftover annotations between assemblies with Liftoff and MiniProt |

## Example Prompts

- "Annotate my bacterial genome assembly with Bakta"
- "Run BRAKER3 on my eukaryotic assembly with RNA-seq evidence"
- "Add functional annotations to my predicted proteins"
- "Find all tRNAs and rRNAs in my genome"
- "Mask repeats in my assembly before gene prediction"
- "Transfer annotations from the reference to my new assembly"
- "Predict genes in my fungal genome using protein evidence"
- "Annotate transposable elements and quantify TE expression"

## Requirements

```bash
# Prokaryotic annotation
conda install -c bioconda bakta prokka

# Eukaryotic gene prediction
conda install -c bioconda braker3 augustus galba

# Functional annotation
conda install -c bioconda eggnog-mapper interproscan

# Non-coding RNA
conda install -c bioconda infernal trnascan-se barrnap

# Repeat annotation
conda install -c bioconda repeatmodeler repeatmasker

# Annotation transfer
conda install -c bioconda liftoff miniprot

# Python utilities
pip install gffutils biopython pandas
```

## Related Skills

- **genome-assembly** - Assemble genomes before annotation
- **genome-intervals** - Work with GFF/GTF annotation files
- **comparative-genomics** - Ortholog and synteny analysis across species
- **rna-quantification** - Quantify expression from annotated genes
