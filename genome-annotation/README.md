# genome-annotation

## Overview

Annotate assembled genomes with repeat masking, gene prediction (prokaryotic and eukaryotic), non-coding RNA detection, functional assignment, completeness QC, and annotation transfer between assemblies. The category encodes the decision-grade reality that gene calling is largely solved while trustworthy *function*, the right masking choice, and honest *completeness* are not.

**Tool type:** cli | **Primary tools:** Bakta, BRAKER3, RepeatModeler2, Infernal, eggNOG-mapper, BUSCO, Liftoff

## Skills

| Skill | Description |
|-------|-------------|
| repeat-annotation | De novo TE library (RepeatModeler2/EDTA/EarlGrey) and soft-masking (RepeatMasker), the prerequisite for gene prediction |
| prokaryotic-annotation | Bacterial/archaeal annotation with Bakta or Prokka, with PGAP/DFAST for submission |
| eukaryotic-gene-prediction | Evidence-first gene prediction with BRAKER3, GALBA, Funannotate, or deep-learning ab initio |
| ncrna-annotation | Structure-aware ncRNA detection with Infernal/Rfam, tRNAscan-SE, and barrnap |
| functional-annotation | GO/KEGG/Pfam/EC assignment with eggNOG-mapper and InterProScan, tiered by evidence |
| annotation-qc | Completeness and sanity with BUSCO, OMArk, CheckM2, and a gene-set sanity panel |
| annotation-transfer | Liftover (liftOver/CrossMap) and projection (Liftoff/miniprot/TOGA) between assemblies |

## Example Prompts

- "Mask repeats in my assembly before gene prediction"
- "Annotate my bacterial genome assembly with Bakta and check coding density"
- "Run BRAKER3 on my eukaryotic assembly with RNA-seq and protein evidence"
- "Find all tRNAs, rRNAs, and other ncRNAs in my genome"
- "Add GO/KEGG/Pfam functional annotations to my predicted proteins"
- "Is my annotation good enough to publish - run BUSCO on the proteome and the genome"
- "Transfer annotations from the reference to my new assembly and validate them"
- "Predict genes in my fungal genome with Funannotate including UTRs and isoforms"

## Requirements

```bash
# Repeat annotation
conda install -c bioconda repeatmodeler repeatmasker edta earlgrey

# Prokaryotic annotation
conda install -c bioconda bakta prokka

# Eukaryotic gene prediction
conda install -c bioconda braker3 galba funannotate augustus

# Non-coding RNA
conda install -c bioconda infernal trnascan-se barrnap aragorn

# Functional annotation
conda install -c bioconda eggnog-mapper=2.1.15 interproscan kofamscan agat

# Annotation QC
conda install -c bioconda busco omark checkm2 compleasm

# Annotation transfer
conda install -c bioconda miniprot ucsc-liftover crossmap
pip install liftoff lifton

# Python utilities
pip install gffutils biopython pandas
```

## Related Skills

- **genome-assembly** - Assemble and QC genomes (purge haplotigs) before annotation
- **genome-intervals** - Work with GFF/GTF annotation files
- **comparative-genomics** - Ortholog, synteny, and pangenome analysis across species
- **pathway-analysis** - GO and KEGG enrichment downstream of functional annotation
- **rna-quantification** - Quantify expression from annotated genes
