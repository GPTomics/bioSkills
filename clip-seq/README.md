# clip-seq

## Overview

Analyze CLIP-seq data (CLIP, PAR-CLIP, iCLIP, eCLIP) to identify protein-RNA binding sites at nucleotide resolution for understanding post-transcriptional regulation.

**Tool type:** mixed | **Primary tools:** CLIPper, Piranha, umi_tools, HOMER

## Skills

| Skill | Description |
|-------|-------------|
| clip-preprocessing | Adapter trimming, UMI extraction, and deduplication |
| clip-alignment | Alignment with crosslink site awareness |
| clip-peak-calling | Call binding site peaks or clusters |
| binding-site-annotation | Annotate peaks to genomic features |
| clip-motif-analysis | De novo and known motif enrichment |

## Example Prompts

- "Process my eCLIP data from FASTQ to peaks"
- "Extract UMIs and deduplicate my CLIP reads"
- "Call binding sites with CLIPper"
- "Annotate peaks to 3'UTR, CDS, introns"
- "Find enriched RBP motifs at binding sites"
- "Compare binding between conditions"

## Requirements

```bash
# UMI handling
pip install umi_tools

# Peak calling
conda install -c bioconda clipper piranha

# Motif analysis
conda install -c bioconda homer meme
```

## Related Skills

- **read-qc** - UMI processing concepts
- **chip-seq** - Peak calling concepts
- **genome-intervals** - Peak annotation
