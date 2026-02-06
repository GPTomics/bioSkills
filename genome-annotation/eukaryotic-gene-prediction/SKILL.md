---
name: bio-genome-annotation-eukaryotic-gene-prediction
description: Predict protein-coding genes in eukaryotic genomes using BRAKER3 for combined RNA-seq and protein evidence, or GALBA for protein-only evidence. Runs Augustus with trained parameters for accurate gene models. Use when annotating a newly assembled eukaryotic genome or improving existing gene predictions.
tool_type: cli
primary_tool: BRAKER3
---

# Eukaryotic Gene Prediction

Predict protein-coding genes in eukaryotic genomes using evidence-based methods. BRAKER3 combines RNA-seq and protein homology evidence for the most accurate predictions. GALBA provides an alternative when only protein evidence is available.

## Prerequisites

**CRITICAL:** The input assembly must be softmasked (repeats in lowercase). Run repeat-annotation first to softmask the genome. Unmasked assemblies produce many false positive gene predictions in repetitive regions.

## BRAKER3 (RNA-seq + Protein Evidence)

BRAKER3 is the preferred pipeline when RNA-seq data exists. It combines GeneMark-ETP with Augustus training for high-accuracy gene models.

### Evidence Preparation

```bash
# Align RNA-seq with HISAT2
hisat2-build assembly_softmasked.fasta genome_index
hisat2 -x genome_index -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz \
    --dta -p 16 | samtools sort -@ 4 -o rnaseq_sorted.bam
samtools index rnaseq_sorted.bam

# Download OrthoDB proteins for the relevant clade
# Available: Metazoa, Viridiplantae, Fungi, Arthropoda, Vertebrata
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz
gunzip Viridiplantae.fa.gz
```

### Run BRAKER3

```bash
braker.pl \
    --genome=assembly_softmasked.fasta \
    --bam=rnaseq_sorted.bam \
    --prot_seq=Viridiplantae.fa \
    --softmasking \
    --threads=16 \
    --species=my_species \
    --workingdir=braker3_out \
    --gff3
```

### Key Options

| Option | Description |
|--------|-------------|
| `--genome` | Softmasked genome assembly (FASTA) |
| `--bam` | RNA-seq alignments (BAM, can specify multiple comma-separated) |
| `--prot_seq` | Protein evidence (FASTA) |
| `--softmasking` | Genome is softmasked (required) |
| `--species` | Species name for Augustus training |
| `--threads` | CPU threads |
| `--gff3` | Output in GFF3 format (default is GTF) |
| `--fungus` | Use fungal-specific parameters |
| `--workingdir` | Output directory |
| `--UTR` | Predict UTRs (requires sufficient RNA-seq coverage) |

### Output Files

```
braker3_out/
├── braker.gtf           # Gene predictions (GTF)
├── braker.gff3          # Gene predictions (GFF3, if --gff3)
├── braker.codingseq     # CDS nucleotide sequences
├── braker.aa            # Protein sequences
├── hintsfile.gff        # All evidence hints
├── Augustus/            # Trained Augustus parameters
└── GeneMark-ETP/        # GeneMark intermediate files
```

## GALBA (Protein-Only Evidence)

GALBA works best with closely related species proteins. Always prefer BRAKER3 when RNA-seq data is available, as intron evidence from RNA-seq substantially improves splice site accuracy.

```bash
galba.pl \
    --genome=assembly_softmasked.fasta \
    --prot_seq=closely_related_proteins.fa \
    --softmasking \
    --threads=16 \
    --species=my_species \
    --workingdir=galba_out \
    --gff3
```

### When to Use GALBA vs BRAKER3

| Scenario | Tool |
|----------|------|
| RNA-seq + proteins available | BRAKER3 |
| Proteins only, closely related species | GALBA |
| Proteins only, distant species | BRAKER3 --prot_seq only (uses ProtHint) |
| No external evidence | Augustus with pre-trained species (last resort) |

## Augustus Standalone

For genomes with pre-trained parameters or when rerunning predictions with different hints.

```bash
# List available pre-trained species
augustus --species=help

# Run with pre-trained species
augustus \
    --species=arabidopsis \
    --softmasking=1 \
    --gff3=on \
    --UTR=off \
    assembly_softmasked.fasta > augustus_predictions.gff3

# Run with hints (evidence)
augustus \
    --species=my_species \
    --softmasking=1 \
    --gff3=on \
    --hintsfile=hints.gff \
    --extrinsicCfgFile=extrinsic.cfg \
    assembly_softmasked.fasta > augustus_hints.gff3
```

## Evaluation with BUSCO

```bash
# Evaluate predicted proteins
busco -i braker3_out/braker.aa -m proteins -l embryophyta_odb10 -o busco_proteins -c 8

# Compare to genome-mode BUSCO
busco -i assembly_softmasked.fasta -m genome -l embryophyta_odb10 -o busco_genome -c 8
```

### Interpretation

| Metric | Meaning |
|--------|---------|
| Protein BUSCO > Genome BUSCO | Annotation captures most genes |
| Protein BUSCO < Genome BUSCO | Missing gene models, re-examine evidence |
| High duplication | Check for retained duplicates vs artifacts |

## Python: Parse Gene Models

```python
import gffutils
import pandas as pd

def load_gene_models(gff_file):
    '''Load eukaryotic gene models from GFF3/GTF.'''
    db = gffutils.create_db(str(gff_file), ':memory:', merge_strategy='merge')
    return db

def gene_model_stats(db):
    '''Compute statistics for predicted gene models.'''
    genes = list(db.features_of_type('gene'))
    transcripts = list(db.features_of_type(['mRNA', 'transcript']))

    exon_counts, intron_lengths, gene_lengths = [], [], []
    for gene in genes:
        gene_lengths.append(gene.end - gene.start + 1)
        for tx in db.children(gene, featuretype=['mRNA', 'transcript']):
            exons = list(db.children(tx, featuretype='exon'))
            exon_counts.append(len(exons))
            sorted_exons = sorted(exons, key=lambda e: e.start)
            for i in range(len(sorted_exons) - 1):
                intron_len = sorted_exons[i + 1].start - sorted_exons[i].end - 1
                intron_lengths.append(intron_len)

    stats = {
        'total_genes': len(genes),
        'total_transcripts': len(transcripts),
        'transcripts_per_gene': len(transcripts) / len(genes) if genes else 0,
        'median_exons_per_transcript': pd.Series(exon_counts).median() if exon_counts else 0,
        'median_gene_length': pd.Series(gene_lengths).median() if gene_lengths else 0,
        'median_intron_length': pd.Series(intron_lengths).median() if intron_lengths else 0,
        'single_exon_genes': sum(1 for e in exon_counts if e == 1),
    }
    return stats

db = load_gene_models('braker3_out/braker.gff3')
stats = gene_model_stats(db)
for key, val in stats.items():
    print(f'{key}: {val}')
```

## Troubleshooting

### Low Gene Count
- Verify softmasking (lowercase repeats present in assembly)
- Check RNA-seq alignment rate (>70% expected)
- Ensure OrthoDB proteins match the correct clade

### Many Single-Exon Genes
- May indicate bacterial contamination
- Check if assembly was properly softmasked
- Consider filtering predictions by evidence support

### BRAKER3 Fails
- Ensure all dependencies are in PATH (GeneMark, Augustus, DIAMOND)
- Check that genome headers have no special characters (pipes, spaces)
- Simplify FASTA headers: `sed 's/ .*//' assembly.fasta > clean.fasta`

## Related Skills

- repeat-annotation - PREREQUISITE: softmask repeats before gene prediction
- functional-annotation - Add GO/KEGG/Pfam to predicted proteins
- read-alignment/star-alignment - Alternative RNA-seq aligner for evidence
- genome-assembly/assembly-qc - Verify assembly quality before prediction
