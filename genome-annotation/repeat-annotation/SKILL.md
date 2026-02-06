---
name: bio-genome-annotation-repeat-annotation
description: Identify and classify repetitive elements and transposable elements using RepeatModeler for de novo repeat library construction and RepeatMasker for genome-wide repeat annotation. Quantify TE expression from RNA-seq with TEtranscripts. Use when masking repeats before gene prediction or analyzing transposable element activity.
tool_type: cli
primary_tool: RepeatMasker
---

# Repeat and Transposable Element Annotation

Identify, classify, and mask repetitive elements using RepeatModeler (de novo library construction) and RepeatMasker (genome-wide annotation). Softmasked output is a prerequisite for eukaryotic gene prediction.

## RepeatModeler (De Novo Library)

RepeatModeler builds a species-specific repeat library by detecting repetitive elements de novo from the assembly.

### Build Database and Run

```bash
# Build RepeatModeler database
BuildDatabase -name my_genome -engine ncbi assembly.fasta

# Run RepeatModeler (this takes hours to days depending on genome size)
# -LTRStruct enables LTR structural detection (recommended)
RepeatModeler -database my_genome -pa 16 -LTRStruct
```

### Key Options

| Option | Description |
|--------|-------------|
| `-database` | Database name from BuildDatabase |
| `-pa` | Parallel processes |
| `-LTRStruct` | Enable LTR structural detection pipeline |
| `-engine` | Search engine: ncbi (RMBLAST) or abblast |

### Output

```
my_genome-families.fa     # Consensus repeat library
my_genome-families.stk    # Stockholm alignments
RM_*/                     # Working directory with intermediate files
```

The output `*-families.fa` is the repeat library used by RepeatMasker.

## RepeatMasker (Genome-Wide Annotation)

### With De Novo Library

```bash
# Use species-specific de novo library (recommended)
RepeatMasker \
    -lib my_genome-families.fa \
    -pa 16 \
    -xsmall \
    -gff \
    -dir repeatmasker_out \
    assembly.fasta
```

### With Dfam/RepBase Library

```bash
# Use Dfam curated library for a known species
RepeatMasker \
    -species "Homo sapiens" \
    -pa 16 \
    -xsmall \
    -gff \
    -dir repeatmasker_out \
    assembly.fasta
```

### Combined Library (De Novo + Known)

```bash
# Combine de novo and curated libraries for best results
cat my_genome-families.fa known_repeats.fa > combined_lib.fa

RepeatMasker \
    -lib combined_lib.fa \
    -pa 16 \
    -xsmall \
    -gff \
    -dir repeatmasker_out \
    assembly.fasta
```

### Key Options

| Option | Description |
|--------|-------------|
| `-lib` | Custom repeat library FASTA |
| `-species` | Species name (uses Dfam database) |
| `-pa` | Parallel processes |
| `-xsmall` | Softmask output (lowercase repeats, required for gene prediction) |
| `-gff` | Generate GFF output |
| `-dir` | Output directory |
| `-nolow` | Skip low-complexity masking |
| `-noint` | Skip interspersed repeats |
| `-e` | Search engine: crossmatch, ncbi, hmmer, abblast |
| `-s` | Slow/sensitive search |
| `-q` | Quick search (5-10% less sensitive) |

### Output Files

```
repeatmasker_out/
├── assembly.fasta.masked    # Hardmasked genome (N's replace repeats)
├── assembly.fasta.out       # Detailed repeat annotation table
├── assembly.fasta.tbl       # Summary statistics table
├── assembly.fasta.out.gff   # GFF annotation of repeats
└── assembly.fasta.cat.gz    # Search result alignments
```

### Softmasking for Gene Prediction

The `-xsmall` flag produces softmasked output where repeats are lowercase. This is the required input format for BRAKER3 and most gene prediction tools.

```bash
# The softmasked genome is written in place of the input
# Copy original first
cp assembly.fasta assembly_original.fasta

RepeatMasker -lib my_genome-families.fa -pa 16 -xsmall assembly.fasta

# assembly.fasta.masked is the softmasked output
mv assembly.fasta.masked assembly_softmasked.fasta
```

## TEtranscripts (TE Expression)

Quantify transposable element expression from RNA-seq data using TEtranscripts, which works with DESeq2 for differential TE expression.

```bash
# Requires STAR alignment with multi-mapping reads retained
STAR --runThreadN 16 \
    --genomeDir star_index \
    --readFilesIn reads_R1.fq.gz reads_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --winAnchorMultimapNmax 100 \
    --outFilterMultimapNmax 100 \
    --outFileNamePrefix sample_

# Run TEtranscripts for differential expression
TEtranscripts \
    --treatment sample1.bam sample2.bam sample3.bam \
    --control ctrl1.bam ctrl2.bam ctrl3.bam \
    --GTF genes.gtf \
    --TE te_annotation.gtf \
    --mode multi \
    --sortByPos
```

### Key TEtranscripts Options

| Option | Description |
|--------|-------------|
| `--treatment` | Treatment BAM files |
| `--control` | Control BAM files |
| `--GTF` | Gene annotation GTF |
| `--TE` | TE annotation GTF (from RepeatMasker) |
| `--mode` | multi (recommended) or uniq |
| `--sortByPos` | Input sorted by position |
| `--stranded` | Strand-specific protocol (yes, no, reverse) |

## Python: Repeat Statistics

```python
import pandas as pd
import re

def parse_repeatmasker_out(out_file):
    '''Parse RepeatMasker .out file into a DataFrame.'''
    records = []
    with open(out_file) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            records.append({
                'score': int(parts[0]),
                'perc_div': float(parts[1]),
                'perc_del': float(parts[2]),
                'perc_ins': float(parts[3]),
                'seqid': parts[4],
                'start': int(parts[5]),
                'end': int(parts[6]),
                'strand': '+' if parts[8] == '+' else '-',
                'repeat_name': parts[9],
                'repeat_class': parts[10],
                'length': int(parts[6]) - int(parts[5]) + 1,
            })
    return pd.DataFrame(records)

def repeat_summary(rm_df, genome_size):
    '''Summarize repeat content by class.'''
    class_summary = rm_df.groupby('repeat_class').agg(
        count=('repeat_name', 'count'),
        total_bp=('length', 'sum'),
    ).sort_values('total_bp', ascending=False)

    class_summary['pct_genome'] = class_summary['total_bp'] / genome_size * 100
    total_masked = rm_df['length'].sum()

    print(f'=== Repeat Summary ===')
    print(f'Total masked: {total_masked:,} bp ({total_masked/genome_size:.1%} of genome)')
    print(f'\nBy class:')
    for cls, row in class_summary.iterrows():
        print(f'  {cls}: {row["count"]:,} elements, {row["total_bp"]:,} bp ({row["pct_genome"]:.1f}%)')

    return class_summary

def repeat_landscape(rm_df, output_file='repeat_landscape.png'):
    '''Plot repeat divergence landscape (Kimura distance).'''
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 6))

    major_classes = ['LINE', 'SINE', 'LTR', 'DNA']
    colors = {'LINE': '#1f77b4', 'SINE': '#ff7f0e', 'LTR': '#2ca02c', 'DNA': '#d62728'}

    for cls in major_classes:
        subset = rm_df[rm_df['repeat_class'].str.contains(cls, case=False, na=False)]
        if len(subset) > 0:
            ax.hist(subset['perc_div'], bins=50, range=(0, 50), alpha=0.6, label=cls, color=colors.get(cls))

    ax.set_xlabel('Kimura Divergence (%)')
    ax.set_ylabel('Count')
    ax.set_title('Repeat Landscape')
    ax.legend()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
```

## Expected Repeat Content

| Organism | Repeat Content | Notes |
|----------|---------------|-------|
| Bacteria | 1-5% | Mostly IS elements |
| Yeast | 3-5% | Ty elements |
| Drosophila | 15-25% | LTR-rich |
| Zebrafish | 45-55% | DNA transposon-rich |
| Human | 45-50% | LINE/SINE-rich |
| Maize | 80-85% | LTR-rich |

## Troubleshooting

### RepeatModeler Runs Very Slowly
- Normal for large genomes (days for mammalian-size)
- Use `-pa` for parallelization
- Consider EDTA as alternative for plant genomes

### Low Masking Percentage
- May indicate novel repeats not in database
- Always run RepeatModeler before RepeatMasker
- Combine de novo + known libraries

### Gene Prediction Finds Too Many Genes After Masking
- Verify softmasking with: `grep -v '^>' assembly.fasta | tr -cd 'a-z' | wc -c`
- Ensure using `-xsmall` not default hardmasking

## Related Skills

- eukaryotic-gene-prediction - Requires softmasked genome from repeat annotation
- genome-assembly/assembly-qc - Assess assembly quality including repeat content
- differential-expression/deseq2-basics - Differential TE expression analysis
