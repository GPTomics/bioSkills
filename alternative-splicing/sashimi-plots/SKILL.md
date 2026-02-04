---
name: bio-sashimi-plots
description: Creates sashimi plots showing RNA-seq read coverage and splice junction counts using ggsashimi or rmats2sashimiplot. Visualizes differential splicing events with grouped samples and junction read support. Use when visualizing specific splicing events or validating differential splicing results.
tool_type: python
primary_tool: ggsashimi
---

# Sashimi Plot Visualization

Create sashimi plots to visualize splicing events with read coverage and junction counts.

## ggsashimi Usage

```python
import subprocess
import pandas as pd

# Create sample grouping file (TSV: path, group, color)
groups = pd.DataFrame({
    'bam': ['sample1.bam', 'sample2.bam', 'sample3.bam', 'sample4.bam'],
    'group': ['control', 'control', 'treatment', 'treatment'],
    'color': ['#1f77b4', '#1f77b4', '#ff7f0e', '#ff7f0e']
})
groups.to_csv('sashimi_groups.tsv', sep='\t', index=False, header=False)

# Basic sashimi plot for a region
subprocess.run([
    'ggsashimi.py',
    '-b', 'sashimi_groups.tsv',
    '-c', 'chr1:1000000-1010000',  # Genomic coordinates
    '-o', 'sashimi_output',
    '-M', '10',  # Minimum junction reads to show
    '--alpha', '0.25',  # Coverage transparency
    '--height', '3',
    '--width', '8',
    '-g', 'annotation.gtf'
], check=True)
```

## Batch Plotting Significant Events

```python
import subprocess
import pandas as pd

# Load differential splicing results
diff_results = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')
significant = diff_results[
    (diff_results['FDR'] < 0.05) &
    (diff_results['IncLevelDifference'].abs() > 0.1)
]

# Generate plots for top events
for idx, event in significant.head(20).iterrows():
    chrom = event['chr']
    # Extend region around the exon
    start = event['upstreamES'] - 500
    end = event['downstreamEE'] + 500
    region = f'{chrom}:{start}-{end}'
    gene = event['geneSymbol']

    subprocess.run([
        'ggsashimi.py',
        '-b', 'sashimi_groups.tsv',
        '-c', region,
        '-o', f'sashimi_plots/{gene}_{chrom}_{start}',
        '-M', '5',
        '--shrink',  # Shrink introns for better visualization
        '-g', 'annotation.gtf',
        '--fix-y-scale'  # Same y-axis across groups
    ], check=True)
```

## rmats2sashimiplot

```bash
# For rMATS output specifically
rmats2sashimiplot \
    --b1 sample1.bam,sample2.bam \
    --b2 sample3.bam,sample4.bam \
    -t SE \
    -e rmats_output/SE.MATS.JC.txt \
    --l1 Control \
    --l2 Treatment \
    -o sashimi_rmats \
    --exon_s 1 \
    --intron_s 5
```

## Customization Options

```python
# Advanced ggsashimi options
subprocess.run([
    'ggsashimi.py',
    '-b', 'sashimi_groups.tsv',
    '-c', 'chr1:1000000-1010000',
    '-o', 'custom_sashimi',
    '-g', 'annotation.gtf',

    # Visual options
    '-M', '10',           # Min junction reads
    '--alpha', '0.25',    # Coverage alpha
    '--height', '3',      # Plot height per track
    '--width', '10',      # Plot width
    '--base-size', '14',  # Font size

    # Layout options
    '--shrink',           # Shrink introns
    '--fix-y-scale',      # Same y-axis
    '-A', 'mean',         # Aggregate: mean, median, or none

    # Annotation options
    '--gtf-filter', 'protein_coding',  # Filter GTF features

    # Output format
    '-F', 'pdf'           # pdf, png, svg, eps
], check=True)
```

## Best Practices

| Tip | Rationale |
|-----|-----------|
| Use `--shrink` for large introns | Keeps exons visible |
| Set `--fix-y-scale` for comparisons | Fair visual comparison |
| Aggregate replicates with `-A mean` | Reduces clutter |
| Limit to 3-4 groups | More groups become hard to read |
| Include flanking exons | Show full splicing context |

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No junctions shown | Lower `-M` threshold |
| Plot too crowded | Use `--shrink`, reduce samples |
| Annotation missing | Check GTF format, gene name field |
| Memory issues | Plot smaller regions |

## Related Skills

- differential-splicing - Identify events to plot
- splicing-quantification - Context for PSI values
- data-visualization/ggplot2-fundamentals - Further customization
