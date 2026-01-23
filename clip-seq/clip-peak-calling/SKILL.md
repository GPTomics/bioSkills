---
name: bio-clip-seq-clip-peak-calling
description: Call protein-RNA binding site peaks from CLIP-seq data using CLIPper, Piranha, or other peak callers. Use when identifying RBP binding sites from aligned CLIP reads.
tool_type: cli
primary_tool: CLIPper
---

# CLIP-seq Peak Calling

## CLIPper

```bash
clipper \
    -b deduped.bam \
    -s hg38 \
    -o peaks.bed \
    --save-pickle

# Output: BED file with binding clusters
```

## Piranha

```bash
piranha -s deduped.bam \
    -o peaks.bed \
    -p 0.01

# -p: p-value threshold
```

## Parse Peaks

```python
import pandas as pd

def load_clip_peaks(bed_path):
    peaks = pd.read_csv(bed_path, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    return peaks
```

## Related Skills

- **clip-alignment** - Generate aligned BAM
- **binding-site-annotation** - Annotate peaks
