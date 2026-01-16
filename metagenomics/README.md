# Metagenomics

Taxonomic profiling and analysis of metagenomic data using Kraken2 and MetaPhlAn.

## Overview

This category covers taxonomic classification and community profiling from metagenomic sequencing data. Kraken2 provides fast k-mer based classification, while MetaPhlAn uses clade-specific marker genes for accurate species-level profiling. Both tools are essential for understanding microbial community composition.

**Tool type:** `cli`
**Primary tools:** Kraken2, MetaPhlAn, Bracken

## Skills

| Skill | Description |
|-------|-------------|
| [kraken-classification](kraken-classification/) | Taxonomic classification with Kraken2 |
| [metaphlan-profiling](metaphlan-profiling/) | Marker gene profiling with MetaPhlAn |
| [abundance-estimation](abundance-estimation/) | Species abundance with Bracken |
| [metagenome-visualization](metagenome-visualization/) | Visualize taxonomic profiles |

## Workflow

```
Raw Metagenomic Reads (FASTQ)
    |
    v
[Quality Control] -----------> Trim adapters, filter low quality
    |
    +------------------+
    |                  |
    v                  v
[kraken-classification]    [metaphlan-profiling]
    |                  |
    v                  v
[abundance-estimation]     Relative abundances
    |                  |
    v                  v
[metagenome-visualization] -> Stacked bars, heatmaps, diversity
```

## Tool Comparison

| Feature | Kraken2 | MetaPhlAn |
|---------|---------|-----------|
| Method | K-mer matching | Marker genes |
| Speed | Very fast | Moderate |
| Database | RefSeq/custom | ChocoPhlAn |
| Output | Read counts | Relative abundance |
| Best for | Screening, viruses | Accurate profiling |

## Example Prompts

### Kraken2 Classification
- "Classify my metagenomic reads with Kraken2"
- "Run Kraken2 with the standard database"
- "What species are in my microbiome sample?"

### MetaPhlAn Profiling
- "Profile my metagenome with MetaPhlAn"
- "Get species-level abundances from my reads"
- "Run MetaPhlAn on paired-end reads"

### Abundance Estimation
- "Estimate species abundances with Bracken"
- "Convert Kraken2 output to abundance table"
- "Get genus-level abundance estimates"

### Visualization
- "Create a stacked bar chart of my samples"
- "Make a heatmap of species abundances"
- "Calculate alpha diversity metrics"

## Requirements

### Kraken2 and Bracken

```bash
# Conda installation
conda install -c bioconda kraken2 bracken

# Download standard database (8GB)
kraken2-build --standard --db kraken2_standard_db

# Or download pre-built
# https://benlangmead.github.io/aws-indexes/k2

# Verify
kraken2 --version
bracken -v
```

### MetaPhlAn

```bash
# Conda installation
conda install -c bioconda metaphlan

# Download database (on first run)
metaphlan --install

# Verify
metaphlan --version
```

## Database Options

### Kraken2 Databases

| Database | Size | Content |
|----------|------|---------|
| Standard | ~50GB | RefSeq bacteria, archaea, viral, human |
| Standard-8 | ~8GB | Reduced standard |
| PlusPF | ~70GB | Standard + protozoa, fungi |
| MinusB | ~8GB | Standard without bacteria |
| Viral | ~1GB | Viral only |

### MetaPhlAn Database

MetaPhlAn uses the ChocoPhlAn database with clade-specific marker genes. Database is downloaded automatically on first run (~1GB).

## Key Functions

| Tool | Command | Purpose |
|------|---------|---------|
| kraken2 | kraken2 | Classify reads |
| kraken2-build | Build database |
| kraken2-inspect | View database |
| bracken | Estimate abundance |
| metaphlan | Profile metagenome |
| merge_metaphlan_tables.py | Merge samples |

## Notes

- **Database size** - Kraken2 standard database is large (~50GB); consider pre-built indexes
- **Memory requirements** - Kraken2 loads database into RAM; ~50GB for standard
- **Paired-end reads** - Both tools support paired-end; concatenate or use specific flags
- **Relative vs absolute** - MetaPhlAn gives relative abundance; Kraken2/Bracken give counts
- **Unknown reads** - Both tools have unclassified reads; this is normal

## Related Skills

- **sequence-io** - FASTQ handling
- **database-access** - Download reference sequences
- **pathway-analysis** - Functional profiling (with HUMAnN)

## References

- [Kraken2 Manual](https://github.com/DerrickWood/kraken2/wiki)
- [MetaPhlAn Documentation](https://github.com/biobakery/MetaPhlAn)
- [Bracken Documentation](https://github.com/jenniferlu717/Bracken)
