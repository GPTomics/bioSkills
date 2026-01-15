# Variant Calling

Skills for variant calling and VCF/BCF file manipulation using bcftools and cyvcf2.

## Overview

This category covers the variant calling workflow: calling SNPs/indels from alignments, filtering, normalization, annotation, and downstream analysis.

**Tool type:** `cli` (with Python alternatives via cyvcf2)
**Primary tools:** bcftools, cyvcf2

## Skills

| Skill | Description |
|-------|-------------|
| [vcf-basics](vcf-basics/) | View, query, understand VCF/BCF format structure |
| [variant-calling](variant-calling/) | Call SNPs/indels from BAM files using mpileup/call |
| [vcf-filtering](vcf-filtering/) | Filter variants by quality, depth, type, expressions |
| [vcf-manipulation](vcf-manipulation/) | Merge, concat, sort, intersect VCF files |
| [variant-normalization](variant-normalization/) | Left-align indels, split multiallelic sites |
| [variant-annotation](variant-annotation/) | Add annotations, predict functional consequences |
| [vcf-statistics](vcf-statistics/) | Generate quality metrics, Ti/Tv ratio, concordance |
| [consensus-sequences](consensus-sequences/) | Apply variants to reference FASTA |

## Workflow

```
BAM file (sorted, indexed)
    |
    v
[variant-calling] - mpileup + call
    |
    v
Raw VCF
    |
    +---> [vcf-basics] - View, inspect, query
    |
    v
[variant-normalization] - Left-align, split multiallelics
    |
    v
[vcf-filtering] - Quality/depth filters
    |
    +---> [vcf-statistics] - QC metrics, Ti/Tv ratio
    |
    v
[variant-annotation] - rsIDs, AF, consequences
    |
    +---> [vcf-manipulation] - Merge, compare call sets
    |
    v
[consensus-sequences] - Apply to reference
```

## Example Prompts

### Calling and Viewing
- "Call variants from my aligned BAM file"
- "View the first 20 variants in my VCF"
- "Extract chromosome 1 variants to a new file"
- "List all sample names in this VCF"

### Filtering
- "Filter variants with QUAL < 30"
- "Keep only SNPs with depth >= 10"
- "Extract PASS variants only"
- "Get rare variants with AF < 0.01"

### Manipulation
- "Merge VCF files from different samples"
- "Compare variants between two callers"
- "Find variants shared between two files"
- "Concatenate per-chromosome VCFs"

### Normalization and Annotation
- "Normalize indels to left-aligned representation"
- "Split multiallelic sites to biallelic"
- "Add rsIDs from dbSNP"
- "Predict functional consequences"

### Statistics and Consensus
- "Calculate transition/transversion ratio"
- "Check sample identity between files"
- "Generate consensus sequence from variants"

## Installation

### bcftools
```bash
# Conda
conda install -c bioconda bcftools

# macOS
brew install bcftools

# Ubuntu
sudo apt install bcftools
```

### cyvcf2
```bash
pip install cyvcf2
```

## Related Skills

- **alignment-files** - Prepare BAM files for variant calling
- **database-access** - Download reference databases (dbSNP, gnomAD)
- **sequence-io** - Work with consensus FASTA output

## References

- [bcftools documentation](http://www.htslib.org/doc/bcftools.html)
- [cyvcf2 documentation](https://brentp.github.io/cyvcf2/)
- [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
