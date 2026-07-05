# VCF/BCF Basics - Usage Guide

## Overview

View, query, and understand VCF/BCF variant files using bcftools and cyvcf2. Essential for working with variant calling results.

## Prerequisites

```bash
conda install -c bioconda bcftools htslib
pip install cyvcf2
```

## Quick Start

Tell your AI agent what you want to do:

- "View the first 20 variants in my VCF file"
- "List all samples in this VCF"
- "Extract chromosome, position, and genotypes to a TSV"
- "Convert my VCF to BCF format"
- "Query variants in the region chr1:1000000-2000000"

## Example Prompts

### Viewing Files

> "Show me the header of my VCF file"

> "View variants in chromosome 1 between positions 1M and 2M"

> "How many SNPs and indels are in this VCF?"

### Extracting Data

> "Extract variant positions and allele frequencies to a table"

> "Get genotypes for all samples as a TSV file"

> "List all variant IDs (rsIDs) in the file"

### Format Conversion

> "Compress my VCF with bgzip and index it"

> "Convert this VCF to BCF format for faster processing"

### Python Analysis

> "Iterate through variants and count heterozygous calls per sample"

> "Filter variants with QUAL > 30 and write to a new file"

### Interpreting Fields

> "Explain the difference between QUAL and GQ for this site"

> "Compute allele balance from AD for the heterozygous calls"

> "Why does this record have a <NON_REF> allele -- is this a gVCF?"

> "Check whether my missing genotypes are being treated as reference"

## What the Agent Will Do

1. Identify the VCF/BCF file location and check if it's indexed
2. Select appropriate tool (bcftools for CLI, cyvcf2 for Python)
3. Construct the query or view command with proper flags
4. Execute the operation and format output as requested
5. For region queries, ensure the file is indexed first

## Tips

- Always use bgzip (not gzip) for VCF compression -- bcftools requires it
- Index files (.tbi or .csi) enable fast region queries
- BCF format is faster to process than VCF for large files
- Use `-H` flag with bcftools query to skip the header line
- cyvcf2 is faster than PyVCF for large files
- Filter status of `None` in cyvcf2 means the variant passed all filters
- QUAL measures site-level variant existence; GQ measures per-sample genotype confidence -- they answer different questions and are not interchangeable, so filter on both
- PL is phred-scaled and rebased so the called genotype is 0; GQ is the difference between the two smallest PL values (GL is the same info as raw log10 likelihoods)
- The sum of AD values may be less than DP because DP includes uninformative reads -- this is expected, not an error; allele balance for a het is derived as alt_AD / (ref_AD + alt_AD)
- QD (quality by depth) is generally preferred over raw QUAL for filtering because it normalizes for coverage
- A field's header Number (A/R/G/.) controls how it re-subsets on a multiallelic split; a field mis-declared Number=. silently keeps the wrong allele's value
- `.` (missing) is never `0/0` (hom-ref) -- treating `./.` as reference biases allele frequencies and burden tests
- `|` (phased) is only meaningful within its PS phase-set block; a bare `*` ALT is a spanning-deletion placeholder, not a real allele
- A raw gVCF (with `<NON_REF>` records and reference blocks) is an intermediate, not a filtered callset -- joint-genotype it before filtering, annotating, or counting variants
- At multiallelic sites, splitting into biallelic records loses compound heterozygosity information -- consider whether downstream tools need multiallelic or biallelic input

## Related Skills

- variant-calling/variant-calling - Generate VCF from alignments
- variant-calling/variant-normalization - Split multiallelics, left-align, Number-code reapportionment
- variant-calling/filtering-best-practices - Filter variants by site and genotype quality
- variant-calling/joint-calling - gVCF reference-confidence model and joint genotyping
- variant-calling/vcf-manipulation - Merge, concat, intersect VCFs
- alignment-files/pileup-generation - Generate pileup for calling
