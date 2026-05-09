---
name: bio-reference-operations
description: Generate consensus sequences and manage reference files using samtools. Use when creating consensus from alignments, indexing references, or creating sequence dictionaries.
tool_type: cli
primary_tool: samtools
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, bcftools 1.19+, pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Reference Operations

Generate consensus sequences and manage reference files using samtools.

**"Prepare a reference genome"** → Index the FASTA and create a sequence dictionary for downstream tools.
- CLI: `samtools faidx ref.fa` + `samtools dict ref.fa -o ref.dict`
- Python: `pysam.FastaFile('ref.fa')` (auto-uses .fai index)

**"Build a consensus from BAM"** → Derive the most-supported base at each position from aligned reads.
- CLI: `samtools consensus input.bam -o consensus.fa`
- Python: iterate pileup columns and take majority base (pysam)

## samtools faidx - Index Reference FASTA

Create index for random access to reference sequences.

### Create Index
```bash
samtools faidx reference.fa
# Creates reference.fa.fai
```

### Fetch Region from Reference
```bash
samtools faidx reference.fa chr1:1000-2000
```

### Fetch Multiple Regions
```bash
samtools faidx reference.fa chr1:1000-2000 chr2:3000-4000
```

### Fetch Entire Chromosome
```bash
samtools faidx reference.fa chr1
```

### Output to File
```bash
samtools faidx reference.fa chr1:1000-2000 > region.fa
```

### Reverse Complement
```bash
samtools faidx -i reference.fa chr1:1000-2000
```

### FAI File Format
```
chr1    248956422    6    60    61
chr2    242193529    253404903    60    61
```
Columns: name, length, offset, line bases, line width

## samtools dict - Create Sequence Dictionary

Create SAM header dictionary for reference (used by GATK, Picard).

### Create Dictionary
```bash
samtools dict reference.fa -o reference.dict
```

### With Assembly Info
```bash
samtools dict -a GRCh38 -s "Homo sapiens" reference.fa -o reference.dict
```

### Dictionary Format
```
@HD VN:1.6 SO:unsorted
@SQ SN:chr1 LN:248956422 M5:6aef897c3d6ff0c78aff06ac189178dd UR:file:reference.fa
@SQ SN:chr2 LN:242193529 M5:f98db672eb0993dcfdabafe2a882905c UR:file:reference.fa
```

The `M5:` (MD5) tag is the only definitive reference-identity check -- two references named "GRCh38" with different decoy/alt content have different M5s. CRAM enforces M5 match on read-back. See alignment-validation for BAM-vs-reference M5 cross-check.

### GRCh38 Is Not One Reference

| Reference flavor | ALT | Decoy | EBV | HLA | Use case |
|------------------|-----|-------|-----|-----|----------|
| GRCh38 no-alt | no | no | no | no | Conservative analyses |
| GRCh38 + decoy + EBV (1000G analysis set) | no | yes | yes | no | Cohort projects |
| GRCh38 ALT + decoy + EBV + HLA (Broad / hs38DH) | yes | yes | yes | yes | GATK Best Practices |
| T2T-CHM13 v2.0 | n/a | n/a | n/a | n/a | Distinct coordinates -- NOT interchangeable |

Mixing no-alt and ALT-aware BAMs in one cohort produces inconsistent multi-mapping behavior at HLA, KIR, and segmental-duplication regions. Standardize before joint calling.

### Contig Naming: The Silent Killer

| Convention | Source | chr1 | mitochondrion |
|-----------|--------|------|---------------|
| UCSC (hg19, hg38) | UCSC Genome Browser | chr1 | chrM |
| Ensembl (GRCh37, GRCh38) | Ensembl, ENA | 1 | MT |
| NCBI RefSeq (recent) | NCBI | chr1 | chrM |
| 1000G analysis sets | 1000G Phase II/III | chr1 | chrM |

A BAM with `@SQ SN:chr1` cannot be analyzed against a `1`-named reference (and vice versa). Detect:
```bash
samtools view -H sample.bam | grep '^@SQ' | head -3
samtools dict ref.fa | head -3
```

Convert: `bcftools annotate --rename-chrs` for VCF; for BAM there is no clean conversion -- re-align.

## samtools consensus - Generate Consensus

Create consensus sequence from alignments.

### Basic Consensus
```bash
samtools consensus input.bam -o consensus.fa
```

### From Specific Region
```bash
samtools consensus -r chr1:1000-2000 input.bam -o region_consensus.fa
```

### Output Formats
```bash
# FASTA (default)
samtools consensus -f fasta input.bam -o consensus.fa

# FASTQ (includes quality)
samtools consensus -f fastq input.bam -o consensus.fq
```

### Quality Options
```bash
# Minimum depth to call base
samtools consensus -d 5 input.bam -o consensus.fa

# Call all positions (including low coverage)
samtools consensus -a input.bam -o consensus.fa
```

### IUPAC Ambiguity for Heterozygotes
```bash
# Emit IUPAC codes (R, Y, S, W, K, M, B, D, H, V, N) for heterozygous columns
# --ambig is REQUIRED -- without it, output is restricted to A,C,G,T,N,*
samtools consensus --ambig --het-fract 0.2 --call-fract 0.5 input.bam -o consensus.fa
```

`--het-fract` controls the fraction of the second-most-common base relative to the most common required to call a heterozygote (default ~0.15). Without `--ambig`, columns where the second base passes `--het-fract` resolve to `N` rather than the IUPAC code. `--show-ins` / `--show-del` control insertion / deletion display, not ambiguity.

### Platform-Aware Consensus
```bash
# Default: Bayesian Illumina profile
samtools consensus -f fasta input.bam -o consensus.fa

# Platform-specific profiles (samtools 1.21+; verify via samtools consensus --help for installed version)
samtools consensus --config hifi    input.bam -o consensus.fa   # PacBio HiFi
samtools consensus --config ont     input.bam -o consensus.fa   # ONT R10.4+
samtools consensus --config ultima  input.bam -o consensus.fa   # Ultima Genomics
samtools consensus --config illumina input.bam -o consensus.fa  # default

# Report ref base where consensus unavailable (low coverage)
samtools consensus -T ref.fa input.bam -o consensus.fa
```

### samtools consensus vs bcftools consensus

Different operations -- conflating them produces nonsense:

| Tool | Input | Output | Use case |
|------|-------|--------|----------|
| `samtools consensus` | BAM | Consensus FASTA derived from reads (Bayesian) | Viral, de novo / amplicon, low-coverage species |
| `bcftools consensus` | reference + VCF | Reference with VCF variants applied | Apply called variants (haplotype reconstruction, custom ref for re-mapping) |

For viral consensus from BAM:
```bash
# Modern: samtools consensus
samtools consensus --config illumina -d 10 --het-fract 0.5 \
    --show-ins yes --show-del yes input.bam -o consensus.fa

# Apply called variants to reference (different question)
bcftools consensus -f reference.fa variants.vcf.gz -o sample_consensus.fa
bcftools consensus -f reference.fa -H 1 phased.vcf.gz -o haplotype1.fa   # phased haplotype 1
```

For bacterial / phage assembly polishing, prefer Pilon (short-read) or medaka (ONT); `samtools consensus` is not iterative.

## pysam Python Alternative

### Fetch from Indexed FASTA
```python
import pysam

with pysam.FastaFile('reference.fa') as ref:
    seq = ref.fetch('chr1', 999, 2000)  # 0-based
    print(seq)
```

### Get Reference Lengths
```python
with pysam.FastaFile('reference.fa') as ref:
    for name in ref.references:
        length = ref.get_reference_length(name)
        print(f'{name}: {length:,} bp')
```

### Fetch All Chromosomes
```python
with pysam.FastaFile('reference.fa') as ref:
    for chrom in ref.references:
        seq = ref.fetch(chrom)
        print(f'>{chrom}')
        print(seq[:100] + '...')
```

### Generate Simple Consensus
```python
import pysam
from collections import Counter

def consensus_at_position(bam, chrom, pos):
    bases = Counter()
    for pileup in bam.pileup(chrom, pos, pos + 1, truncate=True):
        if pileup.pos == pos:
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    bases[read.alignment.query_sequence[read.query_position]] += 1
    if bases:
        return bases.most_common(1)[0][0]
    return 'N'

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    consensus = consensus_at_position(bam, 'chr1', 1000000)
    print(f'Consensus at chr1:1000000 = {consensus}')
```

### Build Consensus Sequence (Pedagogical Only)

The Python majority-vote consensus below is illustrative, NOT production. `samtools consensus` is Bayesian, quality-aware, and platform-aware; majority vote ignores base qualities and produces wrong calls on low-coverage / low-quality regions. Use for teaching pileup iteration mechanics; use `samtools consensus` for any real consensus.

```python
import pysam
from collections import Counter

def build_consensus(bam_path, chrom, start, end, min_depth=3):
    consensus = []

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup in bam.pileup(chrom, start, end, truncate=True):
            bases = Counter()
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base] += 1

            if sum(bases.values()) >= min_depth:
                consensus.append(bases.most_common(1)[0][0])
            else:
                consensus.append('N')

    return ''.join(consensus)
```

### Create Dictionary Header
```python
import pysam

def create_dict_header(fasta_path):
    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': []}

    with pysam.FastaFile(fasta_path) as ref:
        for name in ref.references:
            length = ref.get_reference_length(name)
            header['SQ'].append({'SN': name, 'LN': length})

    return header

header = create_dict_header('reference.fa')
for sq in header['SQ'][:5]:
    print(f'{sq["SN"]}: {sq["LN"]:,} bp')
```

## Reference Preparation Workflow

**Goal:** Set up a reference genome with all indices needed by common analysis tools.

**Approach:** Create FASTA index (.fai), sequence dictionary (.dict), and aligner-specific indices in sequence.

### Prepare Reference for Analysis
```bash
# 1. Index FASTA for samtools/pysam
samtools faidx reference.fa

# 2. Create sequence dictionary for GATK/Picard
samtools dict reference.fa -o reference.dict

# 3. Pre-populate CRAM REF_CACHE (for offline HPC nodes)
seq_cache_populate.pl -root $REF_CACHE_DIR reference.fa
```

For aligner-specific indices (BWA, Bowtie2, STAR, minimap2, Salmon), see read-alignment.

### Check Reference Setup
```bash
# Verify FAI exists
ls -la reference.fa.fai

# Verify dict exists
head reference.dict

# Test fetch
samtools faidx reference.fa chr1:1-100
```

## Common Operations

### Extract Chromosome
```bash
samtools faidx reference.fa chr1 > chr1.fa
samtools faidx chr1.fa  # Index the subset
```

### Get Chromosome Sizes
```bash
cut -f1,2 reference.fa.fai > chrom.sizes
```

### Subset Reference
```bash
samtools faidx reference.fa chr1 chr2 chr3 > subset.fa
samtools faidx subset.fa
```

### Compare Consensus to Reference
```bash
# Generate consensus
samtools consensus input.bam -o consensus.fa

# Align consensus back to reference
minimap2 -a reference.fa consensus.fa > comparison.sam
```

## Quick Reference

| Task | Command |
|------|---------|
| Index FASTA | `samtools faidx ref.fa` |
| Fetch region | `samtools faidx ref.fa chr1:1-1000` |
| Create dict | `samtools dict ref.fa -o ref.dict` |
| Build consensus | `samtools consensus in.bam -o out.fa` |
| Chrom sizes | `cut -f1,2 ref.fa.fai` |

## Related Skills

- sam-bam-basics - CRAM reference resolution (REF_PATH, REF_CACHE)
- alignment-indexing - faidx for reference access
- alignment-validation - BAM-vs-reference M5 cross-validation
- pileup-generation - Pileup for consensus building
- variant-calling/vcf-basics - VCF I/O for `bcftools consensus`
- variant-calling/consensus-sequences - Consensus from VCF (different operation)
- read-alignment/bwa-alignment - BWA index preparation
- sequence-io/read-sequences - Parse FASTA with Biopython
