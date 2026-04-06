---
name: bio-variant-calling-structural-variant-calling
description: Call structural variants (SVs) from sequencing data using Manta, Delly, GRIDSS, and LUMPY. Detects deletions, insertions, inversions, duplications, and translocations too large for standard SNV callers. Use when detecting structural variants from short-read or long-read data and building consensus callsets.
tool_type: cli
primary_tool: manta
---

## Version Compatibility

Reference examples tested with: Manta 1.6+, Delly 1.2+, GRIDSS 2.13+, bcftools 1.19+, samtools 1.19+, SURVIVOR 1.0.7+, Sniffles2 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structural Variant Calling

**"Call structural variants from my WGS data"** -> Detect large genomic rearrangements (deletions, insertions, inversions, duplications, translocations) using split-read, discordant-pair, and assembly-based evidence.
- CLI: `configManta.py` (Manta), `delly call`, `gridss` (GRIDSS), `lumpyexpress`/`smoove call`

## SV Detection Limitations by Platform

Not all SV types are equally detectable across sequencing platforms. This table reflects practical detection performance, not theoretical capability:

| SV Type | Short-read Detection | Long-read Detection | Key Limitation |
|---------|---------------------|---------------------|----------------|
| Deletion | Good (read-pair + split-read) | Excellent | Short reads miss deletions in repetitive regions |
| Duplication | Moderate (read-pair + depth) | Good | Tandem vs dispersed distinction unreliable with short reads |
| Inversion | Moderate (read-pair) | Good | Breakpoints in repeats cause false negatives |
| Insertion | Poor (limited by read length) | Excellent | Short reads cannot resolve insertions >read length |
| Translocation | Moderate (discordant pairs) | Good | High false positive rate near centromeres/telomeres |
| Complex/nested | Poor | Good (with assembly) | Multiple overlapping SVs confound short-read signals |

## Caller Comparison

| Feature | Manta | Delly | GRIDSS | Smoove/LUMPY |
|---------|-------|-------|--------|--------------|
| Method | Read-pair + split-read + local assembly | Read-pair + split-read | Positional de Bruijn graph assembly | Read-pair + split-read |
| Speed | Fastest | Moderate | Slowest (2-5x Manta) | Moderate |
| DEL detection | Good | Good | Best precision | Good |
| INS detection | Good | Limited (small INS only) | Good | Cannot detect |
| Somatic mode | Yes | Yes | Yes (GRIDSS2/GRIPSS) | Limited |
| RNA-seq | Yes | No | No | No |
| Single breakends | No | No | Yes | No |
| Complex SVs | Limited | No | Yes (via LINX) | No |

GRIDSS produces the highest precision for deletions and uniquely detects single breakend events (one side of a breakpoint where the partner cannot be mapped). Manta provides the best speed-to-accuracy ratio for most applications. Delly excels at joint calling across cohorts. LUMPY/Smoove lacks insertion detection entirely.

## Consensus Calling Strategy

Current best practice: run Delly + GRIDSS + Manta + SvABA, require 2/4 caller agreement. This consensus approach yields best sensitivity with minimized false positives. Each caller has distinct algorithmic biases, so union sets are noisy while strict intersection is too conservative.

## Manta

```bash
configManta.py \
    --bam sample.bam \
    --referenceFasta reference.fa \
    --runDir manta_run

manta_run/runWorkflow.py -j 8

# Output: manta_run/results/variants/
# - diploidSV.vcf.gz (germline SVs)
# - candidateSV.vcf.gz (all candidates before scoring)
# - candidateSmallIndels.vcf.gz (50-1000bp indels for Strelka input)
```

## Manta Tumor-Normal Mode

```bash
configManta.py \
    --tumorBam tumor.bam \
    --normalBam normal.bam \
    --referenceFasta reference.fa \
    --runDir manta_somatic

manta_somatic/runWorkflow.py -j 8

# Output includes:
# - somaticSV.vcf.gz (somatic SVs, scored by tumor/normal evidence ratio)
# - diploidSV.vcf.gz (germline SVs)
```

## Manta Options

```bash
# WES mode (adjusts depth filters for uneven exome coverage)
configManta.py \
    --bam sample.bam \
    --referenceFasta reference.fa \
    --exome \
    --callRegions regions.bed.gz \
    --runDir manta_exome

# RNA-seq mode (handles split alignments across splice junctions)
configManta.py \
    --bam rnaseq.bam \
    --referenceFasta reference.fa \
    --rna \
    --runDir manta_rna
```

## Delly

```bash
delly call -g reference.fa -o sv_calls.bcf sample.bam
bcftools view sv_calls.bcf > sv_calls.vcf

# Joint calling across cohort (recommended for population studies)
delly call -g reference.fa -o joint_svs.bcf sample1.bam sample2.bam sample3.bam
```

## Delly Somatic Mode

```bash
delly call -g reference.fa -o svs.bcf tumor.bam normal.bam

echo -e "tumor\ttumor\nnormal\tcontrol" > samples.tsv

delly filter -f somatic -o somatic_svs.bcf -s samples.tsv svs.bcf
```

## Delly SV Types

```bash
delly call -t DEL -g ref.fa -o deletions.bcf sample.bam
delly call -t DUP -g ref.fa -o duplications.bcf sample.bam
delly call -t INV -g ref.fa -o inversions.bcf sample.bam
delly call -t BND -g ref.fa -o translocations.bcf sample.bam
delly call -t INS -g ref.fa -o insertions.bcf sample.bam
```

## GRIDSS

GRIDSS uses positional de Bruijn graph assembly to reconstruct breakpoints, producing the highest precision among short-read callers. It detects single breakend events where only one side of a rearrangement maps to the reference--critical for viral integrations, centromeric breakpoints, and highly rearranged cancer genomes.

```bash
gridss \
    --reference reference.fa \
    --output gridss_svs.vcf \
    --assembly gridss_assembly.bam \
    --threads 8 \
    sample.bam
```

## GRIDSS Somatic Mode (GRIDSS2 + GRIPSS)

```bash
# GRIDSS2 with paired tumor-normal
gridss \
    --reference reference.fa \
    --output gridss_raw.vcf \
    --assembly gridss_assembly.bam \
    --labels normal,tumor \
    --threads 8 \
    normal.bam tumor.bam

# GRIPSS post-filtering (somatic/germline classification)
gripss \
    -ref_genome reference.fa \
    -ref_genome_version 38 \
    -sample tumor \
    -reference normal \
    -vcf gridss_raw.vcf \
    -output_dir gripss_output/
```

Complex rearrangement reconstruction is available via LINX, which interprets GRIDSS breakpoints into higher-order SV events (chromothripsis, breakage-fusion-bridge cycles).

## LUMPY

```bash
samtools view -b -F 1294 sample.bam > discordant.bam
samtools view -h sample.bam | \
    /path/to/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | \
    samtools view -Sb - > splitters.bam

lumpyexpress \
    -B sample.bam \
    -S splitters.bam \
    -D discordant.bam \
    -o lumpy_svs.vcf
```

## Smoove (LUMPY Wrapper)

```bash
smoove call \
    --name sample \
    --fasta reference.fa \
    --outdir smoove_output \
    -p 8 \
    sample.bam

# Output: smoove_output/sample-smoove.genotyped.vcf.gz
```

## Merge Multiple Callers with SURVIVOR

**Goal:** Increase confidence in SV calls by requiring support from multiple callers with distinct algorithmic approaches.

**Approach:** Run 2-4 callers independently, then merge callsets with SURVIVOR requiring agreement on breakpoint proximity and SV type. Using max_dist=1000bp allows for the breakpoint imprecision inherent in short-read callers while min_callers=2 filters false positives unique to any single algorithm.

```bash
ls manta_svs.vcf delly_svs.vcf gridss_svs.vcf smoove_svs.vcf > vcf_list.txt

# max_dist=1000  min_callers=2  type_agree=1  strand_agree=1  estimate_dist=0  min_size=50
SURVIVOR merge vcf_list.txt 1000 2 1 1 0 50 merged_svs.vcf
```

The 1000bp max_dist accounts for breakpoint position uncertainty across callers (Manta and GRIDSS resolve breakpoints more precisely than Delly/LUMPY). Requiring type_agree=1 prevents merging a deletion call with a duplication call at the same locus.

## Filter SV Calls

```bash
bcftools view -i 'QUAL >= 20' svs.vcf > svs.filtered.vcf
bcftools view -i 'ABS(SVLEN) >= 50' svs.vcf > svs.min50.vcf

# Filter by SV type
bcftools view -i 'SVTYPE="DEL"' svs.vcf > deletions.vcf
bcftools view -i 'SVTYPE="INS"' svs.vcf > insertions.vcf
bcftools view -i 'SVTYPE="INV"' svs.vcf > inversions.vcf
bcftools view -i 'SVTYPE="DUP"' svs.vcf > duplications.vcf
bcftools view -i 'SVTYPE="BND"' svs.vcf > translocations.vcf

bcftools view -f PASS svs.vcf > svs.pass.vcf
```

## Annotate SVs

```bash
AnnotSV \
    -SVinputFile svs.vcf \
    -genomeBuild GRCh38 \
    -outputFile annotated_svs

# Output includes: gene overlap, DGV frequency, gnomAD-SV population AF, ClinVar pathogenicity
```

## SV Types

| Type | Code | Description | Typical Size Range |
|------|------|-------------|--------------------|
| Deletion | DEL | Sequence removed | 50bp - 100Mb |
| Insertion | INS | Novel sequence inserted | 50bp - 10kb (short-read); unlimited (long-read) |
| Inversion | INV | Sequence orientation reversed | 1kb - 10Mb |
| Duplication | DUP | Sequence copied (tandem or dispersed) | 1kb - 10Mb |
| Translocation | BND | Breakend connecting different chromosomes | N/A (inter-chromosomal) |

## Coverage Guidelines

| Coverage | Detection Ability | Practical Guidance |
|----------|-------------------|--------------------|
| 10x | Large SVs only (>1kb) | Limited breakpoint accuracy; high false negative rate for SVs <1kb; suitable only for large deletion screening |
| 30x | Most SVs detected | Standard for WGS; good sensitivity for DEL/DUP/INV >300bp; moderate INS detection |
| 50x+ | Small SVs, precise breakpoints | Better sensitivity near repetitive regions; resolves complex SVs; recommended for clinical applications |

Below 30x, split-read evidence becomes sparse and callers rely more heavily on read-pair signals, which have lower breakpoint resolution (~300-500bp uncertainty vs ~10bp for split-reads).

## Short-read vs Long-read Decision Framework

Short reads are sufficient for: deletions >300bp, balanced translocations, large tandem duplications, and population-scale screening where cost per sample matters.

Long reads are necessary for: insertions exceeding read length, complex/nested SVs, SVs in repetitive regions (segmental duplications, LINE/SINE elements), complete breakpoint resolution, and phased SV haplotyping.

Cost consideration: short reads for population-scale SV surveys (hundreds of samples), long reads for clinical-grade SV characterization where completeness matters more than throughput.

## Long-Read SV Callers

| Caller | Best For | Key Strengths |
|--------|----------|---------------|
| Sniffles2 | ONT/HiFi general | 11.8x faster than v1; population merging with `sniffles --merge`; mosaic SV detection; best overall accuracy |
| CuteSV2 | ONT data | Highest recall for ONT; signature-based clustering handles noisy reads |
| pbsv | PacBio HiFi | Official PacBio tool; best paired with PBMM2 aligner; tandem repeat aware |
| Severus | Somatic SVs | Phased breakpoint graph approach; resolves complex somatic rearrangements (Nature Biotechnology 2025) |

### Recommended Aligner-Caller Pairings

- Minimap2 + CuteSV2: ONT general purpose; fastest end-to-end
- Winnowmap + Sniffles2: high accuracy in repetitive regions (Winnowmap downweights repetitive k-mers)
- PBMM2 + pbsv: PacBio HiFi data; PBMM2 produces the CIGAR strings pbsv expects

See long-read-sequencing/structural-variants for long-read SV calling workflows with full pipeline examples.

## Related Skills

- long-read-sequencing/structural-variants - Long-read SV calling with Sniffles2, CuteSV, pbsv
- copy-number/cnvkit-analysis - Copy number variant detection (complements SV calling for dosage changes)
- variant-calling/filtering-best-practices - VCF filtering strategies applicable to SV callsets
- variant-calling/variant-annotation - Functional annotation of variants including SVs
- alignment-files/alignment-filtering - BAM preparation and quality filtering before SV calling
