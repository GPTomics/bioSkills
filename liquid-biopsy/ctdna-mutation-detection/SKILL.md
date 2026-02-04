---
name: bio-ctdna-mutation-detection
description: Detects somatic mutations in circulating tumor DNA using variant callers optimized for low allele fractions with UMI-based error suppression. Reliably detects mutations at VAF above 0.5 percent using consensus-based approaches. Use when identifying tumor mutations from plasma DNA or tracking specific variants.
tool_type: python
primary_tool: VarDict
---

# ctDNA Mutation Detection

Detect somatic mutations in cfDNA at low variant allele fractions.

## Input Requirements

| Requirement | Specification |
|-------------|---------------|
| Data type | Targeted panel or WES (NOT sWGS) |
| Depth | >= 1000x for low VAF detection |
| UMIs | Highly recommended for < 1% VAF |
| Input | Preprocessed BAM (UMI consensus if available) |

## VAF Detection Limits

| VAF Range | Reliability | Notes |
|-----------|-------------|-------|
| > 1% | Reliable | Standard callers work |
| 0.5-1% | Good with UMIs | Requires error suppression |
| 0.1-0.5% | Challenging | Needs deep UMI consensus |
| < 0.1% | Unreliable | Near noise floor |

## VarDict for High Sensitivity

```bash
# VarDict is highly sensitive for low VAF
# Use on UMI-consensus BAM for best results

vardict-java \
    -G reference.fa \
    -f 0.005 \      # Min VAF 0.5%
    -N sample_id \
    -b sample.bam \
    -c 1 -S 2 -E 3 -g 4 \
    regions.bed | \
teststrandbias.R | \
var2vcf_valid.pl \
    -N sample_id \
    -E \
    -f 0.005 \
    > sample.vcf
```

## Python Implementation

```python
import subprocess
import pandas as pd
import pysam


def call_variants_vardict(bam_file, reference, bed_file, output_vcf, min_vaf=0.005, min_depth=100):
    '''
    Call variants with VarDict.

    Args:
        bam_file: UMI-consensus BAM preferred
        reference: Reference FASTA
        bed_file: Target regions BED
        output_vcf: Output VCF path
        min_vaf: Minimum VAF (0.005 = 0.5%)
        min_depth: Minimum read depth
    '''
    sample_id = bam_file.split('/')[-1].replace('.bam', '')

    cmd = f'''
    vardict-java \
        -G {reference} \
        -f {min_vaf} \
        -N {sample_id} \
        -b {bam_file} \
        -c 1 -S 2 -E 3 -g 4 \
        {bed_file} | \
    teststrandbias.R | \
    var2vcf_valid.pl \
        -N {sample_id} \
        -E \
        -f {min_vaf} \
        > {output_vcf}
    '''

    subprocess.run(cmd, shell=True, check=True)
    return output_vcf


def filter_ctdna_variants(vcf_file, chip_genes=None):
    '''
    Filter ctDNA variants, removing CHIP.

    CHIP genes commonly mutated in elderly:
    DNMT3A, TET2, ASXL1, PPM1D, TP53, SF3B1, etc.
    '''
    if chip_genes is None:
        chip_genes = ['DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'JAK2',
                      'SF3B1', 'SRSF2', 'TP53', 'CBL', 'BCOR']

    import vcfpy
    reader = vcfpy.Reader.from_path(vcf_file)

    somatic = []
    chip = []

    for record in reader:
        gene = record.INFO.get('GENE', [''])[0]

        if gene in chip_genes:
            chip.append(record)
        else:
            somatic.append(record)

    print(f'Somatic variants: {len(somatic)}')
    print(f'Potential CHIP variants: {len(chip)}')

    return somatic, chip
```

## UMI-VarCal for Best Specificity

```python
def call_with_umi_varcal(bam_file, reference, bed_file, output_vcf, min_vaf=0.005):
    '''
    UMI-VarCal: Best specificity with UMI data.
    '''
    subprocess.run([
        'umi-varcal',
        '--bam', bam_file,
        '--ref', reference,
        '--bed', bed_file,
        '--out', output_vcf,
        '--min-vaf', str(min_vaf),
        '--min-alt-reads', '3',
        '--min-depth', '100'
    ], check=True)
```

## Variant Annotation

```python
def annotate_ctdna_variants(vcf_file, output_vcf):
    '''Annotate variants with clinically relevant information.'''
    # Use VEP or snpEff for annotation
    subprocess.run([
        'vep',
        '--input_file', vcf_file,
        '--output_file', output_vcf,
        '--format', 'vcf',
        '--vcf',
        '--cache',
        '--canonical',
        '--protein',
        '--sift', 'b',
        '--polyphen', 'b',
        '--af_gnomad'
    ], check=True)
```

## Tracking Known Mutations

```python
def track_specific_mutations(bam_file, mutations, min_depth=100):
    '''
    Track specific known mutations across samples.
    Useful for MRD monitoring.

    Args:
        bam_file: Aligned BAM
        mutations: List of (chrom, pos, ref, alt) tuples
    '''
    import pysam

    bam = pysam.AlignmentFile(bam_file, 'rb')
    results = []

    for chrom, pos, ref, alt in mutations:
        counts = {'ref': 0, 'alt': 0, 'other': 0}

        for pileupcolumn in bam.pileup(chrom, pos-1, pos):
            if pileupcolumn.pos != pos - 1:
                continue

            for read in pileupcolumn.pileups:
                if read.is_del or read.is_refskip:
                    continue
                base = read.alignment.query_sequence[read.query_position]
                if base == ref:
                    counts['ref'] += 1
                elif base == alt:
                    counts['alt'] += 1
                else:
                    counts['other'] += 1

        total = counts['ref'] + counts['alt'] + counts['other']
        vaf = counts['alt'] / total if total > 0 else 0

        results.append({
            'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
            'depth': total, 'alt_count': counts['alt'], 'vaf': vaf
        })

    bam.close()
    return pd.DataFrame(results)
```

## Related Skills

- cfdna-preprocessing - Preprocess with UMI consensus
- tumor-fraction-estimation - Estimate overall tumor burden
- longitudinal-monitoring - Track mutations over time
- variant-calling/variant-calling - General variant calling concepts
