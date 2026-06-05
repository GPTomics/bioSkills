---
name: bio-liquid-biopsy-pipeline
description: Cell-free DNA analysis pipeline from plasma sequencing to tumor monitoring. Preprocesses cfDNA reads, analyzes fragment patterns, estimates tumor fraction from sWGS, and optionally detects mutations from targeted panels. Use when analyzing liquid biopsy samples for cancer detection or monitoring.
tool_type: mixed
primary_tool: ichorCNA
---

## Version Compatibility

Reference examples tested with: BWA 0.7.17+, VarDict 1.8+, fgbio 2.1+, ichorCNA 0.6.0+, FinaleToolkit 0.7+, MethylDackel 0.6+, numpy 1.26+, pandas 2.2+, pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Liquid Biopsy Analysis Pipeline

**"Analyze my liquid biopsy cfDNA data end-to-end"** -> Orchestrate UMI-aware preprocessing (fgbio), ctDNA mutation detection (VarDict), tumor fraction estimation (ichorCNA), fragmentomics analysis, and longitudinal monitoring for treatment response.

Complete workflow for cfDNA analysis from sequencing to clinical interpretation.

## Pipeline Overview

```
Pre-analytical QC -> cfDNA Preprocessing -> Fragment QC
                          ↓
        ┌─────────────────┴─────────────────┐
        ↓                                   ↓
   sWGS Branch                        Panel Branch
        ↓                                   ↓
   ichorCNA                          VarDict/smCounter2
   (Tumor Fraction)                  (Mutation Detection)
        ↓                                   ↓
        └─────────────────┬─────────────────┘
                          ↓
                 Longitudinal Tracking
```

## Step 0: Pre-Analytical QC

```python
def check_preanalytical_quality(sample_metadata):
    '''
    Pre-analytical factors critical for cfDNA quality.

    Requirements:
    - Streck tube: up to 7 days at room temperature
    - EDTA tube: process within 6 hours
    - Avoid hemolysis
    - Store extracted DNA at -80C
    '''
    issues = []

    if sample_metadata['tube_type'] == 'EDTA':
        if sample_metadata['processing_delay_hours'] > 6:
            issues.append('EDTA tube processed > 6 hours - risk of gDNA contamination')

    if sample_metadata['hemolysis_score'] > 1:
        issues.append('Hemolysis detected - expect cellular DNA contamination')

    return issues
```

## Step 1: cfDNA Preprocessing with UMI Consensus

```bash
# For UMI-tagged libraries (targeted panels)
# fgbio pipeline

# Extract UMIs. Read-structure is library-specific; see liquid-biopsy/cfdna-preprocessing.
fgbio ExtractUmisFromBam \
    --input raw.bam \
    --output with_umis.bam \
    --read-structure 3M2S+T 3M2S+T \
    --single-tag RX

# Align
bwa mem -t 8 -Y reference.fa with_umis.bam | \
    samtools view -bS - > aligned.bam

# Group by UMI
fgbio GroupReadsByUmi \
    --input aligned.bam \
    --output grouped.bam \
    --strategy adjacency \
    --edits 1

# Consensus calling: keep the caller permissive (fgbio #1009), apply strictness at the filter
fgbio CallMolecularConsensusReads \
    --input grouped.bam \
    --output consensus.bam \
    --min-reads 1

# Filter: this is the real quality gate
fgbio FilterConsensusReads \
    --input consensus.bam \
    --output final.bam \
    --ref reference.fa \
    --min-reads 2
```

## Step 2: Fragment QC Checkpoint

```python
import pysam
import numpy as np

def verify_cfdna_quality(bam_path):
    '''
    QC Checkpoint: Verify cfDNA fragment profile.
    Expected: peak at ~167bp (mononucleosome)
    '''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    sizes = []

    for read in bam.fetch():
        if read.is_proper_pair and not read.is_secondary and read.template_length > 0:
            sizes.append(read.template_length)

    bam.close()
    sizes = np.array(sizes)

    modal_size = np.bincount(sizes[:400]).argmax()
    mono_frac = np.sum((sizes >= 150) & (sizes <= 180)) / len(sizes)

    qc_pass = 150 <= modal_size <= 180 and mono_frac > 0.3

    return {
        'modal_size': modal_size,
        'mononucleosome_fraction': mono_frac,
        'qc_pass': qc_pass,
        'message': 'Good cfDNA profile' if qc_pass else 'Atypical fragment distribution'
    }
```

## Step 3a: Tumor Fraction Estimation (sWGS)

ichorCNA is a command-line script (`Rscript scripts/runIchorCNA.R`), NOT an importable `runIchorCNA()` function, and it is preceded by HMMcopy `readCounter` to bin the BAM. The ~3% tumor-fraction floor is an analytical limit of detection; below it, route to fragmentomics or methylation rather than trusting a low value (see liquid-biopsy/analytical-validation and liquid-biopsy/tumor-fraction-estimation).

```bash
# For shallow WGS data (0.1-1x coverage); GavinHaLab fork
readCounter --window 1000000 --quality 20 \
    --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" \
    sample.bam > sample.wig

Rscript scripts/runIchorCNA.R \
    --id sample_id --WIG sample.wig \
    --gcWig gc_hg38_1000kb.wig --mapWig map_hg38_1000kb.wig \
    --centromere GRCh38.GCA_000001405.2_centromere_acen.txt \
    --normalPanel HD_ULP_PoN_1Mb_median.rds \
    --normal "c(0.5,0.6,0.7,0.8,0.9)" --ploidy "c(2,3)" --maxCN 7 \
    --estimateNormal TRUE --estimatePloidy TRUE --estimateScPrevalence TRUE \
    --outDir ichor_results/
# Tumor fraction = 1 - n in sample_id.params.txt
```

## Step 3b: Mutation Detection (Targeted Panel)

```bash
# For deep targeted sequencing
# Use UMI-consensus BAM from Step 1

vardict-java \
    -G reference.fa \
    -f 0.005 \
    -N sample_id \
    -b consensus.bam \
    -c 1 -S 2 -E 3 -g 4 \
    panel.bed | \
teststrandbias.R | \
var2vcf_valid.pl \
    -N sample_id \
    -E \
    -f 0.005 \
    > sample.vcf
```

## Step 4: CHIP Filtering

Clonal hematopoiesis (CHIP) is the dominant false-positive source in plasma: ~81.6% of cfDNA variants in controls and ~53.2% in cancer patients trace to white blood cells (Razavi 2019 Nat Med 25:1928). A gene-list filter is a weak fallback; the definitive control is sequencing matched buffy-coat/WBC DNA and subtracting any variant present there. See liquid-biopsy/ctdna-mutation-detection.

```python
CHIP_GENES = ['DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'JAK2', 'SF3B1', 'SRSF2', 'TP53']

def filter_chip(variants_df, wbc_variants=None, chip_genes=CHIP_GENES):
    '''Subtract WBC-matched variants when available; else fall back to a CHIP gene list.'''
    if wbc_variants is not None:
        wbc_keys = set(zip(wbc_variants['chrom'], wbc_variants['pos'], wbc_variants['alt']))
        in_wbc = variants_df.apply(lambda r: (r['chrom'], r['pos'], r['alt']) in wbc_keys, axis=1)
        return variants_df[~in_wbc], variants_df[in_wbc]

    chip = variants_df[variants_df['gene'].isin(chip_genes)]
    somatic = variants_df[~variants_df['gene'].isin(chip_genes)]
    return somatic, chip
```

## Step 5: Fragmentomics Analysis (Optional)

FinaleToolkit (MIT license, not DELFI software) exposes real hyphenated CLI subcommands and an underscored `finaletoolkit.frag` Python API; `delfi` GC-corrects the short/long ratio (raw ratios are dominated by GC and sequencing batch). DELFI is a methodology and a company, not a `pip install`-able tool.

```bash
# GC-corrected genome-wide DELFI profile and end-motif diversity.
# delfi positionals: input chrom_sizes reference bins_file; GC correction is on by default; -R for non-hg19.
finaletoolkit delfi consensus.bam hg38.chrom.sizes hg38.fa bins_100kb.bed -g gaps.bed -R -o sample.delfi.bed
finaletoolkit end-motifs consensus.bam hg38.fa -o sample.end_motifs.tsv
finaletoolkit mds sample.end_motifs.tsv
```

```python
from finaletoolkit.frag import delfi  # see liquid-biopsy/fragment-analysis

def run_fragmentomics(bam_path, chrom_sizes, reference, bins_bed, gap_bed):
    '''GC-corrected DELFI short/long profile (MDS comes from end_motifs().motif_diversity_score()).
    Python positional order is (input, chrom_sizes, bins_file, reference_file) - note this differs
    from the CLI order (input, chrom_sizes, reference, bins), so pass by keyword to be safe.'''
    return delfi(bam_path, chrom_sizes=chrom_sizes, bins_file=bins_bed,
                 reference_file=reference, gap_file=gap_bed)
```

## Step 6: Longitudinal Tracking

```python
import pandas as pd
import numpy as np

def track_longitudinal(samples_df):
    '''
    Track ctDNA over treatment.

    samples_df columns: [sample_id, timepoint, tumor_fraction, mutations...]
    '''
    samples_df = samples_df.sort_values('timepoint')

    baseline = samples_df.iloc[0]['tumor_fraction']
    samples_df['log2_fc'] = np.log2(samples_df['tumor_fraction'] / baseline)

    nadir = samples_df['tumor_fraction'].min()

    response = 'unknown'
    if nadir < 0.001:
        response = 'Complete molecular response'
    elif nadir < baseline * 0.01:
        response = 'Major molecular response (>2 log)'
    elif nadir < baseline * 0.5:
        response = 'Partial molecular response'

    return samples_df, response
```

## Complete Pipeline Script

```python
def run_liquid_biopsy_pipeline(sample_config):
    '''
    Complete liquid biopsy analysis pipeline.

    sample_config: dict with keys:
        - bam_file: Input BAM
        - data_type: 'swgs' or 'panel'
        - reference: Reference FASTA
        - bed_file: Panel BED (for panel data)
        - output_dir: Output directory
    '''
    results = {}

    # Step 1: Preprocess (if UMI data)
    if sample_config.get('has_umis'):
        preprocessed_bam = preprocess_with_fgbio(sample_config['bam_file'])
    else:
        preprocessed_bam = sample_config['bam_file']

    # Step 2: Fragment QC
    frag_qc = verify_cfdna_quality(preprocessed_bam)
    if not frag_qc['qc_pass']:
        print(f"WARNING: {frag_qc['message']}")
    results['fragment_qc'] = frag_qc

    # Step 3: Analysis based on data type
    if sample_config['data_type'] == 'swgs':
        # Tumor fraction estimation
        results['tumor_fraction'] = run_ichorcna(preprocessed_bam)
    elif sample_config['data_type'] == 'panel':
        # Mutation detection
        variants = call_variants(preprocessed_bam, sample_config['bed_file'])
        somatic, chip = filter_chip(variants)
        results['variants'] = somatic
        results['chip_variants'] = chip

    # Step 4: Optional fragmentomics
    if sample_config.get('run_fragmentomics'):
        results['fragmentomics'] = run_fragmentomics(preprocessed_bam)

    return results
```

## Related Skills

- liquid-biopsy/cfdna-preprocessing - UMI/duplex consensus error suppression
- liquid-biopsy/analytical-validation - molecule-counting limits of detection and honest LoD reporting
- liquid-biopsy/ctdna-mutation-detection - low-VAF calling and CHIP subtraction
- liquid-biopsy/tumor-fraction-estimation - ichorCNA tumor fraction from sWGS
- liquid-biopsy/fragment-analysis - fragmentomics features
- liquid-biopsy/methylation-based-detection - methylation detection and tissue-of-origin
- liquid-biopsy/longitudinal-monitoring - serial MRD tracking
