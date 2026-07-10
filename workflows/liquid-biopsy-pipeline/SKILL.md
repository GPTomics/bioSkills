---
name: bio-workflows-liquid-biopsy-pipeline
description: Orchestrates the cell-free DNA / liquid-biopsy pipeline from plasma sequencing to tumor monitoring, forking tumor-naive (screening) vs tumor-informed (MRD), and chaining pre-analytic QC, UMI/duplex error-suppression (fgbio), fragment QC, ichorCNA tumor fraction (sWGS) or VarDict low-VAF calling (panel), CHIP subtraction against matched WBC, optional fragmentomics/methylation, and longitudinal tracking. Use when treating pre-analytics as the irreversible sensitivity ceiling (tube/time-to-plasma/hemolysis), running error-suppression BEFORE calling (single-strand consensus does not remove deamination; only duplex does), reporting a VAF only with input genome-equivalents (TF ~ 2x VAF only for clonal-het-diploid), subtracting CHIP before reporting somatic, or keeping tube/panel/pipeline identical across a longitudinal MRD series. Hands mechanism to the liquid-biopsy component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: ichorCNA
goal_approach_exempt: true
workflow: true
depends_on:
  - liquid-biopsy/cfdna-preprocessing
  - liquid-biopsy/analytical-validation
  - liquid-biopsy/ctdna-mutation-detection
  - liquid-biopsy/tumor-fraction-estimation
  - liquid-biopsy/fragment-analysis
  - liquid-biopsy/methylation-based-detection
  - liquid-biopsy/longitudinal-monitoring
qc_checkpoints:
  - pre_analytics: "Tube type/time-to-plasma in window; double-spin; low hemolysis; gDNA-contamination fragment check"
  - after_consensus: "Unique-molecule coverage / genome-equivalents recovered (NOT raw depth); duplication plateaued -- read both off the GroupReadsByUmi --family-size-histogram output"
  - after_fragment: "Modal insert ~167bp (150-180); mononucleosome fraction >0.3"
  - after_tf: "TF above the ~3% ichorCNA LoD to trust the value (below = below detection, not low burden)"
  - after_mutation: "VarDict -f 0.005 is a reporting floor, NOT a detection threshold: confirm each call against a per-locus background-error model (PoN or smCounter2) before reporting; CHIP-subtracted against matched WBC"
---

## Version Compatibility

Reference examples tested with: BWA 0.7.17+, VarDict 1.8+, fgbio 2.1+, ichorCNA 0.6.0+, FinaleToolkit 0.9+ (the `delfi` param was `autosomes` before 0.9, renamed `chrom_sizes`), MethylDackel 0.6+, numpy 1.26+, pandas 2.2+, pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Liquid Biopsy Analysis Pipeline

**"Analyze my liquid biopsy cfDNA data end-to-end"** -> Orchestrate UMI-aware preprocessing (fgbio), ctDNA mutation detection (VarDict), tumor fraction estimation (ichorCNA), fragmentomics analysis, and longitudinal monitoring for treatment response.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A liquid-biopsy result is decided at four seams, two of them set before the sequencer ever runs.

1. **Pre-analytics is irreversible and is the sensitivity ceiling — set at the blood draw, not the pipeline.** WBC lysis dumps high-quality germline/CHIP DNA into the denominator; a 0.5% tumor fraction diluted to 0.1% by a processing delay is a false negative with NO bioinformatic recovery. Tube chemistry (Streck ~7d RT vs EDTA <6h), double-spin, and hemolysis are the FIRST gate the whole pipeline is conditional on. Mixing tube types within a longitudinal series is a batch confound.
2. **The error-suppression / UMI scheme is committed at library prep and cannot be upgraded after sequencing.** Single-strand UMI consensus floors error ~1e-4 to 1e-5 but does NOT remove deamination (C>T) or oxidation (G>T) — those lesions are on the template, so every PCR copy inherits them and the family votes unanimously for the artifact. Only DUPLEX (both-strand concordance) reaches <1e-7. Call on the consensus BAM, never the raw BAM.
3. **A VAF is undefined without input genome-equivalents, and TF != VAF.** 0.1% on 100 GE is noise; on 30,000 GE it is solid. The LoD is set by input GE and assay design (tumor-informed bespoke integrates 16-50 loci to ppm; tumor-naive fixed panel is error/sampling-limited ~0.1-0.5%; ichorCNA sWGS floors ~3% TF). TF ~ 2x VAF only for a clonal, heterozygous, diploid locus — copy number breaks the factor of 2.
4. **Matched buffy-coat/WBC is the CHIP-subtraction commitment.** ~81.6% of cfDNA variants in controls and ~53.2% in cancer patients are CHIP, not tumor (Razavi 2019). Without matched WBC, a tumor-naive cfDNA call is presumptively CHIP-contaminated; a gene-list filter is a weak fallback. Decide at study design.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Pre-analytics (tube/time-to-plasma/double-spin) | Irreversible TF dilution; the whole pipeline's sensitivity ceiling; consistent across a longitudinal series |
| Error-suppression scheme (single-strand vs DUPLEX) | The achievable error floor; single-strand cannot remove deamination/oxidation; cannot upgrade post-hoc |
| Assay design (tumor-naive vs tumor-informed; panel vs sWGS) + input GE | The LoD and whether it is per-locus or panel-integrated; a VAF is undefined without input GE |
| Matched WBC available? | Whether CHIP is definitively subtracted (WBC) or only gene-list-filtered (weak) |

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

# Align. bwa reads FASTQ, not BAM, and stripping to FASTQ drops the RX/UMI tag -- so emit FASTQ
# carrying RX (samtools fastq -T RX), align, then re-zip the uBAM tags back on with fgbio ZipperBams
# so GroupReadsByUmi still sees RX.
samtools fastq -T RX with_umis.bam | \
    bwa mem -C -p -t 8 -Y reference.fa - | \
    fgbio ZipperBams --unmapped with_umis.bam --ref reference.fa --output aligned.bam

# Group by UMI. --family-size-histogram is not optional here: rule 3 says a VAF is undefined without
# input genome-equivalents, and this histogram is where GE and the duplication plateau are read off.
# Raw depth after consensus collapse is not a sensitivity measure -- unique molecules are.
fgbio GroupReadsByUmi \
    --input aligned.bam \
    --output grouped.bam \
    --strategy adjacency \
    --edits 1 \
    --family-size-histogram qc/family_sizes.txt

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
        # cfDNA fragments are short; cap at 400 bp (fragments beyond are gDNA/noise) so the
        # modal-size bincount is bounded and not skewed by a few long outliers.
        if read.is_proper_pair and not read.is_secondary and 0 < read.template_length <= 400:
            sizes.append(read.template_length)

    bam.close()
    sizes = np.array(sizes)

    modal_size = np.bincount(sizes).argmax()
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
    --normalPanel HD_ULP_PoN_hg38_1Mb_normAutosomes_median.rds \   # GavinHaLab/ichorCNA extdata name; PoN build MUST match the gc/map/centromere build (the non-hg38 1Mb PoN is hg19)
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
        # REF must be in the key. Left-aligned indels share (chrom, pos, alt): a 1bp deletion
        # (REF=AT, ALT=A) and a 2bp deletion (REF=ATT, ALT=A) collide, so a (chrom, pos, alt) key
        # would silently subtract a real somatic indel as CHIP.
        key = ['chrom', 'pos', 'ref', 'alt']
        wbc_keys = set(map(tuple, wbc_variants[key].itertuples(index=False, name=None)))
        in_wbc = variants_df[key].apply(lambda r: tuple(r) in wbc_keys, axis=1)
        return variants_df[~in_wbc], variants_df[in_wbc]

    chip = variants_df[variants_df['gene'].isin(chip_genes)]
    somatic = variants_df[~variants_df['gene'].isin(chip_genes)]
    return somatic, chip
```

## Step 5: Fragmentomics Analysis (Optional)

FinaleToolkit (MIT license, not DELFI software) exposes real hyphenated CLI subcommands and an underscored `finaletoolkit.frag` Python API; `delfi` GC-corrects the short/long ratio (raw ratios are dominated by GC and sequencing batch). DELFI is a methodology and a company, not a `pip install`-able tool.

```bash
# GC-corrected genome-wide DELFI profile and end-motif diversity.
# delfi positionals: input chrom_sizes reference bins_file; GC correction is on by default.
# --no-remove-nocov keeps all bins on non-hg19 references (the default removes two hardcoded hg19 no-coverage regions).
finaletoolkit delfi consensus.bam hg38.chrom.sizes hg38.2bit bins_100kb.bed -g gaps.bed --no-remove-nocov -o sample.delfi.bed
finaletoolkit end-motifs consensus.bam hg38.2bit -o sample.end_motifs.tsv
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
        # Mutation detection. CHIP subtraction is not optional (rule 4): without matched WBC the
        # calls are presumptively CHIP-contaminated and the gene list is only a weak fallback.
        variants = call_variants(preprocessed_bam, sample_config['bed_file'])
        wbc = call_variants(sample_config['wbc_bam'], sample_config['bed_file']) if sample_config.get('wbc_bam') else None
        somatic, chip = filter_chip(variants, wbc_variants=wbc)
        results['variants'] = somatic
        results['chip_variants'] = chip

    # Step 4: Optional fragmentomics
    if sample_config.get('run_fragmentomics'):
        results['fragmentomics'] = run_fragmentomics(preprocessed_bam)

    return results
```

## Choosing the branch (tumor-naive vs tumor-informed)

Pipeline-level selection only; mechanism lives in the component skills.

| Situation | Branch | Hand off to |
|-----------|--------|-------------|
| Screening / no known tumor / MCED | tumor-NAIVE (fixed panel, sWGS, or methylation) | liquid-biopsy/methylation-based-detection |
| MRD/recurrence of a KNOWN tumor | tumor-INFORMED bespoke (16-50 clonal variants from tissue WES) | liquid-biopsy/longitudinal-monitoring |
| Tumor fraction from sWGS (0.1-1x) | ichorCNA (copy-number-based, floor ~3% TF) | liquid-biopsy/tumor-fraction-estimation |
| Low-VAF mutations from a deep panel | VarDict / smCounter2 on the UMI-consensus BAM | liquid-biopsy/ctdna-mutation-detection |
| Below mutation LoD, still need signal | fragmentomics (DELFI/end-motifs) or methylation | liquid-biopsy/fragment-analysis |
| Stating/trusting a sensitivity claim | LoB/LoD/LoD95/LoQ, per-locus vs panel-integrated | liquid-biopsy/analytical-validation |

Tumor-informed reaches ppm by integrating across 16-50 loci (beats the single-locus Poisson floor) but is BLIND by design to new/resistance variants; tumor-naive sees any panel variant but is CHIP-dominated and error-limited.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| False negative on a real low-TF sample | Bad blood draw irreversibly diluted TF | Pre-analytic gate; cannot be rescued downstream |
| Residual C>T / G>T false positives | Single-strand consensus assumed to remove damage | Duplex, or UDG for deamination; flag the substitution spectrum |
| Money spent, no sensitivity gain | Sequenced past molecular saturation | Report unique-molecule coverage; invest in plasma volume / conversion efficiency |
| Blood clones reported as tumor | CHIP not subtracted (esp. TP53/PPM1D, both CHIP and driver) | Matched buffy-coat/WBC sequencing (gene list is a weak fallback) |
| Apparent burden halved or doubled | VAF reported as TF (or vice versa) | State VAF vs TF; TF ~ 2x VAF only for clonal-het-diploid; prefer CNA-based TF when CN is non-neutral |
| A "low" TF trusted below the assay floor | ~3% ichorCNA LoD ignored | Below the floor = below detection; route to fragmentomics/methylation |
| Longitudinal "response" that is a batch shift | Tube/panel/pipeline changed between draws | Same tube, same panel, same pipeline across timepoints |

## References

- Razavi P, Li BT, Brown DN, et al (2019) High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. *Nature Medicine* 25:1928-1937. DOI 10.1038/s41591-019-0652-7. (CHIP dominance of cfDNA variants.)
- Adalsteinsson VA, Ha G, Freeman SS, et al (2017) Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. *Nature Communications* 8:1324. DOI 10.1038/s41467-017-00965-y. (ichorCNA; ~3% TF LoD.)
- Schmitt MW, Kennedy SR, Salk JJ, et al (2012) Detection of ultra-rare mutations by next-generation sequencing. *PNAS* 109:14508-14513. DOI 10.1073/pnas.1208715109. (duplex error floor.)
- Reinert T, Henriksen TV, Christensen E, et al (2019) Analysis of plasma cell-free DNA by ultradeep sequencing in patients with stages I to III colorectal cancer. *JAMA Oncology* 5:1124-1131. DOI 10.1001/jamaoncol.2019.0528. (tumor-informed bespoke MRD.)

## Related Skills

- liquid-biopsy/cfdna-preprocessing - UMI/duplex consensus error suppression
- liquid-biopsy/analytical-validation - molecule-counting limits of detection and honest LoD reporting
- liquid-biopsy/ctdna-mutation-detection - low-VAF calling and CHIP subtraction
- liquid-biopsy/tumor-fraction-estimation - ichorCNA tumor fraction from sWGS
- liquid-biopsy/fragment-analysis - fragmentomics features
- liquid-biopsy/methylation-based-detection - methylation detection and tissue-of-origin
- liquid-biopsy/longitudinal-monitoring - serial MRD tracking
