# ctDNA Mutation Detection - Usage Guide

## Overview
Detect somatic mutations in circulating tumor DNA at low variant allele fractions, distinguishing de novo calling (scanning a panel for unknown variants) from tumor-informed tracking (quantifying a pre-specified variant set). Low-VAF detection is a signal-versus-noise problem set by error suppression and the number of molecules sampled, with clonal hematopoiesis (CHIP) as the dominant false positive requiring matched-WBC subtraction.

## Prerequisites
```bash
# VarDict (de novo calling)
conda install -c bioconda vardict-java

# GATK Mutect2 (de novo calling with orientation-bias model)
conda install -c bioconda gatk4

# UMI-aware caller (optional); not on PyPI, install from source
# git clone https://gitlab.com/vincent-sater/umi-varcal-master

# Python dependencies for known-variant tracking
pip install pysam pandas
```

## Quick Start
Tell your AI agent what you want to do:
- "Run de novo low-VAF calling on my UMI-consensus panel BAM"
- "Call variants with Mutect2 and filter FFPE/OxoG orientation artifacts"
- "Track this list of known tumor mutations across my serial plasma samples"
- "Tell me whether this low-VAF TP53 call is tumor or CHIP"
- "Help me choose a VAF threshold for my consensus depth and input mass"

## Example Prompts

### De Novo Calling
> "Run VarDict at 0.5% VAF on my UMI-consensus targeted-panel BAM and convert to VCF."

> "Call somatic variants with Mutect2 tumor-only, learn the read-orientation model, and filter with a panel of normals."

### Tumor-Informed Tracking
> "I have 30 clonal variants from the patient's tumor WES. Quantify their VAF in this plasma BAM for MRD."

> "Aggregate alt-read support across my patient-specific loci to make a panel-level tumor-present call."

### CHIP and Interpretation
> "Subtract the matched buffy-coat genotype and flag any remaining CHIP-gene variants."

> "Separate germline, CHIP, and candidate tumor variants using gnomAD and matched-WBC presence."

> "Annotate the surviving calls with VEP including gnomAD exome and genome frequencies."

## What the Agent Will Do
1. Decide de novo calling vs tumor-informed tracking from whether a known variant set exists
2. Confirm the input is error-suppressed (UMI/duplex consensus) appropriate to the target VAF regime
3. Run VarDict or Mutect2 for de novo calling, or pileup fixed loci for tracking
4. Subtract matched WBC and flag CHIP-gene variants before reporting anything as somatic-tumor
5. Annotate surviving variants for clinical interpretation
6. Report per-locus vs panel-integrated detection limits with input mass

## Tips
- **Detection is not calling** - tracking a known set integrates across loci to ppm; scanning a panel is bounded by per-locus error and multiple testing.
- **Match the consensus to the VAF** - above 1% any caller works; 0.1-1% needs single-strand UMI consensus; below 0.1% needs duplex plus tumor-informed integration.
- **CHIP is the null** - most non-germline cfDNA variants are not tumor; matched-WBC subtraction is mandatory, gnomAD filtering removes germline only.
- **Watch the VarDict BED flags** - `-c -S -E -g` are column indices, not coordinates; match var2vcf_valid.pl `-f` to VarDict's `-f`; add `-P 0` for amplicons.
- **Do not call below the error floor** - lowering `-f` past the demonstrated input error rate just manufactures calls that scale with depth.
- **Quote LoD honestly** - a bare VAF without input mass and replicate detection rate is not a sensitivity spec.

## Related Skills
- cfdna-preprocessing - UMI/duplex consensus input that sets the error floor
- analytical-validation - LoD/LoB and the panel-integration math behind detection
- longitudinal-monitoring - track detected variants across serial samples
- tumor-fraction-estimation - orthogonal burden estimate to cross-check
- variant-calling/variant-calling - general somatic calling principles
- clinical-databases/variant-prioritization - clinical annotation and interpretation
