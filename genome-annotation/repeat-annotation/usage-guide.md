# Repeat Annotation - Usage Guide

## Overview

Discover, classify, and mask repetitive and transposable elements in genome assemblies with the RepeatModeler2 (de novo library) + RepeatMasker (masking) pair, EarlGrey (auto-curating wrapper, the non-model default), or EDTA (plant/structural TEs). Soft-masked output is a prerequisite for eukaryotic gene prediction. The library is the experiment: an uncurated de novo library still masks roughly the right real estate (so % looks fine) while silently corrupting classification, age landscapes, and - worst - masking real gene families. Assembly quality caps everything: short reads collapse the youngest, most active TEs.

## Prerequisites

```bash
# Classic pair
conda install -c bioconda repeatmodeler repeatmasker

# Non-model wrapper / plant pipeline
conda install -c bioconda earlgrey edta

# TE expression (EM multimapper handling)
pip install TEtranscripts

# Python utilities
pip install pandas matplotlib biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Soft-mask repeats in my genome assembly before gene prediction"
- "Build a de novo TE library for my non-model insect genome with EarlGrey"
- "Annotate transposable elements in my plant genome with EDTA"
- "Quantify transposable element expression from RNA-seq without unique-only mapping"

## Example Prompts

### Repeat Masking

> "Run RepeatModeler2 to build a de novo library, decontaminate it against a protein DB, then RepeatMasker to soft-mask my genome with -xsmall for gene prediction."

> "Soft-mask my vertebrate assembly using the Dfam library, interspersed repeats only (skip low-complexity)."

### TE Analysis

> "Classify the transposable elements in my genome, report content by class, and tell me what fraction is Unknown."

> "Generate a Kimura repeat-divergence landscape and tell me whether there's a recent TE burst."

### TE Expression

> "Quantify TE expression in treated vs control with TEtranscripts (subfamily level) and flag any read-through confounds."

## What the Agent Will Do

1. Check the assembly type (short-read vs HiFi/T2T) and flag that short reads under-count TEs and bias age toward old
2. Build a de novo library (RepeatModeler2 `-LTRStruct`, or EarlGrey/EDTA)
3. Decontaminate the library against a protein DB so real gene families (NLR/ZNF/OR, domesticated TEs) are not masked away
4. Soft-mask (`-xsmall`) genome-wide, union de novo + Dfam, interspersed repeats only before gene prediction
5. Report content by class, the Unknown fraction, and a Kimura landscape (a relative, right-censored age readout)
6. For TE expression, use EM tools and separate autonomous transcription from read-through

## Tips

- **Soft-mask, never hard-mask before gene prediction** - `-xsmall` keeps sequence; hard-masking (N) silently deletes genes overlapping repeats with no error or log.
- **Decontaminate the library first** - A de novo library contains real multicopy gene families (KRAB-ZNF, NLR, RAG1, CENP-B, syncytins) that look repetitive; masking them deletes the genome's most interesting genes. BLAST the library against a protein DB and drop host-gene hits with no TE domain.
- **Curation gates biological claims** - Masking and gross % can use an automated library; any per-family claim (young/active/novel) needs curation (Goubert protocol, TE-Aid). RepeatModeler2 output is a draft.
- **The library, not the masker, sets quality** - and masking % is robust to a bad library, so the headline number hides the damage.
- **Assembly is destiny** - Short-read assemblies collapse and drop young copies; always report what assembly a % repeat came from. LAI (<10 draft / 10-20 reference / >20 gold) gauges LTR-RT resolution but only for LTR-rich genomes.
- **% repeat is not comparable across studies** - It depends on library + engine (rmblast vs nhmmer) + assembly; report all three.
- **TE expression needs EM tools** - Unique-only mapping discards most TE signal; a TE in an intron is not "expressed" (read-through). Use TEtranscripts/SQuIRE and separate autonomous from passenger signal.

## Related Skills

- genome-annotation/eukaryotic-gene-prediction - Requires soft-masked genome
- genome-annotation/annotation-qc - LAI and repeat-content sanity checks
- genome-assembly/assembly-qc - Assembly quality sets the TE-annotation ceiling
- differential-expression/deseq2-basics - Differential TE expression from EM counts
