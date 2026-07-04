# MiXCR Analysis - Usage Guide

## Overview

MiXCR aligns immune-repertoire reads to V/D/J/C germline and assembles them into clonotypes (CDR3 + V + J). In MiXCR 4.x correctness is dominated by one decision: the preset, which encodes the library chemistry (RNA vs gDNA, 5'RACE vs multiplex primers, bulk vs 10x single-cell, UMI vs no-UMI). The wrong preset does not error -- it silently produces plausible-but-wrong clonotypes (truncated CDR3, mis-called V, inflated diversity from uncollapsed UMIs), so the analytical work is choosing and auditing the preset, then reporting the correct quantitation unit (reads for non-UMI bulk, molecules for UMI, cells for single-cell). MiXCR 4.x also requires an activated license (free for academics) or it refuses to run, and the D-segment call in TRB/IGH is near-unassignable and must never key a clonotype.

## Prerequisites

```bash
conda install -c bioconda mixcr
# or download the release jar: https://github.com/milaboratory/mixcr/releases
# MiXCR 4.x needs Java 17 (mixcr --version prints the JVM)

# License is mandatory for 4.x (academic use is free):
#   get a key at https://platforma.bio/getlicense
mixcr activate-license          # paste the key, or:
export MI_LICENSE_FILE=/path/to/mi.license
```

## Quick Start

Tell your AI agent what you want to do:
- "Pick the right MiXCR preset for my Takara SMARTer TCR kit and run it"
- "Process my 10x VDJ FASTQ with MiXCR and pair the chains"
- "Assemble clonotypes from bulk TCR-seq and export CDR3 with V/J usage"
- "Mine TCR clonotypes from my bulk RNA-seq BAM/FASTQ"
- "Export my clonotypes as AIRR TSV for Immcantation"
- "Check the alignment rate and chain composition QC on my run"

## Example Prompts

### Preset selection and running

> "This is a 5'RACE template-switch TCR library with a 12nt UMI -- which MiXCR preset, and run it on human data."

> "Run MiXCR on my multiplex-primer BCR amplicon and set the floating boundary on the primer side."

> "Process my 10x Genomics single-cell VDJ FASTQ and assemble paired-chain cells."

### Export and handoff

> "Export TRB clonotypes with CDR3 nucleotide and amino-acid sequences, V and J genes, and the UMI count."

> "Give me AIRR rearrangement TSV from these clones so I can run Immcantation lineage analysis."

> "Restrict the export to productive clonotypes only."

### Quality control

> "Show the MiXCR alignment-rate QC and tell me if the preset matches the chemistry."

> "Check the chain-usage QC for cross-contamination between my TRB samples."

> "Why did my RNA-seq run produce so few clonotypes?"

## What the Agent Will Do

1. Confirm a MiXCR license is activated (4.x refuses to run otherwise).
2. Match a preset to the exact library chemistry (material, boundary model, UMI/cell barcodes, bulk vs single-cell) and audit it with `mixcr exportPreset`.
3. Run `mixcr analyze <preset>` (or the hand-run align -> refineTagsAndSort -> assemblePartial/extend -> assemble -> assembleCells chain), adding `--assemble-clonotypes-by` when the preset lacks an intrinsic assembling feature.
4. Export clonotypes as native MiXCR TSV (`exportClones`) or AIRR (`exportAirr`), selecting the correct quantitation field (reads vs UMIs vs cells).
5. Run `mixcr qc` and `exportQc align`/`chainUsage` and interpret alignment rate and chain composition against the chemistry.

## Tips

- The preset is the analysis. Match it to the exact kit/chemistry; a mismatch is silent, not an error. Audit with `mixcr exportPreset --preset-name <name>`.
- Activate the license before anything else, and whitelist the phone-home IPs on firewalled clusters or a job stalls waiting on egress.
- Report the right denominator: reads for non-UMI bulk, `uniqueMoleculeCount` for UMI libraries, cells for single-cell. Reporting reads on a UMI library re-adds PCR bias.
- Never analyze a UMI/10x library without a tag pattern and `refineTagsAndSort` -- uncollapsed barcode errors become fake clonotypes and inflate diversity.
- For RNA-seq/fragmented data, align with `--keep-non-CDR3-alignments`, run `assemblePartial` twice, then `extend`; judge success by absolute clonotype yield, not percent aligned.
- Never trust or key on the D-gene call in TRB/IGH -- it is near-unassignable.
- From MiXCR 4.7, add `--assemble-clonotypes-by CDR3` (or `VDJRegion`) if a preset has no intrinsic assembling feature.
- `extend` on BCR can fabricate germline over hypermutated ends; use cautiously and prefer longer reads for SHM work.

## Related Skills

- vdjtools-analysis - Downstream diversity and overlap on bulk clonotypes
- immcantation-analysis - BCR clonal clustering, SHM and lineage from AIRR output
- scirpy-analysis - Single-cell VDJ integration with gene expression
- repertoire-visualization - Plot V/J usage and clonal structure
- specificity-annotation - Antigen-specificity clustering and database lookup
- read-qc/adapter-trimming - Upstream read QC and adapter handling
- workflows/tcr-pipeline - End-to-end orchestration
