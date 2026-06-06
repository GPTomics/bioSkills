---
name: bio-long-read-sequencing-basecalling
description: Basecalls raw Oxford Nanopore signal (POD5/FAST5) into reads with Dorado, choosing the chemistry-matched model and accuracy tier (fast/hac/sup), requesting modified bases (5mCG_5hmCG, 6mA, m6A) at basecall time, and handling duplex, demultiplexing, trimming, and HERRO read correction. Covers why the model+version is an irreversible analysis decision, why methylation cannot be recovered later, and why downstream polish/variant models must match the basecaller. Use when converting POD5/FAST5 to reads, picking a Dorado model for R9/R10 or RNA004, enabling methylation calling, basecalling duplex, demultiplexing barcoded runs, or correcting reads for assembly.
tool_type: cli
primary_tool: dorado
---

## Version Compatibility

Reference examples tested with: Dorado 1.0+, pod5 0.3+, samtools 1.19+, chopper 0.7+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Results depend on inputs that outlive the binary version - record them:
- The basecaller MODEL string (e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`) sets the entire error profile and must be propagated to every downstream tool. Pin it.
- Modified-base models carry a SECOND version (`..._sup@v5.0.0_5mCG_5hmCG@v3`); the mod version can lag the simplex version - check `dorado download --list`.
- R9.4.1 and RNA002 models were removed from Dorado v1.0 defaults; legacy data needs an archived model path.

If code throws an error, introspect the installed tool (`dorado --help`, `dorado basecaller --help`) and adapt the example to the actual API rather than retrying.

# Nanopore Basecalling

**"Basecall my Nanopore data"** -> Convert raw signal (POD5) into reads with Dorado using the chemistry-matched model, deciding the accuracy tier and whether to call modifications now - because the model choice is baked irreversibly into the output.
- CLI: `dorado basecaller sup pod5s/ > calls.bam` (simplex), `dorado basecaller sup,5mCG_5hmCG pod5s/ > calls.bam` (with methylation), `dorado duplex sup pod5s/ > duplex.bam` (duplex)

PacBio note: PacBio "basecalling" (CCS -> HiFi reads) runs on-instrument/in SMRT Link; users receive HiFi BAMs already at Q20-Q30+. This skill is Oxford Nanopore / Dorado. HiFi assembly lives in genome-assembly/hifi-assembly.

## The Single Most Important Modern Insight -- There Is No "The Reads," Only "The Reads As Called By This Model"

Basecalling is not fixed preprocessing that yields a neutral FASTQ. The model and version chosen are an analysis decision written permanently into the BAM, with three consequences a naive user misses:

1. **Methylation is a basecalling decision, not a later analysis step.** Modified bases are inferred from raw signal at basecall time by Remora models and emitted as MM/ML tags. A plain BAM/FASTQ with no MM/ML tags has thrown the signal away - mods CANNOT be recovered without re-basecalling from POD5. If methylation might ever matter, request it now (`sup,5mCG_5hmCG`) and KEEP the POD5. See nanopore-methylation.
2. **Downstream polish/variant models must match the basecaller model+version.** medaka and Clair3 ship per-model weights (Clair3 `r1041_e82_400bps_sup_v500`; medaka the dotted `r1041_e82_400bps_sup_v5.2.0`). A mismatched model silently degrades accuracy with no error. Propagate the basecaller model name to every downstream step.
3. **Mixing model versions across a cohort is a batch effect.** Different model versions have different identity and homopolymer-indel error profiles. Re-basecall the WHOLE cohort with ONE current model before joint or differential analysis.

## Dorado Subcommand Taxonomy

Dorado (one GPU-first executable) replaced Guppy, which is end-of-life. Bonito is ONT's research/training basecaller (not production); Rerio hosts research-release models (niche mods, bacterial methylation).

| Subcommand | Purpose | Canonical invocation |
|------------|---------|----------------------|
| `basecaller` | simplex basecalling | `dorado basecaller hac pod5s/ > calls.bam` |
| `duplex` | template+complement duplex | `dorado duplex sup pod5s/ > duplex.bam` |
| `demux` | barcode classification/split | `dorado demux --kit-name SQK-NBD114-24 --output-dir out/ calls.bam` |
| `trim` | standalone adapter/primer trim | `dorado trim reads.bam > trimmed.bam` |
| `aligner` | minimap2 alignment (carries MM/ML) | `dorado aligner ref.mmi reads.bam > aln.bam` |
| `correct` | HERRO single-read correction | `dorado correct reads.fastq > corrected.fasta` |
| `summary` | sequencing-summary TSV from BAM | `dorado summary calls.bam > summary.tsv` |
| `download` | model management | `dorado download --model <name>` / `--list` |

## Model Naming Scheme (load-bearing)

Format `{analyte}_{pore}_{chemistry}_{speed}@v{ver}` + optional mod suffix, e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v3`.

| Token | Meaning | Examples |
|-------|---------|----------|
| analyte | molecule | `dna`, `rna004` |
| pore | flow-cell generation | `r10.4.1` (current), `r9.4.1` (legacy) |
| chemistry | kit chemistry | `e8.2` (Kit 14) |
| speed | translocation speed -> sampling rate | `400bps` (5 kHz DNA), `130bps` (RNA004, 4 kHz) |
| tier | model size/accuracy | `fast`, `hac`, `sup` |
| version | model version | `@v4.3.0`, `@v5.2.0`, `@v6.0.0` |

Passing the bare tier (`sup`) lets Dorado auto-detect chemistry from POD5 metadata and fetch the matching latest model; pin a version (`sup@v5.2.0`) or a full path for reproducibility. Append mods comma-separated (`sup,5mCG_5hmCG,6mA`); only one mod model per canonical base may be active.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Any analysis (variant/assembly/methylation) | `sup` + matched model, pinned version | `fast`/`hac` error profile leaks into calls |
| Live run / adaptive sampling / quick QC only | `fast` | speed; never for downstream analysis |
| Routine work, compute-limited | `hac` | strong accuracy/compute balance (v5.2 closed much of the gap to sup) |
| Methylation wanted now or maybe later | `sup,5mCG_5hmCG` (DNA), keep POD5 | mods are unrecoverable from a plain BAM -> nanopore-methylation |
| Per-molecule accuracy, low input, phasing | `dorado duplex sup` | ~Q30 reads, but expect <10% duplex yield |
| Diploid/phased T2T assembly from simplex | `dorado correct` (HERRO) before assembler | haplotype-aware Q22->Q40 -> genome-assembly/long-read-assembly |
| Barcoded multiplexed run | basecall `--no-trim`, then `dorado demux` | trimming first strips barcodes before demux sees them |
| Legacy R9.4.1 / RNA002 data | explicit archived model path | removed from Dorado v1.0 default downloads |
| PacBio data | already HiFi; no Dorado step | CCS runs on-instrument -> genome-assembly/hifi-assembly |

## Core Commands

```bash
# Simplex, super-accuracy, auto-detected chemistry-matched model (BAM is the default output)
dorado basecaller sup pod5s/ > calls.bam

# Pin the model version for reproducibility
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.2.0 pod5s/ > calls.bam

# Call methylation AT basecall time (CpG 5mC + 5hmC); KEEP pod5s/ - mods are unrecoverable later
dorado basecaller sup,5mCG_5hmCG pod5s/ > calls.bam
dorado basecaller sup,6mA pod5s/ > calls.bam               # all-context 6mA
# RNA004 direct RNA (cDNA CANNOT call mods - PCR erases the signal):
dorado basecaller rna004_130bps_sup@v5.1.0,m6A_DRACH pod5s/ > rna_mods.bam

# FASTQ output and a per-read quality floor (relative filter, not a calibrated accuracy)
dorado basecaller sup pod5s/ --emit-fastq --min-qscore 10 > calls.fastq

# Duplex (needs raw POD5; cannot be recovered from simplex FASTQ); dx tag marks read types
dorado duplex sup pod5s/ > duplex.bam

# Demultiplex: basecall WITHOUT trimming, then demux (demux trims barcodes itself)
dorado basecaller sup pod5s/ --no-trim > calls.bam
dorado demux --kit-name SQK-NBD114-24 --output-dir demux/ calls.bam
dorado demux --kit-name SQK-NBD114-24 --barcode-both-ends --output-dir demux/ calls.bam  # stringent

# HERRO read correction for diploid/phased assembly (input FASTQ of HAC/SUP R10 reads >=10kb -> FASTA)
dorado download --model herro-v1
dorado correct reads.fastq > corrected.fasta
```

POD5 is ONT's default raw format (faster random access than FAST5). Convert FAST5 first:

```bash
pod5 convert fast5 raw/*.fast5 --output pod5s/    # FAST5 is legacy; basecalling it directly is slow
pod5 view pod5s/                                   # summary table (replaces deprecated `pod5 inspect reads`)
pod5 merge pod5s/*.pod5 --output merged.pod5
```

## Per-Method Failure Modes

### Methylation gone forever
**Trigger:** basecalling without a mod model, then wanting 5mC later. **Mechanism:** Remora infers mods from raw signal at basecall time; a plain BAM has only bases. **Symptom:** no MM/ML tags; modkit pileup returns nothing. **Fix:** re-basecall from POD5 with `sup,5mCG_5hmCG`; keep POD5 archives.

### Barcodes land in unclassified
**Trigger:** default `--trim all` basecall, then a separate `dorado demux`. **Mechanism:** trimming removes the barcode before demux can read it. **Symptom:** most reads in `unclassified.bam`, low classification rate. **Fix:** basecall `--no-trim`, then demux (it trims barcodes itself).

### Silent accuracy loss downstream
**Trigger:** polishing/calling with a medaka/Clair3 model that doesn't match the basecaller model+version. **Mechanism:** per-model neural weights expect a specific error profile. **Symptom:** no error, just quietly worse consensus/calls. **Fix:** propagate the basecaller model name; use `medaka tools resolve_model --auto_model`; pick the matching Clair3 model dir.

### Duplex double-counting
**Trigger:** treating every read in a duplex BAM as an independent molecule. **Mechanism:** a simplex parent and its duplex offspring both appear. **Symptom:** inflated coverage/allele counts. **Fix:** the `dx:i:-1` tag marks simplex parents of duplex reads - filter them when counting molecules (`dx:i:1` = duplex, `dx:i:0` = simplex-only).

### Cohort batch effect
**Trigger:** runs basecalled with different model versions joined for analysis. **Mechanism:** version-specific identity/indel error profiles confound a technical batch with biology. **Symptom:** spurious between-run differences. **Fix:** re-basecall the whole cohort with one model version.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `sup` for any analysis | ONT model guidance | `fast`/`hac` error profiles contaminate variant/assembly/methylation calls |
| R10.4.1 SUP modal accuracy ~Q20 (99%) | Sereika 2022 | dual-reader head fixes homopolymers; enables nanopore-only near-finished genomes |
| Duplex read ~Q30; yield typically <10% of reads | community benchmarks | duplex is library-prep/loading-limited, not free accuracy |
| A "Q20" base errs at ~Q12.5 empirically | Delahaye 2021 | nanopore qscores >Q10 are overconfident posteriors; use for relative filtering only |
| HERRO input reads >=10 kbp, HAC/SUP R10 | Dorado correct docs | HERRO operates on 4096-bp chunks; shorter reads dropped |
| `--min-qscore 10` as a permissive QC floor | convention | Q10 ~ 90% nominal; a starting filter, not a hard rule |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "Failed to determine sequencing chemistry from data" | R9/RNA002 or non-standard kit; bare tier can't auto-resolve | pass an explicit model path; for legacy chemistry use an archived model |
| No MM/ML tags in BAM | basecalled without a mod model | re-basecall from POD5 with `sup,5mCG_5hmCG` |
| Most reads `unclassified` after demux | trimmed before demux | basecall `--no-trim`, then demux |
| `--model sup` errors | model is the positional arg, not a flag | `dorado basecaller sup pod5s/` |
| `dorado correct reads.bam` fails | input is FASTQ(.gz), output FASTA | `dorado correct reads.fastq > corrected.fasta` |
| Out of GPU memory | batch too large for VRAM (sup is heaviest) | lower `--batchsize`; or drop to `hac` |
| cDNA m6A calling returns nothing | PCR erased native modifications | use direct RNA (RNA004), not cDNA |

## References

- Sereika M, Kirkegaard RH, Karst SM, et al. 2022. Oxford Nanopore R10.4 long-read sequencing enables the generation of near-finished bacterial genomes from pure cultures and metagenomes without short-read or reference polishing. *Nat Methods* 19:823-826.
- Stanojević D, Lin D, Nurk S, Florez de Sessions P, Šikić M. 2026. Telomere-to-telomere assembly using HERRO-corrected Nanopore simplex reads. *Nature* (online ahead of print). DOI 10.1038/s41586-026-10563-y.
- Wick RR, Judd LM, Holt KE. 2019. Performance of neural network basecalling tools for Oxford Nanopore sequencing. *Genome Biol* 20:129.
- Pagès-Gallego M, de Ridder J. 2023. Comprehensive benchmark and architectural analysis of deep learning models for nanopore sequencing basecalling. *Genome Biol* 24:71.
- Delahaye C, Nicolas J. 2021. Sequencing DNA with nanopores: troubles and biases. *PLoS ONE* 16(10):e0257521.
- Gamaarachchi H, Samarakoon H, et al. 2025. The enduring advantages of the SLOW5 file format for raw nanopore sequencing data. *GigaScience* giaf118.

## Related Skills

- long-read-qc - Assess read length/quality and run health after basecalling
- nanopore-methylation - Pile up the MM/ML tags this skill must request at basecall time
- long-read-alignment - Map the reads; use `-y` to carry MM/ML tags through alignment
- medaka-polishing - Consensus model that must match this basecaller model+version
- clair3-variants - Variant model that must match this basecaller model+version
- genome-assembly/long-read-assembly - Assemble the reads (HERRO-corrected for diploid/T2T)
- genome-assembly/hifi-assembly - PacBio HiFi (basecalled on-instrument, not here)
- epitranscriptomics/m6anet-analysis - ONT direct-RNA m6A from signal
- workflows/longread-sv-pipeline - End-to-end basecall -> align -> SV call
