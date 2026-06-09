---
name: bio-metagenomics-contamination-controls
description: Cleans a shotgun metagenome of everything that is not the target community before profiling - host-read depletion (Hostile, bowtie2/T2T-CHM13), reagent/kitome contamination control with blanks and decontam, mock-community validation, and depth-adequacy checks (Nonpareil). Covers why a metagenomic result is a position in a choice-chain rather than a direct observation, why extraction is the experiment, why a low-biomass community can be entirely kitome, why absence means not-detectable-by-this-chain, and why a confident classifier call can still be wrong when the reference is contaminated. Use when designing controls, removing host reads, identifying reagent contaminants, validating with mocks, or judging whether a low-biomass result is real. For adapter/quality trimming see read-qc; for MAG-level decontamination see genome-assembly/metagenome-assembly.
tool_type: mixed
primary_tool: decontam
---

## Version Compatibility

Reference examples tested with: decontam 1.22+, Hostile 1.1+, Bowtie2 2.5+, Nonpareil 3.4+, pandas 2.2+, R 4.3+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('decontam')` then `?isContaminant` to verify parameters
- CLI: `hostile --version`, `nonpareil -h` to confirm flags and indexes
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The controls define the result: an extraction blank defines the kitome for its batch/lot, a mock community defines the limit of detection and the extraction-lysis bias, and the host reference (prefer T2T-CHM13 over GRCh38) defines what host is removed. Record the extraction kit and lot, the host index, the blanks, the mock version (whole-cell vs DNA), and the reads removed at each step.

# Contamination Controls

**"Is this signal real, or did my pipeline create it?"** -> Remove host reads, define the kitome with blanks, validate with a mock, and confirm depth - because a low-biomass community can be entirely reagent contamination.
- R: `decontam::isContaminant(seqtab, conc=, neg=, method='combined')` on the classifier output table
- CLI: `hostile clean --fastq1 R1.fq.gz --fastq2 R2.fq.gz --index human-t2t-hla`

Scope: sample-level pre-analysis cleanup and controls - host depletion, kitome/blank/mock controls, decontam, depth adequacy. Adapter/quality trimming mechanics -> read-qc/adapter-trimming, read-qc/quality-filtering. MAG-level decontamination (CheckM2/GUNC chimerism, FCS-GX foreign sequence) -> genome-assembly/metagenome-assembly - a different, genome-level problem. Classification -> kraken-classification, metaphlan-profiling.

## The Single Most Important Modern Insight -- A Metagenomic Result Is a Position in a Choice-Chain

A metagenomic profile is the product of a chain of choices - extraction, host/contaminant depletion, depth, read-vs-assembly, classifier, database, normalization - and each link silently sets what is observable. The community is never observed directly; the report is the community as refracted by this pipeline. Three consequences a newcomer misses:

1. **There is no raw truth to recover.** Even the input DNA is already a biased sample of the cells (lysis bias). Cleaning the data does not get closer to the community; it changes which lens dominates. The honest framing is relative-within-a-consistent-pipeline.
2. **Absence means not-detectable-by-this-chain.** A zero is below the depth detection limit, OR not in the database, OR lost in extraction, OR removed by depletion - almost never simple biological absence. Force the question "which link is responsible?" before interpreting absence.
3. **The pipeline can manufacture the result.** In low biomass the entire community can BE the kitome; with the wrong database the profile is the database's bias; over-aggressive depletion deletes real taxa. Controls plus a consistent chain plus explicit reporting are what convert an uninterpretable number into a defensible measurement.

## Extraction Is the Experiment

Lysis efficiency is taxon-dependent: tough-walled Gram-positives (Firmicutes, *Staphylococcus*, *Enterococcus*), endospores (*Bacillus*, *Clostridium*), acid-fast *Mycobacterium*, and fungi/archaea resist lysis and are under-represented unless bead-beating is used. Gentle/enzymatic kits inflate easy-to-lyse Gram-negatives - so a Firmicutes:Bacteroidetes shift can be an extraction artifact. Extraction had the largest effect on observed composition across 21 protocols (Costea 2017 *Nat Biotechnol* 35:1069). Use bead-beating, hold one method constant across a study, validate lysis with a whole-cell mock, and report kit and lot. Note the tradeoff: aggressive bead-beating shears DNA and hurts long-read assembly, so the best extraction depends on the read-vs-assembly choice.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Host-associated sample (gut, oral, skin, tissue) | host-deplete first (Hostile, T2T-CHM13) | host reads waste depth and leak false calls; also a data-sharing/ethics requirement |
| Low-biomass sample (skin, BAL, CSF, tissue, blood) | blanks + DNA quantification + decontam mandatory | the lower the biomass, the larger the kitome fraction |
| Novel taxa claimed in low biomass | treat as kitome until proven (canonical genera) | placenta/tumor "microbiomes" were largely kitome |
| Need limit of detection / lysis check | run a mock (ZymoBIOMICS whole-cell) | the only sample with a known answer |
| Is my depth enough for this question? | Nonpareil coverage curve | depth sets the detection limit; host depletion halves usable depth |
| Confident classifier call, odd taxon | suspect a contaminated reference | confidence is not correctness if the reference is mislabeled |
| Adapter/quality trimming | -> read-qc | this skill owns metagenomics-specific cleanup, not generic trimming |
| MAG chimerism / foreign sequence in a bin | -> genome-assembly/metagenome-assembly | genome-level decontamination is a different problem |

## Host-Read Depletion

```bash
# Hostile removes >99.5% of human reads while discarding far fewer microbial reads than naive mapping.
# Prefer the T2T-CHM13-based index over GRCh38; high-sensitivity Bowtie2 drives removal more than the reference.
hostile clean --fastq1 sample_R1.fq.gz --fastq2 sample_R2.fq.gz \
    --index human-t2t-hla --aligner bowtie2
# Report the reads removed - it is a QC metric, not a footnote. For long reads use --aligner minimap2.
```

Remove host for two reasons: analytical (depth, false positives, runtime) and ethical (raw human-associated reads carry identifiable host genotype; depleting before deposit is increasingly required). Wet-lab depletion (saponin/DNase, methyl-CpG capture) saves sequencing but adds its own bias - a genuine tradeoff.

## Identify Reagent Contaminants with decontam

**Goal:** Separate real low-abundance taxa from the kitome using blanks and DNA concentration.

**Approach:** Run decontam on the classifier output table (taxa x samples) using the frequency signal (contaminants scale inversely with input DNA) and the prevalence signal (contaminants are enriched in blanks); raise the prevalence threshold for low biomass.

```r
library(decontam)
# seqtab: samples x taxa from the Bracken/MetaPhlAn table; conc: per-sample DNA concentration; neg: TRUE for blanks.
contam <- isContaminant(seqtab, conc = dna_conc, neg = is_blank, method = 'combined', threshold = 0.1)
# Low-biomass studies: use the prevalence method at the more aggressive threshold 0.5, and inspect the calls.
contam_lowbio <- isContaminant(seqtab, neg = is_blank, method = 'prevalence', threshold = 0.5, batch = batch_id)
seqtab_clean <- seqtab[, !contam$contaminant]
```

decontam runs per batch (`batch=`) because the kitome differs by lot/run. Always inspect the called contaminants against the canonical kitome genera (*Bradyrhizobium*, *Ralstonia*, *Burkholderia*, *Pseudomonas*, *Acinetobacter*, *Sphingomonas*, *Methylobacterium*, *Stenotrophomonas*) rather than applying blindly - over-aggressive removal deletes real taxa.

## Depth Adequacy

```bash
# Nonpareil estimates how much of the community's sequence space you have sampled, without assembly or a DB.
nonpareil -s reads.fasta -T kmer -f fasta -b sample_np
# Plot the coverage-vs-effort curve in R (Nonpareil.curve); a non-detection below the implied limit is meaningless.
```

Depth is set by the question: dominant taxa need a few million reads, rare-pathogen detection sets a limit of detection, and strain SNVs need high per-genome coverage. Host depletion can silently halve usable depth - budget for it.

## Per-Method Failure Modes

### Extraction bias reported as biology
**Trigger:** a Firmicutes:Bacteroidetes shift or "low Gram-positive" community from a gentle-lysis kit. **Mechanism:** taxon-dependent lysis under-represents tough-walled organisms. **Symptom:** composition differences tracking the kit, not the sample. **Fix:** bead-beating, one method held constant, a whole-cell mock to prove hard taxa are lysed.

### Kitome called as novel taxa in low biomass
**Trigger:** reporting novel low-abundance taxa from skin/BAL/tissue/blood without blanks. **Mechanism:** reagent DNA is a fixed dose; at low biomass it dominates the signal. **Symptom:** canonical kitome genera presented as discovery. **Fix:** extraction blanks through the full workflow, DNA quantification, decontam (prevalence + frequency), skepticism toward the kitome genera.

### Absence read as biological absence
**Trigger:** "taxon/function not present" or "low diversity." **Mechanism:** detection is bounded by depth, database, extraction, and depletion. **Symptom:** a negative interpreted as biology. **Fix:** report the classified fraction and the limit of detection; state which link is responsible before interpreting absence.

### Confident call from a contaminated reference
**Trigger:** trusting a high-confidence classifier assignment. **Mechanism:** >2 million GenBank/RefSeq entries carry mislabeled or chimeric sequence (Steinegger & Salzberg 2020 *Genome Biol* 21:115). **Symptom:** a confident, systematic wrong assignment (the classic stray human/vector in a microbial genome). **Fix:** treat confidence as not equal to correctness; cross-check surprising calls against a cleaner database.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| decontam threshold 0.1 default; 0.5 prevalence for low biomass | Davis 2018 *Microbiome* 6:226 | aggressive prevalence call needed when the kitome dominates |
| >= 1 extraction blank per batch | Salter 2014 *BMC Biol* 12:87 | blanks define the kitome for that lot; more for very low biomass |
| Hostile removes > 99.5% host | Constantinides 2023 *Bioinformatics* 39:btad728 | high host removal with low microbial loss |
| Prefer T2T-CHM13 over GRCh38 + mask rDNA | host-removal practice | GRCh38 gaps let host reads escape; rDNA masking spares microbial reads |
| Whole-cell vs DNA mock | mock-standard practice | whole-cell tests extraction/lysis; DNA tests classifier/library only |
| Nonpareil coverage before interpreting absence | Rodriguez-R 2018 *mSystems* 3:e00039-18 | a non-detection below the limit of detection is uninformative |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| decontam finds nothing useful | no blanks or DNA concentration supplied | add blanks (`neg=`) and/or DNA quant (`conc=`); run per batch |
| Host removal leaves human reads | GRCh38 with gaps, low-sensitivity aligner | use a T2T-CHM13 index and high-sensitivity Bowtie2 |
| Real microbes deleted in host removal | rDNA / conserved regions not masked | mask host rDNA; check microbial reads removed |
| Low-biomass "novel taxon" not reproducible | kitome | blanks + decontam; check canonical kitome genera |
| Cross-study profiles disagree | different extraction/depth/DB chains | hold the chain constant; do not meta-analyze across links |

## References

- Salter SJ, Cox MJ, Turek EM, et al. 2014. Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biol* 12:87.
- Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. 2018. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. *Microbiome* 6:226.
- Costea PI, Zeller G, Sunagawa S, et al. 2017. Towards standards for human fecal sample processing in metagenomic studies. *Nat Biotechnol* 35:1069-1076.
- Steinegger M, Salzberg SL. 2020. Terminating contamination: large-scale search identifies more than 2,000,000 contaminated entries in GenBank. *Genome Biol* 21:115.
- Constantinides B, Hunt M, Crook DW. 2023. Hostile: accurate decontamination of microbial host sequences. *Bioinformatics* 39:btad728.
- Rodriguez-R LM, Gunturu S, Tiedje JM, Cole JR, Konstantinidis KT. 2018. Nonpareil 3: fast estimation of metagenomic coverage and sequence diversity. *mSystems* 3:e00039-18.

## Related Skills

- kraken-classification - Classification after host removal; database bias and contaminated references
- metaphlan-profiling - Marker-gene profiling after cleanup
- abundance-estimation - decontam runs on the classifier output abundance table
- metagenome-visualization - Plot blanks alongside samples; depth-adequacy curves
- read-qc/adapter-trimming - Generic adapter/quality trimming before this step
- genome-assembly/metagenome-assembly - MAG-level decontamination (a different, genome-level problem)
- workflows/metagenomics-pipeline - End-to-end pipeline with a controls/depletion stage up front
