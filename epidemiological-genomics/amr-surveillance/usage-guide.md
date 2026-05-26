# AMR Surveillance Usage Guide

## Overview

This skill drives postdoc-grade antimicrobial-resistance surveillance from bacterial WGS. It covers acquired AMR gene detection (AMRFinderPlus, ResFinder 4.0, CARD-RGI, abricate), chromosomal point-mutation calling (PointFinder, AMRFinderPlus species mode, TB-Profiler, Mykrobe), mobile-genetic-element context (MOB-suite, PlasmidFinder, MobileElementFinder), and phenotype prediction against EUCAST / CLSI breakpoints. Cross-tool reporting is harmonised through hAMRonization to the PHA4GE schema, and TB drug-resistance prediction is interpreted against the WHO 2nd-edition catalogue. The skill is decision-grade: it tells the agent when to switch from AMRFinderPlus to TB-Profiler, how to handle WHO Group 3 (Uncertain) calls without misleading clinicians, when heteroresistance demands deep read-based variant calling, and when to insist on long-read assembly for plasmid resolution.

Distinct from `metagenomics/amr-detection` (community-level ARG quantification, not isolate-focused) and from `variant-calling/variant-calling` (which owns SNP-calling mechanics that this skill consumes).

## Prerequisites

```bash
conda install -c bioconda ncbi-amrfinderplus resfinder rgi abricate hamronization abritamr staramr tb-profiler mykrobe mob_suite plasmidfinder mefinder lofreq snippy

amrfinder -u
tb-profiler update_tbdb

pip install pandas biopython
```

For long-read or hybrid assembly upstream:

```bash
conda install -c bioconda flye unicycler medaka
```

## Quick Start

Tell the AI agent what surveillance question to answer:

- "Detect acquired AMR genes and chromosomal point mutations in this Klebsiella pneumoniae assembly"
- "Predict drug-susceptibility phenotype against EUCAST 2024 breakpoints for these E. coli bloodstream isolates"
- "Run TB-Profiler against the WHO 2nd-edition catalogue and report per-drug resistance with Group 1/2/3/4/5 grading explicit"
- "Identify carbapenemase variants in these 40 hospital isolates and tell me which are plasmid-borne via MOB-suite"
- "Harmonise AMRFinderPlus + ResFinder + RGI output for 200 isolates via hAMRonization in PHA4GE format for WHO GLASS submission"
- "Screen these isolates for heteroresistance to colistin using deep read-based variant calling at 1% MAF"
- "Cross-check TB-Profiler and Mykrobe for any XDR-TB call before clinical handover"

## Example Prompts

### WGS-to-WHO-GLASS clinical surveillance report

> "We have 250 Illumina paired-end runs of clinical *E. coli* and *Klebsiella pneumoniae* from a tertiary hospital ICU surveillance project spanning Q1-Q4 2025. Assemble with Unicycler, run AMRFinderPlus with `--organism` set per isolate plus `--plus`, run ResFinder 4.0 with `--acquired --point -s 'species'`, harmonise via hAMRonization, and produce a long-format CSV mapping isolate -> antibiotic class -> highest-confidence determinant -> predicted phenotype, in WHO GLASS reporting categories (S / I / R / Insufficient evidence) against EUCAST 2025 breakpoints. Document every tool version and database date."

### Mobile-genetic-element context for an NDM-1 outbreak

> "Twenty *Klebsiella pneumoniae* isolates from a suspected NDM-1 outbreak across two wards. For each, identify the carbapenemase via AMRFinderPlus, reconstruct plasmid contigs via MOB-suite `mob_recon`, type each plasmid via `mob_typer`, check for IS*26* / IS*Kpn7* mobilisation via MobileElementFinder, and tell me whether the NDM-1-bearing plasmid is the same replicon and MOB-suite cluster across isolates. Recommend a long-read confirmation isolate if plasmid context appears fragmented in any short-read assembly (`Method=PARTIAL_CONTIG_END` in AMRFinderPlus output)."

### Mycobacterium tuberculosis drug-resistance prediction

> "Sixteen MTB clinical assemblies from a high-burden setting. Run TB-Profiler with the current WHO 2nd-edition catalogue; report per-drug resistance call (R / R-interim / Uncertain / S / Not predicted) for the 13 anti-TB drugs in the catalogue, flag Coll/Napier lineage assignment, and call out high-confidence rpoB / katG / pncA / gyrA / inhA mutations explicitly. Run Mykrobe as a cross-check on any R or XDR call. Do NOT collapse Group 3 (Uncertain) mutations to 'S' -- report them as 'phenotypic DST recommended'."

### Heteroresistance screen for colistin

> "Five *Klebsiella pneumoniae* clinical isolates with treatment failure on colistin. Run deep read-based variant calling with lofreq at MAF >= 1% across mgrB / pmrAB / phoPQ, report any non-reference allele frequency between 1% and 20% as heteroresistance candidates, and compare against the assembly-based AMRFinderPlus output. For any heteroresistance candidate flag for targeted deep amplicon sequencing."

### mcr surveillance with phenotypic context

> "Hundred isolates from a one-year colistin-resistance surveillance project across food-animal and clinical sources. Detect mcr-1 through mcr-10 acquired genes via AMRFinderPlus; for each mcr-9 and mcr-10 hit, flag that phenotypic colistin MIC confirmation is required before reporting as colistin-R. Cross-reference with MOB-suite to determine plasmid host-range. Produce a per-sample table: mcr variant, plasmid replicon, predicted host range, phenotypic confirmation status."

### Cross-tool reconciliation for a single isolate

> "I have AMRFinderPlus calling `bla_NDM-5`, ResFinder calling `bla_NDM-1`, and RGI calling `bla_NDM-1`. Pass everything through hAMRonization, reconcile against NCBI ReferenceGeneCatalog nomenclature, and tell me whether the discordance is real allele difference or database-naming drift. If real, explain the MIC implications."

## What the Agent Will Do

1. Identify the surveillance question (presence vs phenotype prediction vs MGE context vs heteroresistance vs cross-lab reporting) and the species.
2. Select the primary tool by the decision tree: AMRFinderPlus with `--organism` for *E. coli* / *Salmonella* / *Klebsiella* / *Staphylococcus* / etc.; TB-Profiler for *M. tuberculosis*; ResFinder 4.0 or abritAMR for phenotype prediction.
3. Pin tool versions AND database versions; record AMRFinderPlus DB date, TB-Profiler WHO catalogue edition, MOB-suite version.
4. For Mtb, run TB-Profiler primary + Mykrobe cross-check; interpret against WHO 2nd-edition catalogue; report Group 3 mutations as "uncertain significance".
5. For mobility claims, run MOB-suite `mob_recon` + `mob_typer` and recommend long-read confirmation if short-read assembly produces fragmented plasmid contigs (`Method=PARTIAL_CONTIG_END`).
6. For heteroresistance, supplement assembly-based calls with deep read-based variant calling (lofreq at MAF >= 1%) for clinically critical drugs.
7. Harmonise all per-tool output via `hamronize` with mandatory PHA4GE provenance metadata, then `hamronize summarize` for the cohort table.
8. Translate gene-level calls to drug-class -> EUCAST/CLSI breakpoint-year-anchored S/I/R; populate WHO GLASS categories.
9. Flag any reconciliation ambiguity (AMRFinderPlus vs ResFinder allele assignment, RGI tier semantics) for analyst review rather than silent resolution.

## Tips

- AMRFinderPlus has no Mtb organism mode in any v4.x release. Calling `amrfinder -n mtb.fa` returns an empty AMR table on an XDR-TB genome. Use TB-Profiler.
- WHO 2nd-edition catalogue (2023) supersedes the 2021 1st edition. The 3rd edition is anticipated late 2025-2026; check `tb-profiler list_db` for the bundled edition before clinical handover.
- Never report OXA-48-like at family level. The allele identity (OXA-48 vs OXA-181 vs OXA-232 vs OXA-244) determines whether ertapenem alone is elevated, whether ceftazidime-avibactam will work, and whether the strain looks phenotypically susceptible.
- `Method=PARTIAL_CONTIG_END` in AMRFinderPlus output is the smoking gun for plasmid-context fragmentation in short-read assembly. Re-assemble with long-read (Nanopore R10.4.1 Q20+ or PacBio HiFi) or hybrid before reporting MGE context.
- mcr-9 and mcr-10 frequently report without elevated MIC. Always require phenotypic confirmation before triggering infection-control as if mcr-1.
- AMRFinderPlus `--organism Klebsiella_pneumoniae` suppresses intrinsic-gene noise (fosA, oqxAB). Cross-organism comparisons that mix `--organism`-on with `--organism`-off output are not pooled-comparable.
- hAMRonization is mandatory for any multi-lab or multi-tool surveillance report. Raw AMRFinderPlus vs ResFinder vs RGI comparison is uninterpretable because of CARD ARO vs ResFinder vs NCBI ReferenceGeneCatalog ontology differences.
- For Tn4401 promoter variants (Tn4401a vs Tn4401b vs Tn4401d -- expression differences of ~5-fold), short-read pipelines cannot resolve the variant. Long-read confirmation required for any expression-level claim.
- EUCAST and CLSI breakpoint tables update annually. Document the breakpoint year alongside the S/I/R call -- a 2023-breakpoint "S" may be a 2025-breakpoint "I".
- For heteroresistance screening, default variant-caller MAF thresholds (10% lofreq, 20% GATK) hide the clinically relevant 0.1-1% subpopulations documented by Andersson et al. 2019. Use targeted deep amplicon sequencing or lofreq with `-q 13 -a 0.01`.
- MOB-suite v2 and v3 cluster codes are non-interoperable. Document version explicitly; v3.1+ adds MGE reporting. Longitudinal surveillance crossing the v2 -> v3 boundary requires re-running.

## Related Skills

- pathogen-typing - cgMLST / Kleborate / MLST for clonal context of AMR dissemination
- transmission-inference - SNP-cluster / outbreak transmission inference for resistant clones
- variant-surveillance - Lineage-level antiviral-drug-resistance tracking
- metagenomics/amr-detection - Community ARG quantification (not isolate AMR)
- variant-calling/variant-calling - Per-isolate SNP calling for point-mutation panels
- variant-calling/filtering-best-practices - MAF threshold discipline for heteroresistance
- long-read-sequencing/long-read-alignment - Plasmid context resolution
- comparative-genomics/whole-genome-alignment - Reference-based coordinate handling
- workflows/somatic-variant-pipeline - End-to-end orchestration patterns
