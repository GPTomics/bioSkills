# Pathogen Typing Usage Guide

## Overview

This skill assigns isolate identity at the resolution that matches the surveillance question. It covers species triage (Mash / skani / fastANI), historical 7-locus MLST, cgMLST / wgMLST for outbreak investigation (chewBBACA + EnteroBase HierCC for stable nomenclature; Ridom SeqSphere and BIGSdb as alternatives), in-silico serotyping per organism (SISTR / SeqSero2 for *Salmonella*; SerotypeFinder for *E. coli*; Kaptive K/O via Kleborate for *Klebsiella*; SeroBA for *Streptococcus pneumoniae*; spa + SCCmec for *Staphylococcus aureus*), MTBC lineage via the Coll 2014 + Napier 2020 90-SNP barcode (TB-Profiler / Mykrobe), and SARS-CoV-2 Pango lineage + Nextclade clade with explicit pangolin-data and dataset version pinning. The skill is decision-grade: it routes the agent to cgMLST vs core-SNP based on the cluster question, refuses to apply Walker 2013 UK TB thresholds in high-transmission settings without caveat, flags chewBBACA missing-locus codes that must be excluded from distance counts, and insists on Pango / Nextclade version pinning before any cross-time-point comparison.

Distinct from `variant-surveillance` (which owns ongoing lineage-frequency tracking) and from `transmission-inference` (which owns the who-infected-whom inference downstream of cluster definition).

## Prerequisites

```bash
conda install -c bioconda mlst chewbbaca sistr_cmd seqsero2 serotypefinder kleborate kaptive seroba poppunk pangolin nextclade tb-profiler mykrobe mash skani snippy snp-dists gubbins iqtree

amrfinder -u
tb-profiler update_tbdb

pip install pandas biopython
```

## Quick Start

Tell the AI agent what typing question to answer:

- "Run 7-locus MLST and cgMLST on these 40 Salmonella assemblies and define outbreak clusters at the EFSA <=5 allele threshold"
- "Serotype these Salmonella isolates with SISTR and flag any monophasic Typhimurium 1,4,[5],12:i:- variants"
- "Run Kleborate on these Klebsiella pneumoniae assemblies and tell me which are hypervirulent vs classical"
- "Type these MRSA isolates with spa, SCCmec, and MLST CC; flag CC8 USA300, CC22 EMRSA-15"
- "Assign Pango lineages and Nextclade clades for these 1200 SARS-CoV-2 consensus sequences using UShER mode; pin pangolin-data version"
- "Call MTBC lineage with TB-Profiler against the Napier 2020 90-SNP barcode and report sublineage"
- "Compute core-SNP distances after Gubbins recombination masking and apply the Walker 2013 12-SNP TB threshold (caveat: UK setting)"

## Example Prompts

### Salmonella surveillance cohort with cgMLST + serotyping

> "Three hundred Illumina paired-end runs of *Salmonella enterica* from a one-year national surveillance project. Assemble with SKESA, run SISTR for serovar prediction (flag monophasic 1,4,[5],12:i:- variants), run 7-locus MLST via Seemann's mlst, run chewBBACA AlleleCall against the EnteroBase SalmonellaEnterica cgMLST schema fetched 2025-09, extract cgMLST with 0.95 completeness threshold, compute pairwise-complete allele distances, and define outbreak clusters at the EFSA harmonised <=5 allele threshold. Report serovar -- ST -- cgMLST HC5 / HC10 / HC50 (EnteroBase HierCC nomenclature) -- cluster membership per isolate."

### Klebsiella hypervirulence vs classical surveillance

> "Eighty *Klebsiella pneumoniae* clinical isolates from a regional surveillance project. Run Kleborate to get integrated MLST + K/O typing (via Kaptive) + ICEKp / yersiniabactin / colibactin / aerobactin virulence loci + AMR. Flag hypervirulent isolates explicitly. Document Kleborate version and Kaptive K-locus DB version. Cross-reference K-locus call against MOB-suite plasmid context for any KL with elevated virulence score."

### Mtb outbreak investigation with the Napier barcode

> "Forty *Mycobacterium tuberculosis* isolates from a hospital cluster investigation. Run TB-Profiler to call lineage (verify Napier 2020 90-SNP barcode is loaded, not the older Coll 2014 62-SNP version) and drug resistance against the WHO 2nd-edition catalogue. Run Mykrobe as a cross-check. Compute core-SNP distances after Snippy + Gubbins, define clusters at the Walker 2013 <=12 SNP threshold (note: this is a UK low-transmission threshold; flag the caveat). Distinguish Beijing sublineages explicitly (modern Beijing 2.2.1.1 vs ancestral 2.2.1.2)."

### SARS-CoV-2 lineage assignment with version pinning

> "Twelve hundred SARS-CoV-2 consensus FASTAs spanning January through July 2025. Run Pangolin with `--analysis-mode usher` and Nextclade with the current sars-cov-2 dataset. Record `pangolin --all-versions` and `nextclade dataset list --tag latest` for every batch. Reconcile any UShER vs pangoLEARN disagreement (UShER takes precedence). For any sample where Pango and Nextclade disagree on recombinant assignment (X-prefix vs parent lineage), flag for manual review. Pin pangolin-data version before producing the final report."

### S. pneumoniae serotype + GPSC for vaccine surveillance

> "Two hundred *Streptococcus pneumoniae* invasive disease isolates from a post-PCV20 surveillance period. Run SeroBA for serotyping from raw reads (>=15x coverage). Run PopPUNK in query mode against the up-to-date GPS reference DB to assign GPSC. Distinguish PCV13 / PCV15 / PCV20 / PCV21 serotype coverage explicitly. For vaccine-replacement analysis, compute serotype prevalence by year and flag emerging non-vaccine GPSCs."

### Multi-country outbreak with harmonised schema

> "We have *Listeria monocytogenes* assemblies from four countries; the source labs each used their own cgMLST schema. Re-run chewBBACA against the BIGSdb-Lm 1748-locus schema for all assemblies; compute pairwise-complete allele distances; define outbreak clusters at the PulseNet <=4 allele convention; reconcile against the originating-lab cluster definitions and explain any discordance."

## What the Agent Will Do

1. Identify the typing question and select resolution (species triage / MLST / cgMLST / SNP / lineage caller / serotype) to match.
2. For bacteria, pin the schema source AND fetch date (chewBBACA / Ridom / EnteroBase); for SARS-CoV-2, record `pangolin --all-versions` and Nextclade dataset tag; for Mtb, verify TB-Profiler bundled barcode is Napier 2020 or later.
3. Run the appropriate caller(s); for *Klebsiella* always use Kleborate (integrated MLST + Kaptive + virulence + AMR); for *Salmonella* use SISTR + SeqSero2; for *S. aureus* use spa + SCCmec + MLST + CC.
4. For outbreak cluster definition, compute pairwise-complete cgMLST allele distance OR core-SNP distance (recombination-masked via Gubbins for bacteria); apply the published pathogen-specific threshold and cite its source population.
5. For lineage callers, never collapse Group 3 / Uncertain calls to a default; report explicitly.
6. Cross-check at decision points: Pangolin UShER vs Nextclade clade; TB-Profiler vs Mykrobe for any R/XDR call; SISTR vs slide agglutination for any monophasic-variant ambiguity.
7. Document missing-locus handling explicitly for cgMLST distance calculation -- pairwise-complete intersection, never count missing-vs-anything.
8. For multi-country comparisons, agree on a single schema (typically EnteroBase HierCC for *Salmonella* / *E. coli* / *Listeria*) before running.

## Tips

- Apply Walker 2013 UK TB SNP thresholds outside a UK / low-transmission context with explicit caveat. In high-burden settings (e.g., Cape Town, Mumbai), the same threshold inflates apparent recent-transmission rates 2-5x.
- chewBBACA LNF / PLOT / NIPH / ASM / ALM are MISSING DATA, not informative differences. Compute pairwise allele distances on the intersection of called loci.
- pangolin-data is updated weekly. Cross-lab lineage comparisons over time are uninterpretable without a pinned pangolin-data version.
- pangoLEARN was deprecated mid-2023. Pangolin `--analysis-mode usher` is the modern default.
- For pre-Napier 2020 TB-Profiler bundled barcodes, lineages 7-9 (Ethiopia, Rwanda, East Africa) get miscalled as "unknown" or partial L4 sublineages. Update before running.
- "Beijing" in pre-2014 Mtb literature is a paraphyletic spoligotype-defined grouping. Modern WGS surveillance must specify sublineage (modern 2.2.1.1 vs ancestral 2.2.1.2 vs proto 2.1).
- Mash distances are non-metric below ANI ~80%. For cross-genus or distant comparisons, switch to skani or pyANI.
- Cross-schema cgMLST distances are NOT comparable. Ridom SeqSphere 5 differs from chewBBACA 5 differs from EnteroBase 5.
- SISTR's monophasic Typhimurium 1,4,[5],12:i:- call differs from slide agglutination because the fljB gene is deleted; trust the genome.
- Kleborate has had multiple K-locus DB version flips between Kaptive v1 / v2 / v3. Document both versions for any longitudinal K-locus prevalence claim.
- For SARS-CoV-2 recombinants (X-prefix), Pango designation lags emergence; tools assign to a parent lineage until manual designation. Flag for review.
- Snippy + Gubbins: `run_gubbins.py core.full.aln` (full alignment with reference), NOT `core.aln` (variable-only). The latter cannot estimate background SNP density.
- For *S. pneumoniae*, serotype IS the vaccine-actionable surveillance unit. Reporting MLST ST without serotype misses the vaccine-replacement signal.

## Related Skills

- amr-surveillance - Strain context for AMR dissemination; Kleborate integrates typing + AMR for Klebsiella
- transmission-inference - SNP-cluster definition feeds outbreak transmission inference
- phylodynamics - Time-scaled tree from typed isolates for R_e estimation
- variant-surveillance - SARS-CoV-2 longitudinal lineage frequency tracking
- comparative-genomics/pangenome-analysis - Core / accessory genome partitioning
- comparative-genomics/whole-genome-alignment - Core-genome alignment for SNP-typing
- variant-calling/variant-calling - Per-isolate variant calling
- read-alignment/bwa-alignment - Read mapping upstream of variant calling
- database-access/entrez-fetch - Reference genome retrieval for Snippy
- metagenomics/strain-tracking - Community strain tracking (NOT isolate-focused)
