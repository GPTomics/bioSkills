---
name: bio-metagenomics-amr-detection
description: Profiles the antimicrobial-resistance gene content (resistome) of shotgun metagenomes - read-based quantification with RGI bwt, AMR++/MEGARes, ARGs-OAP/SARG, deepARG, or GROOT, and presence calling with AMRFinderPlus/ABRicate on assembled contigs or MAGs. Covers why an ARG hit is a sequence match not a phenotype, why a metagenomic ARG has no host and no genomic context until assembly (and assembly breaks at ARGs), per-gene curated thresholds vs a flat 80/80, gene-fraction false-positive control, and cross-study normalization pitfalls. Use when quantifying a community resistome, normalizing ARG abundance, or calling ARGs from metagenome contigs. For pure-culture isolate AMR, point mutations, and phenotype/MIC prediction see epidemiological-genomics/amr-surveillance.
tool_type: cli
primary_tool: AMRFinderPlus
---

## Version Compatibility

Reference examples tested with: AMRFinderPlus 3.12+, RGI 6+ (CARD 3.2+), ABRicate 1.0+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `amrfinder -V` (reports software AND database version), `rgi main --version`, `abricate --list` to confirm DB snapshots
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The AMR reference DATABASE is versioned and updated roughly monthly; `amrfinder -V` reports both software and database version, and ABRicate ships pinned database snapshots so two labs on different versions get different calls. Record the tool version, the database version, and (for read-based work) the normalization unit and sequencing depth - none are recoverable later.

# AMR Detection (Community Resistome)

**"What resistance genes are in my community, and how abundant?"** -> Match reads or contigs to a curated ARG database - reporting that ARG sequences are present at some relative abundance, never that the sample is resistant, because a metagenomic hit has no host and no expression.
- CLI (reads): `rgi bwt -1 R1.fq.gz -2 R2.fq.gz -a kma -n 16 -o sample --local`
- CLI (contigs/MAGs): `amrfinder -n contigs.fasta --plus -o amr.tsv`

Scope: community/metagenomic resistome - read-based quantification and contig/MAG presence calling. Pure-culture isolate AMR, point-mutation resistance, in-silico antibiogram, MLST/clone/outbreak context, and GLASS reporting -> epidemiological-genomics/amr-surveillance. Assembly/binning and ARG-host linkage mechanics -> genome-assembly/metagenome-assembly. General gene-family/pathway abundance -> functional-profiling.

## The Single Most Important Modern Insight -- An ARG Hit Is a Sequence Match, Not a Phenotype

An ARG hit is a match against a reference database, not a measured resistance phenotype - and in a metagenome it is a match with no host and no genomic context until assembly. Three inferences a naive pipeline silently makes, all wrong:

1. **match -> gene** - defeated by partial hits (a 30 bp conserved-domain hit to a 1 kb ARG is not a gene); guarded by gene-fraction / breadth-of-coverage.
2. **gene -> resistance** - defeated by expression and regulation: silent sul2 below ECOFF, ampC/blaOXA driven only when an IS lands in the promoter, efflux that needs overexpression, truncations that still score a partial hit, point mutations where only the SNP matters.
3. **gene -> who carries it / is it mobile** - defeated by read shortness: a short read cannot see its neighbors, so host and plasmid context need assembly, long reads, or Hi-C.

The honest deliverable is "these ARG sequences are present at this relative abundance in this community," never "this sample is resistant to drug X." The moment a report says resistant, it has smuggled in a host, an expression assumption, and a clinical breakpoint the data never contained. On a pure culture the organism can be grown and an MIC measured - that is a different skill (epidemiological-genomics/amr-surveillance).

## Read-Based vs Assembly-Based: the Core Metagenomic Tradeoff

| Axis | Read-based (RGI bwt, AMR++, ARGs-OAP, deepARG, GROOT) | Assembly-based (AMRFinderPlus/RGI main/ABRicate on contigs/MAGs) |
|------|---------------------------------|----------------------------|
| Low-abundance sensitivity | high - every read counts, below assembly coverage | low - ARGs at low coverage do not assemble |
| Quantification | yes - abundance + normalization | presence/absence per contig |
| Host / MGE context | none without binning | possible via contig taxonomy / MAG |
| Point-mutation resistance | weak/unreliable (RGI bwt cannot screen the SNP) | yes, with organism/model |
| False positives | partial hits unless gene-fraction filtered | chimeric contigs, but vettable |

The assembly paradox: metagenomic assemblies preferentially break exactly at ARG/MGE boundaries, recovering only a small fraction of true ARG genomic contexts and underestimating the resistome (Abramova 2024 *BMC Genomics* 25:959). So "assemble to get host" is necessary but not sufficient - long reads (Nanopore/PacBio) and Hi-C metagenomics are the real remedy for ARG-host/MGE linkage.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| AMRFinderPlus | Feldgarden 2021 *Sci Rep* 11:12728 | NCBI Reference Gene Catalog; per-gene curated cutoffs + HMMs | contig/MAG presence calling; the default contig caller |
| RGI bwt | Alcock 2023 *Nucleic Acids Res* 51:D690 | CARD homolog-model read mapping (KMA/bowtie2/bwa) | read-based resistome with coverage/depth per allele |
| AMR++ / MEGARes 3.0 | Bonin & Doster 2023 *Nucleic Acids Res* 51:D744 | BWA-MEM + gene-fraction filter + rarefaction | quantitative resistome with built-in partial-hit control |
| ARGs-OAP / SARG | Yin 2023 *Engineering* 27:234 | two-stage read annotation + 16S/cell normalization | copies-ARG-per-16S / per-cell units |
| deepARG | Arango-Argoty 2018 *Microbiome* 6:23 | deep-NN over dissimilarity features | catches divergent ARGs best-hit BLAST misses |
| GROOT | Rowe & Winn 2018 *Bioinformatics* 34:3601 | variation-graph alignment | types SNP-bearing alleles that flat references conflate |
| ABRicate | Seemann (no paper) | flat 80/80 BLASTn, bundled DB snapshots | quick contig screen; acquired genes only, no point mutations |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Quantitative resistome from reads | RGI bwt or AMR++ or ARGs-OAP | abundance + normalization; no host/context |
| Divergent / novel ARGs from reads | deepARG (confirm surprising calls) | dissimilarity features beat top-hit BLAST |
| Type a specific high-similarity allele | GROOT | graph carries SNP-bearing variants |
| Presence per contig / MAG | AMRFinderPlus (`--plus`) on contigs | curated per-gene cutoffs; possible host via binning |
| Quick multi-DB contig screen | ABRicate | fast; but flat 80/80, no point mutations |
| Is the ARG mobile / in a pathogen? | assemble+bin, long read, or Hi-C | short reads cannot link ARG to host |
| Pure culture / phenotype / MIC | -> epidemiological-genomics/amr-surveillance | isolate AMR is a different skill |
| Cross-study abundance comparison | within-study only, same DB+normalization+depth | "total ARG abundance" is rarely comparable |

## Read-Based Resistome Quantification

```bash
# CARD read mapping (homolog models). RGI bwt CANNOT screen point-mutation SNPs, so this is for
# acquired/homolog ARGs only - never report a gyrA read hit as fluoroquinolone resistance.
rgi bwt -1 reads_R1.fq.gz -2 reads_R2.fq.gz \
    -a kma -n 16 \
    -o sample_resistome --local
# Outputs *.gene_mapping_data.txt with percent coverage and depth per gene.

# AMR++/MEGARes applies the gene-fraction filter (default 80%): the minimum proportion of a
# reference covered by >=1 read for "present" - the read-based analog of breadth-of-coverage.
```

ARGs-OAP/SARG normalizes to copies-of-ARG-per-16S or per-cell; report the unit. Gene fraction (breadth) is the single most important false-positive guard - without it a conserved-domain fragment counts as a present gene.

## Contig / MAG Presence Calling

```bash
amrfinder -n contigs.fasta \
    --plus \                  # also report biocide/metal (STRESS) and virulence elements
    --threads 8 -o amr.tsv
# --ident_min default -1 = use the per-gene CURATED cutoffs; overriding with a global value is usually a mistake.
# Point mutations require --organism (a single known species) - inappropriate for a mixed community;
# use it only on a taxonomically resolved MAG, and defer isolate point-mutation work to amr-surveillance.
```

AMRFinderPlus uses manually curated per-gene BLAST cutoffs (plus HMM cutoffs with protein), not a flat 80/80 - catching divergent real variants while rejecting partial housekeeping homologs. ABRicate, by contrast, is flat 80/80 and acquired-genes-only; it will never report a point mutation.

## Per-Method Failure Modes

### ARG presence reported as resistance
**Trigger:** an output column or summary that says "resistant." **Mechanism:** presence is not expression and not a host-linked MIC. **Symptom:** a sewage metagenome described as "resistant to carbapenems." **Fix:** report "ARG detected at abundance X"; reserve phenotype claims for isolates (amr-surveillance).

### `--organism` on a mixed community
**Trigger:** `amrfinder --organism Escherichia` on community contigs. **Mechanism:** organism mode assumes a single known species and calls organism-specific point mutations/intrinsic genes. **Symptom:** spurious point-mutation calls; filtered "intrinsic" genes wrong for the community. **Fix:** run organism mode only on a taxonomically resolved MAG; otherwise omit it.

### Partial hit counted as a gene
**Trigger:** read mapping or BLAST with no breadth filter. **Mechanism:** a short conserved-domain match to a long ARG passes an identity threshold. **Symptom:** inflated ARG counts dominated by fragments. **Fix:** require gene-fraction / breadth-of-coverage (AMR++ default 80%); inspect coverage, not just identity.

### Cross-study abundance comparison
**Trigger:** comparing "total ARG abundance" across papers. **Mechanism:** different databases (CARD/MEGARes/SARG/ResFinder), normalization units, aligners, and depth all change the number. **Symptom:** apparent resistome differences that are pipeline artifacts. **Fix:** compare only within a study with one pipeline; report DB version, tool version, normalization unit, and depth; hAMRonization harmonizes format, not the metric.

### Loose/Discovery hits reported as ARGs
**Trigger:** CARD-RGI `--include_loose` in a surveillance report. **Mechanism:** Loose is below the curated bit-score cutoff - discovery only. **Symptom:** a flood of low-similarity false positives. **Fix:** report Perfect/Strict; reserve Loose for novel-variant discovery with manual curation.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| AMRFinderPlus `--ident_min` -1 (use curated) | Feldgarden 2021 *Sci Rep* 11:12728 | per-gene curated cutoffs beat a global 80/80; overriding is usually wrong |
| AMRFinderPlus `--coverage_min` 0.5 | AMRFinderPlus docs | minimum reference coverage for a call |
| AMR++ gene fraction 80% | Bonin & Doster 2023 *Nucleic Acids Res* 51:D744 | breadth guard against partial-hit false positives |
| ResFinder acquired 0.80 id / 0.60 cov | Bortolaia 2020 *J Antimicrob Chemother* 75:3491 | the documented default (not 0.90) |
| deepARG `--min-prob` 0.8 | Arango-Argoty 2018 *Microbiome* 6:23 | category probability cutoff; confirm surprising calls |
| CARD Loose tier = discovery only | Alcock 2023 *Nucleic Acids Res* 51:D690 | below curated bit-score; not for surveillance |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No point mutations reported | ABRicate or read-based homolog mapper used | those cannot call SNPs; use AMRFinderPlus `--organism` on a MAG / defer to amr-surveillance |
| AMRFinderPlus calls changed silently | stale reference database | `amrfinder -u`; record `amrfinder -V` software + DB version |
| ABRicate results differ between labs | pinned DB snapshot version differs | record `abricate --list` versions; update with abricate-get_db |
| Inflated ARG abundance | no gene-fraction/breadth filter | apply breadth-of-coverage; inspect coverage |
| gyrA "hit" from read-based RGI | `--include_other_models` reports the gene, not the SNP | do not call resistance; SNP screening needs an isolate/organism |

## References

- Feldgarden M, Brover V, Gonzalez-Escalona N, et al. 2021. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep* 11:12728.
- Alcock BP, Huynh W, Chalil R, et al. 2023. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res* 51:D690-D699.
- Bortolaia V, Kaas RS, Ruppe E, et al. 2020. ResFinder 4.0 for predictions of phenotypes from genotypes. *J Antimicrob Chemother* 75:3491-3500.
- Bonin N, Doster E, Worley H, et al. 2023. MEGARes and AMR++, v3.0: an updated comprehensive database of antimicrobial resistance determinants and an improved software pipeline. *Nucleic Acids Res* 51:D744-D752.
- Yin X, Zheng X, Li L, et al. 2023. ARGs-OAP v3.0: antibiotic-resistance gene database curation and analysis pipeline optimization. *Engineering* 27:234-241.
- Arango-Argoty G, Garner E, Pruden A, et al. 2018. DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data. *Microbiome* 6:23.
- Rowe WPM, Winn MD. 2018. Indexed variation graphs for efficient and accurate resistome profiling. *Bioinformatics* 34:3601-3608.
- Abramova A, Karkman A, Bengtsson-Palme J. 2024. Metagenomic assemblies tend to break around antibiotic resistance genes. *BMC Genomics* 25:959.

## Related Skills

- epidemiological-genomics/amr-surveillance - Isolate AMR, point mutations, phenotype/MIC, typing, GLASS
- functional-profiling - General gene-family/pathway abundance (HUMAnN can surface ARG families)
- kraken-classification - Taxonomic context for the community
- genome-assembly/metagenome-assembly - Assembly/binning and ARG-host linkage
- contamination-controls - Host depletion before resistome profiling
- workflows/metagenomics-pipeline - End-to-end shotgun analysis
