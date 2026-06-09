---
name: bio-metagenomics-strain-tracking
description: Resolves and compares bacterial strains below the species level from shotgun metagenomes with inStrain (popANI/conANI microdiversity), StrainPhlAn (marker-SNV consensus phylogeny and nGD), MIDAS2, metaSNV, and StrainGE, plus genome-vs-genome ANI (skani/fastANI/MASH) for isolate/MAG comparison. Covers why a strain is a threshold not a thing, why ANI answers same-genome while popANI/nGD answer same-population-in-situ, the 99.999% popANI and per-species nGD definitions, the coverage detection limit (absence is not absence), why sharing is not transmission direction, and mapping to the dataset's own dRep MAGs. Use when detecting shared strains, tracking transmission, resolving within-host strain dynamics, or deconvoluting co-occurring strains. For pure-culture isolate outbreak SNP trees see epidemiological-genomics; for MAG assembly see genome-assembly/metagenome-assembly.
tool_type: mixed
primary_tool: inStrain
---

## Version Compatibility

Reference examples tested with: inStrain 1.8+, StrainPhlAn/MetaPhlAn 4.1+, dRep 3.4+, skani 0.2+, Bowtie2 2.5+, samtools 1.19+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `inStrain profile -h`, `strainphlan -h`, `skani dist -h` to confirm flags and defaults
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The mapping REFERENCE defines the answer. Map to genomes present in the sample set (dRep-dereplicated MAGs from this dataset), not generic database genomes - a distant reference inflates apparent SNVs and corrupts popANI. Record the reference set, the popANI/nGD threshold, the minimum coverage and breadth, and the co-detection rate; the "strain" is defined by these, not by nature.

# Strain Tracking

**"Is the same strain in samples A and B?"** -> Compare per-position SNV populations (not consensus genomes) at adequate depth - because a strain is defined by the chosen threshold, and ANI cannot resolve the difference that matters.
- CLI: `inStrain profile sample.bam reps.fasta -o sample.IS -s reps.stb -p 8` then `inStrain compare`

Scope: in-situ strain resolution, sharing, and deconvolution from a community. Pure-culture isolate outbreak SNP/cgMLST trees -> epidemiological-genomics. MAG assembly/binning -> genome-assembly/metagenome-assembly. Species presence/abundance -> kraken-classification, metaphlan-profiling.

## The Single Most Important Modern Insight -- A Strain Is a Threshold, Not a Thing

There is no universal definition of a metagenomic strain. A strain is an operational construct fixed by the reference mapped to, the genome fraction comparable at adequate depth, and the cutoff drawn. inStrain's popANI >= 99.999% over >= 50% of the genome IS the strain definition; Valles-Colomer's per-species nGD threshold IS the strain definition. Changing the cutoff changes how many strains exist. Two corollaries:

1. **ANI answers "same genome?"; popANI/nGD answer "same population, in situ?"** Genome-to-genome ANI (MASH/skani/fastANI) saturates - two genuinely distinct, separately transmissible strains routinely share >99.9% ANI, and ANI's resolution floor sits above the difference that matters. Strain sharing cannot be done with an ANI number; it needs microdiversity-aware, position-level concordance.
2. **Detecting a shared strain is a narrow statement:** across the genome fraction both samples covered at >= 5x, their SNV populations were >= 99.999% concordant. That is not proof an organism was transmitted between two people - the threshold certifies genomic identity within a coverage window, nothing more.

## The Three Tasks (Do Not Conflate)

| Task | Question | Tools |
|------|----------|-------|
| Identification | which known reference strain is present? | StrainGST, sourmash gather, MIDAS |
| Tracking / sharing | is the SAME strain in samples A and B? | inStrain compare, StrainPhlAn, MIDAS2, metaSNV, SameStr |
| Deconvolution | how many strains coexist in ONE sample, what are their haplotypes? | DESMAN, Strainberry, strainFlye, Strainy |

Genome-to-genome ANI (MASH/skani/fastANI) is a fourth, orthogonal task - "are these two assembled genomes the same?" - isolate/MAG comparison and dereplication, NOT in-situ strain resolution.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| inStrain | Olm 2021 *Nat Biotechnol* 39:727 | popANI/conANI microdiversity from reads mapped to MAGs | the reference standard for shared-strain detection |
| StrainPhlAn | Truong 2017 *Genome Res* 27:626 | marker-SNV consensus -> phylogeny -> nGD | large cross-sample marker surveys, no assembly needed |
| MIDAS2 | Zhao 2023 *Bioinformatics* 39:btac713 | UHGG pan-genome SNV + gene CNV | accessory-genome strain signal at scale |
| StrainGE | van Dijk 2022 *Genome Biol* 23:74 | k-mer search + low-coverage variant calling | low-abundance strains down to 0.5x coverage |
| metaSNV v2 | Van Rossum 2022 *Bioinformatics* 38:1162 | SNV distances + subspecies clustering | subspecies structure across samples |
| skani | Shaw 2023 *Nat Methods* 20:1661 | sparse-chaining ANI | genome-vs-genome ANI; robust on fragmented MAGs (prefer over fastANI) |
| Strainberry / strainFlye | Vicedomini 2021 *Nat Commun* 12:4485; Fedarko 2022 *Genome Res* 32:2119 | long-read haplotype separation | deconvolute co-occurring strains (-> genome-assembly) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Is a strain shared between two metagenomes? | inStrain compare (popANI) | microdiversity-aware; the field standard |
| Cross-sample transmission survey, many samples | StrainPhlAn (per-species nGD) | marker-based, scalable, no assembly |
| Low-abundance pathogen (< 1% / < 5x) | StrainGE | detects/compares down to 0.5x |
| Accessory-genome / pan-genome strain signal | MIDAS2 | adds gene-content axis SNV tools miss |
| Separate co-occurring strains into haplotypes | DESMAN (many samples) or long-read Strainberry/strainFlye | SNV tools do not partition a mixture |
| Compare two assembled genomes / dereplicate | skani (or fastANI) | genome-vs-genome ANI, not in-situ strains |
| Pure-culture isolate outbreak tree | -> epidemiological-genomics | cgMLST/SNP-distance on one genome per sample |

## inStrain: popANI vs conANI

**Goal:** Decide whether two metagenomes share a strain without being fooled by which allele happens to be the majority.

**Approach:** dRep the dataset's MAGs into representative genomes, map reads to the concatenated references, profile each sample, then compare on popANI (microdiversity-aware) over the co-covered genome fraction.

```bash
# 1. dRep -> representative genomes (97-99% ANI); concatenate; build scaffold-to-bin (.stb).
# 2. Map reads to the concatenated reps - your OWN MAGs, not database genomes.
bowtie2 -x reps -1 r1.fq.gz -2 r2.fq.gz | samtools sort -o sampleA.bam
inStrain profile sampleA.bam reps.fasta -o sampleA.IS -s reps.stb -g genes.fna -p 8
inStrain profile sampleB.bam reps.fasta -o sampleB.IS -s reps.stb -g genes.fna -p 8
inStrain compare -i sampleA.IS sampleB.IS -o compare.out -s reps.stb -p 8
```

conANI calls a difference whenever the consensus base differs - confounded by within-sample microdiversity (a minor-allele flip fakes a difference). popANI calls a difference only if the two samples share NO alleles at all, including minor ones, so popANI >= conANI always and is what detects shared strains consensus tools miss. Read genome-level calls from `genomeWide_compare.tsv` (breadth column `percent_compared`); the per-scaffold `comparisonsTable.tsv` uses `percent_genome_compared`.

## StrainPhlAn: Marker SNVs and nGD

```bash
metaphlan sample.fq.gz --input_type fastq -s sample.sam.bz2 --bowtie2out sample.bz2 -o profile.tsv  # need the SAM (-s)
sample2markers.py -i sams/*.sam.bz2 -o consensus_markers -n 8
extract_markers.py -c t__SGB1877 -o clade_markers/
strainphlan -s consensus_markers/*.json -m clade_markers/t__SGB1877.fna \
    -r reference_genomes/*.fna.bz2 -o output -c t__SGB1877 \
    --marker_in_n_samples_perc 80 --sample_with_n_markers 20 --nproc 8  # 4.0 named this --marker_in_n_samples
```

The output tree gives a pairwise nGD (normalized genetic distance). There is no universal nGD strain cutoff - derive a per-species threshold from the data (same-individual-different-timepoint pairs fall below it, unrelated pairs above), as in Valles-Colomer 2023. Low coverage means too few markers pass the filters and the sample is dropped from the species tree silently - so a missing shared-strain call is not evidence of no shared strain.

## Genome-vs-Genome ANI (Isolate/MAG Comparison, NOT In-Situ Strains)

```bash
skani dist genomeA.fasta genomeB.fasta   # prefer skani over fastANI: robust on fragmented MAGs
```

~95% ANI is the species boundary (Jain 2018 *Nat Commun* 9:5114). ANI saturates above that and cannot resolve same-vs-different strain - use it to compare isolates/MAGs and to dereplicate, never to call transmission.

## Per-Method Failure Modes

### ANI distance reported as strain resolution
**Trigger:** "MASH distance < 0.001 = same strain" or "fastANI > 99% = same strain." **Mechanism:** ANI operates on consensus genomes, saturates above 99.9%, and ignores microdiversity. **Symptom:** distinct transmissible strains called identical; transmission inferred from an ANI number. **Fix:** use ANI for isolate/MAG comparison; use inStrain popANI / StrainPhlAn nGD for strain sharing.

### Coverage detection limit (absence is not absence)
**Trigger:** concluding "no transmission" or "strain turnover." **Mechanism:** a shared strain can only be called for a species detected at adequate depth in BOTH samples (inStrain >= 5x and >= 50% breadth; StrainPhlAn enough markers). **Symptom:** a coverage dropout misread as biological absence; sharing rates biased to abundant taxa. **Fix:** report co-detection rates alongside sharing rates; use StrainGE for low-abundance targets.

### Sharing read as transmission direction
**Trigger:** narrating "A infected B." **Mechanism:** a shared strain is an undirected edge. **Symptom:** directionality claimed from one cross-sectional comparison. **Fix:** direction comes from timepoints, contact metadata, or a known index case - the published landscapes infer it from study design, not the genomic comparison.

### Wrong reference genome
**Trigger:** mapping to a generic database genome. **Mechanism:** a distant reference inflates apparent SNVs. **Symptom:** corrupted popANI; spurious differences. **Fix:** map to dRep-dereplicated MAGs from the sample set.

### Asking a SNV tool to deconvolute a mixture
**Trigger:** "what are the two strains here?" from inStrain. **Mechanism:** SNV/marker tools characterize population diversity; they do not partition it into haplotypes. **Symptom:** a category error. **Fix:** use DESMAN (many samples) or long-read Strainberry/strainFlye/Strainy for haplotype separation.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| popANI >= 99.999% same strain | Olm 2021 *Nat Biotechnol* 39:727 | empirical shared-strain cutoff; IS the operational definition |
| percent_compared >= 50% (genome-level breadth) | Olm 2021 *Nat Biotechnol* 39:727 | a genome below 50% breadth is not confidently present |
| min_cov 5x | Olm 2021 *Nat Biotechnol* 39:727 | lowest coverage at which sub-50% minor alleles are reliable |
| StrainGE detection ~0.5x | van Dijk 2022 *Genome Biol* 23:74 | tracks low-abundance strains below the inStrain floor |
| Per-species nGD threshold (derive it) | Valles-Colomer 2023 *Nature* 614:125 | no universal cutoff; separate within-host timepoints from unrelated |
| ~95% ANI species boundary | Jain 2018 *Nat Commun* 9:5114 | ANI saturates above this; cannot resolve strains |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Everything looks like one strain | ANI/MASH used for strain calls | switch to inStrain popANI / StrainPhlAn nGD |
| Sample missing from the StrainPhlAn tree | too few markers passed filters at low coverage | report co-detection; do not read absence as no-sharing |
| popANI implausibly low across the board | mapped to a distant database reference | map to dRep MAGs from the dataset |
| inStrain compare gives no genomes | < 50% breadth or < 5x in one sample | deepen sequencing or use StrainGE for that taxon |
| "Who infected whom" asked of one timepoint | sharing is undirected | need longitudinal/epi design for direction |

## References

- Olm MR, Crits-Christoph A, Bouma-Gregson K, et al. 2021. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. *Nat Biotechnol* 39:727-736.
- Truong DT, Tett A, Pasolli E, Huttenhower C, Segata N. 2017. Microbial strain-level population structure and genetic diversity from metagenomes. *Genome Res* 27:626-638.
- Zhao C, Dimitrov B, Goldman M, Nayfach S, Pollard KS. 2023. MIDAS2: Metagenomic Intra-species Diversity Analysis System. *Bioinformatics* 39:btac713.
- Van Rossum T, Costea PI, Paoli L, et al. 2022. metaSNV v2: detection of SNVs and subspecies in prokaryotic metagenomes. *Bioinformatics* 38:1162-1164.
- van Dijk LR, Walker BJ, Straub TJ, et al. 2022. StrainGE: a toolkit to track and characterize low-abundance strains in complex microbial communities. *Genome Biol* 23:74.
- Vicedomini R, Quince C, Darling AE, Chikhi R. 2021. Strainberry: automated strain separation in low-complexity metagenomes using long reads. *Nat Commun* 12:4485.
- Valles-Colomer M, Blanco-Miguez A, Manghi P, et al. 2023. The person-to-person transmission landscape of the gut and oral microbiomes. *Nature* 614:125-135.
- Shaw J, Yu YW. 2023. Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nat Methods* 20:1661-1665.
- Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. 2018. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat Commun* 9:5114.

## Related Skills

- metaphlan-profiling - StrainPhlAn builds on MetaPhlAn markers; profile species first
- kraken-classification - Species presence before strain resolution
- genome-assembly/metagenome-assembly - dRep MAGs to map against; long-read deconvolution
- epidemiological-genomics/amr-surveillance - Isolate outbreak SNP/cgMLST trees from pure cultures
- workflows/metagenomics-pipeline - End-to-end shotgun analysis
