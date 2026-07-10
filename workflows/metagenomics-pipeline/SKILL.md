---
name: bio-workflows-metagenomics-pipeline
description: End-to-end shotgun metagenomics workflow from FASTQ to taxonomic and functional profiles, orchestrating controls/host depletion, Kraken2+Bracken classification, MetaPhlAn marker profiling, and HUMAnN functional profiling. Covers the controls-first ordering, why Kraken2 read counts are not abundances and MetaPhlAn cell fractions do not equal Bracken read fractions, and the consistent-pipeline framing. Use when profiling shotgun metagenomic samples end to end, or chaining classification, abundance, and function. For resistome see metagenomics/amr-detection; for strains see metagenomics/strain-tracking; for assembly see genome-assembly/metagenome-assembly.
tool_type: cli
primary_tool: Kraken2
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - metagenomics/contamination-controls
  - metagenomics/kraken-classification
  - metagenomics/metaphlan-profiling
  - metagenomics/abundance-estimation
  - metagenomics/functional-profiling
  - metagenomics/metagenome-visualization
qc_checkpoints:
  - after_qc: "Q30 >80%, host reads removed"
  - after_classification: "Classification rate >60%, known taxa dominant"
  - after_functional: "Pathway coverage reasonable, unmapped <50%"
---

## Version Compatibility

Reference examples tested with: Bowtie2 2.5.3+, Bracken 2.9+, HUMAnN 3.8+, Kraken2 2.1+, MetaPhlAn 4.1+, fastp 0.23+, samtools 1.19+, matplotlib 3.8+, pandas 2.2+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metagenomics Pipeline

**"Analyze my metagenomic samples from FASTQ to taxonomic and functional profiles"** -> Orchestrate controls and host depletion, Kraken2/Bracken taxonomic classification, MetaPhlAn profiling, and HUMAnN3 functional analysis - reporting results relative to a consistent pipeline, never as a direct observation of the community.

Complete workflow from metagenomic FASTQ to taxonomic and functional profiles. This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

A taxonomic/functional profile is a position in a choice-chain (extraction -> depletion -> depth -> classifier -> DB -> normalization), never a direct observation; the trustworthy result is decided at these seams.

1. **The reference DB + version is THE inherited commitment.** The Kraken2/GTDB (or standard/RefSeq/UHGG) DB chosen at classification fixes what is detectable — a zero means below-detection OR not-in-DB OR lost-in-extraction OR removed-by-depletion, almost never biological absence. Pin the DB build (version alone moves species/genus calls); match DB to habitat (UHGG for gut); report the CLASSIFIED FRACTION (a low fraction is the tell that the DB is wrong).
2. **Host removal against T2T-CHM13 is a made-once commitment done BEFORE profiling.** Prefer the complete T2T over gapped GRCh38 (which lets human reads masquerade as novel microbes); mask rDNA; discard both mates if either maps host. It is a privacy obligation (leaked human reads are identifiable), not just QC.
3. **Controls-first is a design commitment, not a step added later.** Extraction blanks + a whole-cell mock carried through the WHOLE workflow. NO low-biomass result (skin, BAL, CSF, blood, tissue) is interpretable without blanks + DNA-concentration + decontam; at near-zero biomass the signal IS the kitome (Salter 2014). A blank cannot be retrofitted.
4. **Read-fraction is not cell-fraction, and tools/DBs are not comparable.** Kraken2 read-fraction (genome-size/copy-number biased) and MetaPhlAn cell-fraction must never be merged into one table. Holding tool+DB constant within a study is the only way a comparison measures biology and not the tool.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Reference DB + version + habitat match | What is detectable; a zero is below-detection/not-in-DB, not absence; report classified fraction |
| Host-removal reference (T2T-CHM13) | Human reads masquerading as microbes; a privacy obligation; removed before profiling |
| Controls (blanks + mock) carried through | Whether any low-biomass result is interpretable; the signal IS the kitome without them |
| Extraction method held constant | Extraction bias outweighs much biological signal (Costea 2017); interacts with read-vs-assembly choice |

## Workflow Overview

```
FASTQ files (+ extraction blanks, mock)
    |
    v
[0. QC, Host Removal & Controls] --> fastp + Hostile/Bowtie2(T2T) + blanks/decontam + Nonpareil depth check
    |
    v
[1. Taxonomic Classification]
    |
    +---> Kraken2 (+confidence, +hit-groups) + Bracken -> read fraction
    |
    +---> MetaPhlAn 4 (marker-based, pinned --index) -> cell fraction (NOT comparable to Bracken %)
    |
    v
[2. Functional Profiling] --> HUMAnN (potential, not activity; keep UNMAPPED)
    |
    v
Taxonomic profiles + Pathway abundances (+ AMR/strain via their own skills)
```

## Primary Path: Kraken2 + Bracken + HUMAnN

### Step 0: Quality Control, Host Removal, and Controls

Carry extraction blanks and a mock through the whole workflow; host-deplete against T2T-CHM13; confirm depth with Nonpareil. See metagenomics/contamination-controls for the controls/decontam detail.

```bash
# QC with fastp (trimming mechanics: read-qc/fastp-workflow)
for sample in sample1 sample2 sample3; do
    fastp -i ${sample}_R1.fastq.gz -I ${sample}_R2.fastq.gz \
        -o trimmed/${sample}_R1.fq.gz -O trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --html qc/${sample}_fastp.html
done

# Remove host reads - Hostile with a T2T-CHM13 index removes >99.5% host with low microbial loss.
# Report the reads removed; host depletion can halve usable depth.
for sample in sample1 sample2 sample3; do
    hostile clean --fastq1 trimmed/${sample}_R1.fq.gz --fastq2 trimmed/${sample}_R2.fq.gz \
        --index human-t2t-hla --aligner bowtie2 --output host_removed/
    # hostile names each paired output from its OWN input basename: fastq1 -> {R1}.clean_1.fastq.gz,
    # fastq2 -> {R2}.clean_2.fastq.gz. Rename to the ${sample}_R1/_R2.fq.gz the steps below consume.
    mv host_removed/${sample}_R1.clean_1.fastq.gz host_removed/${sample}_R1.fq.gz
    mv host_removed/${sample}_R2.clean_2.fastq.gz host_removed/${sample}_R2.fq.gz
done
# Then run decontam on the classifier output table using the blanks (contamination-controls),
# and confirm depth adequacy with Nonpareil before interpreting any non-detection.
```

### Step 1A: Kraken2 Classification

```bash
# Classify reads. Raise --confidence above the default 0 to suppress single-k-mer false positives,
# and require >=2 hit groups. The database defines what can be detected.
for sample in sample1 sample2 sample3; do
    kraken2 --db kraken2_db \
        --threads 8 \
        --paired \
        --confidence 0.1 \
        --minimum-hit-groups 2 \
        --report kraken/${sample}.report \
        --output kraken/${sample}.output \
        host_removed/${sample}_R1.fq.gz \
        host_removed/${sample}_R2.fq.gz
done
```

### Step 1B: Bracken Abundance Estimation

```bash
# Estimate species abundance
for sample in sample1 sample2 sample3; do
    bracken -d kraken2_db \
        -i kraken/${sample}.report \
        -o bracken/${sample}.species.txt \
        -r 150 \
        -l S \
        -t 10
done

# Combine samples into abundance matrix
combine_bracken_outputs.py \
    --files bracken/*.species.txt \
    -o bracken/combined_species.txt
```

### Step 1C: Alternative - MetaPhlAn Profiling

```bash
# Profile with MetaPhlAn 4. Pin --index (DB version is a batch variable). MetaPhlAn % is a cell
# fraction - do NOT merge it with Bracken read fractions. In 4.2 --bowtie2out is renamed --mapout.
for sample in sample1 sample2 sample3; do
    metaphlan host_removed/${sample}_R1.fq.gz,host_removed/${sample}_R2.fq.gz \
        --bowtie2out metaphlan/${sample}.bowtie2.bz2 \
        --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
        --input_type fastq \
        --nproc 8 \
        -o metaphlan/${sample}_profile.txt
done

# Merge profiles
merge_metaphlan_tables.py metaphlan/*_profile.txt > metaphlan/merged_abundance.txt
```

### Step 2: Functional Profiling with HUMAnN

```bash
# Run HUMAnN
for sample in sample1 sample2 sample3; do
    # Concatenate paired reads
    cat host_removed/${sample}_R1.fq.gz host_removed/${sample}_R2.fq.gz > \
        host_removed/${sample}_concat.fq.gz

    humann --input host_removed/${sample}_concat.fq.gz \
        --output humann/${sample} \
        --threads 8 \
        --metaphlan-options "--bowtie2db metaphlan_db"
done

# Normalize and join tables. HUMAnN names outputs from the input STEM, so the ${sample}_concat.fq.gz
# input above yields sample1_concat_pathabundance.tsv (not sample1_pathabundance.tsv).
humann_renorm_table --input humann/sample1/sample1_concat_pathabundance.tsv \
    --output humann/sample1/sample1_concat_pathabundance_cpm.tsv \
    --units cpm

# --search-subdirectories: per-sample outputs live in humann/<sample>/ subdirs; the join is
# non-recursive by default and would otherwise find zero files.
humann_join_tables --input humann \
    --search-subdirectories \
    --output humann/merged_pathabundance.tsv \
    --file_name pathabundance
```

### Visualization

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load Bracken species table. combine_bracken_outputs.py emits name, taxonomy_id, taxonomy_lvl,
# then per-sample {sample}_num / {sample}_frac columns. Keep ONLY the fractions: taxonomy_lvl is a
# string ('S') and summing it raises TypeError; summing taxonomy_id would be meaningless anyway.
species = pd.read_csv('bracken/combined_species.txt', sep='\t', index_col=0)
species = species.filter(regex='_frac$').rename(columns=lambda c: c.replace('_frac', ''))

# Top 20 species heatmap
top20 = species.sum(axis=1).nlargest(20).index
plt.figure(figsize=(12, 8))
sns.heatmap(species.loc[top20], cmap='viridis', annot=False)
plt.title('Top 20 Species Abundance')
plt.tight_layout()
plt.savefig('top20_species_heatmap.pdf')

# Stacked bar plot
species_norm = species.div(species.sum()) * 100
top10 = species_norm.sum(axis=1).nlargest(10).index
other = species_norm.loc[~species_norm.index.isin(top10)].sum()

plot_data = species_norm.loc[top10].T
plot_data['Other'] = other
plot_data.plot(kind='bar', stacked=True, figsize=(10, 6))
plt.ylabel('Relative Abundance (%)')
plt.legend(bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.savefig('species_barplot.pdf')
```

## Parameter Recommendations

| Step | Parameter | Value |
|------|-----------|-------|
| fastp | --length_required | 50 (metagenomic reads) |
| Kraken2 | --confidence | 0.1-0.4 (default 0.0 over-classifies; see metagenomics/kraken-classification) |
| Kraken2 | --minimum-hit-groups | 2 (cut single-region false positives) |
| Bracken | -r | Read length (e.g., 150; must match the DB build) |
| Bracken | -l | S (species) or G (genus) |
| Bracken | -t | 10 (min reads threshold) |
| MetaPhlAn | --min_cu_len | 2000 (default) |
| HUMAnN | --threads | 8+ |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| A "novel community" that is the kitome | Negative controls skipped / contamination bleed | Blanks + mock through the full workflow; decontam; skepticism toward canonical kitome genera |
| Same organism appears twice under different names | Merged GTDB-Tk-labelled and NCBI-labelled tables (Firmicutes vs Bacillota) | State which taxonomy each table uses; never name-merge without a crosswalk |
| Cross-study comparison is really tool differences | Compared across tools/DBs | Hold tool + DB constant within a study; benchmark on a mock with OPAL |
| A zero read as biological absence | Confused detection-limit with biology | Name which link (depth/DB/extraction/depletion) is responsible before any biological reading; Nonpareil depth check |
| Read-fraction and cell-fraction merged | Kraken2 % joined to MetaPhlAn % | Keep separate tables; they use different absence semantics |
| Low classification rate | Database mismatch / novel organisms | Match DB to habitat; report classified fraction; a low fraction = wrong/incomplete DB |
| High host reads | Incomplete host removal | Use the complete T2T host reference; mask rDNA |

## References

- Salter SJ, Cox MJ, Turek EM, et al (2014) Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biology* 12:87. DOI 10.1186/s12915-014-0087-z. (the kitome.)
- Costea PI, Zeller G, Sunagawa S, et al (2017) Towards standards for human fecal sample processing in metagenomic studies. *Nature Biotechnology* 35:1069-1076. DOI 10.1038/nbt.3960. (extraction dominates.)
- Meyer F, Bremges A, Belmann P, et al (2019) Assessing taxonomic metagenome profilers with OPAL. *Genome Biology* 20:51. DOI 10.1186/s13059-019-1646-y.
- Sczyrba A, Hofmann P, Belmann P, et al (2017) Critical Assessment of Metagenome Interpretation (CAMI). *Nature Methods* 14:1063-1071. DOI 10.1038/nmeth.4458.

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

THREADS=8
KRAKEN_DB="kraken2_standard_db"
HOST_INDEX="human_bt2_index"   # MUST be built from T2T-CHM13 (CHM13v2, +HLA) per the made-once host-removal commitment, not a legacy GRCh38 index
SAMPLES="sample1 sample2 sample3"
OUTDIR="metagenomics_results"

mkdir -p ${OUTDIR}/{trimmed,host_removed,kraken,bracken,metaphlan,humann,qc}

# Step 1: QC
echo "=== QC ==="
for sample in $SAMPLES; do
    fastp -i ${sample}_R1.fastq.gz -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --length_required 50 \
        --html ${OUTDIR}/qc/${sample}_fastp.html -w ${THREADS}
done

# Host removal
echo "=== Host Removal ==="
for sample in $SAMPLES; do
    # -f 12 (read unmapped AND mate unmapped) keeps only pairs where NEITHER mate hit the host.
    # --un-conc-gz would instead keep every non-CONCORDANT pair, retaining pairs whose mate mapped
    # human -- a privacy leak, not just a QC lapse. -F 256 drops secondary alignments.
    bowtie2 -p ${THREADS} -x ${HOST_INDEX} --very-sensitive \
        -1 ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -2 ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        2> ${OUTDIR}/qc/${sample}_host.log \
      | samtools view -b -f 12 -F 256 - \
      | samtools sort -n -@ ${THREADS} - \
      | samtools fastq -1 ${OUTDIR}/host_removed/${sample}_R1.fq.gz \
                       -2 ${OUTDIR}/host_removed/${sample}_R2.fq.gz \
                       -0 /dev/null -s /dev/null -n -
done

# Step 2: Kraken2
echo "=== Kraken2 ==="
for sample in $SAMPLES; do
    kraken2 --db ${KRAKEN_DB} --threads ${THREADS} --paired \
        --confidence 0.1 --minimum-hit-groups 2 \
        --report ${OUTDIR}/kraken/${sample}.report \
        --output ${OUTDIR}/kraken/${sample}.output \
        ${OUTDIR}/host_removed/${sample}_R1.fq.gz \
        ${OUTDIR}/host_removed/${sample}_R2.fq.gz
done

# Bracken
echo "=== Bracken ==="
for sample in $SAMPLES; do
    bracken -d ${KRAKEN_DB} \
        -i ${OUTDIR}/kraken/${sample}.report \
        -o ${OUTDIR}/bracken/${sample}.species.txt \
        -r 150 -l S -t 10
done

echo "=== Pipeline Complete ==="
echo "Kraken reports: ${OUTDIR}/kraken/"
echo "Bracken abundances: ${OUTDIR}/bracken/"
```

## Related Skills

- database-access/sra-data - Pull metagenomic FASTQ from SRA / ENA (16S amplicon or shotgun)
- database-access/ncbi-datasets-cli - Bulk-pull reference genomes for read mapping
- database-access/remote-homology - DIAMOND --ultra-sensitive for predicted-ORF annotation
- metagenomics/contamination-controls - Host depletion, blanks/decontam, depth checks up front
- metagenomics/kraken-classification - Kraken2 details
- metagenomics/metaphlan-profiling - MetaPhlAn parameters
- metagenomics/abundance-estimation - Bracken options and compositional handling
- metagenomics/functional-profiling - HUMAnN workflow
- metagenomics/amr-detection - Community resistome from the same reads
- metagenomics/strain-tracking - Strain resolution from the same reads
- metagenomics/metagenome-visualization - Plotting and community statistics
