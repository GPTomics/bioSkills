# Kraken Classification - Usage Guide

## Overview
Kraken2 assigns shotgun reads to taxa by matching their k-mers to a reference database via lowest-common-ancestor, using minimizers and a compact probabilistic hash (not exact k-mers). It is fast and sensitive, but its output is a database-conditioned similarity ledger: what each read most resembles in the chosen database, not a verified inventory and not an abundance. Useful results come from tuning three things - the database, the confidence threshold, and false-positive control - and from always re-estimating abundance with Bracken.

## Prerequisites
```bash
conda install -c bioconda kraken2 bracken krakentools

# Download a prebuilt database (sizes range from ~8 GB capped to ~70 GB Standard to ~480 GB nt)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz
mkdir kraken2_db && tar -xzf k2_standard_20240605.tar.gz -C kraken2_db
```

Conceptual prerequisites:
- The database defines what can be detected. Standard = bacteria/archaea/viral/plasmid/human/UniVec; PlusPF adds protozoa/fungi; PlusPFP adds plant; nt is broadest but slowest and most contamination-prone; GTDB is curated. Capped 8/16 GB databases are random k-mer downsamples, not curated subsets.
- Remove host reads before or during classification (include `human` in the database, or pre-deplete). Host reads dominate low-biomass samples.
- Prebuilt databases ship matching Bracken distributions (typically 100/150/200mers). Bracken `-r` must match both the database and the actual read length.

## Quick Start
Tell your AI agent what you want to do:
- "Classify my shotgun reads with Kraken2 and the Standard database, with a confidence threshold"
- "Run Kraken2 on paired-end reads and re-estimate species abundance with Bracken at read length 150"
- "Check whether my low-abundance hits are real using unique k-mer counts"
- "Remove human reads before classifying my stool metagenome"

## Example Prompts

### Choosing a database and confidence
> "I have paired-end gut metagenome reads at 150 bp. Run Kraken2 with the PlusPF database, set confidence to suppress single-k-mer false positives, require at least two hit groups, then re-estimate species abundance with Bracken at read length 150."

### False-positive control for low biomass
> "This is a low-biomass clinical sample. Classify with high confidence, report distinct minimizer counts, and flag any species that has many reads but few distinct minimizers as a likely false positive."

### Host removal first
> "These are human stool shotgun reads. Remove human reads against a T2T reference before classifying, then profile with Kraken2 and the Standard database."

### Diagnosing a low classification rate
> "Only 25% of my reads classified. Help me decide whether that means novel taxa, host contamination, or the wrong database."

## What the Agent Will Do
1. Confirm the Kraken2 and Bracken versions, the database build, and the matching Bracken read length.
2. Recommend (or apply) host removal if the reads are not host-depleted.
3. Run classification with a non-default confidence and `--minimum-hit-groups`, producing a `.kreport`.
4. Apply false-positive control (unique minimizers / KrakenUniq) on tail taxa when biomass is low.
5. Re-estimate abundance with Bracken at the correct read length and level.
6. State explicitly that Bracken `fraction_total_reads` is a read fraction, not a cell fraction, and hand off genome-size caveats to abundance-estimation.

## Tips
- Default `--confidence 0` over-classifies by design; raise it. On a capped database, high confidence collapses classification toward zero - use a larger database instead.
- `--minimum-hit-groups 2` is a cheap precision lever the defaults rarely advertise; raise to 3 for clinical samples.
- A taxon at the bottom of the report is a hypothesis. Confirm it with distinct minimizers or KrakenUniq before believing it.
- Bracken never removes false positives; run confidence, hit-groups, and unique-k-mer filtering first.
- A low classification rate can mean novel taxa, host contamination, or the wrong database - diagnose which rather than assuming novelty.

## Output Files

| File | Description |
|------|-------------|
| out.kraken | Per-read classifications (large; can be dropped to /dev/null) |
| out.kreport | Taxonomic report that feeds Bracken |
| out.bracken | Bracken abundance table (read fractions) |

## Related Skills

- abundance-estimation - Bracken command mechanics and read-count-to-abundance conversion
- metaphlan-profiling - Marker-gene alternative; FP-conservative
- metagenome-visualization - Plot and run community stats on profiles
- contamination-controls - Host depletion, blanks, and decontam before classification
- genome-assembly/metagenome-assembly - Assembly/MAG recovery; this category is read-based
- read-qc/contamination-screening - Host/vector read screening before classification
- workflows/metagenomics-pipeline - End-to-end shotgun profiling

## Resources
- [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)
- [Prebuilt databases](https://benlangmead.github.io/aws-indexes/k2)
- [Bracken GitHub](https://github.com/jenniferlu717/Bracken)
