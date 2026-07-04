---
name: bio-comparative-genomics-genome-distance-and-species-delineation
description: Compute genome-to-genome distances (ANI, AAI, dDDH, k-mer Mash) and assign taxonomic classifications using skani (Shaw 2023), FastANI (Jain 2018), pyani / pyANI ANIb / ANIm, OrthoANI (Lee 2016), AAI (amino-acid identity), dDDH via TYGS / GGDC, GTDB-Tk (Chaumeil 2020 standard prokaryote taxonomy), and Mash MinHash (Ondov 2016). Use when delineating prokaryote species (95% ANI threshold; Jain 2018 Nat Commun 9:5114), assigning genomes to GTDB taxonomy with ANI radius, computing genome similarity matrices for clustering, classifying archaea, evaluating MAG (metagenome-assembled genome) species assignment, applying skani for fast metagenomic ANI screening, or reconciling 16S rRNA-based taxonomy with whole-genome ANI.
tool_type: cli
primary_tool: skani
---

## Version Compatibility

Reference examples tested with: skani 0.2.5+ (Shaw & Yu 2023 Nat Methods 20:1661; bluenote-1577/skani), FastANI 1.34+ (Jain 2018 Nat Commun 9:5114), pyani 0.3.0+ (Pritchard 2016 Anal Methods 8:12), pyskani 0.1+ (Larralde 2025), OrthoANI 1.40+ (Lee 2016 Int J Syst Evol Microbiol 66:1100), OrthoANIu 1.2+, GTDB-Tk 2.7.1+ (Chaumeil 2022 Bioinformatics 38:5315), GTDB release 220 (2024-Q3+), TYGS web (Meier-Kolthoff & Goker 2019 Nat Commun 10:2182), GGDC v3.0 (web), Mash 2.3+ (Ondov 2016 Genome Biol 17:132), Dashing 2 (Baker & Langmead 2023 Genome Res 33:1218), CompareM 0.1.2+ for AAI (Parks/Cherubini), pyANI 0.3.1+, BLAT 36+, DIAMOND 2.1+. JSpeciesWS web (Richter & Rossello-Mora 2009 PNAS 106:19126).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `skani --version`; `fastANI --version`; `gtdbtk --version`; `mash --version`; `pyani --version`
- Python: `pip show gtdbtk pyani`

If code throws `GTDB-Tk database not found`, `skani sketch incompatible`, `Mash sketch version`, these tools have database-version coupling: GTDB-Tk requires the GTDB release matched to the binary version; Mash/skani sketches are forward-compatible but not always backward. Check `gtdbtk check_install` for database completeness.

# Genome Distance and Species Delineation

**"Are these genomes the same species, and what species are they?"** -> Prokaryote species delineation has shifted from 16S rRNA identity (now considered insufficient at < 98.7%) to **whole-genome ANI** at a 95% threshold (Jain 2018 Nat Commun 9:5114; corroborating Goris 2007 and Konstantinidis 2005). The modern operational standard for taxonomy is **GTDB-Tk** (Chaumeil 2020/2022 Bioinformatics 38:5315), which assigns genomes to the Genome Taxonomy Database (GTDB) using ANI radius + marker-gene placement. **skani** (Shaw & Yu 2023 Nat Methods 20:1661) has replaced FastANI as the default ANI tool in GTDB-Tk 2.4+ for being 20-30x faster while maintaining accuracy. The 95% ANI threshold is robust but not absolute -- the species circumscription radius varies by genus (Parks 2018 Nat Biotech 36:996).

- CLI: `skani dist genomes1.fa genomes2.fa -t 16` -- fast ANI computation
- CLI: `fastANI -q query.fa -r reference.fa -o output.txt` -- standard ANI
- CLI: `gtdbtk classify_wf --genome_dir genomes/ --out_dir gtdbtk_out --cpus 32` -- GTDB classification
- CLI: `mash dist *.fa` -- k-mer MinHash distance
- Web: TYGS (https://tygs.dsmz.de/) and GGDC (https://ggdc.dsmz.de/) for dDDH

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| skani (Shaw & Yu 2023 Nat Methods 20:1661) | Sparse chaining on minimizers; ANI estimation | ANI percent + alignment fraction | 20-30x faster than FastANI; default in GTDB-Tk 2.4+; supports MAGs | Currently no AAI; not for cross-domain (archaea vs bacteria) |
| FastANI (Jain 2018 Nat Comm 9:5114) | Mashmap-based fragment alignment | ANI percent + orthologous fraction | Standard ANI tool 2018-2023; well-validated | Slower than skani; designed for >=80% identity |
| pyani / pyANI (Pritchard 2016 Anal Methods 8:12) | Multi-method ANI: ANIb (BLASTN), ANIm (MUMmer), TETRA | ANI matrix + visualization | Multiple algorithm consensus; reproducible | Slower than skani / FastANI; legacy for many studies |
| OrthoANI / OrthoANIu (Lee 2016 IJSEM 66:1100) | Reciprocal-best-orthologs ANI | ANI percent (more robust than blast-based) | Considered more precise than ANIb | Slower; less integrated |
| GTDB-Tk (Chaumeil 2020/2022 Bioinformatics 38:5315) | Marker-gene phylogeny + ANI radius (skani; was FastANI v2.3.x) | Taxonomic classification at all ranks (GTDB nomenclature) | Modern prokaryote taxonomy standard | Specific to GTDB; some classifications differ from NCBI |
| TYGS (Meier-Kolthoff & Goker 2019 Nat Commun 10:2182) | dDDH (digital DNA-DNA hybridization) | Pairwise dDDH + species delineation | Most rigorous species delineation (vs traditional DDH) | Web-only; rate-limited; specific platform |
| GGDC (Auch 2010 Stand Genomic Sci 2:117; v3 web) | Digital DDH calculation | dDDH percent + thresholds | Validated against laboratory DDH | Web-only; computational cost |
| Mash (Ondov 2016 Genome Biol 17:132) | MinHash k-mer sketches | Approximate distance (1 - similarity) | Extremely fast for large-scale clustering | k-mer-based; loses biological interpretation |
| Dashing 2 (Baker & Langmead 2023 GR 33:1218) | Sketching with Bloom filter optimization | Same as Mash but faster | 5-10x faster than Mash | Newer; less broadly used |
| pyskani (Larralde 2025 bioRxiv) | Python wrapper around skani | ANI in Python | Programmatic access; CI/CD friendly | Newer; ecosystem still developing |
| CompareM | All-vs-all AAI (amino acid identity) | AAI percent | Cross-genus comparison via protein | Slow; needs all proteomes |
| JSpeciesWS (web; Richter & Rossello-Mora 2009) | ANIb / ANIm | Web ANI + species delineation | Standard for clinical microbiology | Web rate limits; slow |
| ANI Calculator (CGB Korea) | Web ANI | Web ANI | Quick check | Web-only |

Methodology evolves; verify GTDB release (currently r220 / 2024-Q3) and GTDB-Tk version compatibility. The 95% ANI species threshold has been confirmed across 90,000+ prokaryote genomes (Jain 2018; Parks 2018 demonstrating clear bimodality).

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Classify a bacterial genome to species | GTDB-Tk classify_wf | Standard prokaryote taxonomy; ANI + marker-gene placement |
| Compute ANI between two genomes | skani | Fast (vs FastANI); accurate |
| Compute ANI for 1000+ genome pairs | skani all-vs-all | Scales; preferred for screening |
| Verify species delineation publication | TYGS + GGDC dDDH | Gold standard for novel species |
| MAG species assignment | GTDB-Tk + skani; CheckM2 first for completeness | MAGs need quality assessment + taxonomy |
| Bacterial strain typing | ANIb or dDDH; 99-99.99% for same strain | Strain resolution requires ANI > 99% |
| Sub-species / serotype level | ANI > 99.5% + epidemiological context | Sub-species requires biology + ANI |
| Across deep prokaryote divergence | AAI (CompareM); ANI saturates below 75% | AAI better for cross-genus |
| Cross-archaeal vs bacterial taxonomy | Skani-archaea-aware OR separate analysis | Default skani doesn't differentiate domains explicitly |
| Fast metagenomic taxonomy screen | Mash or Dashing 2 | k-mer-based; sketches reusable |
| Distance for genome clustering | skani matrix -> hierarchical clustering or NJ | Standard workflow |
| Reconcile 16S vs genome taxonomy | Run both; check for inter-genus conflicts | 16S < 98.7% typically inadequate |
| Build a reference database for ANI lookup | skani sketches indexed; query against | Pre-sketched reference for repeated queries |
| Phylogenetic placement | GTDB-Tk + IQ-TREE on extracted markers | Place new genome in known tree |
| AAI for genera-level comparison | CompareM or aai.rb (KBase) | AAI < ANI signal at deep divergence |
| Subspecies-level pathogenicity | ANI > 99% + virulence-gene annotation | ANI alone insufficient |
| Type-strain comparison | TYGS automatic type-strain matching | Built-in type-strain database |

## Per-Tool Failure Modes

### skani / FastANI ANI saturating below 75%

**Trigger:** Computing ANI between two genomes at < 75% nucleotide identity.

**Mechanism:** ANI is computed only on alignable regions; at < 75% identity, alignment fraction drops dramatically (< 50%); the few alignable regions are biased toward conserved regions, inflating apparent ANI.

**Symptom:** Reported ANI 75-80% with alignment fraction < 50%; meaningless biologically.

**Fix:** Below 75% ANI, switch to AAI (amino-acid identity from translated proteins); AAI is more meaningful at deep divergence (~50% AAI between distant genera). For ANI matrix at long range, use Mash distance (k-mer-based, no alignment).

### GTDB-Tk database version mismatch

**Trigger:** Using GTDB-Tk binary with mismatched GTDB reference data version.

**Mechanism:** GTDB releases (r207, r214, r220, ...) include reference trees, marker-gene HMMs, and ANI sketch files; GTDB-Tk versions are aligned to specific releases. Mismatch causes silent or loud failures.

**Symptom:** GTDB-Tk fails with "marker gene HMM not found" or runs but produces inconsistent classifications.

**Fix:** Check release compatibility: `gtdbtk check_install` shows the expected vs found versions. Pin via `conda env`. Download the matching release from data.gtdb.ecogenomic.org.

### Below-95% ANI but same species (genus-specific radius)

**Trigger:** Strict 95% ANI threshold; rejecting closely related genomes as different species.

**Mechanism:** Parks 2018 demonstrated species-circumscription radius varies (typically 95-99% but tighter for some clonal lineages). For some genera (e.g., Pseudomonas), the species threshold is 94% per genus-specific analysis.

**Symptom:** Two clearly biologically-related strains (epidemiologically connected outbreak) have ANI 94%; using strict 95% calls them different species.

**Fix:** Use GTDB-Tk's species radius approach which uses genus-specific cutoffs. For novel-species naming, consider 95% as primary + ecology + biology. Report ANI alongside additional context (gene content, phenotype, ecology).

### High alignment fraction required (AF >= 0.5)

**Trigger:** ANI reported without alignment fraction; or AF < 0.5.

**Mechanism:** A 95% ANI with 20% AF (only 20% of genome alignable) is biologically meaningless; the comparison covers a small fraction of the genomes. Standard species delineation requires AF >= 0.5 + ANI >= 95%.

**Symptom:** ANI calls "same species" but only 20-30% of genome aligned.

**Fix:** Require AF >= 0.5 + ANI >= 95% for species call. Below this, classification is ambiguous; consider AAI or HGT (e.g. a phage-rich genome will have low AF to its actual species).

### MAG contamination / incompleteness

**Trigger:** Running GTDB-Tk on MAGs without CheckM2 pre-screening.

**Mechanism:** Low-completeness MAGs (< 50%) may have marker-gene gaps that break GTDB-Tk's placement; high-contamination MAGs (> 5%) place anomalously.

**Symptom:** GTDB-Tk reports "no genome retained" or makes implausible classifications.

**Fix:** Pre-filter MAGs with CheckM2 (Chklovski 2023 Nat Methods 20:1203): require >= 70% completeness, < 5% contamination, < 10% strain heterogeneity for species-level. For lower-quality MAGs, report genus-level only or exclude.

### Mash distance vs ANI inconsistency

**Trigger:** Using Mash distance for species delineation directly.

**Mechanism:** Mash distance = 1 - (Mash similarity) is correlated with 1 - ANI but is not the same metric. The 0.05 Mash distance threshold (sometimes cited as "5% Mash = 95% ANI") is approximate; exact ANI varies +-0.5%.

**Symptom:** Genomes with Mash distance 0.04-0.06 inconsistently called "same species" or not.

**Fix:** Use Mash for fast screening + clustering, but verify species delineation with ANI on candidate pairs. Mash 0.05 ~ ANI 95% is a rough heuristic; not a publishable threshold.

### Ortholog-based ANI conservative vs alignment-based

**Trigger:** Comparing OrthoANI / OrthoANIu to FastANI / skani.

**Mechanism:** OrthoANI uses reciprocal-best-orthologs; FastANI uses fragment alignment. They differ by 0.5-2% systematically; OrthoANI is more conservative.

**Symptom:** OrthoANI = 94.5%, FastANI = 96% on same pair; different "same species" calls.

**Fix:** Document tool used; for taxonomy, use the tool aligned to the reference database (GTDB-Tk uses skani; NCBI uses FastANI). Cross-validation for unclear cases.

### Tetranucleotide bias inflating Mash for low-GC genomes

**Trigger:** Comparing genomes with extreme GC content (Streptomyces ~70% vs Mycoplasma ~25%) via Mash.

**Mechanism:** k-mer frequency distributions are GC-dependent; Mash distance between extreme-GC genomes is inflated by background composition rather than biology.

**Symptom:** Two unrelated extreme-GC genomes (e.g. Mycoplasma + Mycoplasma) appear closer than they biologically are.

**Fix:** Use ANI / AAI for cross-GC comparisons; Mash is reliable only within a GC-comparable range. Document GC range when reporting Mash distances.

### Type strain conflicts in TYGS

**Trigger:** Submitting genome to TYGS where type strain is missing or misclassified.

**Mechanism:** TYGS depends on type strain database; if the type strain is genome-sequenced incompletely or misclassified, the placement may be incorrect.

**Symptom:** TYGS reports unexpected nearest type strain.

**Fix:** Cross-validate with GTDB-Tk; check type-strain genome quality. For novel-species naming, TYGS report should be supplemented with manual taxonomic check.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Species delineation ANI | >= 95% (Jain 2018) | Standard; based on 90,000 prokaryote genomes |
| Genus-specific species ANI radius | 94-99% (varies by clade) | Parks 2018 Nat Biotech 36:996 |
| Alignment fraction (AF) for ANI | >= 0.5 (Jain 2018) | Below this, comparison too small |
| Strain delineation ANI | >= 99% (typical); >= 99.5% strict | Operational |
| Sub-species delineation | >= 99% ANI + epidemiology | Manual |
| AAI species delineation | >= 70% | Konstantinidis 2005; varies clade |
| AAI genus delineation | >= 60% | Standard |
| dDDH species delineation | >= 70% | Goris 2007; matches ANI 95% |
| dDDH genus delineation | >= 50% | Auch 2010 |
| 16S rRNA species threshold (deprecated) | >= 98.7% | Stackebrandt 2006; superseded by ANI |
| Mash distance ~ ANI heuristic | 0.05 ~ 95% ANI | Ondov 2016; rough |
| GTDB-Tk completeness for placement | >= 50% (CheckM2 OBLIGATE >= 70%) | Chaumeil 2022 |
| GTDB-Tk contamination | < 10% (< 5% for species-level) | Chaumeil 2022 |
| MAG quality for taxonomy | CheckM2 >= 70% comp, < 5% cont, < 10% strain het | Chklovski 2023 |
| skani CLI threads | up to 64; scales linearly | skani docs |
| skani sketch size | default 1000 minimizers; tunable | skani docs |
| FastANI default fragment | 3000 bp; varies | Jain 2018 default |
| GTDB-Tk classify_wf time per genome | 2-30 min on 16 CPUs | Empirical |
| Type strain ANI uncertainty | +- 1% | Operational |

## skani Standard Workflow

**Goal:** Compute ANI between query and reference set; classify species.

**Approach:** Build skani sketch -> compute distances -> apply species delineation logic.

```bash
# 1. Pre-sketch reference set
skani sketch reference_genomes/*.fa -o reference_sketches

# 2. Compute ANI for a single query
skani dist query.fa reference_sketches/* -t 16 \
    --robust --slow > query_distances.tsv

# 3. Compute all-vs-all matrix (large set)
skani triangle genomes/*.fa -t 32 --robust --sparse -o ani_matrix.tsv
# --sparse emits tabular Ref_file Query_file ANI Align_fraction_ref Align_fraction_query;
# without --sparse, `skani triangle` emits a Phylip-style square matrix.

# 4. Filter and visualize
awk '$3 >= 95 && $5 >= 50' query_distances.tsv > species_matches.tsv
```

```python
'''Apply 95% ANI species delineation with AF >= 0.5 constraint.'''
import pandas as pd
import numpy as np


def parse_skani(path):
    '''skani output: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query'''
    df = pd.read_csv(path, sep='\t')
    df.columns = ['ref_file', 'query_file', 'ani', 'af_ref', 'af_query']
    return df


def species_delineate(df, ani_threshold=95.0, af_threshold=0.5):
    '''Return genome pairs called same species.'''
    df['min_af'] = df[['af_ref', 'af_query']].min(axis=1)
    same_species = df[(df['ani'] >= ani_threshold) & (df['min_af'] >= af_threshold)]
    return same_species
```

## GTDB-Tk Classification Workflow

**Goal:** Assign GTDB taxonomy (kingdom -> species) to a set of bacterial/archaeal genomes.

**Approach:** GTDB-Tk classify_wf identifies markers, builds tree placement, calculates ANI radius, returns full taxonomy.

```bash
# Set environment
export GTDBTK_DATA_PATH=/path/to/release220_data

# Verify install
gtdbtk check_install

# Run classify_wf
gtdbtk classify_wf \
    --genome_dir genomes/ \
    --out_dir gtdbtk_out \
    --cpus 32 \
    --extension fa \
    --skip_ani_screen   # Skip if want phylogeny-based only

# Output:
#   gtdbtk_out/classify/gtdbtk.bac120.summary.tsv   bacterial classifications
#   gtdbtk_out/classify/gtdbtk.ar53.summary.tsv     archaeal classifications
#   gtdbtk_out/identify/                              marker gene tables
#   gtdbtk_out/align/                                 multiple sequence alignments
```

```python
'''Parse GTDB-Tk summary for ranked classification.'''
import pandas as pd


def parse_gtdbtk_summary(path):
    '''GTDB-Tk summary columns include: user_genome, classification, classification_method, ani, msa_percent'''
    df = pd.read_csv(path, sep='\t')
    df['classification_split'] = df['classification'].str.split(';')
    df['kingdom'] = df['classification_split'].str[0].str.replace('d__', '')
    df['phylum'] = df['classification_split'].str[1].str.replace('p__', '')
    df['class'] = df['classification_split'].str[2].str.replace('c__', '')
    df['order'] = df['classification_split'].str[3].str.replace('o__', '')
    df['family'] = df['classification_split'].str[4].str.replace('f__', '')
    df['genus'] = df['classification_split'].str[5].str.replace('g__', '')
    df['species'] = df['classification_split'].str[6].str.replace('s__', '')
    return df


df = parse_gtdbtk_summary('gtdbtk_out/classify/gtdbtk.bac120.summary.tsv')
df_species_level = df[df['species'] != '']  # Species-level classification
df_high_quality = df[df['msa_percent'] > 80]
```

## TYGS / GGDC for dDDH

For publication-grade novel-species naming, dDDH is required. Use the TYGS web service (https://tygs.dsmz.de/) or GGDC (https://ggdc.dsmz.de/). Submit FASTA -> receive dDDH against type strains + phylogeny.

## Mash for Fast Clustering

```bash
# Sketch all genomes (do once)
mash sketch -p 16 -o all_sketches genomes/*.fa

# All-vs-all distance
mash dist -p 16 all_sketches.msh all_sketches.msh > mash_distances.tsv

# Cluster via NJ or hierarchical clustering
python -c "
import pandas as pd
from scipy.cluster import hierarchy
df = pd.read_csv('mash_distances.tsv', sep='\t', header=None,
                 names=['ref', 'query', 'distance', 'pvalue', 'shared_hashes'])
matrix = df.pivot('ref', 'query', 'distance').fillna(0)
linkage = hierarchy.linkage(matrix.values, method='average')
hierarchy.dendrogram(linkage, labels=matrix.columns)
"
```

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| skani says 96%, FastANI says 95.5% | Tool difference; both > 95% threshold | Consistent species call |
| skani says 95.2%, FastANI says 94.8% | Borderline case, tool difference | Cross-validate with OrthoANI; check GTDB-Tk |
| ANI = 94%, dDDH = 70% | Different metrics; both borderline same-species | Standard concordant ambiguous; report both |
| 16S rRNA same; ANI < 95% | 16S insufficient for species (Stackebrandt 2006) | Trust ANI; reclassify |
| GTDB says species X, NCBI says species Y | GTDB taxonomy differs from NCBI for ~10% of species | Use GTDB for genome-based; cite both if external comparison |
| ANI 96%, AF 0.3 | Insufficient alignment fraction | Same-species call invalid; report cautious |
| ANI 99%, AF 0.9 | Same species, high confidence | Robust species call |
| Mash 0.04 distance, FastANI 94.5% ANI | Different units; both around species threshold | Use ANI for definitive call |
| OrthoANI 94%, FastANI 96% | Method variation | Cross-check with skani; report range |
| GTDB-Tk classifies but MSA % < 50 | Low-quality placement | Report at higher rank (genus) only |
| TYGS gives unexpected nearest type strain | Type-strain database issue | Verify type-strain quality; cross-validate with GTDB |
| MAG GTDB-Tk classify fails | Incomplete MAG; missing markers | Improve assembly; report at family or order only |

**Operational rule for publication:** GTDB-Tk classify_wf as primary classification + ANI to nearest type strain (skani or FastANI) + AF >= 0.5 + report Tettelin partition where relevant + cross-validate with dDDH (TYGS) for novel-species claims.

## Cohort Gotchas

- **Archaea:** GTDB has separate ar53 marker set; specify `--archaea` or let GTDB-Tk auto-detect
- **Cyanobacteria:** large genomes, sometimes split GTDB-Tk markers; cross-validate
- **MAGs from metagenomes:** require CheckM2 quality filter first
- **Strain-level resolution:** ANI > 99% needed; epidemiological context useful
- **Endosymbionts:** small / reduced genomes; ANI may be unreliable
- **High-GC genomes (Streptomyces, Mycobacterium):** GTDB-Tk specific markers handle these
- **Type strains:** authoritative anchor for taxonomy; TYGS automatic
- **Pre-2017 ANI publications:** likely used 30-fold lower-precision tools; verify current ANI
- **Sub-species naming:** ANI 99% + biology + epidemiology; not standardized
- **Cross-domain comparisons (Bacteria vs Archaea):** rarely meaningful; AAI better

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why 95% ANI?" | Jain 2018 Nat Comm 9:5114; 90,000 prokaryote genomes demonstrate clear bimodality at 95% |
| "Genus-specific radius?" | GTDB-Tk uses genus-specific ANI radius; ranges 94-99% per genus |
| "AF >= 0.5 reported?" | Yes; standard convention; sub-0.5 invalidates species call |
| "Tool choice?" | skani 2.4+ default in GTDB-Tk; FastANI for verification; OrthoANI as third check for borderline cases |
| "GTDB vs NCBI?" | GTDB taxonomy is genome-based; NCBI is heritage; cite GTDB primary, NCBI secondary |
| "Type strain comparison?" | TYGS automatic type-strain matching; reported alongside ANI |
| "MAG quality?" | CheckM2 >= 70% completeness, < 5% contamination required; reported per MAG |
| "dDDH for novel species?" | TYGS performed; dDDH >= 70% threshold for same species |
| "Multiple methods agree?" | skani + FastANI + OrthoANI converged within 1% on candidate pairs |
| "GTDB release version?" | r220 (2024-Q3); database version pinned |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| skani "sketch incompatible" | Old sketch + new skani | Rebuild sketches; use same skani version |
| GTDB-Tk "marker not found" | Wrong GTDB version | Update GTDB-Tk + database to matching release |
| FastANI "no alignment" | < 75% ANI | Switch to AAI or Mash |
| pyani memory exhaustion | > 50 genomes | Use skani or pre-cluster with Mash |
| Mash distance pvalue uninformative | Default | Use NA filter; trust distance directly |
| TYGS web rate limit | Submitting too many genomes | Batch submissions; use offline GGDC if available |
| GTDB-Tk classify_wf classification empty | Genomes < 50% complete | Filter MAGs with CheckM2 first |
| OrthoANI very slow | All-vs-all on > 50 genomes | Use skani for screening |
| skani sketch GTDB databases | GTDB-Tk auto-sketches | Don't re-sketch separately |
| Custom species ID needed | TYGS supports type-strain comparison + GGDC for dDDH | Use both |
| ANI vs 16S contradict | 16S insufficient | Trust ANI |
| Subspecies confusion | ANI > 99% + biology required | Standardize naming |

## Tool Installation Notes

```bash
# skani
conda install -c bioconda skani

# FastANI
conda install -c bioconda fastani

# GTDB-Tk
conda install -c bioconda gtdbtk
# Download database
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_r220_data.tar.gz
tar xf gtdbtk_r220_data.tar.gz
export GTDBTK_DATA_PATH=$PWD/release220

# pyani / pyANI
pip install pyani-plus

# Mash + Dashing 2
conda install -c bioconda mash
conda install -c bioconda dashing2

# CheckM2 for MAG QC
conda install -c bioconda checkm2

# OrthoANI / OrthoANIu
git clone https://github.com/EzbioCloud-Bioinformatics-Team/OrthoANIu
```

For 1000+ genome scans, use cluster with 64+ cores; skani all-vs-all on 1000 genomes runs in ~10-30 minutes; GTDB-Tk classify_wf runs in ~20-60 min per 100 genomes.

## References

- Jain C et al 2018 Nat Commun 9:5114 (FastANI; 95% ANI threshold)
- Shaw J & Yu YW 2023 Nat Methods 20:1661 (skani)
- Chaumeil P-A et al 2020 Bioinformatics 36:1925 (GTDB-Tk v1)
- Chaumeil P-A et al 2022 Bioinformatics 38:5315 (GTDB-Tk v2)
- Parks DH et al 2018 Nat Biotech 36:996 (GTDB establishment)
- Parks DH et al 2022 Nucleic Acids Res 50:D785 (GTDB r207)
- Meier-Kolthoff JP & Goker M 2019 Nat Commun 10:2182 (TYGS)
- Auch AF et al 2010 Stand Genomic Sci 2:117 (GGDC)
- Ondov BD et al 2016 Genome Biol 17:132 (Mash MinHash)
- Baker DN & Langmead B 2023 Genome Res 33:1218 (Dashing 2)
- Pritchard L et al 2016 Anal Methods 8:12 (pyani)
- Lee I et al 2016 IJSEM 66:1100 (OrthoANIu)
- Richter M & Rossello-Mora R 2009 PNAS 106:19126 (JSpeciesWS)
- Goris J et al 2007 IJSEM 57:81 (ANI threshold validation)
- Konstantinidis KT & Tiedje JM 2005 PNAS 102:2567 (AAI concept)
- Stackebrandt E & Ebers J 2006 Microbiol Today 33:152 (16S thresholds; superseded)
- Chklovski A et al 2023 Nat Methods 20:1203 (CheckM2)
- Olm MR et al 2017 ISME J 11:2864 (dRep; ANI clustering)
- Rodriguez-R LM & Konstantinidis KT 2014 Microbe 9:111 (ANI reference)
- Yoon S-H et al 2017 IJSEM 67:1613 (ANI species delimitation)
- Larralde M et al 2025 bioRxiv (pyskani, pyfastani, pyorthoani)

## Related Skills

- comparative-genomics/pangenome-analysis - ANI-based clustering precedes pangenome construction
- comparative-genomics/ortholog-inference - Cross-species ANI as orthology benchmark
- comparative-genomics/hgt-detection - High-ANI same-species genomes for HGT context
- comparative-genomics/gene-tree-species-tree-reconciliation - Species-tree construction precedes ANI species delineation
- phylogenetics/species-trees - Marker-gene tree alongside ANI
- metagenomics/kraken-classification - Metagenomic classification different problem
- metagenomics/metaphlan-profiling - Profile vs taxonomic placement
- genome-assembly/assembly-qc - Quality before classification
- read-qc/quality-reports - CheckM2 on MAGs
- variant-calling/clinical-interpretation - Pathogen typing context
