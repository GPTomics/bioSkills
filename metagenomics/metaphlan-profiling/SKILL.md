---
name: bio-metagenomics-metaphlan
description: Profiles shotgun metagenomes to species/SGB relative abundance with MetaPhlAn 4's clade-specific marker genes (bowtie2 short reads, minimap2 long reads). Covers why a MetaPhlAn percentage is a cell fraction (genome-size-normalized taxonomic abundance) and must never be merged with Kraken/Bracken read fractions, kSGB vs uSGB units for quantifying database-absent taxa, the unknown-fraction rescaling and its version-default flip, --index pinning as a batch variable, and when mOTUs3 or sourmash gather beat marker profiling. Use when profiling who-is-there with high precision, needing HMP-comparable species abundances, quantifying novel taxa, or deciding marker-gene vs k-mer profiling. For k-mer classification see kraken-classification; for strains see strain-tracking; for 16S amplicon see the microbiome category.
tool_type: cli
primary_tool: MetaPhlAn
---

## Version Compatibility

Reference examples tested with: MetaPhlAn 4.1+, Bowtie2 2.5.3+, minimap2 2.26+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `metaphlan --version` then `metaphlan --help` to confirm flag names and defaults
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The marker DATABASE version is the experimental variable. Results track the index (e.g. `mpa_vJun23_CHOCOPhlAnSGB_202403` vs the live `vJan25` build); MetaPhlAn 3 and MetaPhlAn 4 databases are not interchangeable. Pin `--index` and report it like a reagent lot. Two flags were renamed in 4.2: `--bowtie2out` -> `--mapout` and `--bowtie2db` -> `--db_dir` (the `--input_type` value `bowtie2out` likewise becomes `mapout`); unknown-fraction estimation flipped from opt-in (`--unclassified_estimation`) to on-by-default (`--skip_unclassified_estimation` to disable). Confirm against `metaphlan --help`.

# MetaPhlAn Profiling

**"Who is in my metagenome, by cell fraction?"** -> Detect which clades' private marker genes are present, average their per-marker coverage, and normalize to a genome-size-aware relative abundance - so the percentage is a fraction of cells, not of reads.
- CLI: `metaphlan reads_1.fq.gz,reads_2.fq.gz --input_type fastq --index mpa_vJun23_CHOCOPhlAnSGB_202403 -o profile.txt --mapout sample.bz2`

Scope: marker-gene species/SGB profiling and its alternatives (mOTUs3, sourmash gather). K-mer read classification -> kraken-classification. Strain-resolved SNV haplotypes -> strain-tracking. Functional profiling -> functional-profiling. Compositional stats and plotting -> metagenome-visualization. 16S amplicon -> the microbiome category.

## The Single Most Important Modern Insight -- A MetaPhlAn Percentage Is a Cell Fraction, Not a Read Fraction

A MetaPhlAn percentage estimates what fraction of the CELLS in the community belong to a clade - a genome-size-normalized taxonomic abundance. A Kraken/Bracken percentage estimates what fraction of the READS came from a clade - a sequence abundance. There is no sample-independent conversion between them, because sequence abundance under-estimates small-genome microbes and over-estimates large-genome ones by a factor that depends on the whole community's genome-size distribution (Sun 2021 *Nat Methods* 18:618). Therefore:

- Never merge MetaPhlAn percentages with Kraken/Bracken percentages into one table, correlate them, or benchmark one against the other. Disagreement between them is expected even when both are correct.
- Marker profiling is not "classify every read." It detects which clades' PRIVATE markers are present (default presence gate: reads cover roughly 20% of a clade's markers) and averages their per-marker coverage. Most reads are never assigned - by design, not failure.

Mnemonic: markers measure WHO is there (cells); k-mers measure HOW MUCH DNA is there (reads).

## SGBs: the Unit Is Species-Level, and uSGBs Quantify the Unnamed

MetaPhlAn 4's atomic taxon is the SGB (species-level genome bin, a ~95% ANI cluster), not an NCBI species. A kSGB contains a cultured reference genome and gets a Latin name; a uSGB is defined only from MAGs (>=5 required) and is reported with a placeholder ID and no name. Quantifying uSGBs - taxa with no reference genome - is MetaPhlAn 4's headline advance over MetaPhlAn 3 and explains ~20% more gut reads, >40% more in under-characterized environments (Blanco-Miguez 2023 *Nat Biotechnol* 41:1633). Consequences: an unnamed `t__SGB...` row is a real quantified taxon - do not drop it; one named species can split into several SGBs; MetaPhlAn 3 species profiles and MetaPhlAn 4 SGB profiles are not row-compatible (use `sgb_to_gtdb_profile.py` for GTDB names). The `t__` tier is the SGB, NOT a strain - strain resolution is StrainPhlAn (-> strain-tracking).

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| MetaPhlAn 4 | Blanco-Miguez 2023 *Nat Biotechnol* 41:1633 | ~189 clade-specific markers/SGB; robust coverage average | high-precision species/SGB %, HMP-comparable, characterized communities |
| mOTUs3 | Ruscheweyh 2022 *Microbiome* 10:212 | 10 universal single-copy marker genes | higher recall of novel/divergent taxa; transparent marker-hit confidence |
| sourmash gather | Pierce 2019 *F1000Res* 8:1006 | FracMinHash containment, minimum metagenome cover | genome-resolved hits vs all of GTDB + an honest unknown fraction |
| Kraken2 + Bracken | Wood 2019 *Genome Biol* 20:257 | k-mer LCA + Bayesian reestimation | -> kraken-classification; max recall, willing to filter false positives |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Human gut species %, low false positives, HMP-comparable | MetaPhlAn 4 | curated SGB markers; high precision; huge corpus |
| Quantify novel / database-absent taxa | MetaPhlAn 4 uSGBs OR mOTUs3 ext-mOTUs | reference-independent units |
| Maximize recall in under-characterized environments | mOTUs3 or sourmash gather | universal markers / containment vs everything |
| Genome-resolved + explicit unknown fraction | sourmash gather | minimum metagenome cover reports what is unexplained |
| Max recall of every read, speed | -> kraken-classification | k-mer LCA; filter the false-positive tail |
| Need cell fraction, not read fraction | MetaPhlAn / mOTUs | k-mer tools report read fraction |
| Strain-level resolution | -> strain-tracking | per-SNV haplotypes, not species profiling |
| Composition stats next | -> metagenome-visualization (CLR/ANCOM-BC) | output is closed; naive stats on percentages are invalid |

## Basic Profiling

```bash
# Paired-end reads are passed as ONE comma-separated argument (MetaPhlAn treats them as two
# single-end files - it does not use insert/pairing info). Pin the index for reproducibility.
metaphlan reads_R1.fastq.gz,reads_R2.fastq.gz \
    --input_type fastq \
    --index mpa_vJun23_CHOCOPhlAnSGB_202403 \   # pin it; DB version is a batch variable
    --nproc 8 \
    --mapout sample.map.bz2 \                    # cache the read->marker mapping (pre-4.2: --bowtie2out)
    --output_file profile.txt
```

## Re-Profile from the Mapping Cache (the real operational lever)

**Goal:** Try different analysis types, levels, or estimator settings without realigning.

**Approach:** Save the mapping once with `--mapout`, then re-run from it with `--input_type mapout` (pre-4.2: `bowtie2out`). Realignment is the expensive step; everything downstream is free.

```bash
metaphlan sample.map.bz2 --input_type mapout \
    --tax_lev s \           # k,p,c,o,f,g,s,t (t = SGB tier)
    --stat_q 0.2 \          # quantile-truncated robust mean of per-marker coverages: drop top/bottom 20%, average the middle 60%
    --output_file profile_species.txt
```

`--stat_q` down-weights markers in HGT/mobile and conserved cross-clade regions; the default 0.2 is a sensible robust mean. Changing it changes the reported abundances - report it if it is changed. Long reads (4.1+) route to minimap2 with `--long_reads`.

## The Unknown Fraction Rescales Everything

Relative abundance sums to 100% only over DETECTED clades. With unknown estimation OFF (pre-4.2 default), known taxa absorb 100% and the database-absent community is invisible - overstating every known taxon. With it ON (4.2 default), an `UNCLASSIFIED` row appears and every known abundance shrinks proportionally. In soil/marine/rumen the unknown fraction can be the largest "taxon" in the sample.

```bash
# 4.2 default includes the UNCLASSIFIED row. To force it on pre-4.2: --unclassified_estimation
# For SAM input, pass --nreads <total> or the unknown fraction is wrong.
metaphlan reads.fastq.gz --input_type fastq -o profile.txt   # 4.2: UNCLASSIFIED row present by default
```

Pre-4.2-default and 4.2-default outputs are not comparable abundances - mixing them is a hidden batch effect.

## Merge and Convert

```bash
# All inputs MUST come from the SAME database index or rows mismatch silently.
merge_metaphlan_tables.py profiles/*_profile.txt > merged_abundance.txt
sgb_to_gtdb_profile.py -i merged_abundance.txt -o merged_gtdb.txt   # recover GTDB names for SGBs
```

## Per-Method Failure Modes

### MetaPhlAn percentages merged with Kraken percentages
**Trigger:** putting MetaPhlAn and Bracken abundances in one matrix or correlating them. **Mechanism:** cell fraction vs read fraction - different quantities (Sun 2021). **Symptom:** "tools disagree," spurious scatter, broken ML/differential-abundance features. **Fix:** keep them separate; if harmonizing, convert via genome length (Bracken counts / genome length, renormalize) and accept it is approximate.

### Unknown-fraction default mismatch across samples
**Trigger:** profiles built with different MetaPhlAn versions or `--unclassified_estimation` settings. **Mechanism:** the UNCLASSIFIED row rescales all known abundances. **Symptom:** a batch effect aligned to processing date, not biology. **Fix:** pin one version and one unknown-estimation setting across the whole study; for environmental samples always include the unknown fraction.

### Treating a low mapping rate as a QC failure
**Trigger:** alarm at <1% of reads mapping. **Mechanism:** only clade-specific markers are targeted; low mapping is expected. **Symptom:** unnecessary re-runs. **Fix:** low mapping is normal; a large unknown fraction means database-absent community (consider mOTUs3/sourmash), and a very low rate plus low microbial yield suggests host contamination -> contamination-controls.

### Recall ceiling in under-characterized environments
**Trigger:** profiling soil/marine and reporting only named taxa. **Mechanism:** a marker tool is structurally blind to clades whose markers are not in the database (high precision, low recall; CAMI2 Meyer 2022). **Symptom:** most of the community missing; lowering thresholds does not recover it. **Fix:** use mOTUs3 (universal markers) or sourmash gather (containment vs all of GTDB), or accept Kraken false positives and filter - do not just lower MetaPhlAn thresholds and call it sensitivity.

### Index mismatch on merge
**Trigger:** merging profiles built on different `--index` builds. **Mechanism:** SGB IDs and marker sets differ between releases. **Symptom:** rows silently fail to align; abundances look implausible. **Fix:** rebuild all samples on one pinned index before merging.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Presence gate ~20% of an SGB's markers | Blanco-Miguez 2023 *Nat Biotechnol* 41:1633 | enough markers covered to call a clade present (precision mechanism) |
| `--stat_q` 0.2 default | MetaPhlAn docs | truncated mean drops top/bottom 20% of marker coverages; robust to HGT/conserved outliers |
| uSGB requires >=5 MAGs | Blanco-Miguez 2023 *Nat Biotechnol* 41:1633 | false-positive control for unnamed taxa |
| Pin `--index` | MetaPhlAn docs | DB version changes profiles for identical reads; report like a reagent lot |
| `--min_cu_len` 2000 | MetaPhlAn docs | minimum cumulative marker length to report a clade (low-evidence filter) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "No database found" | DB not installed | `metaphlan --install` (optionally `--index <ver> --db_dir DIR`) |
| Output all zeros | wrong `--input_type` or empty/host-only input | match `--input_type` to the file; check microbial yield |
| `--bowtie2out` not recognized | running MetaPhlAn 4.2+ | use `--mapout` / `--input_type mapout` (4.2 rename) |
| Rows mismatch after merge | profiles from different indices | rebuild on one pinned `--index` |
| SAM input unknown fraction wrong | `--nreads` not supplied | pass total read count with `--nreads` |
| Viral calls look unreliable | `--add_viruses` calls are low-confidence | treat vSGB calls cautiously (CAMI2) |

## References

- Blanco-Miguez A, Beghini F, Cumbo F, et al. 2023. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nat Biotechnol* 41:1633-1644.
- Beghini F, McIver LJ, Blanco-Miguez A, et al. 2021. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife* 10:e65088.
- Sun Z, Huang S, Zhang M, et al. 2021. Challenges in benchmarking metagenomic profilers. *Nat Methods* 18:618-626.
- Meyer F, Fritz A, Deng ZL, et al. 2022. Critical Assessment of Metagenome Interpretation: the second round of challenges. *Nat Methods* 19:429-440.
- Ruscheweyh HJ, Milanese A, Paoli L, et al. 2022. Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments. *Microbiome* 10:212.
- Sunagawa S, Mende DR, Zeller G, et al. 2013. Metagenomic species profiling using universal phylogenetic marker genes. *Nat Methods* 10:1196-1199.
- Pierce NT, Irber L, Reiter T, Brooks P, Brown CT. 2019. Large-scale sequence comparisons with sourmash. *F1000Res* 8:1006.

## Related Skills

- kraken-classification - K-mer read classification; reports read fraction, not cell fraction
- abundance-estimation - Compositional handling and cross-tool abundance comparison
- strain-tracking - StrainPhlAn strain resolution below the SGB level
- functional-profiling - HUMAnN reuses a MetaPhlAn profile for its taxonomic prescreen
- metagenome-visualization - Compositional stats and plotting of profiles
- genome-assembly/metagenome-assembly - Recover the MAGs that define uSGBs; this category is read-based
- workflows/metagenomics-pipeline - End-to-end shotgun profiling
