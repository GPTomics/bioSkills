---
name: bio-epidemiological-genomics-variant-surveillance
description: Assigns pathogen lineages (SARS-CoV-2 Pangolin UShER mode; Nextclade clade + QC; pango-designation alias resolution) and tracks variant frequencies over time using Nextstrain (Augur + Auspice), wastewater deconvolution (Freyja, COJAC, alcov, lineagespot), lineage-fitness modelling (multinomial logistic), and recombinant detection (3SEQ, RDP4, Bolotie). Covers Pangolin pangolin-data and Nextclade dataset version pinning (mandatory; lineage-defining mutations change with dataset), Freyja barcode forward-only date constraint, ARTIC primer scheme churn (V3/V4/V4.1/V5.3.2/Midnight) with dropout regions, and recombinant X-prefix designation lag. Use when assigning Pango lineages and Nextclade clades to viral consensus sequences, building Nextstrain Augur surveillance pipelines, deconvolving wastewater into lineage frequencies with Freyja, tracking lineage frequencies over time, handling ARTIC primer dropouts, or running surveillance for SARS-CoV-2/influenza/Mpox/RSV/H5N1/measles.
tool_type: mixed
primary_tool: Pangolin
---

## Version Compatibility

Reference examples tested with: pangolin 4.3+ (pangolin-data 1.30+), nextclade 3.8+, augur 24.0+, freyja 1.4+, cojac 0.9+, lineagespot 1.6+ (Bioconductor), usher 0.6+, matUtils 0.6+, samtools 1.20+, lofreq 2.1+, ivar 1.4+, ARTIC pipeline 1.3+, snakemake 8.5+, pandas 2.2+, BioPython 1.84+, jq 1.7+.

Before using code patterns, verify installed versions match. If versions differ:
- `pangolin --all-versions` -- prints pangolin + pangolin-data + scorpio + constellations versions
- `nextclade dataset list --tag latest sars-cov-2` -- list current dataset tags
- `freyja --version`; `freyja barcode-build --help` (note: HYPHEN, not underscore; some legacy docs show `barcode_build`)
- `nextclade run --help` -- v3+ syntax replaced v2; old `nextclade` invocation no longer works
- `augur --version`; `augur refine --help` for current root-strategy flags

If `pangolin --inference usher` is rejected, the flag is `--analysis-mode usher` (no `--inference`). If `nextclade --input-dataset DIR` works, the installed version may be v2; v3 accepts both but `--dataset NAME` is the modern form for built-in datasets. Pangolin and Nextclade output column names differ between major releases -- introspect rather than retry.

# Variant Surveillance

**"Which lineages are circulating, and how fast are they growing?"** -> Assign consensus or wastewater samples to a curated lineage / clade nomenclature, then track frequencies over time with explicit version pinning. The lineage assignment is NOT a stable property of the sequence; it is a property of the sequence interpreted by a specific pangolin-data / Nextclade-dataset version. Two labs running the same Pangolin binary with different pangolin-data versions can produce different calls on the same genome. For published or regulatory output, pin BOTH the executable AND the dataset version (`pangolin --all-versions`; `nextclade dataset list --tag latest`), and re-run the whole archive after every dataset update -- comparing today's BA.2.86 call to last month's "Unassigned" call is invalid.

- CLI: `pangolin sequences.fasta --analysis-mode usher --outfile lineage_report.csv` -- UShER mode is the default since v4 (pangoLEARN deprecated mid-2023)
- CLI: `nextclade run --input-dataset nc_dataset/sars-cov-2 --output-tsv nc.tsv sequences.fasta` -- clade + Pango + QC + mutations
- CLI: `freyja variants sample.bam --variants sample.variants.tsv --depths sample.depths.tsv --ref reference.fa` then `freyja demix sample.variants.tsv sample.depths.tsv --output sample.demix.tsv` -- wastewater lineage deconvolution
- CLI: `augur refine --tree tree.nwk --alignment aln.fasta --metadata meta.tsv --output-tree refined.nwk --root oldest --timetree` -- Nextstrain time-scaling

## The Single Most Important Modern Insight -- Lineage assignment is dataset-version-dependent

A SARS-CoV-2 sequence called BA.5 today might be called BA.5.2.1 next week and KP.3 a month after that. pangolin-data and nextclade-dataset are updated weekly; lineage definitions evolve through pango-designation GitHub issues, often days-to-weeks before pangolin-data releases include the lineage. During the lag window, the same genome submitted in lab A (older pangolin-data) and lab B (current) gets different calls. The cross-lab "different lineage" result is then misread as biology. For any report, pin BOTH the executable AND the dataset version with `pangolin --all-versions` and `nextclade dataset list --tag latest` recorded alongside the call. For longitudinal studies, re-run the WHOLE archive after every dataset update -- comparing today's BA.2.86 call against last month's "Unassigned" call is invalid. Second-order insight: Pangolin's pangoLEARN mode was officially deprecated mid-2023 in favour of UShER mode (Pongmoragot 2024 *Virus Evol* 10:vead085); cross-study comparison of XBB sub-lineage prevalence from 2022 - mid-2023 is contaminated by the pangoLEARN -> UShER mode switch even when the same pangolin-data version is used.

## Algorithmic Taxonomy

| Tool | Mechanism | Inputs | Output | Strength | Fails when |
|------|-----------|--------|--------|----------|------------|
| Pangolin UShER mode (O'Toole 2021 *Virus Evol* 7:veab064; Pongmoragot 2024 *Virus Evol* 10:vead085) | Parsimony placement on daily-updated UShER mutation-annotated tree | SARS-CoV-2 consensus | Pango lineage call | UShER is the default since v4; more accurate than pangoLEARN for recent / divergent lineages | Designation lag for emerging lineages; recombinants require manual Pango-X designation |
| Pangolin pangoLEARN mode | Random-forest classifier trained on pangolin-data | SARS-CoV-2 consensus | Pango lineage call | Fast | DEPRECATED mid-2023; less accurate than UShER for novel sub-lineages |
| Nextclade (Aksamentov 2021 *JOSS* 6:3773) | Reference-tree placement + clade assignment + mutation calling + QC | Viral consensus (multi-pathogen) | Clade + Pango + QC + mutations | Integrated alignment QC; mutation outliers; recombination indicators | Dataset version drift changes lineage-defining mutations |
| Nextstrain Augur (Huddleston 2021 *JOSS* 6:2906) | Python CLI for subsampling + alignment + tree + ancestral-trait + time-tree | Genomes + metadata + sampling config | Auspice JSON for visualization | End-to-end pipeline for curated surveillance builds | Subsampling configuration drives results more than data; nextstrain.org subsamples ~3000-5000 of millions |
| UShER + matUtils + matOptimize + RIPPLES (Turakhia 2021 *Nat Genet* 53:809) | Parsimony placement on daily MAT; SPR refinement; recombination detection | New consensus + existing MAT | Updated MAT, subtrees, recombinant calls | Pandemic-scale (millions of genomes) | Parsimony branch lengths systematically shorter than ML; re-estimate before downstream R_e |
| Freyja (Karthikeyan 2022 *Nature* 609:101) | Depth-weighted LAD regression on barcode-matrix mutation frequencies | Wastewater BAM + barcode | Per-lineage abundance | Recovers expected abundances down to ~5%; quantitative | Lineages absent from barcode invisible; barcode is forward-only -- cannot deconvolve lineages designated after barcode date |
| COJAC (Jahn 2022 *Nat Microbiol* 7:1151) | Co-occurrence of signature mutations on the same read pair | Wastewater BAM | Per-lineage presence / absence | More robust than per-site frequencies; detected Alpha 13 days before clinical | Single-read amplicons (no co-occurrence) cannot resolve; requires paired-end or long-read |
| alcov | Lineage deconvolution similar paradigm to Freyja | Wastewater BAM | Per-lineage abundance | Alternative to Freyja | Less benchmarked |
| lineagespot (Pechlivanis 2022 *Sci Rep* 12:2659) | R/Bioconductor lineage deconvolution from VCF + signature mutations | VCF + reference lineage mutations | Per-lineage abundance | R / Bioconductor integration | Less ML-driven; smaller community |
| Wenseleers / Bedford-Figgins multinomial logistic (Abousamra, Figgins, Bedford 2024 *PLoS Comput Biol* 20:e1012443) | Multinomial logistic regression on lineage frequencies over time | Lineage frequencies + dates | Growth advantage per lineage with 95% CI | Standard for outbreak.info / cov-lineages.org | Marginal CI for one lineage hides covariance with all others; early estimates inflated |
| 3SEQ (Boni 2007 *Genetics* 176:1035) | Triplet-based recombination detection | Aligned sequences | Recombinant candidates | General-purpose | High false-positive rate at low divergence |
| RDP4 / RDP5 (Martin 2015 *Virus Evol* 1:vev003) | Multiple-method recombination detection | Aligned sequences | Recombinant candidates | Multi-method consensus | Slow; parameter-sensitive |
| Bolotie (Varabyou 2021 *Bioinformatics* 37:2298) | SARS-CoV-2-specific recombination detection | SARS-CoV-2 consensus | Recombinant candidates | Tuned for SARS-CoV-2 sub-lineage divergence | Specialist tool |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| Assign lineage to a SARS-CoV-2 consensus | Pangolin with `--analysis-mode usher` (UShER default since v4) + Nextclade cross-check; pin pangolin-data and Nextclade dataset versions | pangoLEARN alone (deprecated); not pinning version (cross-lab calls diverge) |
| Detect emerging variants in wastewater | COJAC for early detection (co-occurrence on amplicon) + Freyja for quantitative tracking; pin Freyja barcode version | Naive site-frequency aggregation; comparing across barcode versions |
| Track lineage frequencies over time | Multinomial logistic regression (Wenseleers / Bedford-Figgins) OR Bayesian renewal equation; report covariance among lineages | Plotting raw counts without CI; reporting single-lineage growth advantage without covariance |
| Build a regional surveillance phylogeny | Nextstrain Augur pipeline; subsample to manageable size; TreeTime for dates; document subsampling | BEAST on raw 10k+ samples (intractable); not documenting subsampling |
| Compare wastewater results across labs | Same primer scheme + same Freyja barcode + same Pangolin / Nextclade version | Mixing primer schemes; mixing barcode versions; comparing across pangolin-data versions |
| QC a new SARS-CoV-2 genome | Nextclade (alignment QC; mutation outliers; recombination indicators) | Pangolin alone (passes confidently on bad genomes) |
| Detect recombinant lineages | Trust Pango-designation X-prefix assignments; for novel candidates use RDP5 / 3SEQ / Bolotie + manual review | Manual eyeballing of mutation patterns; ignoring designation lag |
| Phylogenetic context for outbreak | UShER + matUtils subtree extraction; re-estimate branch lengths via TreeTime for downstream R_e | Re-treeing from scratch every time |
| Estimate vaccine-escape risk | Lab assays (neutralisation, escape mutants) + structural prediction; genomic surveillance flags candidates | Pure genomic prediction without lab validation |
| Wastewater-to-cases conversion | Variant-specific shedding rate (Omicron BA.1 shed less per case than Delta); pin barcode + report uncertainty | Fixed RNA-to-cases ratio across variants is wrong; variant-specific shedding has been documented in the wastewater literature |

Methodology evolves; before any high-stakes lineage report, verify Pangolin's current default analysis-mode and Nextclade's bundled dataset against pango-designation issues for any emerging lineage.

## Pangolin Lineage Assignment With Version Pinning

**Goal:** Assign Pango lineages to SARS-CoV-2 consensus sequences using UShER mode (the default since v4; pangoLEARN deprecated mid-2023), with full pangolin-data version provenance preserved for reproducibility.

**Approach:** Always pass `--analysis-mode usher`; record `pangolin --all-versions` output alongside every lineage call; for published or regulatory output, pin pangolin-data to a specific release tag and re-run the whole archive whenever the version is updated.

```bash
pangolin sequences.fasta --analysis-mode usher --outfile lineage_report.csv
pangolin --all-versions > pangolin_versions.txt
```

`pangolin --all-versions` prints: pangolin executable version, pangolin-data version (weekly updated; mandatory pin for reproducibility), scorpio version, and constellations version. All four are version-sensitive; in published surveillance reports, pin all four.

## Nextclade With Dataset Pinning

**Goal:** Assign Nextstrain clade, Pango lineage, mutations, and QC flags to SARS-CoV-2 consensus sequences with explicit dataset version provenance.

**Approach:** Fetch the current dataset with `nextclade dataset get --name sars-cov-2 --output-dir nc_dataset/sars-cov-2`; record the `pathogen.json` tag / commit hash; run `nextclade run --input-dataset` on the pre-downloaded folder so the dataset version is locked in for the analysis.

```bash
nextclade dataset get --name sars-cov-2 --output-dir nc_dataset/sars-cov-2
NC_DATASET_TAG=$(jq -r '.tag // .version' nc_dataset/sars-cov-2/pathogen.json)

nextclade run \
    --input-dataset nc_dataset/sars-cov-2 \
    --output-tsv nextclade.tsv \
    --output-json nextclade.json \
    sequences.fasta

echo "nextclade_dataset_tag: ${NC_DATASET_TAG}" > nextclade.metadata
```

Different dataset versions assign different mutations as "lineage-defining" because internal-node placement can shift as the tree grows. Cross-version comparison of mutation reports is therefore method-dependent.

## Wastewater Lineage Deconvolution With Freyja

**Goal:** Estimate per-lineage abundance in a wastewater pooled sample with explicit handling of the barcode forward-only date constraint, primer-scheme awareness, and residual mass interpretation.

**Approach:** Confirm barcode date postdates sample collection; if not, `freyja barcode-build` from the current UShER tree; variant call with `freyja variants` then deconvolve with `freyja demix`; inspect the `resid` column (residual mass NOT assigned to known lineages; high resid indicates a novel lineage is invisible); apply primer-scheme-aware coverage masking; report variant-specific uncertainty.

```bash
freyja update

freyja variants \
    sample.bam \
    --variants sample.variants.tsv \
    --depths sample.depths.tsv \
    --ref reference.fa

freyja demix \
    sample.variants.tsv \
    sample.depths.tsv \
    --output sample.demix.tsv
```

The Freyja `--barcodes` (or default bundled) date MUST postdate the sample collection date. Lineages designated after the barcode date cannot be detected -- the demixing silently fails and presents as elevated abundance of the closest parent lineage. For samples potentially containing emerging lineages, regenerate barcodes:

```bash
freyja barcode-build \
    --pb-and-meta usher_tree.pb \
    --output-dir custom_barcodes
```

Subsequent methodological extensions to Karthikeyan have appeared in the wastewater literature, and recent benchmarks comparing the major deconvolution tools (Freyja, COJAC, alcov, lineagespot, LCS) confirm Freyja and COJAC consistently perform well, with performance degrading at low coverage and for divergent lineages.

## COJAC for Co-Occurrence Detection

**Goal:** Detect emerging variants in wastewater earlier than per-site frequency methods by requiring co-occurrence of two signature mutations on the same amplicon (read pair).

**Approach:** COJAC checks read pairs for joint occurrence of variant-defining mutations; the inferential leap is more robust because a single site can have shared mutations across lineages, but two signature mutations on the same read pair strongly imply a single lineage. Detected Alpha 13 days before clinical samples in Swiss data (Jahn 2022 *Nat Microbiol* 7:1151).

```bash
cojac cooc-mutbamscan \
    -a primer_scheme.bed \
    -m variants_definitions.yaml \
    -b sample.bam \
    -o sample.cooc.tsv
```

## Nextstrain Augur Pipeline

**Goal:** Build a curated regional surveillance phylogeny with subsampling, alignment, tree, ancestral-trait inference, and time-scaling -- in Auspice-visualisable format. The Nextstrain platform was introduced by Hadfield 2018 *Bioinformatics* 34:4121; Augur is the Python CLI (Huddleston 2021 *JOSS* 6:2906).

**Approach:** Pull the latest official pathogen build from github.com/nextstrain/<pathogen>; subsample to manageable size (typically 3000-5000 genomes per global build; smaller regional); document subsampling configuration explicitly (it drives the result more than the underlying data per Hodcroft 2021).

```bash
augur align --sequences seqs.fasta --reference-sequence ref.gb --output aligned.fasta
augur tree --alignment aligned.fasta --output tree.nwk
augur refine \
    --tree tree.nwk \
    --alignment aligned.fasta \
    --metadata meta.tsv \
    --output-tree refined.nwk \
    --output-node-data branch_lengths.json \
    --timetree \
    --root oldest \
    --coalescent skyline
augur ancestral --tree refined.nwk --alignment aligned.fasta --output-node-data nt_muts.json
augur translate --tree refined.nwk --ancestral-sequences nt_muts.json --reference-sequence ref.gb --output-node-data aa_muts.json
augur traits --tree refined.nwk --metadata meta.tsv --columns country region --output-node-data traits.json
augur export v2 \
    --tree refined.nwk \
    --metadata meta.tsv \
    --node-data branch_lengths.json nt_muts.json aa_muts.json traits.json \
    --output auspice.json
```

Hodcroft 2021 *Nature* 591:30 documented that Nextstrain subsampling configurations drive lineage-time estimates more than the underlying data. Two researchers using the official pipeline with different subsampling can get different MRCA dates and migration patterns from the same raw genomes.

## Per-Method Failure Modes

### pangolin-data version skew between labs

**Trigger:** Two labs submit the same consensus genome to Pangolin with different pangolin-data versions; the lineage call differs.

**Mechanism:** Lineage designation happens through pango-designation GitHub issues -- days-to-weeks before pangolin-data releases include the lineage. During the lag, the same genome is callable as the parent (older pangolin-data) or the child (current). pangolin-data is updated weekly.

**Symptom:** Cross-lab lineage prevalence comparisons over time show implausible jumps coinciding with pangolin-data release dates rather than biology.

**Fix:** Pin pangolin-data version explicitly with `pangolin --all-versions` recorded alongside every call. For published or regulatory output, re-run the WHOLE archive against a single pangolin-data version before reporting.

### Freyja barcode predates the sample collection date

**Trigger:** Wastewater sample collected after a new lineage was designated; Freyja barcode built before that designation.

**Mechanism:** Freyja barcodes are built from the UShER tree at a specific date; lineages designated AFTER the barcode date cannot be detected. The demixing silently fails -- the new lineage's signal is misassigned to its closest parent.

**Symptom:** Wastewater sample shows implausibly high abundance of a single parent lineage; new lineage that should be present is reported as 0%.

**Fix:** Run `freyja update` regularly; for samples potentially containing emerging lineages, regenerate barcodes with `freyja barcode-build` from the current UShER tree. Report `resid` (residual mass not assigned to known lineages); high resid indicates a novel lineage is being missed.

### ARTIC primer dropout misread as deletion

**Trigger:** SARS-CoV-2 surveillance using ARTIC V4.1 amplicons; new variant has mutation at primer site; amplicons 64 / 76 / 88-90 silently drop out.

**Mechanism:** When a primer fails to bind, the amplicon doesn't amplify; consensus calling produces N's or reference-derived calls in that region. This LOOKS LIKE a deletion in downstream analysis but is actually missing data. Itokawa 2020 *PLoS ONE* 15:e0239403 documented primer interactions specifically.

**Symptom:** "Deletion" calls cluster in known dropout amplicons; Pangolin / Nextclade lineage call shifts when masked positions are filled with reference.

**Fix:** Inspect per-amplicon coverage with `samtools depth -aa`; mask consensus positions in dropped amplicons (use Ns -- Pangolin and Nextclade handle Ns gracefully). Document primer scheme version (V3 / V4 / V4.1 / V5.3.2 / Midnight) per isolate.

### Recombinant assigned to one parent lineage

**Trigger:** A SARS-CoV-2 recombinant (e.g., XEC = KS.1.1 x KP.3.3) emerges; pango-designation has not yet issued the X-prefix designation; Pangolin assigns to one of the parents.

**Mechanism:** Pangolin in either mode assigns a recombinant to one parent lineage if no Pango-X designation exists yet. Identifying recombinants requires breakpoint detection (3SEQ, Bolotie, RDP4) and manual designation through pango-designation; the designation can lag emergence by weeks-to-months for novel recombinants.

**Symptom:** Outbreak interpretation conflates a recombinant lineage with its parent; transmissibility / immune-escape claims are wrong.

**Fix:** For any candidate emerging lineage with unusual mutations, run Bolotie or 3SEQ for recombination detection; cross-check Pangolin vs Nextclade lineage call; submit candidate recombinants to cov-lineages issue tracker if novel.

### pangoLEARN result reported as authoritative

**Trigger:** Pangolin run with `--analysis-mode pangolearn` (or via legacy Docker image that defaults to pangoLEARN); user reports the call.

**Mechanism:** Pongmoragot 2024 *Virus Evol* 10:vead085 demonstrated UShER mode is significantly more accurate for recent / divergent lineages. pangoLEARN was officially deprecated mid-2023.

**Symptom:** Cross-lab comparison reveals one lab using pangoLEARN (legacy) and another using UShER; calls differ at borderline lineages.

**Fix:** Switch to `--analysis-mode usher` (default since v4). For longitudinal datasets crossing the mid-2023 mode-switch, re-run the historical archive against UShER mode.

### Freyja barcode-build vs barcode_build flag

**Trigger:** Script written from older Freyja documentation using `freyja barcode_build` (underscore).

**Mechanism:** Current Freyja versions use `barcode-build` (hyphen); the underscore form may not be recognised.

**Symptom:** Subprocess fails with "unrecognized command".

**Fix:** Use `freyja barcode-build` (hyphen). Verify with `freyja --help`.

### Wenseleers / Bedford lineage-growth CI hides covariance

**Trigger:** Reporting a single lineage's growth advantage 95% CI from a multinomial logistic regression.

**Mechanism:** The CI for any single lineage is conditional on all other lineages being held at their estimated growth rates; the marginal CI hides covariance among lineages. Early growth-advantage estimates are systematically too large; they shrink as more time passes (alternative explanations become identifiable).

**Symptom:** Initial published growth advantage > later refined estimate; "outlier-fast" lineages later moderated.

**Fix:** Report the full multinomial covariance matrix or at minimum the rank-ordered growth advantages with simultaneous CIs. Cite Abousamra 2024 *PLoS Comput Biol* 20:e1012443.

### Nextstrain subsampling drives the result

**Trigger:** Nextstrain Augur build with default subsampling at 3000-5000 genomes from millions; user interprets the tree topology as authoritative.

**Mechanism:** Hodcroft 2021 *Nature* 591:30 commented that subsampling decisions drive lineage-time estimates more than the underlying data. Two researchers using the official Nextstrain pipeline with different subsampling configurations get different MRCA dates and migration patterns from the same raw genomes.

**Symptom:** Published Nextstrain tree differs from another analysis on the same raw data; conclusions sensitive to subsampling.

**Fix:** Document subsampling configuration explicitly in any Nextstrain build; run sensitivity analysis with alternative subsampling; treat MRCA dates and migration calls with appropriate uncertainty.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Pangolin "BA.2.86", Nextclade clade "23I" | Equivalent at different resolutions -- BA.2.86 is within 23I | Report both; Pango lineage for sub-clade resolution |
| Pangolin "BA.5.2", Nextclade "Unassigned" | Nextclade dataset older than pangolin-data; OR Nextclade QC failed | Update Nextclade dataset; re-run; inspect QC fields |
| Pangolin UShER and pangoLEARN disagree | pangoLEARN is the deprecated decision-tree classifier | Trust UShER call |
| Freyja shows 0% of expected lineage | Lineage absent from current barcode (post-dates barcode) | Rebuild barcodes (`freyja barcode-build`); confirm lineage is in the UShER tree the barcode is built from |
| Freyja confidence < 0.7 on dominant lineage | Sub-100x coverage OR amplicon dropout | Inspect per-amplicon coverage; consider re-sequencing; report as indeterminate |
| Nextclade and Pangolin disagree on recombinant | Recombinants inherently ambiguous; depends on which parent's SNPs dominate | Report as recombinant candidate; submit to cov-lineages if novel |
| Two consecutive pangolin-data releases call the same consensus differently | Lineage definitions revised between releases | Pin pangolin-data; record version + date alongside lineage |
| COJAC detects a variant Freyja does not | COJAC's co-occurrence requirement more sensitive at low abundance | Trust COJAC for early detection; Freyja for quantitative tracking |
| Wastewater Freyja result conflicts with clinical lineage prevalence | Barcode staleness; primer dropout in wastewater; faecal shedding rate varies by variant | Update barcode; check per-amplicon coverage; flag variant-specific shedding |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| Pangolin min coverage for lineage call | >=50% genome coverage (~14kb) | Pangolin convention |
| Nextclade QC stop-codon threshold | Per dataset; check `pathogen.json` | Nextclade dataset-specific |
| Freyja minimum coverage per site | >=10x typical | Freyja convention; per-site weighting accounts for variance |
| Freyja `resid` flag threshold | Project-specific; >0.1 typically indicates novel lineage missed | Freyja documentation |
| COJAC early-detection lead time vs clinical | Up to 13 days in Swiss data | Jahn 2022 *Nat Microbiol* 7:1151 |
| Karthikeyan 2022 wastewater Omicron lead | 11 days before clinical detection (San Diego) | Most-favourable configuration; subsequent retrospective analyses produced detection lags ranging from -5 to +3 days |
| ARTIC V4.1 known chronic dropouts | Amplicons 64, 76, 88-90 | Itokawa 2020 / community documentation |
| Augur subsampling typical | 3000-5000 genomes per global build | Nextstrain convention; document explicitly |
| Multinomial logistic growth-advantage early estimate inflation | Systematically too large; shrinks over time | Abousamra 2024 *PLoS Comput Biol* 20:e1012443 |
| GISAID 2024-2025 weekly submission rate | ~5,000-20,000/week (down from ~500,000/week peak early 2022) | Community-documented; emerging-lineage detection lag increased |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `pangolin --inference usher` rejected | Flag is `--analysis-mode usher` | `--analysis-mode usher` |
| `nextclade run --input-dataset DIR` rejected on v2 | v2 used `--input-dataset` differently | Verify `nextclade --version`; v3+ accepts pre-downloaded dataset folder |
| `freyja barcode_build` rejected | Current is `barcode-build` (hyphen) | Use hyphen form |
| Pangolin output column not present | Column names changed between major releases | Introspect output schema; `pangolin --all-versions` |
| Freyja silently misassigns new lineage | Barcode predates lineage designation | Rebuild barcodes; check `resid` |
| Nextclade and Pangolin disagree | Different versions; recombinant; QC | Update both; reconcile per table |
| ARTIC consensus has Ns clustered in one region | Primer dropout in that amplicon | Mask the amplicon; document scheme version |
| Lineage frequency shows implausible jump | pangolin-data version drift | Pin version; re-run archive |
| Augur tree topology changes between runs | Subsampling randomness | Pin random seed; document subsampling |
| Augur `refine` requires `--root` | Shallow tree without explicit root strategy | `--root best` / `oldest` / `residual` |
| Wastewater Freyja result differs from clinical | Barcode staleness or primer-scheme mismatch | Update barcode; document scheme; check coverage |
| COJAC misses a known variant | Single-read amplicons (no co-occurrence) | Re-sequence with paired-end or long-read |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Pangolin version?" | `pangolin --all-versions` recorded; pinned for the analysis; archive re-run on dataset update |
| "Nextclade dataset version?" | Dataset tag from `pathogen.json` recorded; pre-downloaded folder used to lock the version |
| "Why UShER not pangoLEARN?" | pangoLEARN deprecated mid-2023 (Pongmoragot 2024); UShER default since v4 |
| "How were ARTIC dropouts handled?" | Per-amplicon coverage checked; failed amplicons masked; primer scheme documented per isolate |
| "Were recombinants checked for?" | Bolotie / 3SEQ run on candidates; cross-checked Pangolin vs Nextclade; submitted to cov-lineages for novel candidates |
| "Wastewater barcode date?" | Barcode date postdates sample collection; `freyja barcode-build` from current UShER tree if needed |
| "Wastewater-to-cases conversion?" | Variant-specific shedding rate flagged in the wastewater literature; not assumed constant |
| "Lineage growth-advantage CI?" | Multinomial covariance reported; early estimates noted as inflated (Abousamra 2024) |
| "Nextstrain subsampling?" | Configuration explicit; sensitivity analysis run; MRCA / migration treated with uncertainty (Hodcroft 2021) |

## References

- O'Toole Á, Scher E, Underwood A et al (2021) Assignment of epidemiological lineages in an emerging pandemic using the pangolin tool. *Virus Evol* 7(2):veab064. doi:10.1093/ve/veab064
- Aksamentov I, Roemer C, Hodcroft EB, Neher RA (2021) Nextclade: clade assignment, mutation calling and quality control for viral genomes. *J Open Source Softw* 6(67):3773. doi:10.21105/joss.03773
- Karthikeyan S, Levy JI, De Hoff P et al (2022) Wastewater sequencing reveals early cryptic SARS-CoV-2 variant transmission. *Nature* 609(7925):101-108. doi:10.1038/s41586-022-05049-6
- Hadfield J, Megill C, Bell SM et al (2018) Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics* 34(23):4121-4123. doi:10.1093/bioinformatics/bty407
- Huddleston J, Hadfield J, Sibley TR et al (2021) Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens. *J Open Source Softw* 6(57):2906. doi:10.21105/joss.02906
- Turakhia Y, Thornlow B, Hinrichs AS et al (2021) Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic. *Nat Genet* 53(6):809-816. doi:10.1038/s41588-021-00862-7
- Pongmoragot J, Pearson C, Borg ML et al (2024) Comparison of UShER-based and pangoLEARN-based Pangolin lineage assignments for SARS-CoV-2 sequences. *Virus Evol* 10(1):vead085. doi:10.1093/ve/vead085
- Jahn K, Dreifuss D, Topolsky I et al (2022) Early detection and surveillance of SARS-CoV-2 genomic variants in wastewater using COJAC. *Nat Microbiol* 7(8):1151-1160. doi:10.1038/s41564-022-01185-x
- Pechlivanis N, Tsagiopoulou M, Maniou MC et al (2022) Detecting SARS-CoV-2 lineages and mutational load in municipal wastewater and a use-case in the metropolitan area of Thessaloniki, Greece. *Sci Rep* 12:2659. doi:10.1038/s41598-022-06625-6
- Itokawa K, Sekizuka T, Hashino M, Tanaka R, Kuroda M (2020) Disentangling primer interactions improves SARS-CoV-2 genome sequencing by multiplex tiling PCR. *PLoS ONE* 15(9):e0239403. doi:10.1371/journal.pone.0239403
- Hodcroft EB, De Maio N, Lanfear R et al (2021) Want to track pandemic variants faster? Fix the bioinformatics bottleneck. *Nature* 591(7848):30-33. doi:10.1038/d41586-021-00525-x
- Abousamra E, Figgins M, Bedford T (2024) Fitness models provide accurate short-term forecasts of SARS-CoV-2 variant frequency. *PLoS Comput Biol* 20(9):e1012443. doi:10.1371/journal.pcbi.1012443
- Boni MF, Posada D, Feldman MW (2007) An exact nonparametric method for inferring mosaic structure in sequence triplets. *Genetics* 176(2):1035-1047. doi:10.1534/genetics.106.068874
- Martin DP, Murrell B, Golden M, Khoosal A, Muhire B (2015) RDP4: detection and analysis of recombination patterns in virus genomes. *Virus Evol* 1(1):vev003. doi:10.1093/ve/vev003
- Varabyou A, Pockrandt C, Salzberg SL, Pertea M (2021) Rapid detection of inter-clade recombination in SARS-CoV-2 with Bolotie. *Bioinformatics* 37(15):2298-2300. doi:10.1093/bioinformatics/btab080

## Related Skills

- pathogen-typing - Lineage assignment overlaps with typing; this skill owns longitudinal frequency tracking and wastewater deconvolution
- phylodynamics - Lineage-stratified BDSKY / BICEPS R_e estimation runs downstream of lineage assignment
- transmission-inference - SARS-CoV-2 cluster definition combines lineage + 0-2 SNPs + epi link
- amr-surveillance - Antiviral drug-resistance mutation tracking is the variant-surveillance analogue for AMR
- phylogenetics/modern-tree-inference - IQ-TREE / RAxML for non-UShER topology
- phylogenetics/tree-io - Tree parsing and format conversion for Augur output
- comparative-genomics/whole-genome-alignment - Reference-based alignment for SNP calling
- variant-calling/vcf-basics - VCF for lineage-defining mutations
- variant-calling/variant-calling - Variant calling for wastewater (lofreq, ivar)
- variant-calling/filtering-best-practices - Per-amplicon coverage filtering for ARTIC
- read-alignment/bwa-alignment - Read mapping upstream
- read-alignment/minimap2-alignment - Long-read alignment for ARTIC-Midnight 1200
- read-qc/quality-reports - Sequencing QC upstream
- database-access/sra-data - SRA / INSDC retrieval; GISAID is a separate restricted-access source
- data-visualization/multipanel-figures - Lineage frequency / wastewater plotting
- workflows/somatic-variant-pipeline - End-to-end orchestration patterns
