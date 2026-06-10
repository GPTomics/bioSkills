---
name: bio-microbiome-qiime2-workflow
description: Operates the QIIME2 framework as the glue for an amplicon analysis - the .qza/.qzv artifact model, semantic types (FeatureTable[Frequency], SampleData[PairedEndSequencesWithQuality], Phylogeny[Rooted], FeatureData[Taxonomy]), embedded provenance plus provenance replay, import (Casava/manifest/EMP/BIOM), export, the Metadata object, and the q2cli vs Artifact API interfaces. Covers why a .qza is data-plus-executable-history not a file, why export drops provenance, why a .qzv is terminal, why classifier .qza are version-pinned, and the 2026 distribution/rachis rename. Use when importing reads, choosing a manifest/Casava/EMP/BIOM path, reading or replaying provenance, exporting to BIOM/phyloseq, fixing semantic-type or Phred or sklearn-version errors, or orchestrating the pipeline. Denoising -> amplicon-processing; classifier/DB -> taxonomy-assignment; diversity metric/depth -> diversity-analysis; DA tool -> differential-abundance; PICRUSt2 -> functional-prediction; shotgun moshpit -> metagenomics.
tool_type: cli
primary_tool: QIIME2
---

## Version Compatibility

Reference examples tested with: QIIME2 2026.1+ (amplicon distribution; framework now `rachis`), provenance-lib 2024.10+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `qiime --version`, `qiime info`, then `qiime <plugin> <action> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The QIIME2 release tag (calendar-versioned `YYYY.RELEASE`, e.g. `2026.1`) defines the plugin API AND the `.qza` artifact format; an sklearn taxonomy classifier `.qza` trained under one release may need retraining under another (the classifier is pinned to its scikit-learn version). The conda env name encodes the release and distribution (`qiime2-amplicon-2026.1`; renamed toward `rachis-qiime2-<release>` in 2026.4). Names and the install YAML URL are moving targets - verify the current release and distribution names against the live install page before pinning anything.

# QIIME2 Amplicon Workflow

**"Run my amplicon study through QIIME2"** -> Move data through the framework as typed, provenance-carrying artifacts and defer every scientific choice to the owning skill - because a `.qza` is not a file, it is data plus its entire executable history, and that history is the deliverable.
- CLI: `qiime tools import`, `qiime <plugin> <action> --i-* --p-* --m-* --o-*`, `qiime tools peek/export`

Scope: the artifact/provenance/type machinery and the import/export/metadata/interface mechanics - this is the GLUE skill. Denoising params (trunc/trim/maxEE, DADA2 vs Deblur) -> amplicon-processing. Classifier/DB choice + training -> taxonomy-assignment. Diversity metric/sampling-depth/rarefaction + PERMANOVA-vs-dispersion -> diversity-analysis. DA tool choice/consensus -> differential-abundance. PICRUSt2 -> functional-prediction. Shotgun reads (moshpit distribution, Kraken2/MetaPhlAn/HUMAnN) -> metagenomics. This skill shows each scientific action and routes the decision out; it does not re-teach the method.

## The Single Most Important Modern Insight -- A .qza Is Data Plus Its Executable History, Not a File Format

A `.qza` carries the data AND the complete computational graph that produced it. The semantic-type system plus the embedded provenance ARE the reproducibility guarantee - the whole reason to work inside the framework instead of passing loose BIOM/FASTA/Newick files. The cost is exact and unavoidable: there is no `cat`-ing the data. Three corollaries each common misuse violates:

1. **Working THROUGH the framework keeps the chain; exporting breaks it.** `qiime tools export` writes native data and silently drops the QIIME2 wrapper AND the provenance. Export early and go ad-hoc, and the final figure has no history back to the raw reads - the framework overhead was paid and the deliverable thrown away. Export at the LAST step, or use `qiime2R::qza_to_phyloseq` so the chain survives as far as possible.
2. **Semantic types are a type system for biology.** `core-metrics-phylogenetic` refuses a `FeatureData[Taxonomy]` where a `FeatureTable[Frequency]` belongs, BEFORE running. A type error is the guard WORKING - fix the upstream action that made the wrong type, do not launder it by re-importing.
3. **A `.qzv` is terminal and a classifier is version-pinned.** A Visualizer's output can never be another action's input (keep the `.qza` it was made from). An sklearn classifier `.qza` trained under 2024.x raises a version-mismatch under 2026.x - the training version is part of the method.

Organize the analysis around protecting the provenance chain and the type contract, not around listing flags.

## What Is Inside a .qza

A `.qza` (QIIME Zipped Artifact) and `.qzv` (Visualization) are ZIP archives keyed at top level by a UUID. Every artifact carries four things:

1. **UUID** - identifies THIS computation (provenance references inputs/outputs by UUID), not just a file.
2. **Semantic TYPE** - what the data MEANS: `FeatureTable[Frequency]`, `SampleData[PairedEndSequencesWithQuality]`, `Phylogeny[Rooted]`, `FeatureData[Taxonomy]`, `FeatureData[Sequence]`, `DistanceMatrix`, `SampleData[AlphaDiversity]`. Types can carry Properties (`SampleData[AlphaDiversity] % Properties('phylogenetic')`).
3. **FORMAT** - the on-disk layout the bytes live in (e.g. `BIOMV210DirFmt`, a Newick file). Type is the meaning; format is the bytes.
4. **PROVENANCE** - in a `provenance/` subtree: for every upstream action the plugin/action name, every parameter value, input/output UUIDs, plugin + framework versions, execution environment, timestamp, and BibTeX citations. The references form a DAG of the whole analysis.

```bash
qiime tools peek table.qza                    # UUID + Type + Format, without unzipping
qiime tools validate table.qza --level max    # archive integrity + payload conforms to its format
qiime tools extract --input-path table.qza --output-path extracted/   # FULL archive incl provenance (read by hand)
qiime tools export  --input-path table.qza --output-path exported/     # ONLY the native data - DROPS provenance
```

`qiime tools extract` keeps the QIIME2 structure (data + `provenance/`); `qiime tools export` is the one-way door out. A plain `unzip table.qza` works too (it is a standard ZIP).

## Tool / Interface Taxonomy

| Interface / tool | Role | When |
|------------------|------|------|
| q2cli (`qiime ...`) | the command-line interface; `--i-*` inputs, `--p-*` params, `--m-*` metadata, `--o-*`/`--output-dir` outputs | default, most-documented, scriptable; what tutorials/forum answers use |
| Artifact API (`from qiime2 import Artifact, Metadata`) | the Python 3 interface; `Artifact.load`/`.save`/`.view`, actions importable as functions returning `Results` | notebooks, embedding QIIME2 in a larger Python pipeline (no temp files) |
| view.qiime2.org | renders any `.qzv` viz AND the `.qza`/`.qzv` provenance DAG client-side, NO install | sharing results and inspecting provenance without QIIME2 installed |
| provenance-lib (`qiime tools replay-provenance`) | parses an artifact's provenance DAG and regenerates executable code (Keefe 2023) | recovering the commands that made an artifact; reproducing a shared `.qza` |

Neither q2cli nor the Artifact API is "more reproducible" - provenance is identical; pick by host environment. An Action is a Method (Artifacts in -> Artifacts out), a Visualizer (-> exactly one terminal `.qzv`), or a Pipeline (-> many Artifacts and/or Visualizations, e.g. `core-metrics-phylogenetic`). The Method/Visualizer distinction is WHY a `.qzv` is a dead end. The legacy `q2studio` desktop GUI is dead (last release 2022.8); the no-CLI answers are Galaxy + view.qiime2.org.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Demultiplexed per-sample FASTQ, filenames are Casava 1.8 | `--type SampleData[PairedEndSequencesWithQuality] --input-format CasavaOneEightSingleLanePerSampleDirFmt` | sample IDs parsed from filenames; no manifest needed |
| Demultiplexed FASTQ, arbitrary paths | V2 manifest (`PairedEndFastqManifestPhred33V2`) | TSV of absolute paths; the most general/explicit on-ramp |
| Still multiplexed (one big FASTQ + barcodes) | import `EMPPairedEndSequences`, then `qiime demux emp-paired` | demultiplexing is a QIIME2 step, not the import |
| A feature table built elsewhere | `--input-format BIOMV210Format --type FeatureTable[Frequency]` | BIOM v2.1 (HDF5); attach metadata separately |
| Scripting a notebook / larger Python pipeline | Artifact API | returns Artifacts directly, no temp files; same provenance |
| Need to share a result with no-QIIME2 collaborators | upload `.qzv` to view.qiime2.org | renders viz + provenance client-side |
| Handed a single `.qza`, need the commands that made it | `qiime tools replay-provenance` | regenerates executable code from the provenance DAG |
| One-off custom R analysis, fighting the framework | `qiime2R::qza_to_phyloseq` / export, own the provenance loss | the overhead is not worth it; be honest about the exit point |
| Shotgun / WGS reads | -> metagenomics (moshpit distribution) | different distribution and toolchain; cross-link, do not merge |

## Import

**Goal:** Turn raw demultiplexed reads into a typed, provenance-rooted artifact with the correct Phred offset.

**Approach:** Write a V2 manifest (TSV, absolute paths), declare the semantic type and the format whose name encodes the Phred offset, then immediately summarize to confirm the reads decoded sanely.

```bash
# manifest.tsv (TAB-separated, V2; absolute paths):
#   sample-id<TAB>forward-absolute-filepath<TAB>reverse-absolute-filepath
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.tsv \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path demux.qza
# Phred offset is BAKED INTO the format name: Phred33V2 (modern Illumina) vs Phred64V2 (legacy).
# V1 was CSV with a `direction` column; V2 is TSV with separate forward/reverse columns - prefer V2.

qiime demux summarize --i-data demux.qza --o-visualization demux.qzv   # per-base quality (drives trunc choices)
```

For EMP-multiplexed data: import `--type 'EMPPairedEndSequences'`, then `qiime demux emp-paired --i-seqs emp.qza --m-barcodes-file metadata.tsv --m-barcodes-column barcode-sequence --o-per-sample-sequences demux.qza --o-error-correction-details ec.qza`. The per-base quality plot in `demux.qzv` is read by amplicon-processing to pick truncation - not here.

## The Orchestration Skeleton (each science step DEFERS)

The pipeline shape, with every method choice routed to its owning skill:

```bash
# Denoise -> ASV table + rep-seqs.  PARAM CHOICE (trunc/trim/maxEE, DADA2 vs Deblur) -> amplicon-processing
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza \
    --p-trunc-len-f 0 --p-trunc-len-r 0 \
    --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats stats.qza

qiime tools peek table.qza    # confirm Type is FeatureTable[Frequency] before wiring downstream

# Taxonomy.  CLASSIFIER + DB choice and training -> taxonomy-assignment
# Use a classifier .qza trained for THIS release (data.qiime2.org/<release>/common/...); old ones break.
qiime feature-classifier classify-sklearn \
    --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

# Phylogeny (Pipeline) -> rooted tree for UniFrac/Faith PD
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza \
    --o-alignment aln.qza --o-masked-alignment masked-aln.qza \
    --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

# Diversity (Pipeline).  SAMPLING DEPTH + metric + rarefy-or-not -> diversity-analysis (pick depth from alpha-rarefaction)
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza \
    --p-sampling-depth 10000 --m-metadata-file metadata.tsv --output-dir core-metrics/
# PERMANOVA via diversity beta-group-significance; the location-vs-dispersion (betadisper) confound -> diversity-analysis

# Differential abundance - MODERN q2-composition (NOT add-pseudocount+ancom).  Tool choice/consensus -> differential-abundance
qiime composition ancombc --i-table table.qza --m-metadata-file metadata.tsv \
    --p-formula 'group' --o-differentials ancombc.qza
qiime composition da-barplot --i-data ancombc.qza --o-visualization ancombc-barplot.qzv
```

`core-metrics-phylogenetic` and `align-to-tree-mafft-fasttree` are Pipelines (one call, a directory of artifacts + Emperor `.qzv`s out). `--p-formula` takes column names from the Metadata; annotate integer ID/batch columns `categorical` (below) or they enter the model as continuous covariates.

## Metadata

The Metadata TSV is the spine - the same `--m-metadata-file` drives demux barcodes, group-significance, taxa barplots, ANCOM-BC grouping, and Emperor coloring. First column header is the ID column (`sample-id`, `id`, `#SampleID`, ...). An optional second row `#q2:types` overrides type inference per column (`categorical` / `numeric`):

```
sample-id	subject	group
#q2:types	categorical	categorical
s1	101	treatment
s2	102	control
```

Without the `#q2:types` row, a column of only integers is inferred **numeric** - so a subject/batch/timepoint ID silently becomes a continuous covariate. Annotate ID-like integer columns `categorical`. Validate the sheet with Keemei (Rideout 2016 *GigaScience* 5:27) before running - a malformed metadata file is a top cause of cryptic action failures. `qiime metadata tabulate --m-input-file metadata.tsv --o-visualization metadata.qzv` renders any metadata (including an artifact viewed as metadata, e.g. taxonomy or denoising stats) as a table.

## Provenance Replay

**Goal:** Recover the executable commands that produced an artifact, from the artifact alone.

**Approach:** Parse the embedded provenance DAG and regenerate a q2cli (or Artifact-API) script plus a citations BibTeX.

```bash
qiime tools replay-provenance --in-fp core-metrics/ --out-fp replay.sh --usage-driver cli
qiime tools replay-citations  --in-fp core-metrics/ --out-fp citations.bib
# --usage-driver selects cli vs python3/artifact-api output. Verify flag spelling with
# `qiime tools replay-provenance --help` on the installed build (the interface is still maturing).
```

Replay recovers the commands; it is not a guaranteed bit-identical rerun across very different releases (plugin versions are part of the record). The aggregated DAG citations are how a methods section's references come straight from provenance, also via the `.qzv` Citations tab on view.qiime2.org.

## Export (the one-way door)

**Goal:** Hand the data to R/Python when the analysis is no longer expressible in QIIME2 - while losing as little provenance as possible.

**Approach:** Stay in artifacts as long as the work is QIIME2-expressible; export (or read into phyloseq) only at the last step, and keep the upstream `.qza`s so the chain survives up to the exit.

```bash
qiime tools export --input-path table.qza --output-path exported/      # -> exported/feature-table.biom
biom convert -i exported/feature-table.biom -o feature-table.tsv --to-tsv
# FeatureData[Sequence] -> dna-sequences.fasta; FeatureData[Taxonomy] -> taxonomy.tsv; Phylogeny[Rooted] -> tree.nwk
```

Export DROPS the QIIME2 wrapper and the provenance - the exported TSV has no history back to the reads. For R, prefer `qiime2R::qza_to_phyloseq('table.qza', 'taxonomy.qza', 'rooted-tree.qza', 'metadata.tsv')` (Bisanz), which reads artifacts directly and assembles a phyloseq object without manual export. Record where the chain ends.

## Per-Method Failure Modes

### Export-early-loses-provenance
**Trigger:** `qiime tools export` to TSV at step three, then everything else in a notebook. **Mechanism:** export writes native data only and drops the `provenance/` subtree. **Symptom:** the final figure has no provenance back to the raw reads - the framework overhead bought nothing. **Fix:** export at the LAST step; save upstream `.qza`s; or use `qiime2R::qza_to_phyloseq` so the chain survives to the exit.

### Semantic-type mismatch treated as a bug
**Trigger:** feeding a `Phylogeny[Unrooted]` or a `FeatureData[Taxonomy]` where a `FeatureTable[Frequency]` is required. **Mechanism:** the type system refuses incompatible inputs at the interface boundary before running. **Symptom:** "expected an artifact of type ..." error. **Fix:** this is the guard WORKING; `qiime tools peek` to read the actual Type, then fix the UPSTREAM action that produced the wrong type - do not re-import to coerce it.

### Classifier / artifact version break across releases
**Trigger:** a `silva-138-99-nb-classifier.qza` from 2024.x used under 2026.x. **Mechanism:** the sklearn naive-Bayes classifier is pinned to its scikit-learn version; provenance replay assumes recorded plugin versions. **Symptom:** scikit-learn version-mismatch warning/error, or refusal to load. **Fix:** download/train the classifier for YOUR release (the `data.qiime2.org/<release>/common/...` URLs are release-namespaced); retrain or pin the whole env if reusing an old one.

### Manifest Phred / format error on import
**Trigger:** `Phred64V2` on modern Illumina, V1-vs-V2 manifest confusion, relative paths, or `SampleData[...]` for still-multiplexed EMP data. **Mechanism:** the Phred offset is baked into the format name and is applied without checking. **Symptom:** silently mis-decoded quality scores, or an import that "works" but `demux summarize` shows garbage qualities. **Fix:** modern Illumina = Phred33; V2 TSV manifests with absolute paths; `qiime demux summarize` immediately after import; EMP data needs `EMPPairedEndSequences` + `qiime demux`.

### A .qzv treated as data
**Trigger:** trying to feed a `.qzv` into the next action. **Mechanism:** a Visualizer's output is terminal by the framework's type contract. **Symptom:** the action will not accept it as an input. **Fix:** keep and feed the `.qza` the Visualizer was MADE from; a `.qzv` is for viewing only (browser or view.qiime2.org).

### Metadata numeric cast
**Trigger:** an integer subject/batch/timepoint column with no `#q2:types` row. **Mechanism:** inference casts an all-integer column to numeric. **Symptom:** an ID enters a model as a continuous covariate; nonsensical group results. **Fix:** add a `#q2:types` row annotating ID-like columns `categorical`; validate with Keemei.

### Mixing distributions
**Trigger:** expecting amplicon plugins in `moshpit`, or shotgun assembly in `amplicon`. **Mechanism:** distributions are curated, partially-disjoint plugin sets. **Symptom:** "plugin not found." **Fix:** amplicon/marker-gene -> `amplicon` distribution (renamed `qiime2` in 2026.4); shotgun -> `moshpit` (and -> metagenomics); pin both distribution and release.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `--sampling-depth` (rarefaction depth) | -> diversity-analysis | required by `core-metrics`; pick from `alpha-rarefaction`, not a default - the 10000 in examples is a placeholder |
| `--p-formula` integer columns annotated `categorical` | use.qiime2.org metadata reference | otherwise inferred numeric and used as a continuous covariate |
| Phred offset = 33 (modern Illumina) | Illumina format history | Phred64 only for pre-2011 pipelines; wrong choice silently mis-decodes quality |
| Classifier release-match | Bokulich 2018 *Microbiome* 6:90 | the classifier is pinned to its scikit-learn version; cross-release reuse breaks |
| denoise / taxonomy / DA tuning | -> the owning sibling skill | this skill owns no scientific thresholds by design |

Most scientific magic numbers live in the five sibling skills, not here - this skill owns the machinery.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "The scikit-learn version ... could not be found" / classifier won't load | classifier `.qza` trained under a different release | use the release-namespaced classifier or retrain under the current release |
| "Argument ... is not a subtype of ..." / type error | wrong semantic type wired into an action | `qiime tools peek`; fix the upstream action, do not re-import |
| `demux summarize` shows nonsense quality scores | wrong Phred offset in the import format name | re-import with `...Phred33V2`; modern Illumina is Phred33 |
| Import fails on the manifest | V1/V2 confusion, relative paths, wrong delimiter | V2 TSV, absolute paths, tab-separated header `sample-id` |
| A `.qzv` rejected as an action input | Visualizations are terminal | feed the `.qza` it was made from |
| Action treats an ID column as continuous | no `#q2:types` row | annotate the column `categorical`; validate with Keemei |
| Plugin not found | wrong distribution installed | install the `amplicon` (a.k.a. `qiime2` in 2026.4) distribution |

## References

- Bolyen E, Rideout JR, Dillon MR, Bokulich NA, ..., Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nat Biotechnol* 37:852-857.
- Keefe CR, Dillon MR, Gehret E, Herman C, Jewell M, Wood CV, Bolyen E, Caporaso JG. 2023. Facilitating bioinformatics reproducibility with QIIME 2 Provenance Replay. *PLoS Comput Biol* 19(11):e1011676.
- Lin H, Peddada SD. 2020. Analysis of compositions of microbiomes with bias correction. *Nat Commun* 11:3514.
- Rideout JR, Chase JH, Bolyen E, Ackermann G, Gonzalez A, Knight R, Caporaso JG. 2016. Keemei: cloud-based validation of tabular bioinformatics file formats in Google Sheets. *GigaScience* 5:27.
- Bokulich NA, Kaehler BD, Rideout JR, Dillon M, Bolyen E, Knight R, Huttley GA, Caporaso JG. 2018. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. *Microbiome* 6:90.

## Related Skills

- amplicon-processing - DADA2 denoising parameters this skill defers (trunc/trim/maxEE, DADA2 vs Deblur)
- taxonomy-assignment - Classifier and reference-database choice and training behind classify-sklearn
- diversity-analysis - Sampling depth, diversity metric, rarefaction, and the PERMANOVA-vs-dispersion confound
- differential-abundance - DA tool choice and consensus behind composition ancombc
- functional-prediction - PICRUSt2 functional prediction from the feature table
- metagenomics/kraken-classification - Shotgun (moshpit distribution) read classification, not amplicon
- phylogenetics/tree-io - Phylogenetic tree I/O for UniFrac / Faith PD
- read-qc/adapter-trimming - cutadapt primer removal before import/denoising
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
