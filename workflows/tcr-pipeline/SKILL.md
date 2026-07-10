---
name: bio-workflows-tcr-pipeline
description: Orchestrates an end-to-end immune-repertoire pipeline from FASTQ to clonotypes, diversity, overlap, somatic hypermutation and lineages, routing on two forks. Use when deciding bulk vs single-cell (bulk amplicon/RNA-seq -> MiXCR analyze preset -> VDJtools/immunarch depth-normalized diversity and overlap -> figures; 10x paired VDJ -> MiXCR 10x preset or Cell Ranger -> scirpy gene-expression integration, chain QC, clonotype clusters); and TCR vs BCR (TCR -> exact CDR3-nt+V/J clonotypes, VDJtools diversity is fine; BCR -> somatic hypermutation makes exact clonotypes wrong -> Immcantation distToNearest/findThreshold clonal clustering, germline reconstruction, SHM, Dowser lineages); selecting the MiXCR 4.x preset by chemistry and activating its license; downsampling to equal depth before comparing diversity or overlap; and optionally annotating antigen specificity.
tool_type: cli
primary_tool: MiXCR
workflow: true
depends_on:
  - tcr-bcr-analysis/mixcr-analysis
  - tcr-bcr-analysis/vdjtools-analysis
  - tcr-bcr-analysis/immcantation-analysis
  - tcr-bcr-analysis/scirpy-analysis
  - tcr-bcr-analysis/repertoire-visualization
  - tcr-bcr-analysis/specificity-annotation
qc_checkpoints:
  - after_align: "Amplicon alignment >80-90%; RNA-seq legitimately low, judge by clonotype yield; chain usage on-target"
  - after_assemble: "Clonotype count plausible; report UMI/molecule counts, not reads, on UMI libraries"
  - before_diversity: "All samples downsampled to a common depth, else diversity/overlap are confounded by library size"
---

## Version Compatibility

Reference examples tested with: MiXCR 4.7+, VDJtools 1.2.1+, Immcantation suite 4.x, scirpy 0.24+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: MiXCR 4.x replaced the 3.x hand-built `mixcr align -s hsa -p rna-seq` chain with the preset-driven `mixcr analyze <preset>` system, and 4.x refuses to run any command without an activated license (`mixcr activate-license`, or `MI_LICENSE_FILE` on HPC/Docker). A copied 3.x recipe fails on both counts.

# TCR/BCR Repertoire Pipeline

**"Analyze my immune-repertoire sequencing data end-to-end"** -> Route the data by chemistry and receptor, assemble clonotypes with MiXCR, then hand off to depth-normalized diversity (bulk TCR), clonal clustering plus somatic-hypermutation and lineages (BCR), or single-cell gene-expression integration (10x), and finally figures.

This workflow is a router, not a fixed line. Two forks decide everything downstream; pick both before running anything.

## The governing principle: the pipeline forks on two axes

A repertoire measurement is a depth- and chemistry-confounded sample of an unevenly-expanded clonal population, so the correct pipeline depends on how the library was made and which receptor was sequenced. Choosing the wrong branch silently produces plausible-but-wrong numbers.

### Fork A -- bulk vs single-cell

| Axis | Bulk (amplicon or RNA-seq) | Single-cell (10x VDJ) |
|------|----------------------------|-----------------------|
| Assembler | `mixcr analyze <bulk preset>` | `mixcr analyze 10x-sc-xcr-vdj` OR Cell Ranger `vdj` |
| Chain pairing | UNPAIRED (TRB or IGH alone) | Native pairing (TRA+TRB, IGH+IGK/L) |
| Depth vs breadth | Deep repertoire, no cell state | Shallower, links receptor to transcriptome |
| Downstream | VDJtools / immunarch diversity + overlap -> figures | scirpy: chain QC, clonotype clusters, GEX integration |
| Best when | Diversity, overlap, tracking, deep clonotype capture | Antigen-specific cell state, alpha-beta / heavy-light pairing |

### Fork B -- TCR vs BCR

| Axis | TCR (TRA/TRB/TRG/TRD) | BCR (IGH/IGK/IGL) |
|------|-----------------------|-------------------|
| Somatic hypermutation | None | Yes -- clone members are NOT identical |
| Clonotype definition | Exact CDR3-nt + V + J (after UMI/error correction) | NEVER exact CDR3; cluster same-V/J/junction-length by distance |
| Diversity path | VDJtools `CalcDiversityStats` on exact clonotypes is fine | Cluster clones FIRST, then diversity on `clone_id` |
| Extra stages | none | germline reconstruction, SHM, selection, Dowser lineage trees |
| Tool | VDJtools / immunarch | Immcantation (Change-O, SHazaM, SCOPer, Dowser, TIGGER) |

The most common pipeline mistake is running BCR through exact-CDR3 VDJtools diversity. SHM shatters one clone into hundreds of near-identical variants, so exact clonotypes over-count diversity and destroy lineage structure. BCR must route to Immcantation clonal clustering (distToNearest -> findThreshold) before any diversity, SHM, or lineage step. TCR has no SHM, so exact CDR3-nt+V/J is the correct, defensible clonotype and VDJtools diversity is appropriate.

## Pipeline overview

```
FASTQ (+ chemistry, species, receptor known)
    |
    v
[0. License + preset choice] --- mixcr activate-license ; pick preset by kit
    |
    v
[1. MiXCR analyze] ---------- <preset> R1 R2 out_prefix  ->  clones.clns + reports
    |
    +-- QC: mixcr qc / exportQc align + chainUsage
    |
    v
  FORK on data type
    |
    |-- bulk --> [2b. Export] exportClones (VDJtools) / exportAirr
    |               |
    |               v
    |            FORK on receptor
    |               |-- TCR --> [3t. DownSample to equal depth] --> CalcDiversityStats + overlap
    |               |-- BCR --> [3b. Immcantation] distToNearest->findThreshold->clones
    |               |                                -> CreateGermlines --cloned -> SHM -> Dowser trees
    |               v
    |            [4. Visualization] VDJtools / immunarch / ggplot
    |
    |-- single-cell --> [2s. exportAirr / Cell Ranger] --> [3s. scirpy]
                            chain_qc -> ir_dist -> define_clonotypes -> GEX integration
    |
    v
[5. Optional] specificity annotation (VDJdb / GLIPH2 / TCRdist) -- hypothesis, not label
```

## Stage 0: License and preset selection

MiXCR 4.x will not run unlicensed. Activate once (academic license is free), then choose the preset by the exact kit -- the preset encodes species, RNA vs DNA, 5' boundary model (floating for multiplex primers, rigid for 5'-RACE), tag pattern, and assembling feature. The wrong preset does not error; it silently mis-calls V and truncates CDR3.

```bash
mixcr activate-license                 # or: export MI_LICENSE_FILE=/path/mi.license
mixcr exportPreset --preset-name generic-amplicon   # audit what a preset actually does
```

Preset by chemistry (verify against `mixcr analyze --help` and the built-in preset list; MiLaboratories renames occasionally):

| Data | Preset |
|------|--------|
| Generic multiplex/RACE amplicon | `generic-amplicon`, `generic-amplicon-with-umi` (+ `--species hsa`, `--rna`/`--dna`, boundary mixins) |
| Bulk RNA-seq mining | `rna-seq` (judge by clonotype yield, not alignment %) |
| 10x single-cell V(D)J | `10x-sc-xcr-vdj` |
| Takara SMARTer | `takara-human-rna-tcr-umi-smarter-v2`, `takara-human-rna-bcr-umi-smarter` |
| BD Rhapsody | `bd-human-sc-xcr-rhapsody-cdr3` |
| Full component-skill preset table | tcr-bcr-analysis/mixcr-analysis |

## Stage 1: MiXCR assembly (all branches)

```bash
# One command runs align -> refineTagsAndSort -> (assemblePartial) -> assemble -> export.
# From MiXCR 4.7, presets without an intrinsic assembling feature require --assemble-clonotypes-by.
# generic-amplicon REQUIRES material type + both alignment-boundary mixins (it errors without them).
# Multiplex primers on both ends -> floating boundaries; 5'RACE -> --rigid-left-alignment-boundary.
mixcr analyze generic-amplicon \
    --species hsa \
    --rna \
    --floating-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    results/sample

# QC every sample -- low alignment or off-target chains means wrong preset/species/contamination
mixcr qc results/sample.clns
mixcr exportQc align results/*.clns results/qc_align.pdf
mixcr exportQc chainUsage results/*.clns results/qc_chains.pdf
```

**QC checkpoint 1 (after align):** amplicon libraries should align at high rate (often >80-90%); a low rate signals wrong species, wrong boundary model, or untrimmed primers. RNA-seq mining legitimately aligns a tiny fraction -- judge it by absolute clonotype yield. chainUsage catches cross-contamination and index hopping (a TRB library showing appreciable IGH).

Detailed alignment, UMI/cell-barcode handling, and export flags: tcr-bcr-analysis/mixcr-analysis.

## Stage 2: Export (branch-specific handoff)

```bash
# Bulk -> VDJtools-readable clonotype table (per chain)
mixcr exportClones -c TRB results/sample.clns results/sample.clones_TRB.tsv

# BCR or single-cell -> AIRR Rearrangement TSV (the Immcantation / scirpy interchange)
mixcr exportAirr results/sample.clns results/sample.airr.tsv
```

**QC checkpoint 2 (after assemble):** a large reads-to-clonotypes drop-off is normal (millions of reads -> thousands of clones), especially after UMI collapse. Report the right denominator: `uniqueMoleculeCount` on UMI libraries (reporting reads re-introduces the PCR bias the UMIs removed), cells on single-cell, reads only on non-UMI bulk.

## Stage 3t: Bulk TCR diversity and overlap -- downsample FIRST

Diversity (richness, Shannon, clonality) and set-based overlap (Jaccard, shared-clonotype counts) are strictly increasing functions of sequencing depth. Comparing raw values across samples of unequal depth measures depth, not biology -- the single most common error in the field. `DownSample` every sample to a common depth (or read rarefaction curves at a common x) before comparing.

```bash
# 1. Equalize depth: set the target near the cohort lower quartile, and EXCLUDE (do not drag
#    everyone down to) any sample far below it -- an under-sampled library cannot support a claim.
vdjtools DownSample -x 50000 -m metadata.txt ds/

# 2. Diversity on depth-normalized samples; report the resampled table for cross-sample claims
vdjtools CalcDiversityStats -m ds/metadata.txt diversity/

# 3. Overlap with a depth-robust, abundance-weighted metric (F2 / Morisita-Horn), not Jaccard
vdjtools CalcPairwiseDistances -m ds/metadata.txt overlap/
```

**Version caveat (MiXCR 4.x -> VDJtools):** VDJtools is unmaintained for MiXCR 4.x and its parser breaks on raw `exportClones` output (`Unable to parse clonotype string`; 4.x injects commas and renames/moves columns). For a MiXCR 4.7+ cohort, prefer MiXCR's own `mixcr postanalysis individual` / `mixcr postanalysis overlap` (its native 4.x replacement for VDJtools diversity/overlap, with the same downsample-first semantics) over the `vdjtools Convert -S mixcr` route; if VDJtools is required, strip the added contig/target-sequence columns before `Convert`. immunarch (R) is the other modern alternative.

**QC checkpoint 3 (before diversity):** confirm all samples share one depth, and drop any sample whose rarefaction curve is still climbing steeply below that depth (under-sampled -- exclude rather than normalize the cohort down to it). Hold the clonotype match key (nt vs aa, +/-V, +/-J) constant study-wide; aa-level matching inflates apparent sharing via convergent recombination. Report clonality alongside a q=2 Hill number (inverse Simpson) and a rarefaction curve, not alone. immunarch is the modern R alternative with the same normalization semantics: tcr-bcr-analysis/vdjtools-analysis.

## Stage 3b: BCR clonal clustering, SHM and lineages (Immcantation)

BCR cannot use exact clonotypes. Feed the AIRR TSV to Immcantation and follow the mandatory order: annotate -> (TIGGER genotype) -> per-sequence germline -> data-derived clonal threshold -> cluster -> per-clone germline -> SHM/selection -> lineage trees. The threshold from the bimodal distance-to-nearest distribution drives every downstream number; a wrong threshold merges or splits clones.

```r
library(shazam); library(scoper); library(dowser)
db <- airr::read_rearrangement('results/sample.airr.tsv')
# 1. Clonal threshold: valley between the intra-clone and inter-clone modes (per dataset, never reused)
dtn <- distToNearest(db, model = 'ham', normalize = 'len')
thr <- findThreshold(dtn$dist_nearest, method = 'density')@threshold
# 2. Cluster within same-V/J/junction-length partitions at that threshold
cl <- hierarchicalClones(dtn, threshold = thr)
# 3. Reconstruct per-clone germline (CreateGermlines.py --cloned), then observedMutations, then Dowser getTrees
```

Cluster clones within an individual only (genotypes and thresholds are private). Diversity for BCR runs on `clone_id`, not exact CDR3. Full germline reconstruction, SHM quantification, BASELINe selection, and IgPhyML/Dowser trees: tcr-bcr-analysis/immcantation-analysis.

## Stage 3s: Single-cell integration (scirpy)

10x paired VDJ carries native chain pairing and links receptor to cell state. The clonotype definition is a choice, not a default, and multichain cells are likely doublets.

```python
import scirpy as ir
import mudata as mu
airr = ir.io.read_airr('results/sample.airr.tsv')           # or read_10x_vdj on Cell Ranger output
mdata = mu.MuData({'gex': gex_adata, 'airr': airr})         # scirpy 0.13+ stores AIRR as an awkward array
ir.pp.index_chains(mdata)                                    # REQUIRED before chain_qc / clonotyping
ir.tl.chain_qc(mdata)                                        # flag multichain (doublet) / orphan cells
# TCR path: exact-identity clonotypes on CDR3-nt + V/J
ir.pp.ir_dist(mdata)                                         # default metric='identity'
ir.tl.define_clonotypes(mdata)
# BCR path: SHM breaks exact identity, so cluster. ir_dist MUST be recomputed with the SAME
# metric+sequence the clustering call uses (scirpy keys the matrix as ir_dist_{sequence}_{metric};
# the identity matrix above will not satisfy a normalized_hamming call).
# ir.pp.ir_dist(mdata, metric='normalized_hamming', sequence='nt')
# ir.tl.define_clonotype_clusters(mdata, sequence='nt', metric='normalized_hamming', same_v_gene=True)
# integrate with the scanpy GEX modality; measure expansion vs cell state
```

Filtering multichain/orphan cells before expansion analysis preferentially deletes small clones and inflates apparent expansion -- state the trade-off, do not blindly drop them. CellRanger BCR contigs are not IMGT-numbered and include partial/nonproductive contigs; reannotate with IgBLAST (dandelion/airrflow) before rigorous BCR clustering. GEX side (clustering, annotation): single-cell/preprocessing and single-cell/clustering. Full clonotype-definition decisions: tcr-bcr-analysis/scirpy-analysis.

## Stage 4: Visualization

Spectratype (CDR3-length), V-J usage circos, clonal-space bars, rarefaction curves, and clonal tracking across timepoints. Every figure inherits the depth caveat -- plot rarefaction at a common x, and track clones only after downsampling timepoints to equal depth. Recipes: tcr-bcr-analysis/repertoire-visualization.

## Stage 5 (optional): Specificity annotation

Annotate or cluster clonotypes by likely antigen (VDJdb/McPAS lookup, or GLIPH2/TCRdist clustering). A database hit is a sequence match to a published antigen-specific receptor, not proof the clone binds that antigen. "Public" clonotypes are enriched for high generation-probability (Pgen), short, low-insertion CDR3s produced independently in many donors by convergent recombination (Venturi 2006 *PNAS* 103:18691-18696) -- publicity is not antigen selection. Treat every specificity call as a hypothesis, condition on Pgen, and validate. Handoff: tcr-bcr-analysis/specificity-annotation.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| MiXCR exits immediately, "no license" | 4.x needs an activated license | `mixcr activate-license`, or set `MI_LICENSE_FILE` on HPC/Docker; whitelist the phone-home IPs on firewalled clusters |
| `mixcr align -s hsa -p rna-seq` unrecognized | 3.x syntax removed in 4.x | Use `mixcr analyze <preset> R1 R2 out_prefix` |
| Analysis runs but clonotypes look wrong (truncated CDR3, inflated diversity) | Wrong preset -- RNA/DNA, boundary model, or missing tag pattern | Match preset to the exact kit; `mixcr exportPreset` to audit; set `--species` on generic presets |
| MiXCR 4.7 errors on `analyze` needing an assembling feature | Preset lacks an intrinsic assembling feature | Add `--assemble-clonotypes-by CDR3` |
| Diversity/clonality differ across samples but tracks read count | Comparing raw diversity at unequal depth | `DownSample` to a common depth first; report the resampled table / rarefaction at common x |
| BCR clones fragmented, diversity absurdly high, no lineages | Exact-CDR3 clonotypes applied to a hypermutating receptor | Route BCR to Immcantation distToNearest -> findThreshold clustering before any diversity/SHM |
| SHM counts inflated, spurious mutations in junction | No germline reconstruction, or junction not masked | `CreateGermlines.py -g dmask` then `--cloned`; restrict `observedMutations` to `IMGT_V` |
| scirpy expansion inflated by doublets | multichain cells not filtered | Run `chain_qc` and drop multichain cells before clonotype/expansion analysis |
| Overlap dominated by the shallower sample | Jaccard / shared-count on unequal depth | Downsample, then use abundance-weighted F2 or Morisita-Horn |

## Related Skills

- tcr-bcr-analysis/mixcr-analysis - V(D)J alignment and clonotype assembly
- tcr-bcr-analysis/vdjtools-analysis - Depth-normalized diversity and overlap
- tcr-bcr-analysis/immcantation-analysis - BCR clonal clustering, SHM and lineages
- tcr-bcr-analysis/scirpy-analysis - Single-cell VDJ + gene-expression integration
- tcr-bcr-analysis/repertoire-visualization - Figures for the pipeline outputs
- tcr-bcr-analysis/specificity-annotation - Optional antigen-specificity annotation

## References

- Bolotin DA, et al. MiXCR: software for comprehensive adaptive immunity profiling. *Nat Methods* 2015; 12:380-381.
- Shugay M, et al. VDJtools: unifying post-analysis of T cell receptor repertoires. *PLoS Comput Biol* 2015; 11:e1004503.
- Gupta NT, et al. Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. *Bioinformatics* 2015; 31:3356-3358.
- Sturm G, et al. Scirpy: a Scanpy extension for analyzing single-cell T-cell receptor-sequencing data. *Bioinformatics* 2020; 36:4817-4818.
- Chao A, et al. Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. *Ecol Monogr* 2014; 84:45-67.
- Venturi V, et al. Sharing of T cell receptors in antigen-specific responses is driven by convergent recombination. *PNAS* 2006; 103:18691-18696.
