---
name: bio-copy-number-focal-amplification-ecdna
description: Resolve the architecture of focal oncogene amplifications — extrachromosomal DNA (ecDNA), breakage-fusion-bridge (BFB) cycles, homogeneously staining regions (HSR), and linear amplification — from whole-genome sequencing with AmpliconArchitect, the AmpliconSuite pipeline, and AmpliconClassifier. Covers copy-number seed selection, breakpoint-graph reconstruction, balanced-flow optimization, ecDNA classification, and the limits of depth-only amplification calls. Use when a focal amplification needs structural characterization, when distinguishing ecDNA from chromosomal amplification, suspecting ecDNA-driven oncogene amplification or therapy resistance, or selecting copy-number seeds for amplicon reconstruction.
tool_type: cli
primary_tool: AmpliconArchitect
---

## Version Compatibility

Reference examples tested with: AmpliconSuite-pipeline 1.3+, AmpliconArchitect 1.3+, AmpliconClassifier 1.2+, CNVkit 0.9.10+, Python 3.10+, samtools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `AmpliconSuite-pipeline.py --help`, `amplicon_classifier.py --help`
- AmpliconArchitect needs a `$AA_DATA_REPO` reference download and a Mosek license (free for academic use); confirm both are configured before running

Verify the reference build — AmpliconArchitect was historically hg19-centric; GRCh38 support and data repos exist but the build must be set explicitly and consistently.

# Focal Amplification and ecDNA

**"This oncogene is amplified — but how, structurally"** -> A depth caller reports "high focal amplification" and stops. The biology depends entirely on the *architecture*: extrachromosomal DNA (ecDNA) behaves utterly differently from a chromosomal homogeneously staining region. Resolving architecture needs the breakpoint graph, not depth.

- CLI: `AmpliconSuite-pipeline.py` (end-to-end), `AmpliconArchitect` (graph reconstruction), `AmpliconClassifier` (architecture call)
- Input: WGS BAM plus copy-number seeds (high-CN focal regions)

## Why Architecture Matters — Four Amplicon Classes

| Class | Structure | Behavior | Why it matters |
|-------|-----------|----------|----------------|
| ecDNA | Circular, episomal, no centromere | Hundreds of copies; unequal mitotic segregation; rapid CN adaptation | Drives oncogene overexpression, intratumor heterogeneity, therapy resistance; ~14% of cancers |
| BFB | Chromosomal, fold-back inversions | Stepwise CN gradient toward telomere | Distinct breakpoint signature; bounded amplification |
| HSR | Linear, integrated chromosomally | Stable inheritance | Chromosomal — segregates evenly, unlike ecDNA |
| Linear/simple | Tandem or simple amplification | Modest copy gain | Often passenger-scale; lowest oncogenic concern |

ecDNA is the highest-stakes call: because it lacks a centromere it segregates unequally, so copy number can surge under selection — a structural basis for resistance. Depth alone cannot distinguish ecDNA from an HSR; both look like a high-amplitude focal gain.

## When to Suspect ecDNA

| Signal | Interpretation |
|--------|----------------|
| Very high focal copy number (CN >> 10) at an oncogene | Consistent with ecDNA; not specific |
| Amplicon spanning multiple non-contiguous genomic segments | Suggestive — ecDNA often fuses distal regions |
| Breakpoint graph forms a closed cycle with balanced flow | AmpliconArchitect signature of circular structure |
| Highly variable per-cell copy number (single-cell / FISH) | Hallmark of unequal ecDNA segregation |
| Co-amplified enhancers distal to the oncogene | ecDNA can hijack regulatory elements |

## The AmpliconSuite Workflow

AmpliconArchitect does not call amplifications from scratch — it reconstructs the architecture of the regions it is seeded with. The pipeline is: (1) call copy number and select high-CN focal seeds, (2) AmpliconArchitect builds the breakpoint graph and optimizes a balanced flow, (3) AmpliconClassifier labels each amplicon ecDNA / BFB / HSR / linear.

```bash
# End-to-end: AmpliconSuite-pipeline runs CNVkit seeding, AmpliconArchitect, and
# AmpliconClassifier in sequence.
AmpliconSuite-pipeline.py \
    -s sample_id \
    -t 8 \
    --bam tumor.bam \
    --ref GRCh38 \
    --run_AA --run_AC

# Output: per-amplicon breakpoint graphs, cycles files, and an AmpliconClassifier
# table assigning each amplicon an architecture class.
```

Supplying explicit seeds (recommended when a vetted CNV callset exists):

```bash
# Seeds: a BED of high-copy focal regions (e.g. from cnvkit-analysis), filtered to
# CN above the seed threshold and to focal (not arm-level) size.
AmpliconSuite-pipeline.py -s sample_id -t 8 --bam tumor.bam --ref GRCh38 \
    --cnv_bed focal_seeds.bed --run_AA --run_AC
```

## Failure Modes

### Garbage copy-number seeds produce garbage amplicons

**Trigger:** Seeding AmpliconArchitect with a noisy CNV callset, a flat-reference tumor-only callset, or arm-level segments.

**Mechanism:** AmpliconArchitect reconstructs the architecture of exactly the regions it is seeded with; false high-CN seeds generate spurious amplicons, and arm-level seeds dilute the focal signal.

**Symptom:** Implausible amplicons at no known oncogene; amplicons spanning whole arms; classifier output dominated by low-confidence calls.

**Fix:** Seed only vetted, focal, high-CN regions. Build the CNV callset from a proper panel of normals (see cnvkit-analysis); filter to focal size and CN above the seed threshold before passing to AA.

### Calling ecDNA from depth alone

**Trigger:** Labeling a high-amplitude focal gain "ecDNA" without breakpoint-graph evidence.

**Mechanism:** ecDNA and a chromosomal HSR both present as high focal copy number; only the breakpoint graph (a closed cycle with balanced flow) distinguishes them.

**Symptom:** ecDNA claimed from a CNVkit/GATK profile; no graph, no cycle.

**Fix:** Require AmpliconArchitect graph reconstruction and an AmpliconClassifier ecDNA call. Where feasible, confirm with orthogonal evidence — FISH, single-cell copy number (variable per-cell CN), or optical mapping.

### Genome-build mismatch

**Trigger:** BAM aligned to one build, `--ref` or `$AA_DATA_REPO` set to another.

**Mechanism:** Coordinates and the bundled annotation diverge; breakpoints and genes are mis-assigned.

**Symptom:** Amplicons at wrong loci; AA errors on contig names.

**Fix:** Set `--ref` to match the BAM's build and confirm the corresponding `$AA_DATA_REPO` is installed; AA was historically hg19-centric, so GRCh38 must be explicit.

### Short-read limits on complex amplicon resolution

**Trigger:** Expecting a fully resolved amplicon structure from short-read WGS on a highly rearranged amplicon.

**Mechanism:** Short reads cannot phase long-range structure or traverse repeats; complex amplicons (many junctions, segmental duplications) are only partially reconstructed.

**Symptom:** Fragmented breakpoint graph; ambiguous or "unknown" classifier calls on a clearly amplified locus.

**Fix:** Treat short-read amplicon structure as a hypothesis for the most complex cases; confirm with optical mapping (AmpliconReconstructor) or long-read sequencing.

### Inadequate coverage or FFPE input

**Trigger:** Low-coverage WGS or degraded FFPE DNA.

**Mechanism:** Breakpoint detection needs sufficient discordant/split-read support; FFPE artifacts add false junctions.

**Symptom:** Missing junctions; noisy graph; unstable classification.

**Fix:** Use adequate-coverage WGS (AmpliconArchitect is designed for WGS, not panels/WES); apply FFPE-aware filtering; corroborate junctions across read-pair and split-read evidence.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Depth caller: "amplification"; AA: ecDNA | Architecture only visible in the graph | Trust AA for architecture; depth gives amplitude only |
| AA ecDNA vs FISH negative | Subclonal ecDNA, or false-positive cycle | Check cell fraction; review graph balanced flow |
| AA "unknown" on a clear amplicon | Complex structure beyond short-read resolution | Escalate to optical mapping / long-read |
| BFB vs ecDNA ambiguous | Fold-back and circular signatures overlap | Inspect CN gradient (BFB) vs closed cycle (ecDNA) |

**Operational rule:** A depth caller establishes *that* a region is amplified and *how much*; it never establishes the architecture. An ecDNA call requires an AmpliconArchitect breakpoint graph with a closed cycle and an AmpliconClassifier ecDNA label, and ideally orthogonal confirmation (FISH, single-cell, optical mapping). Seeds must be vetted focal high-CN regions, not raw or arm-level calls.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| ecDNA prevalence | ~14% of cancers | Kim et al 2020; baseline expectation |
| CN seed threshold | CN >= ~4-5 focal | AmpliconSuite seeding; amplicons, not single-copy gains |
| Seed size | focal (sub-arm), not whole-arm | Arm-level seeds dilute focal amplicon signal |
| Assay | whole-genome sequencing | AmpliconArchitect needs genome-wide breakpoint coverage |
| Confirmation for ecDNA | graph cycle + classifier + orthogonal evidence | Depth alone is insufficient |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Amplicons at no known oncogene | Noisy or arm-level seeds | Seed vetted focal high-CN regions only |
| ecDNA "called" from a CNVkit profile | Depth-only claim, no graph | Run AmpliconArchitect + AmpliconClassifier |
| AA errors on contig names | Build mismatch | Match `--ref` and `$AA_DATA_REPO` to the BAM |
| AA fails to start | Missing Mosek license / data repo | Configure the academic Mosek license and `$AA_DATA_REPO` |
| Fragmented graph on a clear amplicon | Short-read limits / low coverage | Confirm with optical mapping or long reads |
| Classifier output all low-confidence | Coverage too low or FFPE artifacts | Use adequate-coverage WGS; FFPE-aware filtering |

## References

- Turner KM et al 2017. Extrachromosomal oncogene amplification drives tumour evolution and genetic heterogeneity. Nature 543:122
- Deshpande V et al 2019. Exploring the landscape of focal amplifications in cancer using AmpliconArchitect. Nat Commun 10:392
- Kim H et al 2020. Extrachromosomal DNA is associated with oncogene amplification and poor outcome across multiple cancers. Nat Genet 52:891
- Luebeck J et al 2024. AmpliconSuite: an end-to-end workflow for analyzing focal amplifications in cancer genomes. bioRxiv (AmpliconSuite-pipeline)
- Luebeck J et al 2020. AmpliconReconstructor integrates NGS and optical mapping to resolve focal amplifications. Nat Commun 11:4374

## Related Skills

- copy-number/cnvkit-analysis - Generates the copy-number seeds for amplicon reconstruction
- copy-number/recurrent-cnv - Cohort-level recurrent focal amplification (GISTIC2)
- copy-number/allele-specific-copy-number - Absolute copy number of amplified loci
- copy-number/cnv-annotation - Oncogene annotation of amplified regions
- copy-number/subclonal-copy-number - Subclonal dynamics of ecDNA copy number
- long-read-sequencing/structural-variants - Long-read resolution of complex amplicons
