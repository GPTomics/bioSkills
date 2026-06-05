---
name: bio-fragment-analysis
description: Extracts cfDNA fragmentomics features (DELFI genome-wide short/long ratios, WPS nucleosome positioning, Griffin GC-corrected accessibility profiles, end-motifs/MDS, OCF) for cancer detection and tissue-of-origin from plasma WGS. Centers on the nuclease-footprint reframe (every feature re-reads one nucleosome object), the mandatory GC correction, and the cross-protocol non-comparability that breaks naive classifiers. Runs FinaleToolkit (real CLI/Python, MIT) and the Griffin Snakemake pipeline; DELFI is a method, not a package. Use when deriving fragment-based signal from cfDNA, choosing a feature family for detection vs subtyping, or diagnosing why a fragmentomic model failed validation.
tool_type: python
primary_tool: FinaleToolkit
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.2+, pysam 0.22+, finaletoolkit 0.7+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: FinaleToolkit CLI subcommands are hyphenated (`frag-length-bins`, `end-motifs`, `delfi-gc-correct`); the Python functions are underscored in `finaletoolkit.frag`. The filter subcommand is `filter-file`, NOT `filter-bam`. Griffin is a Snakemake pipeline, not an importable function. DELFI is a methodology and a company (DELFI Diagnostics), not a `pip install`-able tool.

# Fragment Analysis

**"Analyze cfDNA fragment patterns for cancer signal"** -> Quantify nucleosome-footprint features (size ratios, protection scores, accessibility profiles, end motifs) from plasma WGS for detection or tissue-of-origin.
- Python/CLI: `finaletoolkit` for DELFI ratios, WPS, end-motifs, MDS, cleavage profiles
- Snakemake: `Griffin` pipeline for GC-corrected nucleosome profiling at TF/accessible sites
- Python: `pysam` for a custom binned short/long ratio when a dependency-light readout is wanted

## The Single Most Important Modern Insight -- Fragmentomics measures nucleosome positioning, not sequence; the biology is real but it can be destroyed in the wet lab or left un-GC-corrected

Plasma cfDNA is the digestion product of chromatin by apoptotic and intracellular nucleases, so every fragmentomic feature is a re-readout of the same physical object: the nucleosome footprint. The ~167 bp mode is the 147 bp histone-protected core plus ~20 bp of linker; the 10.4 bp sawtooth below it is the helical pitch of DNA on the histone surface (the nuclease cuts only where the minor groove faces out, once per turn). DELFI ratios, WPS, Griffin profiles, end-motifs, and OCF are four views of this one object, not four independent measurements -- their correlation inflates apparent multi-feature performance and leaks across train/test splits.

Two consequences dominate practice. First, the single biggest threat to any fragmentomic feature is GC, library, and batch confounding, not biology: an uncorrected genome-wide short/long ratio tracks GC content and library prep far more strongly than tumor fraction, which is why naive fragmentomics works in discovery and dies in validation. Griffin's actual contribution is its fragment-length-specific GC correction, not the nucleosome plot. Second, DELFI-style ratios are partly entangled with copy-number alteration and coverage (a bin's ratio reflects both fragmentation state and how many genomes contributed): deconvolving the fragmentation-specific signal from CNA is a known open problem, so a raw 5 Mb ratio is a hybrid CNA + fragmentation + GC readout, not pure fragmentation.

## Methods / Feature-Family Landscape

| Feature family | Primary method | What it physically measures | Citation |
|---|---|---|---|
| Genome-wide short/long ratio | DELFI: ~5 Mb bins, GC-corrected short(100-150)/long(151-220), boosted classifier | Coarse fragmentation state across the genome (entangled with CNA + GC) | Cristiano 2019 Nature 570:385 |
| Nucleosome positioning | WPS: spanning fragments minus end-containing fragments in a sliding window | Where nucleosomes sit; promoter/gene-body phasing encodes tissue + expression | Snyder 2016 Cell 164:57 |
| GC-corrected accessibility | Griffin: length-specific GC correction then composite coverage around site sets | TF/DHS accessibility in the tissue of origin; robust at low tumor fraction | Doebley 2022 Nat Commun 13:7475 |
| End motifs / MDS | 4-mer at the 5' cut site; MDS = normalized Shannon entropy of the 256 motifs | Nuclease-cleavage signature (DNASE1L3 sculpts the normal CC-ending spectrum) | Jiang 2020 Cancer Discov 10:664 |
| Orientation-aware ends | OCF: phase offset between upstream- and downstream-end peaks in open chromatin | Tissue-of-origin via end orientation, not coverage | Sun 2019 Genome Res 29:418 |

DELFI = DNA EvaLuation of Fragments; the end-motif/MDS biology is anchored in DNASE1L3, whose deletion reorders length and end-motif frequencies (Serpas 2019 PNAS 116:641). Methodology here is still evolving (CNA deconvolution, standalone GC correctors like GCparagon) -- verify current best practice against live FinaleToolkit and Griffin docs before committing to one feature family.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|---|---|---|
| Cancer detection, pan-genome screen | DELFI genome-wide short/long profile (FinaleToolkit `delfi`) | Coarse genome-wide signal is what the boosted classifier was built on; cheap at low-pass WGS |
| Tissue / subtype of origin (e.g. ER status, NEPC) | Griffin nucleosome profiling around TF/accessible site sets | Accessibility composite scales with the contributing tissue; the GC correction makes it portable |
| Low tumor fraction (TF < ~0.03) | Griffin (GC-correction removes dominant technical signal) and/or in-silico size selection (90-150 bp) | Raw ratios are GC-dominated at low TF; Griffin holds AUC ~0.92 vs ~0.99 at TF >= 0.05 |
| Nuclease / cleavage biology, single-scalar comparison | End-motifs + MDS (FinaleToolkit `end-motifs` then `mds`) | MDS is one number per sample; rises in cancer as orderly DNASE1L3 cleavage is lost |
| Nucleosome positions / TF footprints directly | WPS (FinaleToolkit `wps` + `adjust-wps`) | WPS peaks recover nucleosome positions; S-WPS exposes TF footprints |
| Any cross-batch or cross-protocol comparison | GC-correct AND co-process through one pipeline, or do not compare | Uncorrected, cross-protocol fragmentomics is uninterpretable (see Failure Modes) |

## Genome-Wide Short/Long Ratio (custom DELFI-style)

**Goal:** Produce a genome-wide vector of short-to-long fragment ratios in fixed bins as a dependency-light DELFI-style feature, with the explicit caveat that without GC correction it is a GC + CNA readout.

**Approach:** Walk proper-pair fragments per bin from the BAM template length, classify each as short (100-150 bp) or long (151-220 bp), and emit the per-bin ratio. For a publication-grade profile, prefer FinaleToolkit `delfi` (which GC-corrects) over this illustrative version.

```python
import pysam
import numpy as np
import pandas as pd

def binned_short_long_ratio(bam_path, bin_size=5_000_000, chroms=None):
    '''Per-bin short(100-150)/long(151-220) ratio. NOT GC-corrected -- illustrative only.'''
    chroms = chroms or [f'chr{i}' for i in range(1, 23)]
    bam = pysam.AlignmentFile(bam_path, 'rb')
    rows = []
    for chrom in chroms:
        if chrom not in bam.references:
            continue
        n_bins = bam.get_reference_length(chrom) // bin_size + 1
        short = np.zeros(n_bins)
        long = np.zeros(n_bins)
        for read in bam.fetch(chrom):
            if not read.is_proper_pair or read.is_secondary or read.template_length <= 0:
                continue
            size = read.template_length
            b = read.reference_start // bin_size
            if 100 <= size <= 150:
                short[b] += 1
            elif 151 <= size <= 220:
                long[b] += 1
        ratio = np.divide(short, long, out=np.full(n_bins, np.nan), where=long > 0)
        rows.extend({'chrom': chrom, 'bin': i, 'short': short[i], 'long': long[i], 'ratio': ratio[i]} for i in range(n_bins))
    bam.close()
    return pd.DataFrame(rows)
```

## GC-Corrected DELFI Score (FinaleToolkit)

**Goal:** Compute a GC-corrected DELFI score so the genome-wide profile reflects fragmentation rather than base composition.

**Approach:** FinaleToolkit's `delfi` corrects short and long bin counts for GC before forming the ratio; the CLI and Python API are equivalent. Run on a BAM/CRAM or a tabix-indexed `.frag.gz` fragment file.

```bash
# CLI (subcommands are hyphenated). delfi positionals: input chrom_sizes reference bins_file.
# GC correction is ON by default (-G disables it); 100kb bins are merged to 5Mb by default.
# -R keeps no-coverage regions when the genome is not hg19.
finaletoolkit delfi sample.bam hg38.chrom.sizes hg38.fa bins_100kb.bed -g gaps.bed -R -o sample.delfi.bed
finaletoolkit end-motifs sample.bam hg38.fa -o sample.end_motifs.tsv
finaletoolkit mds sample.end_motifs.tsv          # Motif Diversity Score (normalized Shannon entropy)
finaletoolkit wps sample.bam sites.bed -c hg38.chrom.sizes -o sites.wps.bw   # per-site, not a single region
```

```python
from finaletoolkit.frag import delfi, end_motifs, wps  # public finaletoolkit.frag symbols
# delfi() returns GC-corrected short/long per bin; end_motifs() returns an EndMotifFreqs
# object whose .motif_diversity_score() gives the MDS (there is no top-level frag.mds).
```

## Griffin Nucleosome Profiling (Snakemake pipeline)

**Goal:** Obtain GC-corrected composite coverage around a TF/accessible-site set for tissue-of-origin, robust at low tumor fraction and ~0.1x WGS.

**Approach:** Griffin is not an importable function; it is three sequential Snakemake modules. Run them in order against `samples.yaml`, the hg38 reference, and a `sites.yaml` site list.

```bash
# Run each module from Griffin's snakemakes/ dir (config edited per cohort)
snakemake -s griffin_genome_GC_frequency/griffin_genome_GC_frequency.snakefile --cores 8
snakemake -s griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction.snakefile --cores 8
snakemake -s griffin_nucleosome_profiling/griffin_nucleosome_profiling.snakefile --cores 8
# Output: GC-corrected + uncorrected composite coverage profiles around each site set.
```

## Per-Method Failure Modes

### Uncorrected GC dominates the ratio
Trigger: comparing raw short/long ratios across samples without GC correction. Mechanism: PCR and binding-based purification overrepresent GC-balanced fragments, and short vs long fragments have different GC dependence. Symptom: a beautiful discovery-cohort separation that collapses in validation; the profile clusters by sequencing batch. Fix: GC-correct (FinaleToolkit `delfi`/`delfi-gc-correct`, Griffin, or GCparagon) and co-process all samples through one pipeline.

### Cross-protocol non-comparability
Trigger: combining ssDNA and dsDNA libraries, or two end-repair/PCR chemistries, in one analysis. Mechanism: ssDNA prep recovers the sub-100 bp ultrashort population that dsDNA prep loses at the double-strand ligation step, shifting the entire size distribution and every derived feature. Symptom: a model trained on one chemistry mislabels the other systematically. Fix: a fragmentomic model is conditioned on its library chemistry -- match protocols, never cross them, and state the prep as a precondition.

### DELFI / CNA entanglement
Trigger: interpreting a 5 Mb ratio bin as pure fragmentation. Mechanism: copy-number alterations change how many genomes contribute to a bin, moving coverage and therefore the ratio independent of fragmentation. Symptom: ratio "signal" that mirrors the CNA profile. Fix: treat DELFI as a hybrid CNA + fragmentation + GC feature; deconvolve with caution and do not overclaim a pure fragmentation readout.

### Low-coverage WPS noise
Trigger: computing WPS or per-site profiles on too few fragments. Mechanism: WPS is a difference of spanning vs end-containing counts; at low coverage both terms are tiny and the score is dominated by sampling noise. Symptom: no clean nucleosome periodicity, jagged tracks. Fix: aggregate over many copies of a site (composite profiles, Griffin/`multi_wps`), smooth (`adjust-wps`), and require adequate depth before single-locus WPS.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Mononucleosome mode ~167 bp (147 core + ~20 linker) | Snyder 2016 Cell 164:57 | The protected core is 147 bp; the variable ~20 bp linker is the rest of the mode |
| 10.4 bp periodicity below 167 bp | Snyder 2016 | Helical pitch of B-form DNA; the cleanest sanity check that footprints and size estimation are sane |
| Di-/tri-nucleosome ~334 / ~500 bp | Snyder 2016 | Successive nucleosomes add ~167 bp each |
| ctDNA mode ~20-50 bp shorter (toward ~145 bp), enrich 90-150 bp | Mouliere 2018 Sci Transl Med 10:eaat4921 | Tumor chromatin/nuclease processing shifts length down; the lever size selection exploits |
| Short 100-150 bp vs long 151-220 bp | Cristiano 2019 Nature 570:385 | The DELFI ratio numerator/denominator windows |
| ~5 Mb DELFI bins | Cristiano 2019 | Bin scale at which the genome-wide ratio vector was defined and classified |
| In-silico size selection 90-150 bp | Mouliere 2018 | Retaining this window enriches tumor fraction >2x in >95% of cases, >4x in >10% |
| End-motif 4-mer, 256 categories; MDS = normalized Shannon entropy | Jiang 2020 Cancer Discov 10:664 | The categorical end-motif space and its single-scalar diversity summary |

Size selection is a tumor-fraction-vs-depth lever, not a universal win: it discards the 167 bp bulk, so it helps when ctDNA is dilute and short but hurts when already depth-limited.

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `pip install delfi` fails / no DELFI CLI | DELFI is a method + company, not a package | Compute DELFI features via FinaleToolkit `delfi`, or a custom binned ratio |
| `finaletoolkit filter-bam` not found | The subcommand is `filter-file`, not `filter-bam` | Use `finaletoolkit filter-file` for mapq/size/region filtering |
| `import griffin` fails | Griffin is a Snakemake pipeline, not an importable module | Run the three `griffin_*` snakefiles in sequence |
| Profiles cluster by batch, not biology | Uncorrected GC / mixed protocols | GC-correct and co-process one chemistry through one pipeline |
| Confusing end-motif and breakpoint-motif | Different objects (cut-site 4-mer vs k-mer spanning the cut) | FinaleToolkit exposes both: `end-motifs` vs `breakpoint-motifs` |

## References

- Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. 2016. Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. *Cell* 164(1-2):57-68. -- nucleosome footprint, 167 bp, 10.4 bp periodicity, WPS.
- Cristiano S, Leal A, Phallen J, et al. 2019. Genome-wide cell-free DNA fragmentation in patients with cancer. *Nature* 570(7761):385-389. -- DELFI; 5 Mb bins, short/long ratio.
- Mouliere F, Chandrananda D, Piskorz AM, et al. 2018. Enhanced detection of circulating tumor DNA by fragment size analysis. *Sci Transl Med* 10(466):eaat4921. -- size selection, 90-150 bp enrichment.
- Doebley A-L, Ko M, Liao H, et al. 2022. A framework for clinical cancer subtyping from nucleosome profiling of cell-free DNA. *Nat Commun* 13:7475. -- Griffin; length-specific GC correction.
- Jiang P, Sun K, Peng W, et al. 2020. Plasma DNA end-motif profiling as a fragmentomic marker in cancer, pregnancy, and transplantation. *Cancer Discov* 10(5):664-673. -- end motifs and MDS (journal is *Cancer Discovery*, not PNAS).
- Sun K, Jiang P, Wong AIC, et al. 2019. Orientation-aware plasma cell-free DNA fragmentation analysis in open chromatin regions informs tissue of origin. *Genome Res* 29(3):418-427. -- OCF.
- Serpas L, Chan RWY, Jiang P, et al. 2019. Dnase1l3 deletion causes aberrations in length and end-motif frequencies in plasma DNA. *PNAS* 116(2):641-649. -- DNASE1L3 biology behind end motifs.
- Burnham P, Kim MS, Agbor-Enoh S, et al. 2016. Single-stranded DNA library preparation uncovers the origin and diversity of ultrashort cell-free DNA in plasma. *Sci Rep* 6:27859. -- ssDNA prep recovers ultrashort cfDNA.
- FinaleToolkit: accelerating cell-free DNA fragmentation analysis with a high-speed computational toolkit. 2025. *Bioinformatics Advances* 5(1):vbaf236. -- the toolkit; ~50x faster WPS than the original Snyder implementation on BH01 (scope-limited).

## Related Skills

- cfdna-preprocessing - library prep determines which fragments (and features) are recoverable
- tumor-fraction-estimation - fragmentomics enables signal below the CNA-based TF floor
- methylation-based-detection - orthogonal genome-wide cfDNA signal
- atac-seq/nucleosome-positioning - shared nucleosome-footprint biology
