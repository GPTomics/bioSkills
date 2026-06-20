---
name: bio-ribo-seq-orf-detection
description: Detect and quantify translated ORFs from Ribo-seq using 3-nucleotide periodicity, including uORFs, internal ORFs, dORFs, and novel ORFs. Use when finding actively translated regions beyond annotated CDS, classifying ORFs by the 2022 community standard, quantifying ORF-level translation, or choosing between periodicity-based callers.
tool_type: mixed
primary_tool: RiboCode
---

## Version Compatibility

Reference examples tested with: RiboCode 1.2+, ORFquant 1.0+, ORFik 1.22+, DESeq2 1.42+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# ORF Detection

**"Detect translated ORFs from my Ribo-seq data"** -> Identify actively translated open reading frames (uORFs, internal ORFs, dORFs, novel ORFs) using 3-nucleotide periodicity (not mere coverage) as the evidence of translation, then classify and quantify them.
- CLI: `RiboCode` for periodicity-based de novo ORF calling
- R: `ORFquant` for isoform-aware quantification; `ORFik` as the general toolkit

The discriminating signal is PERIODICITY: a translated ORF shows footprint P-sites in frame 0 (F0 >> F1, F2). Coverage alone is not evidence of translation. Every method needs a correct per-read-length P-site offset first (see ribosome-periodicity).

## ORF-type nomenclature (Mudge 2022 standard)

The GENCODE-led standard (Mudge et al 2022) defines the umbrella term "Ribo-seq ORF" and six positional categories. These are positional, not modification-based (N-terminal extensions are not part of the scheme).

| Category | Definition (transcript-relative) |
|----------|----------------------------------|
| uORF | Entirely within the 5' UTR, not overlapping the CDS |
| uoORF | Upstream-overlapping: starts in 5' UTR, overlaps CDS start out-of-frame |
| intORF | Internal/nested: within the CDS in a different frame |
| dORF | Entirely within the 3' UTR, not overlapping the CDS |
| doORF | Downstream-overlapping: overlaps the CDS stop into the 3' UTR |
| lncRNA-ORF | ORF on a transcript annotated as long non-coding RNA |

Related terms: sORF/smORF (<100 codons, product = microprotein), annotated CDS, novel. Catalogs: sORFs.org, OpenProt.

## Near-cognate start codons (the ATG-only blind spot)

Initiation occurs at AUG and at near-cognate codons differing from AUG by one base; the biologically used set is CUG, GUG, ACG, UUG, AUU, AUC, AUA (AAG/AGG also differ by one base but initiate negligibly). CUG is the dominant near-cognate start (~16% of mapped initiation sites; AUG remains >50%). uORFs ESPECIALLY use near-cognate starts, so an ATG-only scanner misses the majority of real uORFs. Periodicity-based callers can be configured with alternative starts (RiboCode `-A`); the manual finder below is ATG-only and is a teaching toy unless extended.

## Tool selection

| Situation | Tool | Why |
|-----------|------|-----|
| De novo discovery (uORFs, novel ORFs), standard Ribo-seq | RiboCode | Periodicity-based, maintained, supports alternative starts |
| Isoform-aware detection + per-ORF quantification | ORFquant | Resolves ORFs across overlapping isoforms; built on Ribo-seQC |
| General R toolkit (uORF finding, P-site shift, TE, plots) | ORFik | Comprehensive; NOT a dedicated de novo caller |
| Initiation-site / non-AUG mapping (needs harringtonine/LTM) | Ribo-TISH, PRICE | TI-seq-aware (see initiation-site-mapping) |
| Assay-agnostic, no TI-seq, short + long ORFs | ribotricer | Phasing-only, species-calibrated cutoff |
| Bacteria/prokaryotes | DeepRibo, smORFer | Prokaryote-trained (eukaryote periodicity tools fit poorly) |
| Differential ORF translation | P-site counts per ORF then DESeq2 | Count-based DE on ORF-level counts |

ORFik and ORFquant are DISTINCT packages (different authors, repos, methods): ORFik (Tjeldnes 2021, general toolkit) is not the same as ORFquant (Calviello 2020, dedicated isoform-aware caller). Do not install ORFik expecting ORFquant.

## Call ORFs de novo with RiboCode

**Goal:** Identify periodicity-significant ORFs, including uORFs at near-cognate starts.

**Approach:** Prepare transcript annotation, run `metaplots` to select periodic read lengths and their P-site offsets, then run `RiboCode` with optional alternative starts.

```bash
# Step 1: annotation
prepare_transcripts -g annotation.gtf -f genome.fa -o ribocode_annot

# Step 2: metaplots picks periodic read lengths + per-length P-site offsets -> config .txt
#   (read lengths come from THIS step, NOT from a -l flag)
metaplots -a ribocode_annot -r transcriptome.bam -o metaplots_out

# Step 3: call ORFs. -A adds near-cognate starts; -l is the longest-ORF toggle (yes/no)
RiboCode -a ribocode_annot -c metaplots_out_pre_config.txt \
    -A CTG,GTG -l no -p 0.05 -o ribocode_result
```

RiboCode works in transcript coordinates, so the `-r` input is the TRANSCRIPTOME-projected BAM (`Aligned.toTranscriptome.out.bam`), not the genome BAM; a genome BAM silently misbehaves. Its core test is a MODIFIED WILCOXON SIGNED-RANK test on the per-codon P-site frame distribution (a separate binomial file is a secondary output). The `-l` flag toggles longest-ORF selection; it is NOT a read-length list.

## Parse and classify RiboCode output

**Goal:** Split called ORFs by the standard categories.

**Approach:** Read the tabular result and group on the `ORF_type` column, whose RiboCode values are `annotated`, `uORF`, `dORF`, `Overlap_uORF`, `Overlap_dORF`, `Internal`, `novel`.

```python
import pandas as pd

def load_ribocode_orfs(path):
    '''Load the RiboCode result table (<output_name>.txt) and group by ORF_type.'''
    df = pd.read_csv(path, sep='\t')
    groups = {t: df[df['ORF_type'] == t] for t in df['ORF_type'].unique()}
    return df, groups
```

RiboCode writes the result as `<output_name>.txt` (plus a `<output_name>_collapsed.txt`), e.g. `ribocode_result.txt` for `-o ribocode_result`; its columns include `ORF_ID`, `ORF_type`, transcript/genome start-stop, `pval_combined`, and `adjusted_pval`.

## Quantify ORFs isoform-aware with ORFquant

**Goal:** Quantify ORF-level translation while resolving footprints across overlapping isoforms.

**Approach:** Prepare annotation once, feed Ribo-seQC-prepared input, then run the master function.

```r
library(ORFquant)

prepare_annotation_files(annotation_directory = "annot/",
                         twobit_file = "genome.2bit",
                         gtf_file = "annotation.gtf")
# Ribo-seQC writes a for_ORFquant object from the Ribo-seq BAM; pass it here
run_ORFquant(for_ORFquant_file = "sample_for_ORFquant",
             annotation_file = "annot/annotation.gtf_Rannot",
             n_cores = 4)
```

For uORF discovery in a general R workflow, ORFik provides `findUORFs()` and the true P-site shift is `detectRibosomeShifts()` then `shiftFootprints()` (there are no `p_offsets`/`lengths` arguments on `fimport`).

## Manual ORF scan (teaching reference, ATG-only)

**Goal:** Illustrate ORF finding mechanics; not a substitute for a periodicity caller.

**Approach:** Scan three frames for start-to-stop pairs. This finds only ATG starts and uses coverage, not periodicity, so it misses near-cognate uORFs and cannot confirm active translation on its own.

```python
def find_orfs(seq, min_codons=10):
    '''Find ATG-to-stop ORFs in all three frames (ATG-only: a teaching toy).'''
    seq = str(seq).upper()
    stops = {'TAA', 'TAG', 'TGA'}
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            if seq[i:i+3] == 'ATG':
                for j in range(i + 3, len(seq) - 2, 3):
                    if seq[j:j+3] in stops:
                        if (j + 3 - i) >= min_codons * 3:
                            orfs.append({'start': i, 'end': j + 3, 'frame': frame})
                        i = j
                        break
            i += 3
    return orfs
```

## Validate called ORFs

**Goal:** Separate genuine translation from coverage artifacts.

**Approach:** Check the in-frame (frame-0) fraction within the ORF; compare the ORF's footprint read-length distribution to annotated CDS with FLOSS; use ORFscore for frame bias; add PhyloCSF/conservation or mass-spec peptide evidence for novel microproteins.

FLOSS (Ingolia 2014) is half the summed absolute difference between an ORF's footprint length-fraction histogram and the CDS reference; a high-coverage region whose length distribution is non-CDS-like is likely not genuine translation. A called ORF with no frame-0 enrichment, a non-ribosomal FLOSS, and no conservation/peptide support should be treated as a candidate, not a finding.

## Differential ORF translation

**Goal:** Compare ORF-level translation across conditions.

**Approach:** Count offset-corrected P-sites per ORF per sample into an integer matrix, then run DESeq2.

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(orf_psite_counts, coldata, ~ condition)
dds <- DESeq(dds)
res <- results(dds)   # adjusted p-value is res$padj
```

Ribosome occupancy is not translation efficiency; a TE change needs an RNA-seq denominator (see translation-efficiency).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| RiboCode runs but uses wrong read lengths | `-l 27,28,29,30` passed as read lengths | `-l` is the longest-ORF toggle; read lengths come from `metaplots` config |
| KeyError on `ORF_type == 'noncoding'` | Wrong category names | RiboCode emits `Overlap_uORF/Overlap_dORF/Internal/novel`, not `noncoding` |
| `library(ORFik)` cannot find ORFquant functions | ORFik and ORFquant conflated | They are different packages; install ORFquant from its own repo |
| `fimport(p_offsets=, lengths=)` errors | Those arguments do not exist | Use `detectRibosomeShifts()` then `shiftFootprints()` |
| `RibORF.py -f -r -g -o` not found | Fabricated single-command CLI | RibORF is a Perl multi-script pipeline (logistic regression) |
| Most uORFs missing | ATG-only scan | Use a periodicity caller with `-A` near-cognate starts |
| High-coverage "ORF" is not real | Coverage used as the translation signal | Require frame-0 enrichment + FLOSS/ORFscore validation |

## Related Skills

- ribosome-periodicity - Calibrate the per-length P-site offsets ORF callers consume
- initiation-site-mapping - Map start codons (including non-AUG) from harringtonine/LTM data
- translation-efficiency - Add an RNA-seq denominator to turn occupancy into TE
- ribosome-stalling - Interpret pause sites within called ORFs
- differential-expression/deseq2-basics - Differential ORF-level translation

## References

- Mudge JM, Ruiz-Orera J, Prensner JR, et al. 2022. Standardized annotation of translated open reading frames. Nat Biotechnol 40(7):994-999. doi:10.1038/s41587-022-01369-0
- Xiao Z, Huang R, Xing X, Chen Y, Deng H, Yang X. 2018. De novo annotation and characterization of the translatome with ribosome profiling data. Nucleic Acids Res 46(10):e61. doi:10.1093/nar/gky179
- Calviello L, Hirsekorn A, Ohler U. 2020. Quantification of translation uncovers the functions of the alternative transcriptome. Nat Struct Mol Biol 27(8):717-725. doi:10.1038/s41594-020-0450-4
- Tjeldnes H, Labun K, Torres Cleuren Y, Chyżyńska K, Świrski M, Valen E. 2021. ORFik: a comprehensive R toolkit for the analysis of translation. BMC Bioinformatics 22:336. doi:10.1186/s12859-021-04254-w
- Choudhary S, Li W, Smith AD. 2020. Accurate detection of short and long active ORFs using Ribo-seq data. Bioinformatics 36(7):2053-2059. doi:10.1093/bioinformatics/btz878
- Ingolia NT, Brar GA, Stern-Ginossar N, et al. 2014. Ribosome profiling reveals pervasive translation outside of annotated protein-coding genes. Cell Rep 8(5):1365-1379. doi:10.1016/j.celrep.2014.07.045
- Bazzini AA, Johnstone TG, Christiano R, et al. 2014. Identification of small ORFs in vertebrates using ribosome footprinting and evolutionary conservation. EMBO J 33(9):981-993. doi:10.1002/embj.201488411
