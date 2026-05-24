---
name: bio-crispr-screens-crispresso-editing
description: Quantifies CRISPR editing outcomes with CRISPResso2 (Clement 2019 Nat Biotechnol) across Cas9-nuclease (indels, HDR), CBE and ABE base editors (target conversion + bystander), and prime editor (pegRNA-templated) modes. Covers single-amplicon (CRISPResso), multi-sample batch (CRISPRessoBatch), pooled-amplicon (CRISPRessoPooled), WGS off-target (CRISPRessoWGS), and sample-comparison (CRISPRessoCompare) workflows; quantification-window math that controls what is called edited; substitution-vs-indel diagnostic to distinguish BE from Cas9 contamination; MMEJ deletion pattern interpretation; allele-frequency tables; and failure modes from amplicon misalignment or contamination. Use when quantifying editing from amplicon sequencing, choosing CRISPResso mode by design, distinguishing intended edits from bystanders and indel byproducts, debugging low-alignment runs, or generating publication-grade editing reports.
tool_type: cli
primary_tool: CRISPResso2
---

## Version Compatibility

Reference examples tested with: CRISPResso2 2.2.14+ (pinellolab/CRISPResso2), pandas 2.2+, numpy 1.26+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `CRISPResso --version`; `CRISPRessoBatch --help`; `CRISPRessoPooled --help`; `CRISPRessoWGS --help`; `CRISPRessoCompare --help`
- Python: `from CRISPResso2 import ...`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## CRISPResso2 Editing Quantification

**"Quantify CRISPR editing from my amplicon sequencing"** -> Align amplicon reads against the reference, classify each read as unmodified / NHEJ / HDR / base-edited / prime-edited within the quantification window, and report per-edit-type frequencies, indel size distributions, allele-frequency tables, and substitution-position profiles.

- CLI: `CRISPResso` -- single amplicon, single sample
- CLI: `CRISPRessoBatch` -- multi-sample with per-sample parameters
- CLI: `CRISPRessoPooled` -- multi-amplicon pooled amplicon sequencing
- CLI: `CRISPRessoWGS` -- off-target quantification from whole-genome BAM
- CLI: `CRISPRessoCompare` -- pairwise outcome comparison (e.g., treated vs untreated)

## Mode Decision Tree

| Experimental design | Mode | Key parameters |
|---------------------|------|----------------|
| Single amplicon, single sample (e.g. pilot edit validation) | `CRISPResso` | `--amplicon_seq`, `--guide_seq` |
| Same amplicon, many samples (e.g. timecourse, dose response) | `CRISPRessoBatch` | `--batch_settings` table |
| Many amplicons, pooled in one library (e.g. arrayed validation pool) | `CRISPRessoPooled` | `--amplicon_file` |
| Off-target survey from whole-genome BAM | `CRISPRessoWGS` | `--bam`, `--reference`, `--regions_file` |
| Comparing two CRISPResso runs (e.g. condition A vs B) | `CRISPRessoCompare` | `--crispresso_output_folder_1`, `_2` |
| HDR / knock-in validation | `CRISPResso` with `--expected_hdr_amplicon_seq` | Same as base CRISPResso |
| Cytosine base editor (C->T) | `CRISPResso --base_editor_output` | `--conversion_nuc_from C --conversion_nuc_to T` |
| Adenine base editor (A->G) | `CRISPResso --base_editor_output` | `--conversion_nuc_from A --conversion_nuc_to G` |
| Prime editor (templated edit) | `CRISPResso` with pegRNA parameters | `--prime_editing_pegRNA_spacer_seq`, `--prime_editing_pegRNA_extension_seq`, `--prime_editing_pegRNA_scaffold_seq` |

**Fails when:**
- Pooled-amplicon mode applied to amplicons that share primer sequences -- reads get misassigned.
- Base editor mode without specifying `--conversion_nuc_from`/`--conversion_nuc_to` -- defaults assume CBE (C->T); ABE runs will misclassify.
- Prime editor mode without `--prime_editing_pegRNA_extension_seq` -- the RTT template is missing, no edit is detectable.

## The Quantification Window

**Why this matters for postdoc-level use:** CRISPResso classifies reads as "edited" or "unmodified" based on whether *modifications fall inside the quantification window* (not the whole amplicon). The window is centered on the predicted cut site (Cas9: 3 bp upstream of PAM; Cas12a: 18 bp downstream of PAM) with a default size of 1 (the cut site itself).

```bash
# Default Cas9 setup
--quantification_window_size 1                  # 1-bp window at cut site
--quantification_window_center -3               # 3 bp upstream of PAM

# Base editor: widen window to cover editing positions 4-8
--quantification_window_size 10                 # extend to capture position 4-8
--quantification_window_center -10              # center on the editing window
```

**Consequences of mis-sized window:**
- Too narrow: misses edits at HDR positions or far bystanders; underestimates editing
- Too wide: includes random sequencing errors; inflates editing rate
- Wrong center: edits at correct position are scored as outside the window

For base editing screens, the convention is `--quantification_window_size 10` to cover positions 4-13 from PAM-distal end. For prime editing with multi-base templated edits, widen to encompass the entire edit region.

## Single-Amplicon Cas9 Editing

**Goal:** Quantify indel frequencies and HDR efficiency from a single target site.

**Approach:** Align FASTQ reads to the reference and (optional) expected-HDR amplicon, classify each read, and report aggregated statistics.

```bash
CRISPResso \
    --fastq_r1 sample_R1.fastq.gz \
    --fastq_r2 sample_R2.fastq.gz \
    --amplicon_seq <amplicon_sequence_ref_genome> \
    --guide_seq <20nt_protospacer_no_PAM> \
    --expected_hdr_amplicon_seq <edited_amplicon_for_HDR> \  # OPTIONAL
    --quantification_window_size 1 \
    --quantification_window_center -3 \
    --min_average_read_quality 30 \                          # Phred quality filter
    --output_folder sample_results \
    --name sample_id

# Outputs:
#   sample_results/<name>/CRISPResso_mapping_statistics.txt
#   sample_results/<name>/CRISPResso_quantification_of_editing_frequency.txt
#   sample_results/<name>/Alleles_frequency_table.zip
#   sample_results/<name>/Indel_size_distribution.png
#   sample_results/<name>/Insertion_deletion_substitution.png
```

**Key outputs:**

| File | Content |
|------|---------|
| `CRISPResso_mapping_statistics.txt` | Reads aligned, reads in quantification window, % alignment, % discarded |
| `CRISPResso_quantification_of_editing_frequency.txt` | % unmodified, % NHEJ, % HDR (if expected), per-edit-class breakdown |
| `Alleles_frequency_table.zip` | Per-allele sequences and frequencies (allele-level resolution) |
| `Nucleotide_percentage_table.txt` | Per-position A/C/G/T/- frequencies (substitutions + deletions) |
| `Quantification_window_nucleotide_percentage_table.txt` | Same, restricted to quantification window (base-editor analysis) |
| `Reference_modified_chr11_115...png` | Sequencing pile-up over amplicon |

## Base Editor Quantification

**Goal:** Distinguish target base conversion from bystander edits and indel byproducts.

**Approach:** Run CRISPResso with `--base_editor_output` flag and specify the conversion direction; widen the quantification window to cover the editing window.

```bash
# Cytosine Base Editor (CBE): C->T conversion
CRISPResso \
    --fastq_r1 cbe_sample.fastq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_protospacer> \
    --base_editor_output \
    --conversion_nuc_from C \
    --conversion_nuc_to T \
    --quantification_window_size 10 \
    --quantification_window_center -10 \
    --output_folder cbe_results \
    --name cbe_sample

# Adenine Base Editor (ABE): A->G conversion
CRISPResso \
    --fastq_r1 abe_sample.fastq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_protospacer> \
    --base_editor_output \
    --conversion_nuc_from A \
    --conversion_nuc_to G \
    --quantification_window_size 10 \
    --quantification_window_center -10 \
    --output_folder abe_results \
    --name abe_sample
```

**Reading the output:**

| Metric | Where | Interpretation |
|--------|-------|----------------|
| Target editing % | `Quantification_window_nucleotide_percentage_table.txt`, target C/A row | Primary endpoint |
| Bystander editing % | Same table, other C/A positions in window | Off-target byproduct in window |
| Indel rate | `CRISPResso_quantification_of_editing_frequency.txt` | Cas9-like cut artifacts; should be <5% for clean BE |
| Substitution-vs-indel ratio | Derived | Ratio >10 indicates clean BE; <3 indicates cut-mediated mutagenesis instead |

**Critical:** Bystander editing is intrinsic to base editors (the deaminase acts across a 5-nt window); it is not noise. Report bystander rates alongside target rates. See [[base-editing-analysis]] for variant-call implications.

## Prime Editor Quantification

**Goal:** Quantify pegRNA-templated edits versus indel byproducts and partial edits.

**Approach:** Provide spacer, extension (PBS + RTT), and scaffold sequences; CRISPResso identifies reads matching the intended edit.

```bash
CRISPResso \
    --fastq_r1 pe_sample.fastq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_protospacer> \
    --prime_editing_pegRNA_spacer_seq <20nt_protospacer> \
    --prime_editing_pegRNA_extension_seq <RTT+PBS_sequence> \
    --prime_editing_pegRNA_scaffold_seq <scaffold_sequence> \
    --output_folder pe_results \
    --name pe_sample

# Output adds:
#   Prime_editing_outcomes.txt - intended-edit vs scaffold-incorporation vs indel vs unmodified
```

**Reading prime-editor output:**

| Metric | Interpretation |
|--------|----------------|
| Intended edit % | The pegRNA-encoded edit was correctly installed |
| Scaffold incorporation % | Reverse transcription read into scaffold instead of stopping at edit; failure mode |
| Indel % | Nick-only editing without templated repair; common at low-PE-activity sites |
| Unmodified % | Read matches the reference exactly |

A high-quality prime-edit run shows intended-edit fraction >5% and scaffold incorporation <2%. See [[prime-editing-screens]] for pegRNA design rules.

## Batch Mode (Multi-Sample, Same Amplicon)

**Goal:** Process tens to hundreds of samples with same amplicon design (e.g., a timecourse, dose response, or replicate panel).

**Approach:** Provide a tab-separated batch settings file with per-sample parameters; CRISPRessoBatch runs all in parallel.

```bash
# batch_settings.txt (tab-separated, headers required)
# name    fastq_r1                fastq_r2                amplicon_seq    guide_seq
# t0      t0_R1.fq.gz             t0_R2.fq.gz             ACGT...         GUIDE
# t6      t6_R1.fq.gz             t6_R2.fq.gz             ACGT...         GUIDE
# t12     t12_R1.fq.gz            t12_R2.fq.gz            ACGT...         GUIDE
# t24     t24_R1.fq.gz            t24_R2.fq.gz            ACGT...         GUIDE

CRISPRessoBatch \
    --batch_settings batch_settings.txt \
    --batch_output_folder batch_run \
    --skip_failed \
    --n_processes 8

# Outputs:
#   batch_run/CRISPRessoBatch_RUNNING_LOG.txt
#   batch_run/CRISPRessoBatch_quantification_of_editing_frequency.txt  (aggregated)
#   batch_run/CRISPResso_on_<name>/ for each sample
```

## Pooled-Amplicon Mode

**Goal:** Process multi-amplicon sequencing libraries (e.g., arrayed validation pools).

**Approach:** Provide an amplicon table with one row per target; CRISPRessoPooled de-multiplexes reads to the correct amplicon.

```bash
# amplicons.txt (tab-separated; header may vary by CRISPResso2 version)
# amplicon_name  amplicon_seq    guide_seq
# BRCA1_exon3    ACGT...         GUIDE1
# TP53_exon7     ACGT...         GUIDE2
# KRAS_codon12   ACGT...         GUIDE3

CRISPRessoPooled \
    --fastq_r1 pooled_R1.fastq.gz \
    --fastq_r2 pooled_R2.fastq.gz \
    --amplicon_file amplicons.txt \
    --output_folder pooled_run \
    --n_processes 8

# Outputs:
#   pooled_run/SAMPLES_QUANTIFICATION_SUMMARY.txt
#   pooled_run/CRISPResso_on_<amplicon>/ for each amplicon
```

**Failure mode:** Amplicons with shared primer regions get reads assigned to whichever amplicon comes first. Design primers with ≥3-bp distinguishing regions or use unique molecular identifiers.

## WGS Off-Target Mode

**Goal:** Quantify off-target editing from whole-genome sequencing.

**Approach:** Provide BAM file + reference + BED file of suspected off-target sites; CRISPResso extracts reads from each region and quantifies edits.

```bash
CRISPRessoWGS \
    --bam aligned.bam \
    --reference genome.fa \
    --regions_file off_targets.bed \
    --output_folder wgs_run \
    --n_processes 8
```

**Use case:** Validate empirically that an in vivo / clinical-grade edit has minimal off-target activity (combine with GUIDE-seq or CIRCLE-seq predicted sites).

## Parse Output in Python

**Goal:** Pull editing metrics into downstream analysis or reports.

**Approach:** Read the tab-separated quantification files and the JSON metadata.

```python
import pandas as pd
import json
from pathlib import Path

def parse_crispresso(output_dir):
    '''Extract key metrics from CRISPResso output directory.'''
    out = {}
    # Mapping statistics
    map_stats = {}
    with open(Path(output_dir) / 'CRISPResso_mapping_statistics.txt') as f:
        for line in f:
            k, v = line.strip().split('\t')
            map_stats[k] = v
    out['mapping_pct'] = float(map_stats.get('READS_ALIGNED_PERCENTAGE', 'nan'))
    out['reads_aligned'] = int(map_stats.get('READS_ALIGNED', '0'))
    # Editing quantification
    quant = pd.read_csv(Path(output_dir) / 'CRISPResso_quantification_of_editing_frequency.txt', sep='\t')
    out['editing_quant'] = quant.set_index('Amplicon').to_dict()
    # JSON metadata
    info_path = Path(output_dir) / 'CRISPResso2_info.json'
    if info_path.exists():
        out['info'] = json.loads(info_path.read_text())
    return out
```

## Failure Modes

### Low alignment rate (<50%)

**Trigger:** Wrong amplicon sequence (off by one nt, wrong strand, primer-trimmed vs untrimmed).
**Mechanism:** CRISPResso fails to align reads beyond the amplicon edges; discards as unmappable.
**Symptom:** `READS_ALIGNED_PERCENTAGE` <50%; per-position coverage drops at amplicon edges.
**Fix:** Re-derive amplicon from genome at primer-trimmed boundaries; verify strand orientation; check that primers are NOT included in `--amplicon_seq`.

### High substitution rate but low indel (Cas9 sample)

**Trigger:** Sample contamination with adjacent amplicon, primer-dimer, or sequencing error inflation.
**Mechanism:** Random substitutions inflate the per-position substitution rate without true indels.
**Symptom:** Substitutions >2% at base positions outside the cut site; alignment metrics look fine.
**Fix:** Increase `--min_average_read_quality` to 30+; filter contaminating amplicons; check primer-dimer in `CRISPResso_RUNNING_LOG.txt`.

### Bystander C/A editing inflates "editing efficiency"

**Trigger:** Base-editor sample with wide quantification window; bystander Cs at adjacent positions counted as edits.
**Mechanism:** Default `--quantification_window_size 10` includes all positions in editing window; bystander edits are real but distinct from target edit.
**Symptom:** Editing efficiency 80%+ but target SNV is 30%; bystander rate is 50%.
**Fix:** Always read the per-position table (`Quantification_window_nucleotide_percentage_table.txt`), not just the aggregate. Report target and bystander rates separately. See [[base-editing-analysis]].

### Prime editor sample with high scaffold incorporation

**Trigger:** RTT is too short relative to PBS, or pegRNA stops short.
**Mechanism:** Reverse transcriptase reads past the edit into scaffold sequence; product is detectable but undesired.
**Symptom:** Scaffold incorporation >5%; intended edit efficiency lower than expected.
**Fix:** Re-design pegRNA with longer RTT; verify with PRIDICT2 (see [[prime-editing-screens]]).

### MMEJ deletion misclassified as NHEJ

**Trigger:** Deletions with microhomology at junction; CRISPResso reports them as indels but doesn't distinguish MMEJ.
**Mechanism:** MMEJ creates predictable deletions using flanking microhomologies; biologically distinct from random NHEJ.
**Symptom:** Recurring same-size deletions in allele table (e.g., -7 bp deletion in 30% of reads).
**Fix:** Examine `Alleles_frequency_table` for over-represented allele patterns; flag MMEJ-mediated deletions for interpretation (these may be inferred from indel hotspots).

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Cas9 editing efficiency (functional KO) | >70% indels | Clement 2019 *Nat Biotechnol*; below this, KO is incomplete |
| Indel rate (clean base editor) | <5% | Clement 2019; >5% = unwanted cut activity |
| Target conversion (CBE) | >30% | Variable by target; below this, screen power is poor |
| Target conversion (ABE) | >30% | ABE typically lower per-base than CBE |
| Bystander rate (BE) | <10% acceptable; <5% ideal | Application-dependent; for variant function studies, must be controlled |
| Intended-edit % (prime editor) | >5% per-edit | Anzalone 2019 baseline; can be 50%+ at favorable sites |
| Scaffold incorporation (PE) | <2% | High-quality pegRNA design |
| Alignment rate | >85% | Below this, amplicon design or contamination issue |
| Minimum read quality | Phred 30 | Joung 2017 sequencing standard |
| Quantification window size (Cas9) | 1 | Clement 2019 default; precise cut-site analysis |
| Quantification window size (BE) | 10 | Cover editing window positions 4-13 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Alignment rate <50% | Wrong amplicon sequence | Re-verify; primers should NOT be in amplicon_seq |
| All reads "modified" | Misaligned reference | Check amplicon strand; reverse-complement test |
| BE shows mostly indels | Cas9 contamination or wrong protein | Re-derive cell line origin; check Cas9 vs nCas9-BE3 |
| Inconsistent batch results | Different amplicon_seq per sample | Use CRISPRessoBatch with consistent amplicon |
| Pooled-amplicon misassignment | Primer overlap between amplicons | Re-design with ≥3-bp distinguishing regions |
| Out-of-window edits ignored | Window too narrow | Increase `--quantification_window_size` |
| Scaffold incorporation high (PE) | RTT too short | Re-design pegRNA |
| Allele frequency dominated by 1 read | Low input / clonal | Verify input cell count; rerun if singleton |

## References

- Clement K et al. 2019. *Nat Biotechnol* 37:224. CRISPResso2 algorithm and modes.
- Pinello L et al. 2016. *Nat Biotechnol* 34:695. Original CRISPResso.
- Anzalone AV et al. 2019. *Nature* 576:149. Prime editing (PE-1/PE-2/PE-3).
- Komor AC et al. 2016. *Nature* 533:420. Base editing (BE3).
- Sanson KR et al. 2020. *Nat Commun* 11:5165. GRACE base-editor screen library.
- Findlay GM et al. 2018. *Nature* 562:217. Saturation genome editing.

## Related Skills

- crispr-screens/base-editing-analysis - Variant-function analysis using CRISPResso2 BE output
- crispr-screens/prime-editing-screens - PRIDICT2 pegRNA design + PE-tiling
- crispr-screens/library-design - sgRNA / pegRNA design for editing screens
- crispr-screens/screen-qc - Editing-efficiency QC for variant interpretation
- variant-calling/variant-annotation - Annotate detected variants downstream
- read-alignment/bwa-alignment - For WGS off-target alignment input
