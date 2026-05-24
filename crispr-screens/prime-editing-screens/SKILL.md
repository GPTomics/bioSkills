---
name: bio-crispr-screens-prime-editing-screens
description: Designs and analyzes pooled prime-editor (PE) screens for installing precise genetic variants without bystander confounding. Covers pegRNA design with PRIDICT and PRIDICT2 (Mathis 2023/2024) for predicting per-pegRNA editing efficiency, pegRNA architecture (spacer + scaffold + PBS + RTT), PE2 / PE3 / PE3b / PEmax / PEAR variants, MOSAIC in situ saturation mutagenesis (Hsu JY et al 2024 bioRxiv), the PRIME pooled-screen methodology (Erwood/Doman 2023 Nat Biotechnol 41:885; ~3,699 ClinVar variant screens), chromatin context as a primary determinant of PE efficiency, scaffold-incorporation and indel byproduct quantification with CRISPResso2, and the cross-modal validation strategy of PE + base-editor screens for variant function. Use when designing a pegRNA library for variant installation, choosing between BE and PE for a specific edit, predicting pegRNA efficiency before library synthesis, analyzing PE screen output, distinguishing intended-edit from scaffold-incorporation, or scaling PE screens to thousands of variants.
tool_type: mixed
primary_tool: PRIDICT2
---

## Version Compatibility

Reference examples tested with: PRIDICT2 v1.0+ (https://github.com/uzh-dqbm-cmi/PRIDICT2), CRISPResso2 2.2.14+, pandas 2.2+, biopython 1.83+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `python pridict2_pegRNA_design.py single --help`; `python pridict2_pegRNA_design.py batch --help`
- Web: PRIDICT2 web interface at https://pridict.it/

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Prime-Editing Screen Analysis

**"Design or analyze a pooled prime-editor screen"** -> Design pegRNAs (spacer + scaffold + PBS + RTT) for intended edits, predict efficiency with PRIDICT2, filter pre-synthesis to efficient candidates, install variants in the screen, quantify intended-edit vs scaffold-incorporation vs indel via CRISPResso2, and aggregate to per-variant fitness scores.

- Python: `PRIDICT2` for pegRNA efficiency prediction
- CLI: `CRISPResso --prime_editing_pegRNA_*` for amplicon-level analysis
- Workflow: pegRNA library design -> PRIDICT2 filtering -> screen execution -> CRISPResso2 quantification -> per-variant scoring

## Prime Editor Chemistry Comparison

| Editor | Year | Mechanism | Indel rate | Use when |
|--------|------|-----------|------------|----------|
| PE2 (Anzalone 2019) | 2019 | nCas9-RT fusion + pegRNA | 1-3% | Standard PE; lowest indel rate |
| PE3 | 2019 | PE2 + nick of opposite strand by additional sgRNA | 2-5% | Higher editing efficiency, slightly more indels |
| PE3b | 2019 | PE3 with edit-blocking ssgRNA | 1-3% | When PE3's added nick risks unwanted indels |
| PEmax (Chen 2021) | 2021 | Engineered RT + nCas9 | 1-2% | Higher editing rate per pegRNA |
| PEAR (Erwood 2023) | 2023 | PE with optimal pegRNA scaffold | 1-2% | Improved PE scaffold |
| PE5max (Chen 2021) | 2021 | PEmax with engineered scaffold variants | 1% | Highest efficiency at favorable sites |
| Dual-pegRNA / PE6 (2024) | 2024 | Twin pegRNA system | Variable | Specific applications |

**Decision rule:** For pooled screens at scale, PE2 or PEmax (less RAM-intensive in cells) is preferred over PE3 (additional sgRNA complicates library architecture). For specific high-efficiency edits, PEmax + PRIDICT2-optimized pegRNA.

## pegRNA Architecture

A pegRNA contains four critical elements that determine efficiency:

```
5'  SPACER (20 nt)  -- standard sgRNA spacer; defines target locus via NGG PAM
    +
    SCAFFOLD (~80 nt) -- canonical or engineered scaffold (Chen 2021 has improved scaffold)
    +
    PBS (Primer Binding Site, 8-15 nt) -- complements protospacer downstream of cut site
    +
    RTT (Reverse Transcription Template, 10-30 nt) -- encodes intended edit; copied by RT
3'
```

**Key design parameters:**
- **PBS length:** 11-13 nt typical; longer for high-GC contexts; PBS GC fraction critical (35-65% target)
- **RTT length:** 10-20 nt typical; longer for distant edits (10+ bp away from cut)
- **RTT-edit position:** intended edit at position 4-30 from cut site
- **Scaffold:** standard sgRNA scaffold OR Chen 2021 engineered scaffold (5-10% higher editing)

## PRIDICT and PRIDICT2 pegRNA Efficiency Prediction

**Mathis N et al 2023 *Nat Biotechnol* 41:1151 (PRIDICT v1) / 2025 *Nat Biotechnol* 43(5):712 (PRIDICT2; published online June 2024)** developed deep-learning predictors of per-pegRNA editing efficiency. PRIDICT2 is the current state of the art.

```bash
# PRIDICT2 is invoked via CLI: pridict2_pegRNA_design.py
# Single sequence input:
python pridict2_pegRNA_design.py single \
    --sequence-name BRCA1_c5135 \
    --sequence "AGCAGCCT(C/T)CTGAATGCCC...60nt_context" \    # parens = intended edit
    --output-dir predictions/ \
    --use_5folds                                              # 5-fold ensemble averaging

# Batch input from CSV:
python pridict2_pegRNA_design.py batch \
    --input-fname variants_to_design.csv \                    # CSV: sequence_name, sequence
    --output-dir predictions/ \
    --cores 4 \
    --summarize                                               # generate summary table

# Output: per-pegRNA predictions in predictions/<sequence_name>/
# Columns: PBS_sequence, PBS_length, RTT_sequence, RTT_length, predicted_editing_efficiency,
#          predicted_indel_rate, deep_ensemble_score, etc.
```

**Loading PRIDICT2 results in Python:**

```python
import pandas as pd
from pathlib import Path

def load_pridict2_predictions(prediction_dir):
    '''Load PRIDICT2 batch outputs from prediction_dir/'''
    summary = pd.read_csv(Path(prediction_dir) / 'pridict2_summary.csv')
    # summary has columns: sequence_name, PBS, RTT, predicted_efficiency, predicted_indel, etc.
    return summary
```

**Key determinants of PE efficiency (Mathis 2024 PRIDICT2):**

| Feature | Effect on efficiency |
|---------|----------------------|
| PBS GC content | 40-55% optimal; high GC slows annealing |
| PBS length | 11-13 nt optimal; longer for high-GC PBS |
| RTT length | 10-20 nt typical; trade-off between coverage and processivity |
| Edit position in RTT | Closest to PBS = highest efficiency |
| Chromatin context | Open chromatin = 2-5x higher efficiency than closed |
| Cell line / Cas9 expression | Variable; piloting required |
| Cell cycle phase | S/G2 = higher efficiency |

**Critical insight from Mathis 2024:** Chromatin context is the dominant determinant. Sequence-based predictions like PRIDICT under-predict at silenced loci and over-predict at open chromatin. For genome-scale screens, validate predictions empirically at representative loci.

## PRIME Pooled Screen Methodology

**Erwood S, Doman JL et al 2023 *Nat Biotechnol* 41:885** established the PRIME pooled-screen methodology (earlier 2022 bioRxiv preprint):

- pegRNA library covering thousands of intended variants
- Filter pegRNAs to PRIDICT2 efficiency >50% (or pilot top 25%)
- Lentiviral delivery at standard MOI 0.3 in PE-expressing cell line
- Selection on integration marker
- Time-course screen for variant function (e.g., drug sensitivity)
- Endpoint amplicon sequencing of each pegRNA target locus
- CRISPResso2 quantification of intended-edit %
- MAGeCK / drugZ-style hit calling on edit-efficient pegRNAs

**Quantified scale:** ~3,699 ClinVar variants installed in a single PRIME screen (Erwood/Doman 2023 *Nat Biotechnol* 41:885), with editing efficiency >5% at >50% of pegRNAs (validation cohort).

## MOSAIC In Situ Saturation Mutagenesis

**MOSAIC (Hsu JY, Lam KC, Shih J, Pinello L, Joung JK 2024 bioRxiv 10.1101/2024.04.25.591078)** is a higher-throughput variant of PRIME with multiplexed read-out:

- Tile pegRNAs across protein domains for systematic mutagenesis
- Saturation: every possible amino acid change in a region
- Identify drug-resistance variants in real-time
- Smaller per-variant cell numbers (more variants total)

**Use case:** Cancer-drug-resistance variant scanning; protein-domain function mapping.

## Run PRIDICT2 on a Custom pegRNA Library

**Goal:** Predict editing efficiency for thousands of pegRNAs before library synthesis.

**Approach:** Build a CSV with one row per intended edit (sequence + edit notation), run PRIDICT2 in batch mode, parse the per-pegRNA efficiency summary, and filter to candidates above the chosen efficiency threshold.

```bash
# Step 1: prepare batch input CSV (sequence_name, sequence with (REF/ALT) edit notation)
cat > variants.csv <<EOF
sequence_name,sequence
BRCA1_R71X,AGCAGCCT(C/T)CTGAATGCCC...
MLH1_c677,GAGCTGAGC(A/G)GAGGCTCTTGAAGC...
EOF

# Step 2: run PRIDICT2 batch
python pridict2_pegRNA_design.py batch \
    --input-fname variants.csv \
    --output-dir predictions/ \
    --cores 8 \
    --summarize
```

```python
# Step 3: parse and filter
import pandas as pd
predictions = pd.read_csv('predictions/pridict2_summary.csv')

# Filter to pegRNAs with predicted efficiency > 50% (Mathis 2024 threshold)
filtered = predictions[predictions['predicted_editing_efficiency'] > 50]
print(f'pegRNAs passing PRIDICT2 >50%: {len(filtered)} / {len(predictions)}')

# Pick top 3 per intended edit
top3 = (filtered.sort_values(['sequence_name', 'predicted_editing_efficiency'],
                              ascending=[True, False])
                 .groupby('sequence_name').head(3))
top3.to_csv('peg_library_filtered.csv', index=False)
```

## Cross-Validate PE with Base Editor Screens

**Goal:** Confirm variant-function calls from PE with orthogonal BE screens.

**Approach:** Design parallel BE library for the same variants; run both screens; intersect hits.

```python
# BE screen output (target conversion + bystander)
be_hits = pd.read_csv('be_screen_hits.tsv', sep='\t')
# PE screen output (intended edit + scaffold-incorp + indel)
pe_hits = pd.read_csv('pe_screen_hits.tsv', sep='\t')

# Intersect on intended variant
concordant = be_hits.merge(pe_hits, on='variant_id', suffixes=('_be', '_pe'))
# Filter to high-confidence: both methods call variant + same direction
concordant['high_confidence'] = (concordant['be_fdr'] < 0.05) & (concordant['pe_fdr'] < 0.05) & \
                                 (np.sign(concordant['be_lfc']) == np.sign(concordant['pe_lfc']))
```

**Critical:** PE-only hits in BE-coverable variants are suspect (BE should detect them). PE-only hits in non-BE-coverable variants (e.g., transversions) are genuinely PE-unique.

## CRISPResso2 for PE Quantification

```bash
CRISPResso \
    --fastq_r1 pe_sample.fq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_spacer> \
    --prime_editing_pegRNA_spacer_seq <spacer> \
    --prime_editing_pegRNA_extension_seq <PBS+RTT> \
    --prime_editing_pegRNA_scaffold_seq <scaffold> \
    --quantification_window_size 25 \              # widen to cover edit
    --output_folder pe_results \
    --name sample_id

# Output: Prime_editing_outcomes.txt
# Columns: intended_edit_pct, scaffold_incorp_pct, indel_pct, unmodified_pct
```

## Failure Modes

### Low pegRNA efficiency despite high PRIDICT prediction

**Trigger:** Sequence-only prediction missed chromatin context.
**Mechanism:** Closed chromatin reduces Cas9 binding and RT activity; PRIDICT2 only sees sequence.
**Symptom:** PRIDICT2 predicts 60% efficiency; observed is 5%.
**Fix:** Cross-reference target with chromatin accessibility data (ATAC-seq) in your cell line; flag pegRNAs at silenced loci; pilot before screen.

### High scaffold incorporation

**Trigger:** RTT too short relative to PBS, or RT processivity issue.
**Mechanism:** RT reads past edit into scaffold; resulting product is detectable but undesired.
**Symptom:** Scaffold incorporation >5%; intended edit efficiency low.
**Fix:** Re-design pegRNA with longer RTT; verify with PRIDICT2 score for scaffold_incorp; pilot at representative loci.

### PE2 cell line lacks RT expression

**Trigger:** PE2 construct expressed at low level; insufficient RT for productive editing.
**Mechanism:** PE2 requires high RT expression; some cell lines down-regulate.
**Symptom:** Library-wide editing <10%; not locus-specific.
**Fix:** Verify PE2 expression by Western blot; consider PEmax (higher activity); use better-validated cell lines (K562, HEK293T, U2OS).

### Multi-base intended edit but only one base installed

**Trigger:** Long RTT designed for multi-base edit; RT prematurely terminates.
**Mechanism:** RT processivity drops with longer RTT; multi-base edits often incomplete.
**Symptom:** Allele table shows partial-edit alleles (some bases installed, not all).
**Fix:** Re-design with shorter RTT covering only the closest edits; or use PE3 to nick opposite strand and force longer RT processivity.

### Library missing intended variant

**Trigger:** No suitable PAM/PBS/RTT combination for the intended edit.
**Mechanism:** PE requires NGG PAM within 30 nt of edit; rare edits cannot be installed.
**Symptom:** Specific variants absent from library.
**Fix:** Use SpRY-PE for relaxed PAM; accept that some variants cannot be PE-installed; consider BE if applicable.

## Cas9 vs BE vs PE for Variant Installation

| Approach | Bystander | Indels | Coverage | When to use |
|----------|-----------|--------|----------|-------------|
| Cas9 + HDR | None | High | Variable (depends on template integration) | Precise edits at scale; high indel byproduct |
| Base editor | YES | Low (<5%) | Limited by editing window | C->T or A->G at editable position |
| Prime editor | NONE | Low (<3%) | NGG-PAM within 30 nt of edit | Precise variants; multi-base; transversions |
| Cas9 (no template) | NONE | 70%+ | Anywhere with NGG | LoF only; not variant-specific |

**Decision tree:**
- C->T or A->G at editing-window position: BE (higher efficiency than PE)
- Multi-base / transversion / out-of-window: PE
- LoF without specifying variant: Cas9
- Random insertions: HDR (lower throughput than PE)

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| PRIDICT2 efficiency for library inclusion | >50% | Mathis 2024 |
| Intended edit % for screen power | >5% (per Anzalone 2019); >20% at favorable sites | Anzalone 2019 |
| Scaffold incorporation | <2% (clean PE); <5% acceptable | Empirical |
| Indel byproduct | <3% (PE2); <5% (PE3) | Anzalone 2019; Chen 2021 |
| PBS GC content | 40-55% | PRIDICT2 |
| PBS length | 11-13 nt | PRIDICT2 |
| RTT length | 10-20 nt | PRIDICT2 |
| Edit position from cut | 1-30 nt | Anzalone 2019 |
| Cell line for PE | K562, HEK293T, U2OS validated | High RT expression |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low editing across library | Cell-line RT inactivity | Verify PE2 expression; switch to validated line |
| Scaffold incorporation >10% | RTT too short | Re-design with longer RTT |
| Partial multi-base edits | RT processivity limit | Shorter RTT or PE3 |
| PRIDICT predicts but observes much lower | Chromatin context | Pilot at chromatin-aware sites |
| Library missing variants | No NGG PAM | SpRY-PE; BE alternative |
| PE concordant with BE on transitions, disagrees on transversions | PE handles transversions BE doesn't | Expected; trust PE |

## References

- Anzalone AV et al. 2019. *Nature* 576:149. Original PE2/PE3 (foundational prime editing paper).
- Mathis N et al. 2023. *Nat Biotechnol* 41:1151. PRIDICT v1 deep-learning pegRNA prediction.
- Mathis N et al. 2025. *Nat Biotechnol* 43(5):712 (published online June 2024). PRIDICT2 + chromatin context (current state-of-the-art).
- Chen PJ et al. 2021. *Cell* 184:5635. PEmax + engineered RT.
- Erwood S et al. 2023. *Nature Biotechnology* 41:885. PEAR pegRNA scaffold design.
- Hsu JY, Lam KC, Shih J, Pinello L, Joung JK 2024 bioRxiv (doi:10.1101/2024.04.25.591078). MOSAIC in situ saturation mutagenesis via prime editing.
- Erwood S, Doman JL et al. 2023. *Nature Biotechnology* 41:885. PRIME pooled-screen methodology (~3,699 ClinVar variants); cross-reference for variant-installation scale.

## Related Skills

- crispr-screens/library-design - pegRNA library design
- crispr-screens/base-editing-analysis - Orthogonal BE for variant attribution
- crispr-screens/crispresso-editing - CRISPResso2 PE mode and quantification
- crispr-screens/hit-calling - Per-variant hit calling
- crispr-screens/screen-qc - Editing-efficiency QC
- variant-calling/variant-annotation - Annotate edited variants
- clinical-databases/clinvar-lookup - Variant pathogenicity
