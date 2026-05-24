---
name: bio-crispr-screens-base-editing-analysis
description: Analyzes base-editing screens for variant function. Covers library design (Sanson 2020 GRACE, Hanna 2021 BRCA1/2 SNV scanning, Cuella-Martin 2021), CBE vs ABE chemistry choice (BE3/BE4 vs ABE7.10/ABE8.20/ABE8e), editing-window math (positions 4-8 from PAM-distal end, wider for ABE8e), bystander-edit quantification and the variant-call ambiguity it creates, sgRNA-efficiency filtering before hit calling, indel byproduct interpretation, the substitution-vs-indel diagnostic, variant annotation against ClinVar / COSMIC, and the Broad be-validation-pipeline. Use when designing a BE variant screen, choosing CBE vs ABE for a specific edit, interpreting bystander-confounded hits, distinguishing functional signal from indel artifact, integrating CRISPResso2 output with screen scoring, or deciding BE vs PE for SNV installation.
tool_type: mixed
primary_tool: CRISPResso2
---

## Version Compatibility

Reference examples tested with: CRISPResso2 2.2.14+, BE-Hive 1.0+ (BE prediction), pandas 2.2+, biopython 1.83+, numpy 1.26+, scipy 1.12+, scikit-learn 1.4+, Broad be-validation-pipeline 1.0+ (Python).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `CRISPResso --version`; `be-validation-pipeline --help`
- Python: `pip show CRISPResso2 be-hive`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Base Editing Screen Analysis

**"Analyze my base-editor variant-function screen"** -> Quantify per-sgRNA target-base conversion, bystander rate, and indel byproducts from amplicon sequencing; filter on editing efficiency; map each sgRNA to its intended SNV (target + bystander pattern); compute per-variant fitness from the screen log-fold change; reconcile target vs bystander variant attribution; annotate against ClinVar / COSMIC.

- CLI: `CRISPResso --base_editor_output` for per-amplicon BE quantification
- CLI: Broad `be-validation-pipeline` for end-to-end pooled-screen analysis with editing-efficiency filtering
- Python: `BE-Hive` (Arbab 2020) for editing-efficiency prediction
- Python: `BE-Designer` (Hwang 2019) for variant-encoding sgRNA design

## Base Editor Chemistry Selection

| Editor | Reaction | Editing window | Indel byproduct rate | When to use |
|--------|----------|----------------|----------------------|-------------|
| BE3 (Komor 2016) | C->T (also G->A on opposite strand) | Pos 4-8 from PAM-distal end | 5-10% | Original; superseded |
| BE4 / BE4max (Koblan 2018) | C->T | Pos 4-8 | <5% | CBE standard |
| eA3A-BE3 | C->T narrow specificity | Pos 5-7 | <5% | Specifically TC contexts (eA3A prefers TC) |
| ABE7.10 (Gaudelli 2017) | A->G (T->C opposite strand) | Pos 4-7 | <2% | First ABE; slow at non-TA contexts |
| ABE8.20 (Richter 2020) | A->G | Pos 4-8 | <2% | Modern ABE; high activity |
| ABE8e (Lapinaite 2020) | A->G | Pos 4-8 | <2% | Highest editing activity; broader window |
| evoCDA-BE | C->T (broader) | Pos 1-9 | 5-10% | Larger editing window; more bystander |
| CGBE1 (Kurt 2021) | C->G | Pos 5-7 | 5-10% | C-to-G transversion; rare use |
| GBE (Zhao 2021) | C->G or C->A | Pos 4-7 | 5-10% | Transversions; less mature |

**Decision rule:** For a target SNV at position 4-8 of a candidate spacer with no bystander Cs/As in the same window, BE3-BE4 or ABE7.10 is sufficient. For high-throughput variant scanning where bystander tolerance must be minimized, use eA3A-BE3 (TC contexts only) or ABE8e (narrower effective window).

## Editing Window Math

**Why this matters for postdoc-level use:** Base editors are tethered to dCas9 (or nCas9) and the deaminase acts on the displaced ssDNA "R-loop" formed when Cas9 binds. The deaminase has a fixed reach -- positions 4-8 from the PAM-distal end of the protospacer for canonical BE3/BE4/ABE7.10. Outside this window, editing efficiency drops by 10-50x.

```
PAM-distal end                                                            PAM-proximal
   |                                                                          |
   1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20    NGG
                  ^^^^^^^^^^^
                  Canonical editing window (positions 4-8)

   For BE4max / ABE7.10: positions 4-8 are 5-50x more efficient than positions 1-3 or 9-13
   For ABE8e: window extends to positions 4-10 due to enhanced TadA8e activity
   For evoCDA-BE: window 1-9 (broader; more bystander)
```

**Critical implication for variant interpretation:** If the intended edit is at position 5 and there is an additional editable C/A at position 7, both will be edited in the same molecule. The screen scores the *combination* of edits, not the intended one alone. This is bystander confounding.

## sgRNA Library Design for BE Screens

**Goal:** Tile editing-window-positioned spacers across a protein region of interest to enable variant scanning.

**Approach:** For each amino acid in the target region, find NGG-adjacent spacers where the SNV-of-interest base falls in editing positions 4-8 with minimal bystander C/A in the same window. Annotate each spacer with the predicted amino acid changes (target + bystander).

```python
import pandas as pd
import re
from Bio.Seq import Seq

def find_be_spacers(cds_sequence, cds_protein_start, target_aa, target_base='C', editor='BE4max'):
    '''Find sgRNAs that place target_base in editor-specific window at target_aa.
    Returns spacers with bystander annotation.

    Args:
        cds_sequence: nucleotide CDS (translated frame 1)
        cds_protein_start: amino acid number of CDS start (usually 1)
        target_aa: amino acid number to install variant (e.g., 130 for residue 130)
        target_base: 'C' (CBE) or 'A' (ABE)
        editor: 'BE3', 'BE4max', 'eA3A-BE3', 'ABE7.10', 'ABE8.20', 'ABE8e', 'evoCDA-BE'

    Returns: DataFrame with spacer, position-in-cds, target-base-position-in-spacer,
             bystander_positions, predicted_aa_changes
    '''
    # Editor-specific editing window (positions from PAM-distal end of spacer)
    window_by_editor = {
        'BE3': (4, 8),       'BE4max': (4, 8),    'eA3A-BE3': (5, 7),
        'ABE7.10': (4, 7),   'ABE8.20': (4, 8),   'ABE8e': (4, 10),    # ABE8e wider!
        'evoCDA-BE': (1, 9),
    }
    window_lo, window_hi = window_by_editor[editor]
    aa_index = target_aa - cds_protein_start  # 0-indexed in protein
    aa_start_nt = aa_index * 3                # nt offset in cds
    candidates = []
    spacer_len = 20
    pam_pattern = re.compile(r'(?=([ACGT]GG))')
    for strand, seq in [('+', cds_sequence), ('-', str(Seq(cds_sequence).reverse_complement()))]:
        for pam_match in pam_pattern.finditer(seq):
            pam_pos = pam_match.start()
            spacer_start = pam_pos - spacer_len
            if spacer_start < 0:
                continue
            spacer = seq[spacer_start:pam_pos]
            # Editor-specific window from PAM-distal end (1-indexed)
            # Find all editable bases in window
            edit_bases_in_window = []
            for i, b in enumerate(spacer[window_lo-1:window_hi], start=window_lo):
                if b == target_base:
                    edit_bases_in_window.append(i)
            if not edit_bases_in_window:
                continue
            # Annotate which edits hit the target_aa codon
            target_codon_start = aa_start_nt
            target_codon_end = target_codon_start + 3
            target_position_in_spacer = []
            for i in edit_bases_in_window:
                genomic_pos = spacer_start + i - 1
                if target_codon_start <= genomic_pos < target_codon_end:
                    target_position_in_spacer.append(i)
            bystander_positions = [i for i in edit_bases_in_window if i not in target_position_in_spacer]
            candidates.append({
                'spacer': spacer,
                'strand': strand,
                'spacer_start': spacer_start,
                'target_positions': target_position_in_spacer,
                'bystander_positions': bystander_positions,
                'n_bystanders': len(bystander_positions),
            })
    return pd.DataFrame(candidates).sort_values('n_bystanders')
```

**Decision rule:** Select spacers with target_positions != empty AND n_bystanders minimized. For variant-by-variant scanning, accept up to 1-2 bystanders if biology of those positions is interpretable; flag for downstream variant attribution.

## Editing Efficiency Filtering (Critical Pre-Hit-Calling)

**Goal:** Drop sgRNAs that do not edit efficiently, since unedited reads represent no biological perturbation.

**Approach:** From CRISPResso2 output, compute target-base-conversion percentage per sgRNA; filter library to sgRNAs with >50% target editing in a pilot or co-screened control.

```python
def filter_by_editing_efficiency(crispresso_outputs_dir, target_pos, target_base, efficiency_threshold=0.5):
    '''Drop sgRNAs that edit <efficiency_threshold of reads at target position.
    crispresso_outputs_dir: directory containing CRISPResso per-sample outputs.'''
    from pathlib import Path
    results = []
    for sample_dir in Path(crispresso_outputs_dir).glob('CRISPResso_on_*'):
        sgrna_id = sample_dir.name.replace('CRISPResso_on_', '')
        quant_file = sample_dir / 'Quantification_window_nucleotide_percentage_table.txt'
        if not quant_file.exists():
            continue
        df = pd.read_csv(quant_file, sep='\t')
        # Find target position in the quantification window
        target_row = df[df['Position'] == target_pos]
        if target_row.empty:
            continue
        # Editing = sum of non-original bases at target position
        original_pct = target_row[target_base].values[0]
        editing_pct = (100 - original_pct) / 100
        results.append({'sgrna_id': sgrna_id, 'editing_pct': editing_pct,
                         'pass_filter': editing_pct >= efficiency_threshold})
    return pd.DataFrame(results)
```

**Convention:** Drop sgRNAs below 50% editing for variant-function screens. Hanna 2021 used a 30% threshold for primary screen and a 50% threshold for confirmed hits. Below 30%, the screen has insufficient power; above 70%, results approach saturation editing.

## Bystander Edit Attribution

**Why this matters:** When a sgRNA's editing window contains the target base AND a bystander base, the screen scores the combination. To attribute screen signal to the target variant alone, you must either (a) include sgRNAs that edit only the target (no bystander) -- often impossible -- or (b) deconvolute via parallel measurements.

**Strategies for variant-by-variant attribution:**

1. **Tile multiple sgRNAs with different bystander patterns:** If 5 different sgRNAs all hit the target base but have different bystanders, common signal across them is target-attributable (Hanna 2021 approach).

2. **Use orthogonal chemistry:** Run the same variant scan with prime editor (no bystanders); cross-validate. See [[prime-editing-screens]].

3. **Bystander stratification:** From CRISPResso2 allele table, partition reads by exact edit pattern (target only, target+bystander_1, target+bystander_2, etc.); separately score each pattern's contribution to the phenotype.

4. **Restrict library:** Use only sgRNAs with zero bystanders in the editing window (rare; may exclude most candidate spacers).

```python
def deconvolute_bystander(allele_table_path, target_pos, bystander_pos_list):
    '''From CRISPResso2 allele table, partition reads by edit pattern at target + bystanders.
    Returns: per-pattern frequency for each combination of target/bystander edits.'''
    alleles = pd.read_csv(allele_table_path, sep='\t', compression='zip')
    # Mark target_edited and per-bystander_edited
    alleles['target_edited'] = alleles['Aligned_Sequence'].str[target_pos-1] != alleles['Reference_Sequence'].str[target_pos-1]
    for bp in bystander_pos_list:
        alleles[f'bystander_{bp}_edited'] = alleles['Aligned_Sequence'].str[bp-1] != alleles['Reference_Sequence'].str[bp-1]
    return alleles.groupby(['target_edited'] + [f'bystander_{bp}_edited' for bp in bystander_pos_list])['Reference_pct'].sum().reset_index()
```

## Hit Calling for Variant-Function Screens

**Goal:** Score per-variant fitness from a base-editor screen.

**Approach:** Filter library to efficiency-passing sgRNAs (>50% editing), then run MAGeCK MLE or drugZ on the sgRNA-level counts; map each significant sgRNA to its predicted variant + bystander pattern; aggregate to per-variant scores.

```python
def aggregate_variant_scores(mageck_sgrna_summary, variant_annotation_df):
    '''Aggregate sgRNA-level scores to per-variant scores.
    variant_annotation_df: per-sgRNA -> predicted variants (target + bystanders).'''
    df = mageck_sgrna_summary.merge(variant_annotation_df, on='sgRNA')
    # Target-only contribution: sgRNAs with no bystanders
    target_only = df[df['n_bystanders'] == 0]
    target_only_scores = target_only.groupby('target_variant')['LFC'].agg(['mean', 'std', 'count'])
    # Mixed signal: sgRNAs with bystanders
    mixed = df[df['n_bystanders'] > 0]
    return target_only_scores, mixed
```

## Hanna 2021 BRCA1/2 Variant-Function Screen Methodology

**Hanna et al 2021 *Cell* 184:1066** established the gold-standard methodology for BE variant scanning:

1. Design CBE library tiling BRCA1 and BRCA2 with 10-15 sgRNAs per amino acid
2. Verify editing efficiency in a control timepoint via amplicon sequencing
3. Drop sgRNAs <30% target editing
4. Run drug-modifier screens (PARPi sensitivity) with vehicle vs drug
5. Use drugZ to identify sensitizing variants (loss of function)
6. Per-variant scoring: aggregate over all sgRNAs that hit that variant; cross-check against bystander-controlled sgRNAs
7. Cross-validate hits via prime-editor scans of the same variants

**Quantified result:** Identified novel loss-of-function variants in BRCA1 RING and BRCT domains, MCL1 and BCL2L1, and PARP1 resistance variants. The approach validated for clinical variant interpretation.

## Cuella-Martin 2021 BRCA1/2 Variant Screening

**Cuella-Martin et al 2021 *Mol Cell* 81:1107** corroborated Hanna's findings:

- Independent CBE design targeting BRCA1/2
- Identified pathogenic/likely-pathogenic variants in critical protein domains
- Combined with proteomic validation
- Demonstrated saturation mutagenesis is feasible at protein-domain scale

**Reconciliation with Hanna 2021:** Both studies converged on similar functional variant calls in BRCA1 RING and BRCT domains; the orthogonal methodology gave confidence in clinical interpretation.

## Cas9 vs Base Editor vs Prime Editor for Variant Installation

| Approach | What it does | Bystander | Indels | When to use |
|----------|--------------|-----------|--------|-------------|
| Cas9 + HDR template | Installs precise edit + template | None | High (NHEJ competition) | When precise edit needed; high indel byproduct |
| Cas9 (no template) | Random indels at cut site | None | 70%+ | Loss-of-function; not variant-specific |
| CBE (BE3/BE4) | C->T at editing window | Yes (multiple Cs) | <5% | C->T variants with manageable bystanders |
| ABE (ABE7.10/ABE8e) | A->G at editing window | Yes (multiple As) | <2% | A->G variants; clean for single-A spacers |
| CGBE / GBE | C->G or C->A | Yes | 5-10% | Transversions; rare use cases |
| Prime editor (PE2/PE3) | Templated edit; any base change | None | 1-3% | Precise variants; lower efficiency |

**Decision:** For C->T or A->G with available editing window: base editor is preferred (higher efficiency than PE). For other transitions/transversions, multi-base edits, or zero-bystander requirements: prime editor.

## Broad be-validation-pipeline

The Broad Institute's `be-validation-pipeline` (https://broadinstitute.github.io/be-validation-pipeline/) is an end-to-end snakemake workflow for BE variant-function screens:

```bash
# Install
git clone https://github.com/broadinstitute/be-validation-pipeline
cd be-validation-pipeline
conda env create -f environment.yaml
conda activate bevalidation

# Configure
edit config.yaml  # Specify library, FASTQ paths, reference, target genes
# Run
snakemake --use-conda --cores 16

# Outputs:
# results/per_sgrna_editing_efficiency.tsv
# results/per_variant_attribution.tsv
# results/hit_calls.tsv
```

The pipeline handles editing-efficiency filtering, bystander attribution, and variant calling in a single end-to-end workflow optimized for the GRACE library (Sanson 2020) but adaptable.

## Failure Modes

### Mostly indels in BE sample

**Trigger:** Cas9 contamination, wrong vector (e.g., used pCas9-BE3 plasmid but selected on Cas9 line), or evoCDA-BE / broader-window chemistry.
**Mechanism:** Cas9 cuts dsDNA; BE relies on nicked-ssDNA deamination. Cas9 expression in the same cell creates indels.
**Symptom:** Substitution-vs-indel ratio <3 in CRISPResso output.
**Fix:** Verify vector (nCas9-BE3 not Cas9-BE3); confirm cell line lacks Cas9 background; restrict to specifically engineered BE-cell lines.

### High editing but no biological signal

**Trigger:** Bystander C/A is dominating; intended variant is not the perturbation driving phenotype.
**Mechanism:** When target is at position 5 and bystander is at position 7, the molecule carries both; phenotype is from the bystander.
**Symptom:** Strong screen signal but variant attribution unclear.
**Fix:** Run orthogonal prime-editor scan of the same intended variants; restrict library to bystander-free spacers when possible; deconvolute via allele-frequency table.

### sgRNA shows perfect editing but no fitness signal

**Trigger:** Intended variant is silent or compensatory; the protein function is unchanged.
**Mechanism:** Variants can be tolerated; not all variants are LoF or GoF.
**Symptom:** High editing efficiency (>70%) but per-sgRNA LFC near zero.
**Fix:** Expected outcome for many variants; flag silent / compensatory variants in the report.

### Low editing across all guides

**Trigger:** Wrong cell line for the BE; cell line has poor BE activity (some lines lack APOBEC or have low expression).
**Mechanism:** BE efficiency depends on cell-line expression of TadA or APOBEC components.
**Symptom:** Median editing <30% across library.
**Fix:** Test in a BE-validated cell line (HEK293T, U2OS, K562 generally work); pilot before full screen.

### Library missing intended-variant sgRNAs

**Trigger:** No NGG-adjacent spacer places target base in editing window for that codon.
**Mechanism:** Editor window is fixed; some codons cannot be targeted with given chemistry.
**Symptom:** Specific variants absent from screen.
**Fix:** Use PAM-relaxed BE variants (SpRY-CBE, SpRY-ABE); use prime editor for variants outside BE accessibility; accept that some variants cannot be installed.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Editing window | Positions 4-8 from PAM-distal end (BE3/BE4); 4-10 (ABE8e) | Komor 2016; Gaudelli 2017; Lapinaite 2020 |
| Editing efficiency for screen power | >30% (primary); >50% (validation) | Hanna 2021 conventions |
| Indel byproduct (clean BE) | <5%; <2% for ABE | Koblan 2018 (BE4max); Gaudelli 2017 (ABE) |
| Substitution-vs-indel ratio | >10 (clean BE); <3 (Cas9-like) | CRISPResso2 diagnostic |
| Target editing % for variant inclusion | >30% Hanna 2021; >50% strict | Empirical |
| Bystander rate (target attribution) | <10% acceptable; <5% ideal for clean attribution | Application-dependent |
| Cell-line BE activity (pilot) | >30% editing at validated target | Below = wrong cell line for BE |
| Per-amino-acid sgRNA density | 10-15 (Hanna 2021); 5-8 (smaller screens) | Tradeoff with library size |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Substitution-vs-indel ratio <3 | Cas9 contamination or wrong BE | Verify vector / cell line; pilot first |
| All edits at bystander positions | Target base outside window | Re-design spacer with target at pos 4-8 |
| Variant attribution unclear | Bystander confounding | Run orthogonal PE; restrict library |
| Library lacks intended variant | No NGG-PAM accessibility | SpRY-CBE; prime editor; accept exclusion |
| Editing <30% library-wide | Cell-line BE inactivity | Re-validate cell line |
| Hit list dominated by single sgRNA | Bystander-driven phenotype | Cross-check with bystander-free sgRNAs |

## References

- Komor AC et al. 2016. *Nature* 533:420. BE3.
- Gaudelli NM et al. 2017. *Nature* 551:464. ABE7.10.
- Koblan LW et al. 2018. *Nat Biotechnol* 36:843. BE4max + improved CBE.
- Lapinaite A et al. 2020. *Science* 369:566. ABE8e mechanism.
- Richter MF et al. 2020. *Nat Biotechnol* 38:883. ABE8.20.
- Sanson KR et al. 2020. *Nat Commun* 11:5165. GRACE library for BE screens.
- Hanna RE et al. 2021. *Cell* 184:1066. Massively parallel BRCA1/2 variant function via CBE.
- Cuella-Martin R et al. 2021. *Mol Cell* 81:1107. BE saturation of BRCA1/2 variants.
- Arbab M et al. 2020. *Cell* 182:463. BE-Hive prediction of editing outcomes.
- Anzalone AV et al. 2019. *Nature* 576:149. Prime editing (PE2/PE3).
- Clement K et al. 2019. *Nat Biotechnol* 37:224. CRISPResso2.
- Kurt IC et al. 2021. *Nat Biotechnol* 39:41. CGBE1.

## Related Skills

- crispr-screens/crispresso-editing - CRISPResso2 BE/PE mode and allele tables
- crispr-screens/library-design - GRACE-style BE library design
- crispr-screens/prime-editing-screens - Orthogonal PE for variant attribution
- crispr-screens/hit-calling - Variant-level hit aggregation
- crispr-screens/screen-qc - Editing-efficiency QC
- crispr-screens/drugz-chemogenomic - drugZ for BE drug-modifier screens
- clinical-databases/clinvar-lookup - Variant pathogenicity annotation
- variant-calling/variant-annotation - VEP for predicted amino acid changes
