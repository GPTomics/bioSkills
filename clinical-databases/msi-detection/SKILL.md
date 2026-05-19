---
name: bio-clinical-databases-msi-detection
description: Calls microsatellite instability from WES/WGS/targeted-panel with MSIsensor, MSIsensor-pro, MSIsensor-ct (panel-aware), MSIngs, MANTIS, MSIPanel, MSIDetect, and ngsMSI for FDA pembrolizumab MSI-H pan-tumor / Lynch syndrome / dMMR ICI biomarker. Use when stratifying ICI eligibility (Le 2015), pairing MSI with TMB-H (Sha 2020 / Salem 2018), screening Lynch syndrome (universal IHC + MSI), or distinguishing MSI-H tumors from POLE-exo hypermutator with overlapping signatures.
tool_type: cli
primary_tool: MSIsensor-pro
---

## Version Compatibility

Reference examples tested with: MSIsensor-pro 1.2+, MSIsensor 0.6+, MSIngs 1.0+, MANTIS 1.0.5+, samtools 1.19+, mSINGS 5.6+, pandas 2.2+, cyvcf2 0.30+. FDA pembrolizumab MSI-H / dMMR pan-tumor approval is from 2017 (Le 2015 *NEJM*; KEYNOTE-016/164/158); approval extended to colorectal first-line in 2020.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. MSIsensor-pro replaces MSIsensor for tumor-only assays; MSIsensor-ct is the bTMB-equivalent for ctDNA panels.

# MSI Detection; The Companion ICI Biomarker to TMB

**'Detect MSI status from this somatic sequencing data'** -> Profile microsatellite instability across canonical loci (Bethesda 5 panel + extended NGS-derived sites); classify MSI-H / MSS / MSI-L per Bethesda / FDA / KEYNOTE convention.

- CLI (recommended tumor-only): `msisensor-pro msi -d microsatellites.list -t tumor.bam -o msi_out -b 16`
- CLI (paired tumor-normal): `msisensor msi -d microsatellites.list -n normal.bam -t tumor.bam -o msi_out`
- CLI (ctDNA / blood MSI): `msisensor-ct ...`
- CLI (older WES standard): `mantis -t tumor.bam -n normal.bam -b targets.bed --threads 8`
- CLI (panel-only): MSIPanel for clinical assays

## The Regulatory and Trial Landscape

| Event | Year | Threshold | Notes |
|-------|------|-----------|-------|
| **Le 2015** *NEJM* | 2015 | MSI-H + ICI in CRC | The seminal paper: pembrolizumab in MSI-H CRC ORR 40% vs 0% MSS |
| **FDA pembrolizumab MSI-H / dMMR pan-tumor** | 2017 | MSI-H | First tissue-agnostic FDA approval (KEYNOTE-016/164/158) |
| **FDA pembrolizumab first-line MSI-H CRC** | 2020 | MSI-H + first-line CRC | KEYNOTE-177 |
| **CheckMate 142** | 2017-2018 | MSI-H + nivolumab/ipilimumab | Pan-tumor MSI-H second-line |
| **ESMO 2024** | 2024 | MSI-H | Maintained pan-tumor MSI-H biomarker |
| **Universal Lynch screening** | -- | IHC + MSI on all CRC <= 70 yr | NCCN / ACG / EGAPP guidelines |

## MSI vs dMMR vs TMB-H: The Conceptual Hierarchy

| Term | Definition | Method | Relationship |
|------|-----------|--------|--------------|
| **dMMR (deficient MMR)** | Loss of MMR protein function | IHC (MLH1, MSH2, MSH6, PMS2) | Causes MSI |
| **MSI-H** | Microsatellite instability high | PCR-based Bethesda or NGS | Consequence of dMMR |
| **Lynch syndrome** | Germline MMR mutation | Germline sequencing | Causes ~50% of MSI-H CRC; rest are sporadic (MLH1 hyper-methylation) |
| **TMB-H** | >= 10 mut/Mb | NGS panel / WES | Statistical correlate of MSI-H |
| **POLE-exo hypermutator** | POLE proofreading defect | Sequencing / signatures | Hypermutator WITHOUT MMR-D; MSI-stable typically |

**MSI-H + TMB-H overlap** (Salem 2018 *Cancer Discov*):
- ~80% of MSI-H tumors in CRC and endometrial are TMB-H.
- ~16% of TMB-H solid tumors are MSI-H.
- Sha 2020 *Cell Rep Med*: MSI-H supersedes TMB-H for ICI biomarker decision; TMB-H not additive.

**POLE-exo vs MMR-D:**
- POLE-exo (SBS10a/10b): hypermutator (100-300 mut/Mb pure); typically MSI-stable.
- MMR-D (SBS6/15/26/44 + ID1/2): 30-50 mut/Mb typical; MSI-H.
- POLE-exo + MMR-D (SBS14 + SBS20): ultra-hypermutator >=500 mut/Mb; MSI-H.

## Tool Taxonomy

| Tool | Paired | Tumor-only | ctDNA | Algorithm | Fails when |
|------|--------|-----------|-------|-----------|-----------|
| **MSIsensor** (Niu 2014 *Bioinformatics*) | Yes | No | No | Bayesian + read-length distribution | Tumor-only data (no baseline); cohort baseline missing |
| **MSIsensor-pro** (Jia 2020 *Genom Proteom Bioinform*) | Optional | **Yes** | No | Distribution comparison to baseline | Baseline cohort not provided; panel < 50 loci |
| **MSIsensor-ct** (Han 2021 *Brief Bioinform*) | -- | -- | **Yes** | cfDNA-aware | Tumor fraction < 3%; low ctDNA shed |
| **MANTIS** (Kautto 2017 *Oncotarget*) | Yes | No | No | Step-wise difference | Tumor-only; low coverage at microsatellites |
| **MSIngs** (Salipante 2014 *Clin Chem*) | Yes | No | No | Number of unstable loci | Tumor-only; outdated vs MSIsensor-pro |
| **mSINGS** (Salipante 2014) | -- | Yes | No | Background panel | Background panel poorly characterized for cohort |
| **MSIPanel** | -- | Yes | -- | Panel-specific calibration | Panel not in vendor calibration table |
| **ngsMSI** | Yes | Yes | -- | Algorithmic variations | Limited benchmarking; rarely first choice |
| **MSIDetect** | Various | -- | -- | DL-based | New tool; reproducibility data still maturing |
| **MMR-MS** | IHC concordance check | -- | -- | -- | Not a direct MSI caller; for orthogonal validation |

**Operational consensus 2024-2026:**
- **Tumor + paired normal WES:** MSIsensor or MANTIS.
- **Tumor-only assay** (commercial panels, often unpaired): MSIsensor-pro with reference baseline.
- **ctDNA / liquid biopsy:** MSIsensor-ct.
- **Lynch screening:** IHC FIRST (rules out 90%+); MSI-PCR / NGS confirmatory.

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Tumor + paired normal WES | MSIsensor (standard) | Reference paired-normal comparison |
| Tumor-only WES/panel | MSIsensor-pro with panel baseline | No matched normal needed |
| ctDNA / liquid biopsy | MSIsensor-ct | cfDNA-aware |
| Lynch syndrome screening | Universal IHC + MSI (NCCN) | IHC catches 90%+; MSI for IHC-equivocal |
| FDA pembrolizumab eligibility | Validate per FoCR PCR + IHC + NGS concordance | Cross-platform required |
| MSI-H + TMB-H concurrence | MSI-H is primary biomarker | Sha 2020; TMB-H not additive |
| POLE+MMR ultra-hypermutator | Sigprofiler signatures (SBS14, SBS20) | Mechanism beyond MSI alone |
| Sporadic MSI-H | Confirm MLH1 hypermethylation; rule out Lynch | Distinguishes sporadic vs germline |
| MSI-stable + TMB-H | Investigate POLE-exo signature (SBS10a/10b) | POLE-exo causes hypermutator without MSI |
| Pan-tumor screening | MSI + IHC + TMB combined | Multiple modalities for ICI eligibility |

## Bethesda Panel and Modern NGS-Derived Loci

The **Bethesda 5-locus panel** (Boland 1998) is the PCR-based reference:
- **BAT-25** (chr4)
- **BAT-26** (chr2)
- **NR-21** (chr14)
- **NR-24** (chr2)
- **MONO-27** (chr2; also called BAT-25/Penta C in some nomenclature)

NGS-based MSI panels use 50-1000+ microsatellite loci. MSI-H requires unstable status at >=40% of tested loci typically (varies by panel calibration).

## Standard Workflow: MSIsensor-pro Tumor-Only

**Goal:** Compute MSI status from tumor-only WES/panel.

**Approach:** Generate baseline from population reference; compare patient tumor.

```bash
# Generate microsatellite list from reference genome (one-time)
msisensor-pro scan -d /reference/GRCh38.fa -o microsatellites.list -p 1 -m 5

# Generate baseline from N normal control samples (one-time per panel)
msisensor-pro baseline -d microsatellites.list -i normal_samples.list -o baseline.list -b 16

# Score tumor sample
msisensor-pro pro \
    -d microsatellites.list \
    -t tumor.bam \
    -i sample_id \
    -o msi_output \
    -b 16 \
    --baseline baseline.list

# Output: msi_output_all (raw); msi_output_unstable (unstable loci); msi_output.txt (summary)
# Critical column: %_unstable. Threshold MSI-H typically >= 20-30% depending on panel.
```

## Paired Tumor-Normal MSIsensor

```bash
msisensor msi \
    -d microsatellites.list \
    -n normal.bam \
    -t tumor.bam \
    -o msi_paired_out \
    -b 16

# Output: %_unstable in paired comparison
# MSI-H threshold: >= 20% by FoCR guidance; varies 10-30% across studies
```

## MANTIS Step-wise Difference

```bash
mantis.py \
    -t tumor.bam \
    -n normal.bam \
    -b microsatellite_targets.bed \
    --threads 8 \
    -o mantis_output

# Output: mantis_output.kmer_counts (raw), mantis_output (status)
# Threshold MSI-H: stepwise difference > 0.4 (default)
```

## MSI-H Classification Logic

```python
import pandas as pd


def classify_msi(unstable_percentage, panel_calibrated_cutoff=20.0):
    '''Classify MSI status from percentage of unstable loci.

    Bethesda PCR: >=2 of 5 unstable -> MSI-H (40% loci)
    NGS: panel-specific cutoffs typically 10-30%
    Concordance: MSI-PCR + IHC + NGS should agree (FoCR)
    '''
    if unstable_percentage >= panel_calibrated_cutoff:
        return 'MSI-H'
    elif unstable_percentage >= panel_calibrated_cutoff / 2:
        return 'MSI-L (intermediate; treat as MSS clinically per FDA)'
    else:
        return 'MSS'


def msi_lynch_workflow(msi_status, ihc_results, mlh1_methylation_status, germline_test):
    '''Standard Lynch syndrome workflow.

    Args:
        msi_status: 'MSI-H' / 'MSS' / 'MSI-L'
        ihc_results: dict {MLH1: 'retained' or 'loss', MSH2, MSH6, PMS2}
        mlh1_methylation_status: 'methylated' (sporadic) / 'unmethylated' (Lynch suspect)
        germline_test: 'positive' / 'negative' / 'not_performed'
    '''
    if msi_status != 'MSI-H':
        return 'No further Lynch screening indicated'

    ihc_loss = [gene for gene, status in ihc_results.items() if status == 'loss']
    if not ihc_loss:
        return 'MSI-H with retained IHC; consider Lynch with germline testing'

    if 'MLH1' in ihc_loss:
        if mlh1_methylation_status == 'methylated':
            return 'Sporadic MSI-H (MLH1 hypermethylation); not Lynch'
        elif mlh1_methylation_status == 'unmethylated':
            return 'Lynch suspect (MLH1 loss without methylation); proceed with germline testing'
        else:
            return 'MLH1 loss; perform methylation test'

    return f'MSH2/6/PMS2 loss ({", ".join(ihc_loss)}); strong Lynch suspect; germline testing'


def msi_tmb_ici_decision(msi_status, tmb_value, tumor_type=None, dmmr_ihc=None):
    '''Integrated ICI eligibility from MSI + TMB.

    Sha 2020: MSI-H is primary biomarker; TMB-H not additive.
    McGrail 2021: TMB-H NOT endorsed for breast/prostate/glioma alone.
    '''
    msi_high = msi_status == 'MSI-H'
    dmmr_positive = dmmr_ihc == 'positive'
    tmb_h = tmb_value >= 10

    if msi_high or dmmr_positive:
        return ('ICI eligible: MSI-H or dMMR (FDA pembrolizumab 2017 pan-tumor; KEYNOTE-016/164/158); '
                'TMB-H is not additive (Sha 2020).')
    if tmb_h and tumor_type and tumor_type.lower() in ('breast', 'prostate', 'glioma'):
        return ('TMB-H but tumor type excluded by ESMO 2024 / McGrail 2021. '
                'Consider tumor-type-specific cutoff.')
    if tmb_h:
        return 'TMB-H pan-tumor (FDA pembrolizumab 2020); ICI eligible.'
    return 'MSS + TMB-low. Standard chemo per tumor type.'
```

## Per-Operation Failure Modes

**1. Tumor-only with paired-normal tool**
- Trigger: Run MSIsensor on tumor-only BAM.
- Mechanism: MSIsensor requires paired normal for baseline comparison.
- Symptom: Tool errors or produces unstable noisy result.
- Fix: Use MSIsensor-pro for tumor-only; or use mSINGS background-panel approach.

**2. Panel size too small**
- Trigger: 30-locus panel called MSI-H based on 20% threshold (= 6 unstable loci).
- Mechanism: Small panel + stochastic unstable rates produce high false-positive rates.
- Symptom: False-positive MSI-H in WES-comparable panels with < 50 microsatellite loci.
- Fix: Validate panel calibration with reference cohort; use panel-specific cutoff; minimum 50 informative loci.

**3. IHC vs MSI discordance not investigated**
- Trigger: IHC retains all four MMR proteins; MSI-H by sequencing.
- Mechanism: IHC may miss subtle loss; MSI may include MSH6-only subtype (more variable); rare germline POLE+MMR ultra-hypermutators show MSI.
- Symptom: Apparent discordance; classification ambiguous.
- Fix: Cross-check with germline MMR sequencing; check for POLE-exo on Sigprofiler.

**4. MSI-H + Lynch syndrome confusion**
- Trigger: Report MSI-H tumor as "Lynch syndrome".
- Mechanism: ~50% of MSI-H CRC is sporadic (MLH1 hypermethylation, not germline Lynch).
- Symptom: Incorrect family counseling; wrong screening.
- Fix: Apply IHC + MLH1 methylation + germline testing workflow.

**5. POLE-exo hypermutator labeled MSI**
- Trigger: Tumor with 200 mut/Mb POLE-exo signature labeled MSI-H.
- Mechanism: Pure POLE-exo causes hypermutator WITHOUT MSI (different repair mechanism); apparent MSI-H call may be a false positive in high-mutation context.
- Symptom: Misclassification; ICI eligibility still positive but for different mechanism.
- Fix: Run Sigprofiler signatures (SBS10a/10b vs SBS6/15/26/44); confirm POLE-exo via SBS10 contribution.

**6. ctDNA MSI without sufficient tumor fraction**
- Trigger: Run MSIsensor-ct on cfDNA with <1% tumor fraction.
- Mechanism: Low ctDNA fraction produces noise-dominated unstable locus counts.
- Symptom: False-negative or unstable MSI call.
- Fix: Estimate tumor fraction first (ichorCNA); require >= 3% for reliable cfDNA MSI.

**7. Universal screening missed**
- Trigger: CRC patient < 70 yr without IHC / MSI.
- Mechanism: NCCN / ACG universal Lynch screening required; without it, Lynch syndrome undiagnosed.
- Symptom: Family loses screening benefit.
- Fix: Universal IHC + MSI on all CRC < 70; institute reflex testing.

**8. MSI-L treated as actionable**
- Trigger: Report MSI-L (intermediate) as ICI-eligible.
- Mechanism: FDA approval specifies MSI-H; MSI-L = MSS clinically.
- Symptom: ICI given on insufficient indication; reimbursement issues.
- Fix: Apply MSI-H threshold strictly per FDA; MSI-L = MSS.

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| PCR Bethesda MSI-H vs NGS MSS | Bethesda panel uses 5 loci only; less sensitive | Trust NGS with >=50 informative loci |
| NGS MSI-H vs IHC retained | Subtle MMR loss; MSH6-only subtype; or POLE-exo | Confirm with germline + POLE-exo signature analysis |
| Paired-normal MSI-H + tumor-only MSS | Sample swap or low tumor purity in tumor-only | Re-validate; check purity (>=20% required) |
| MSIsensor-pro vs MSIsensor (paired) | Different baseline thresholds | Apply panel-specific calibration |
| MSI-H suspected but tools differ | Borderline mutational burden | Use signature analysis (SBS6/15/26/44) as orthogonal evidence |
| ctDNA MSI vs tissue MSI | Tumor fraction low | Trust tissue; estimate ctDNA fraction |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Bethesda MSI-H | >= 2/5 unstable | Boland 1998 |
| NGS MSI-H cutoff | 10-30% unstable loci (panel-specific) | Various |
| MANTIS MSI-H threshold | Step-wise difference > 0.4 | Kautto 2017 |
| MSIsensor MSI-H threshold | >= 20% by FoCR | Friends of Cancer Research |
| Minimum informative loci | >= 50 NGS loci | Panel-design convention |
| ctDNA tumor fraction minimum | >= 3% for reliable cfDNA MSI | Han 2021 |
| Tumor purity minimum | >= 20% | Standard |
| FDA approval | MSI-H or dMMR pan-tumor (2017) | KEYNOTE-016/164/158 |
| First-line MSI-H CRC | KEYNOTE-177 (2020) | -- |
| MSI-H + TMB-H overlap (CRC, endometrial) | ~80% | Salem 2018 |
| TMB-H -> MSI-H rate | ~16% | Salem 2018 |
| Sporadic MSI-H mechanism | ~50% MLH1 hypermethylation | Various |
| Universal screening cutoff | CRC <= 70 yr | NCCN / ACG |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| MSI-H + IHC retained discordance | Subtle loss; MSH6-only; or rare hypermutator | Cross-check germline + signatures |
| Borderline MSI call | Panel too small | Use >= 50 informative loci |
| Tumor-only MSI low confidence | Background subtraction needed | Use MSIsensor-pro with cohort baseline |
| MSI-H + TMB-H reported additive | Tautology per Sha 2020 | MSI-H is primary; TMB-H not additive |
| POLE-exo labeled MMR-D | Different mechanism; mutation count differs | Run Sigprofiler; SBS10a/10b is POLE-exo |
| Sporadic MSI-H mis-labeled Lynch | Need MLH1 methylation test | Confirm MLH1 methylation + germline |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "MSI-H + TMB-H both reported additive" | Sha 2020 *Cell Rep Med*: MSI-H is the primary biomarker; TMB-H is statistical correlate. We report MSI-H first; TMB-H reported but noted not additive. |
| "Why MSIsensor-pro instead of MSIsensor?" | MSIsensor requires paired normal; MSIsensor-pro handles tumor-only via cohort baseline. Most commercial panels are tumor-only. |
| "MSI-PCR vs NGS discordant" | Bethesda 5-locus panel is less sensitive; we use NGS >=50 informative loci for confirmation. |
| "Universal Lynch screening?" | NCCN / ACG recommend reflex IHC + MSI on all CRC <= 70 yr; we implemented universal screening protocol. |
| "POLE-exo hypermutator with MSI-H?" | Sigprofiler signature analysis distinguishes: SBS10a/10b = POLE-exo (typically MSI-stable); SBS6/15/26/44 = MMR-D. POLE+MMR concurrent produces ultra-hypermutator. |
| "MSI-L?" | FDA approval specifies MSI-H; MSI-L = clinically MSS; we apply MSI-H threshold strictly. |
| "ctDNA MSI viability?" | MSIsensor-ct works if tumor fraction >= 3%; we estimate via ichorCNA; below threshold falls back to tissue. |

## References

- Le DT et al. 2015. PD-1 blockade in tumors with mismatch-repair deficiency. *NEJM* 372:2509. (The seminal paper)
- Marabelle A et al. 2020. Efficacy of pembrolizumab in patients with advanced MSI-H or dMMR cancer. *J Clin Oncol* 38:1.
- Niu B et al. 2014. MSIsensor: microsatellite instability detection using paired tumor-normal sequence data. *Bioinformatics* 30:1015.
- Jia P et al. 2020. MSIsensor-pro: fast, accurate, and matched-normal-sample-free detection of microsatellite instability. *Genomics Proteomics Bioinformatics* 18:65.
- Han P et al. 2021. MSIsensor-ct: microsatellite instability detection using cfDNA sequencing data. *Brief Bioinform* 22:bbaa402.
- Kautto EA et al. 2017. Performance evaluation for rapid detection of pan-cancer microsatellite instability with MANTIS. *Oncotarget* 8:7452.
- Salipante SJ et al. 2014. Microsatellite instability detection by NGS. *Clin Chem* 60:1192.
- Boland CR et al. 1998. National Cancer Institute workshop on microsatellite instability for cancer detection and familial predisposition. *Cancer Res* 58:5248.
- Salem ME et al. 2018. Landscape of TMB in different cancers: implications for ICI response. *Cancer Discov* 8:1136.
- Sha D et al. 2020. TMB as a predictive biomarker in solid tumors. *Cell Rep Med* 1:100043.
- Vanderwalde A et al. 2018. MSI and dMMR comparison: the role of NGS. *J Mol Diagn* 20:809.

## Related Skills

- clinical-databases/tumor-mutational-burden - TMB as related ICI biomarker
- clinical-databases/somatic-signatures - SBS6/15/26/44 MMR-D signatures + SBS10a/10b POLE-exo
- clinical-databases/clinvar-lookup - Lynch syndrome variant pathogenicity (MLH1, MSH2, MSH6, PMS2)
- clinical-databases/variant-prioritization - Germline MMR variant prioritization for Lynch
- variant-calling/clinical-interpretation - Clinical reporting
