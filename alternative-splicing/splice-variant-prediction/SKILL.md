---
name: bio-splice-variant-prediction
description: Predicts whether a DNA variant alters mRNA splicing using sequence-based deep-learning tools — SpliceAI (10kb context dilated CNN, clinical default), Pangolin (multi-tissue), MMSplice (modular per-region CNN with calibrated ΔPSI), SpliceTransformer/TrASPr (tissue-aware transformers), SpliceVault (empirical 300K-RNA lookup of likely mis-splicing outcomes), CADD-Splice (composite score). Applies the ClinGen SVI 2023 framework for ACMG/AMP variant interpretation (PVS1, PP3, BP4 evidence codes), HGVS splicing nomenclature (c.123+1G>A, c.123-3T>G, r.spl?), extended-window scoring for deep-intronic pseudoexons, tissue-specific predictions, branchpoint variant detection (BPHunter, LaBranchoR), and splice-switching ASO design. Use when interpreting splice impact of clinical variants, prioritizing VUS, identifying deep-intronic pathogenic variants, or designing ASOs.
tool_type: python
primary_tool: SpliceAI
---

## Version Compatibility

Reference examples tested with: SpliceAI 1.3+, Pangolin 1.0+, MMSplice 2.4+, pyensembl 2.3+, pysam 0.22+, pandas 2.2+, gffutils 0.13+, tensorflow 2.15+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Splice Variant Prediction

Predict whether a DNA variant alters mRNA splicing. **Distinct from "variant pathogenicity" generally**: a variant can be a strong splice disruptor without being pathogenic for the gene's standard mechanism, or pathogenic for reasons orthogonal to splicing. Splice prediction asks specifically: does this variant change splice-site usage?

## Predictor Taxonomy

| Family | Architecture | Output | Fails when |
|--------|--------------|--------|------------|
| Context-aware CNN | 10 kb dilated ResNet | Per-position donor/acceptor probability | Long-range (>5 kb) regulatory effects; tissue-specific events |
| Tissue-aware CNN/transformer | Same arch + multi-tissue training | Per-tissue ΔPSI | Tissue not in training set; novel cell types |
| Modular per-region CNN | Separate sub-models for 5'ss/3'ss/exon/intron | Calibrated quantitative ΔPSI | Atypical events; complex multi-junction effects |
| Foundation transformer | Pretrained on broad genomic context | Splice probability or ΔPSI | New tools; less battle-tested |
| Empirical lookup | Public RNA-seq event database | Top-N most likely mis-splicing outcomes | Variant types not represented in training cohorts |
| Composite score | Blend of multiple predictors | Single scaled score | When component predictors disagree internally |

## Tool Selection Matrix

| Tool | Best for | Output | When to use | Fails when |
|------|----------|--------|-------------|------------|
| SpliceAI | Clinical screening; canonical splice site disruption | Delta score 0-1 | Default for ACMG variant classification | Tissue-specific events; deep-intronic with default 50nt window |
| Pangolin | Tissue-aware predictions | Per-tissue ΔPSI | When disease tissue is known (brain, heart, liver, testis) | Tissue not in 4-tissue training set |
| MMSplice | Quantitative ΔPSI | Δlogit_psi | Research where calibrated effect-size matters | Atypical events outside cassette-exon model |
| SpliceTransformer | 2024+ benchmark improvements | Tissue-specific ΔPSI | When transformer foundation models outperform CNN on benchmark variant sets | New (2024); limited clinical adoption |
| TrASPr | Multi-transformer, 2024-2025 | Tissue-specific PSI/ΔPSI | Strong on tissue-specific test sets | New; verify before clinical use |
| SpliceVault | Empirical mis-splicing outcome | Top-N events at the affected splice site | Predicting consequence (skip vs cryptic) of canonical-disrupting variants | Variants not represented in 300K-RNA training |
| CADD-Splice | Single composite score | Scaled C-score | Clinical pipelines wanting one number | When knowing which sub-component drove the score is needed |

Methodology evolves; verify benchmarks (Smith & Kitzman 2023 *Genome Biol* 24:294; You et al 2024 *Nat Commun*) and ClinGen SVI splicing recommendations before reporting clinical interpretations. Concordance across SpliceAI + Pangolin + MMSplice is gold-standard evidence; discordance flags need RNA validation.

## Decision Tree by Use Case

| Use case | Recommended approach |
|----------|----------------------|
| Clinical variant report (single variant, ACMG classification) | SpliceAI default 50nt + ClinGen SVI 2023 thresholds |
| Tissue-specific clinical question (brain disease, cardiomyopathy) | SpliceAI + Pangolin (tissue-matched) |
| Unsolved Mendelian case (suspect deep-intronic) | SpliceAI extended window (-D 500-2000) + SpliceVault |
| VUS panel screening | SpliceAI + Pangolin + MMSplice concordance scoring |
| Predict consequence of canonical-disrupting variant | SpliceVault top-N empirical events |
| Branchpoint variant suspected | BPHunter (branchpoint screen) — SpliceAI is weak here |
| Splice-switching ASO design (target ESE/ESS occlusion) | SpliceAI on masked sequence + RNAfold accessibility |
| Validate predicted splice change in patient | RNA-seq + FRASER2 (see outlier-splicing-detection) |
| Pseudoexon prediction in deep intron | SpliceAI extended window + CI-SpliceAI; require RNA validation |

## ClinGen SVI 2023 Framework

The ClinGen Sequence Variant Interpretation (SVI) splicing subgroup (Walker 2023 *Am J Hum Genet*) extended the ACMG/AMP 2015 framework with explicit splice-prediction rules.

| Evidence code | Threshold | Notes |
|----------------|-----------|-------|
| **PP3** (supporting pathogenic) | SpliceAI delta >= 0.20 | ClinGen SVI: apply at **supporting** weight (not standalone) |
| **BP4** (supporting benign) | SpliceAI delta <= 0.10 | ClinGen SVI: apply at **supporting** weight |
| **PVS1** (very strong null) | Canonical +/-1, +/-2 site disruption with predicted LoF + NMD | Requires gene where LoF is established mechanism (Abou Tayoun 2018 *Hum Mutat* PVS1 decision tree) |
| **PS3 / BS3** (functional) | RNA evidence (RT-PCR, RNA-seq, minigene) | Supersedes computational evidence |

**Operational rules:** Computational evidence (PP3/BP4) is *supporting*, not standalone. ClinGen SVI 2023 recommends applying predictive splice PP3/BP4 at **supporting** weight only; higher SpliceAI cutoffs (0.5, 0.8) increase precision but are the tool's own tiers (Jaganathan 2019), NOT ClinGen-endorsed evidence-strength upgrades — reaching moderate/strong requires functional/RNA evidence (PS3/BS3), not a higher SpliceAI score alone. Splicing variants benefit from concordance across SpliceAI + Pangolin + MMSplice. RNA validation supersedes prediction. Always log SpliceAI version, distance window, and reference transcript. SpliceAI alone is **not sufficient** for PVS1; canonical site disruption requires gene-level LoF context.

## SpliceAI Workflow

**Goal:** Annotate VCF variants with per-variant delta scores for splice-site change.

**Approach:** Run `spliceai` CLI with reference genome and annotation; parse INFO field for delta scores. **SpliceAI is human-only** (`-A grch37` or `-A grch38`); the model was trained on GENCODE human and does not directly transfer to mouse, fly, or other species. For mouse, retrained variants exist (e.g. mouseSpliceAI); for other species, use Pangolin (4 species: human, mouse, rat, rhesus macaque) or accept that prediction will be unreliable.

```bash
spliceai \
    -I input.vcf \
    -O output.vcf \
    -R GRCh38.primary_assembly.genome.fa \
    -A grch38 \
    -D 50 \
    -M 0
```

`-D 50` = distance window in nt around variant (default 50). For deep-intronic variants suspected of creating pseudoexons, raise to **500-2000**:

```bash
spliceai -I input.vcf -O output_extended.vcf -R genome.fa -A grch38 -D 500 -M 1
```

`-M 0` (default) returns raw scores; `-M 1` masks splice gains at annotated sites and losses at unannotated sites (cleaner for clinical use). Output INFO format: `SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`. Delta score = max(DS_AG, DS_AL, DS_DG, DS_DL).

```python
import pandas as pd
import re

def parse_spliceai_vcf(vcf_path):
    rows = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            info = fields[7]
            m = re.search(r'SpliceAI=([^;]+)', info)
            if not m:
                continue
            for ann in m.group(1).split(','):
                parts = ann.split('|')
                allele, symbol = parts[0], parts[1]
                ds = [float(p) if p != '.' else 0 for p in parts[2:6]]
                dp = parts[6:10]
                rows.append({
                    'chrom': fields[0], 'pos': int(fields[1]),
                    'ref': fields[3], 'alt': allele,
                    'gene': symbol,
                    'DS_AG': ds[0], 'DS_AL': ds[1],
                    'DS_DG': ds[2], 'DS_DL': ds[3],
                    'delta_max': max(ds),
                })
    return pd.DataFrame(rows)

df = parse_spliceai_vcf('output.vcf')
df['acmg_evidence'] = pd.cut(
    df['delta_max'],
    bins=[-0.01, 0.10, 0.20, 0.50, 0.80, 1.01],
    # ClinGen SVI applies splice PP3/BP4 at supporting weight; 0.5/0.8 are SpliceAI
    # precision tiers (Jaganathan 2019), NOT ACMG evidence-strength upgrades
    labels=['BP4', 'inconclusive', 'PP3_supporting', 'PP3_supporting_prec0.5', 'PP3_supporting_prec0.8']
)
```

DS labels: AG = acceptor gain, AL = acceptor loss, DG = donor gain, DL = donor loss.

## Pangolin for Tissue-Specific Prediction

**Goal:** Get tissue-specific splice impact predictions when disease tissue is known.

**Approach:** Run Pangolin CLI with VCF + reference + gffutils annotation database.

```bash
python -c "import gffutils; gffutils.create_db('gencode.v45.annotation.gff3', 'gencode.db', force=True)"

pangolin \
    input.vcf \
    GRCh38.primary_assembly.genome.fa \
    gencode.db \
    pangolin_output \
    -d 500 \
    -m True \
    -s 0.2
```

`-m True` masks splice gains at annotated sites and losses at unannotated sites (recommended for clinical use). `-s 0.2` outputs all sites with predicted change >= cutoff.

Pangolin output is a VCF with per-tissue predictions across the **4 tissues used at training: brain, heart, liver, testis** (Zeng & Li 2022 *Genome Biol*). The model outputs per-species per-tissue predictions but extrapolates poorly to tissues outside this set. Use the tissue closest to disease-relevant context. **For tissues not in the 4-tissue training set, fall back to SpliceAI** — Pangolin extrapolates poorly to unseen tissues.

## SpliceVault for Empirical Mis-Splicing Outcomes

**Goal:** Predict the *type* of mis-splicing (exon skipping vs cryptic site activation) given a canonical-disrupting variant.

**Approach:** Query SpliceVault's database of empirical mis-splicing events from public RNA-seq.

```python
import requests

# Web API: https://kidsneuro.shinyapps.io/splicevault/
# Or use the R/Python package at github.com/kidsneuro-lab/SpliceVault

# Example: NM_000546.6:c.673-2A>G (TP53)
# Returns top-N most likely mis-splicing events: exon skipping, cryptic 3'ss usage, etc.
```

SpliceVault (Dawes 2023 *Nat Genet*) showed that the **Top-4 events** at any splice site predict variant-associated mis-splicing with ~92% sensitivity overall (96% of exon-skipping and 86% of cryptic-activation events) — a striking regularity that makes consequence prediction tractable. Use SpliceVault when the question is not "will splicing change?" but "what specific aberrant splicing will occur?".

## MMSplice for Calibrated ΔPSI

**Goal:** Predict quantitative ΔPSI (not just probability of disruption) for cassette exons.

**Approach:** Score variant impact on each splicing region (5'ss, 3'ss, exon, intron-3'/5') and combine.

```python
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save

dl = SplicingVCFDataloader(
    gtf='gencode.v45.basic.gtf',
    fasta_file='GRCh38.fa',
    vcf_file='input.vcf'
)

model = MMSplice()
predict_save(model, dl, 'mmsplice_predictions.csv', pathogenicity=True)
```

MMSplice (Cheng 2019 *Genome Biol*) reports Δlogit_psi per variant. Useful when calibrated effect sizes matter (research) more than probability of disruption (clinical screening). Companion **MTSplice** (Cheng 2021 *Genome Biol*) adds tissue-specific Δψ predictions.

## HGVS Splicing Nomenclature

Following den Dunnen 2016 *Hum Mutat*:

| Notation | Meaning |
|----------|---------|
| `c.123+1G>A` | +1 of intron downstream of exon ending at cDNA position 123 (canonical 5'ss G) |
| `c.123+5G>A` | +5 position of donor (consensus region) |
| `c.124-1G>A` | -1 of acceptor (canonical AG) |
| `c.124-3T>G` | -3 of acceptor (Py-tract / BPS region) |
| `c.124-50A>G` | Deep-intronic; may activate cryptic site |
| `r.123_456del` | RNA-level deletion (predicted exon skipping) |
| `r.spl?` | Unknown splice consequence |
| `r.0?` | No detectable RNA |
| `p.0?` | Unknown protein consequence |
| `p.(=)` | No predicted protein change (silent) |

Validation tools: VariantValidator (Freeman 2018 *Hum Mutat*), Mutalyzer 2 (Lefter et al 2021 *Bioinformatics* 37:2811-2817).

## Extended-Window Scoring for Deep-Intronic Variants

SpliceAI's default precomputed scores use a **50-nt window**, missing variants that create pseudoexons in deep intronic regions. For unsolved Mendelian cases:

```bash
# Recompute with extended window
spliceai -I input.vcf -O output_2kb.vcf -R genome.fa -A grch38 -D 2000

# Or use CI-SpliceAI (Strauch 2022 PLoS One), SpliceAI retrained on curated GENCODE splice sites
```

| Window | Tradeoff |
|--------|----------|
| -D 50 (default) | Fast; captures canonical-site disruption; misses deep-intronic |
| -D 500 | Captures most pseudoexon-creating deep-intronic variants |
| -D 2000 | Maximum sensitivity; some false positives at large distances |

Pseudoexon creation in deep introns explains a substantial fraction of unsolved Mendelian disease alleles in current cohorts (estimates 5-15% across studies; specific quantitative range will vary by cohort and panel — verify against current literature). Disease examples: CFTR 3849+10kbC>T, USH2A c.7595-2144A>G, CEP290 c.2991+1655A>G (LCA10).

## Concordance Across Predictors

```python
import pandas as pd

merged = (spliceai_df
    .merge(pangolin_df, on=['chrom', 'pos', 'alt'], suffixes=('_sai', '_pang'))
    .merge(mmsplice_df, on=['chrom', 'pos', 'alt'])
)

merged['concordance'] = (
    (merged['delta_max_sai'] >= 0.2).astype(int) +
    (merged['pangolin_score'].abs() >= 0.2).astype(int) +
    (merged['delta_logit_psi'].abs() >= 1.0).astype(int)
)

merged['interpretation'] = merged['concordance'].map({
    0: 'concordant_benign',
    1: 'discordant_low_evidence',
    2: 'concordant_evidence',
    3: 'high_concordance_pathogenic'
})
```

| Concordance | Interpretation | Action |
|-------------|----------------|--------|
| 3/3 above threshold | High confidence | PP3 (supporting); strong candidate for RNA validation (PS3) |
| 2/3 above | Concordant evidence | PP3 (supporting) |
| 1/3 above | Discordant | Report inconclusive; flag for RNA validation |
| 0/3 above | Concordant benign | BP4 (supporting) |

Discordance is the most informative pattern — variants where one model sees impact and others don't are high priority for RNA validation.

## Branchpoint Variant Detection

All current tools are **weak at branchpoint variants** because the BPS motif (yUnAy) has low information content. Specific branchpoint tools:

| Tool | Method | Notes |
|------|--------|-------|
| BPP | Mixture model (BP motif + polypyrimidine tract) | Zhang 2017 *Bioinformatics* 33:3166 |
| LaBranchoR | Bidirectional LSTM | Paggi & Bejerano 2018 *RNA* 24:1647 |
| SVM-BPfinder | SVM on conservation+sequence | Corvelo 2010 *PLoS Comput Biol* |
| BPHunter | Genome-wide branchpoint screen against an aggregated experimental (lariat/RNA-seq) + computational BP database | Zhang 2022 *PNAS* |

Branchpoint variants are under-recognized in clinical pipelines; SpliceAI captures only some because branchpoint motifs have low information content. **Recommendation:** when SpliceAI delta is borderline (0.1-0.3) for a variant in the BPS region (-18 to -40 from 3'ss), run BPHunter as supplement.

## Splice-Switching ASO Design

**Goal:** Design antisense oligonucleotides to modulate splicing therapeutically (e.g. SMA ISS-N1, DMD exon skipping).

**Approach:** Use SpliceAI to predict impact of binding-site occlusion; check accessibility (RNAfold); avoid SR/hnRNP off-target motifs.

```python
# Conceptual workflow - actual design uses ASO synthesis platforms
# 1. Identify target ESE/ESS/ISE/ISS region from MaxEntScan + SpliceAI scan
# 2. Design candidate 18-22 nt ASOs spanning the regulatory element
# 3. For each ASO, simulate splice-site occlusion impact via SpliceAI on the masked sequence
# 4. Filter for RNA accessibility (avoid stable hairpins) using RNAfold
# 5. Whole-transcriptome SpliceAI scan for off-target binding (>=17/20 nt match)
# 6. Avoid TLR9 immunostimulatory CpG motifs

# Chemistry choices:
# - 2'-MOE-PS: nusinersen-like (CNS, intrathecal)
# - PMO: DMD ASOs (systemic IV)
# - GalNAc-conjugated: hepatic targeting
```

Approved precedents: **nusinersen** (SMA ISS-N1 occlusion, exon 7 inclusion); **risdiplam** (small-molecule SMN2 splicing modulator); **eteplirsen/golodirsen/casimersen/viltolarsen** (DMD exon skipping). Design references: Hua 2008 *AJHG*; Roberts et al 2023 *Nat Rev Drug Discov* 22:917 (DMD therapeutic approaches).

## Per-Tool Failure Modes

### SpliceAI: 50nt Window Limitation

**Trigger:** Variant deep in an intron (>50 nt from canonical splice site).

**Mechanism:** Default precomputed scores use ±50 nt window; the model is trained on this context but pre-stored scores limit lookups.

**Symptom:** Known pathogenic deep-intronic variant scores low (<0.2); no pseudoexon detected.

**Fix:** Re-run with `-D 500` or `-D 2000`; or try CI-SpliceAI (SpliceAI retrained on curated GENCODE splice sites) as a second predictor.

### SpliceAI: Tissue Agnosticism

**Trigger:** Variant in a tissue-specific gene (NEFM in neurons, MAPT brain, DMD muscle isoforms).

**Mechanism:** SpliceAI is trained on aggregate GENCODE annotation; tissue-specific events with weak constitutive use score low.

**Symptom:** Tissue-specific pathogenic variant has low SpliceAI delta; functional impact still observed in target tissue.

**Fix:** Use Pangolin for tissue-aware prediction; or SpliceTransformer; require RNA validation in disease-relevant tissue.

### Pangolin: Out-of-Training Tissue

**Trigger:** Disease tissue not represented in Pangolin's 4-species, 4-tissue (Cardoso-Moreira 2019 developmental) training set.

**Mechanism:** Pangolin extrapolates poorly to tissues outside training distribution.

**Symptom:** Pangolin score uncalibrated for queried tissue; doesn't agree with patient RNA-seq from that tissue.

**Fix:** Fall back to SpliceAI for tissues not in Pangolin training; or run patient RNA-seq directly.

### MMSplice: Atypical Events

**Trigger:** Variant affecting a non-cassette event (MXE, complex multi-junction, AFE/ALE).

**Mechanism:** MMSplice modular model is trained primarily on cassette exon events.

**Symptom:** MMSplice ΔPSI doesn't match other predictors or empirical data for non-cassette events.

**Fix:** Use SpliceAI for non-cassette events; restrict MMSplice to cassette exon contexts.

### CADD-Splice: Loss of Component Information

**Trigger:** Wanting to know which sub-component drove a high CADD-Splice score.

**Mechanism:** CADD-Splice combines SpliceAI + MMSplice + CADD into a single C-score; sub-component contributions are abstracted.

**Symptom:** "High CADD-Splice score but unclear why."

**Fix:** Run SpliceAI and MMSplice separately to see which contributed.

### Branchpoint Variants: Low Information Motif

**Trigger:** Variant in the BPS region (-18 to -40 from 3'ss).

**Mechanism:** BPS motif (yUnAy) has low information content; CNNs struggle to learn the consensus.

**Symptom:** Confirmed BPS variant scores SpliceAI delta <0.2 despite functional disruption.

**Fix:** Use BPHunter (Zhang 2022 *PNAS*) for genome-wide branchpoint screening; require RNA validation.

## Population Database Lookup

| Database | Use for |
|----------|---------|
| gnomAD v4 | Allele frequency; SpliceAI annotations integrated |
| ClinVar | Existing classifications; SpliceAI integrated since 2020 |
| SpliceVarDB | Curated splice variants with experimental RNA validation |
| dbNSFP4 | Pre-computed splice scores aggregated |
| Recount3 | Tissue-specific PSI lookups from public RNA-seq |
| GTEx sQTL v8 | Tissue-specific splicing QTLs across 49 tissues |
| MaveDB | Splice MAVE results (e.g. BRCA1 saturation; Findlay 2018 *Nature*) |

Always check ClinVar first for existing classifications; cross-reference with gnomAD for population frequency before committing to PP3/PP4.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `spliceai: tensorflow not found` | TensorFlow not installed | `pip install tensorflow>=2.0` separately |
| `spliceai: chrom not in reference` | VCF chrom name mismatch (chr1 vs 1) | `bcftools annotate --rename-chrs chr_map.txt` |
| `pangolin: no annotations found for variant` | gffutils db doesn't contain queried gene | Rebuild gffutils db with comprehensive GENCODE GFF3 |
| `mmsplice: variant outside any cassette event` | MMSplice model assumes cassette context | Use SpliceAI for non-cassette events |
| `SpliceVault: variant not found` | Variant outside common splice sites in 300K-RNA database | Use SpliceAI for prediction (no empirical baseline available) |
| `VariantValidator: invalid HGVS` | Wrong reference transcript or build | Specify NM_*.* version explicitly |

## Common Pitfalls

- **Using SpliceAI score alone for clinical reporting** — must combine with concordant predictors and ideally RNA validation; ClinGen SVI requires this for non-canonical positions.
- **50nt window for deep intronic variants** — pseudoexon-creating variants 100-2000 nt deep are systematically missed.
- **Tissue-agnostic prediction for tissue-specific genes** — use Pangolin or SpliceTransformer when tissue context matters (NEFM, MAPT, DMD isoforms).
- **Branchpoint variants** — all current predictors are weak here. Use BPHunter for branchpoint screening.
- **Forgetting NMD direction** — confirmed splice disruption needs NMD-status check. Last-exon PTCs escape NMD and can be dominant-negative or gain-of-function.
- **In-silico-only PVS1 application** — PVS1 for non-canonical positions requires functional or strong computational evidence; SpliceAI alone is supporting (PP3), not very strong.
- **Trusting LLMs for variant interpretation** — use as orchestrators on top of SpliceAI/VariantValidator/ClinVar; all clinical-grade calls require human expert sign-off.
- **Skipping HGVS validation** — invalid HGVS leads to silent reference-transcript mismatches; always run VariantValidator first.

## Quality Thresholds

| Metric | Recommendation | Source |
|--------|----------------|--------|
| Default SpliceAI window | -D 50 (clinical screening) | Jaganathan 2019 |
| Deep-intronic SpliceAI window | -D 500-2000 (unsolved Mendelian) | Convention (verify current literature) |
| ACMG PP3 (supporting) | SpliceAI delta >= 0.2 | Walker 2023 *AJHG* (apply at supporting weight) |
| ACMG BP4 (supporting) | SpliceAI delta <= 0.1 | Walker 2023 *AJHG* |
| SpliceAI higher-precision cutoffs | 0.5 / 0.8 raise precision, NOT ACMG strength | Jaganathan 2019 (not ClinGen graded tiers) |
| Off-target ASO match | <=16/20 nt to any non-target transcript | Design convention |
| Concordance for high-confidence | 2/3 predictors above PP3 threshold | Pragmatic |

## Related Skills

- splicing-qc - MaxEntScan + library QC for confirming predicted impact
- splicing-quantification - Empirical PSI from RNA-seq to validate predictions
- outlier-splicing-detection - FRASER2/DROP for RNA-seq confirmation in clinical samples
- variant-calling/clinical-interpretation - Broader ACMG/AMP variant interpretation framework
- variant-calling/variant-annotation - VEP plugin integration for SpliceAI

## References

- Jaganathan et al 2019 *Cell* - SpliceAI
- Zeng & Li 2022 *Genome Biol* - Pangolin
- Cheng et al 2019 *Genome Biol* - MMSplice
- Cheng et al 2021 *Genome Biol* - MTSplice (tissue MMSplice)
- You et al 2024 *Nat Commun* 15:9129 - SpliceTransformer
- Strauch et al 2022 *PLoS One* 17:e0269159 - CI-SpliceAI extended window
- Smith & Kitzman 2023 *Genome Biol* 24:294 - SpliceAI/Pangolin MPSA benchmark
- Rentzsch et al 2021 *Genome Med* - CADD-Splice
- Dawes et al 2023 *Nat Genet* - SpliceVault
- Walker et al 2023 *Am J Hum Genet* - ClinGen SVI splicing recommendations
- Riepe et al 2021 *Hum Mutat* 42:799 - SpliceAI in clinical pipelines (Riepe TV et al)
- Abou Tayoun et al 2018 *Hum Mutat* - PVS1 decision tree
- Richards et al 2015 *Genet Med* - ACMG/AMP framework
- den Dunnen et al 2016 *Hum Mutat* - HGVS standard
- Zhang et al 2022 *PNAS* (PMID 36306325) - BPHunter for branchpoints
- Hua et al 2008 *AJHG* - ISS-N1 / nusinersen mechanism
- Roberts et al 2023 *Nat Rev Drug Discov* 22:917-934 - DMD therapeutic approaches (exon-skipping ASOs)
- Findlay et al 2018 *Nature* - BRCA1 saturation genome editing (MAVE)
