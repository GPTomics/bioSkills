---
name: bio-crispr-screens-library-design
description: Designs pooled sgRNA libraries for CRISPR knockout, interference (CRISPRi), activation (CRISPRa), Cas12a multiplex, base-editor, and prime-editor screens. Covers on-target scoring (Rule Set 2, Azimuth, DeepSpCas9, CRISPRon), off-target scoring (CFD, MIT), TSS-relative positioning for CRISPRi/a (Horlbeck, Dolcetto, Calabrese), PAM-variant chemistries, control-guide composition, oligo cloning architecture, and library QC. Use when choosing a genome-wide library (GeCKOv2 vs Avana vs Brunello vs TKOv3 vs Inzolia), designing a focused or paralog-focused custom library, picking CRISPRi vs CRISPRa TSS windows, deciding control-guide proportions, or diagnosing library skew and dropout in a freshly cloned pool.
tool_type: mixed
primary_tool: CRISPOR
---

## Version Compatibility

Reference examples tested with: CRISPOR 5.01+, BioPython 1.83+, pandas 2.2+, numpy 1.26+, Azimuth 2.0+ (Doench 2016), CRISPRon 1.0+ (Xiang 2021), DeepSpCas9 1.0+ (Kim 2019).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show crispor` then `help(...)` to check signatures
- CLI: `crispor.py --help`, `azimuth --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## sgRNA Library Design

**"Design a CRISPR library for my screen"** -> Pick a chemistry (Cas9 KO, CRISPRi, CRISPRa, Cas12a, base or prime editor), score candidate guides for on-target activity and off-target liability, position them relative to gene/TSS, add appropriate controls, lay out the oligo for synthesis, and validate the cloned pool.

- Python: `crispor.py` (web + CLI) for batch genome-wide guide scoring with CFD+MIT off-target
- Python: `azimuth` (Microsoft Research) for Rule Set 2 on-target predictions (Brunello-style)
- Python: `CRISPRon`, `DeepSpCas9` for modern deep-learning predictors
- R: `crisprDesign` (Bioconductor) for integrated annotation-aware design

## Library Chemistry Decision Tree

| Goal | Chemistry | Canonical library | Guides/gene | TSS / target window |
|------|-----------|-------------------|-------------|---------------------|
| Loss-of-function essentiality, fitness | SpCas9 KO | Brunello, TKOv3, Avana | 4 (Brunello), 4 (TKOv3), 6 (Avana) | Constitutive exons, prefer aa 5-65% from N-terminus |
| Knockdown of non-cuttable genes, dosage-sensitive | dCas9-KRAB (CRISPRi) | Dolcetto, Horlbeck v2 | 6 (Dolcetto), 5 (Horlbeck) | -50 to +300 from FANTOM5 TSS (Dolcetto) or -50 to +300 (Horlbeck) |
| Gain-of-function, gene activation | dCas9-VP64 / SAM / SunTag (CRISPRa) | Calabrese, Horlbeck-CRISPRa | 6 (Calabrese), 5 (Horlbeck) | -150 to -75 from TSS (Calabrese); -400 to -50 (Horlbeck) |
| Paralog buffering, GI screens | enAsCas12a multiplex | Inzolia, in4mer | 4-guide arrays | Constitutive exons |
| Variant function, SNV scanning | CBE / ABE | Sanson 2020 GRACE, custom | Tile editing windows | Editing window pos 4-8 from PAM-distal end |
| Precise edit, indel-free | Prime editor | Custom PRIDICT-designed | Tile pegRNAs | Anywhere with NGG PAM within 30 nt of edit |

**Fails when:**
- CRISPRi/a targeting wrong TSS: any TSS without FANTOM5 CAGE evidence is suspect; Horlbeck-rule guides outside ±300 from the true TSS have ~10x lower knockdown.
- Cas9 KO of essential paralogs: single-KO buffering masks ~half of constitutively-expressed essentials (Dede 2020); switch to Cas12a multiplex.
- Base editor over an exon-intron boundary: editing-window bystanders create splice variants instead of the intended SNV.

## On-Target Scoring: Algorithmic Taxonomy

| Predictor | Year | Training set | Strengths | Fails when |
|-----------|------|--------------|-----------|------------|
| Doench Rule Set 1 | 2014 | Flow-sorted GFP+ knockouts | Simple, interpretable | Limited training data; sub-optimal at >NGG context |
| Doench Rule Set 2 / Azimuth 2.0 | 2016 | Brunello + 1,841 guides | Gold-standard for SpCas9; basis of Brunello | Trained on dropouts; under-predicts efficacy for nuclear-localized targets |
| DeepSpCas9 | 2019 | Endogenous + reporter (12k guides) | Higher Pearson r vs experimental indel rate (~0.78) | Black-box; sensitive to chromatin/context features it wasn't trained on |
| CRISPRon | 2021 | High-throughput indel sequencing | Best for therapeutic-grade target nomination | Slow per-guide; over-fits to its specific cell line |
| DeepHF | 2019 | HF / eSpCas9 high-fidelity variants | Predicts high-fidelity enzyme variants | Not for WT SpCas9 |

**Reconciliation:** When predictors disagree, prefer the model whose training cell line matches the screen line (DeepSpCas9 was trained primarily on HCT116/HEK293T). For Brunello/TKOv3 selection, Azimuth/Rule Set 2 is sufficient because the library was built with it -- introducing a different scorer creates apples-to-oranges ranking with the original library.

## Off-Target Scoring

| Score | Year | Math | Cutoff convention |
|-------|------|------|--------------------|
| MIT (Hsu) | 2013 | Position-weighted mismatch penalty | Sum of off-target scores <50 for libraries |
| CFD (Doench) | 2016 | Position+nucleotide-specific penalty fit on Brunello | Aggregate CFD specificity score >0.2 |
| Elevation | 2018 | ML on CFD + mismatch positions | Tighter than CFD |

CFD remains the default for genome-wide library design. **Critical pitfall:** CFD penalizes only mismatches, not bulges; for ≤1 mismatch + 1-bp bulge off-targets, validate empirically with GUIDE-seq or CIRCLE-seq. CRISPOR aggregates MIT + CFD + Elevation in a single output.

## Score and Rank sgRNAs for a Target Gene

**Goal:** Generate ranked sgRNA candidates for a single gene, jointly scored on on-target activity (Rule Set 2 / Azimuth) and off-target liability (CFD).

**Approach:** Identify all PAM-adjacent 20-nt protospacers in the target gene's coding sequence, retain only those in the first 5-65% of the protein (constitutive-exon convention from Brunello), filter on GC 30-70% and absence of poly-T (≥4 Ts terminates U6), call Azimuth for on-target and CRISPOR for off-target, and select the top N satisfying both criteria.

```python
import re
import pandas as pd
import numpy as np
from Bio.Seq import Seq

def find_sgrna_candidates(cds_sequence, pam='NGG', guide_length=20):
    '''Return all protospacer candidates with PAM coordinates on + strand.
    Caller must filter by exon position and Azimuth/CFD score.'''
    pam_pattern = re.compile(f'(?=([ACGT]{{{guide_length}}}{pam.replace("N", "[ACGT]")}))')
    candidates = []
    for strand, seq in [('+', cds_sequence), ('-', str(Seq(cds_sequence).reverse_complement()))]:
        for m in pam_pattern.finditer(seq):
            spacer = m.group(1)[:guide_length]
            if 'TTTT' in spacer or spacer.count('G') + spacer.count('C') not in range(6, 15):
                continue
            candidates.append({'spacer': spacer, 'strand': strand,
                               'pos_in_cds': m.start() if strand == '+' else len(seq) - m.start() - 23,
                               'gc_frac': (spacer.count('G') + spacer.count('C')) / guide_length})
    return pd.DataFrame(candidates)

def annotate_exon_position(candidates_df, cds_length):
    '''Filter to protospacers within first 5-65% of CDS (Brunello convention).
    Reason: N-terminal indels truncate protein; very-N-terminal hits alt initiation;
    C-terminal hits miss functional domains (Doench 2016 Nat Biotech).'''
    lo, hi = 0.05 * cds_length, 0.65 * cds_length
    return candidates_df[(candidates_df['pos_in_cds'] >= lo) & (candidates_df['pos_in_cds'] <= hi)].copy()
```

## CRISPRi / CRISPRa TSS Targeting

**Goal:** Position guides relative to the empirical TSS for maximum knockdown (CRISPRi) or activation (CRISPRa).

**Approach:** Resolve TSS from FANTOM5 CAGE peaks (highest-ranked peak per gene; fall back to Ensembl/RefSeq if absent), define the modality-specific window, score candidate spacers in that window with Rule Set 2 plus the Horlbeck/Sanson CRISPRi/a-tailored rules, and select 5-6 guides per gene biased toward the window center.

```python
def crispri_window(tss_coord, strand='+'):
    '''Dolcetto convention: -50 to +300 from FANTOM5 highest-rank CAGE peak.
    Reason: dCas9-KRAB knockdown is maximal when the spacer sits 25-75 bp
    downstream of TSS (Horlbeck 2016 eLife); widened to ±300 to cover poorly-
    annotated TSSs and broad CpG-island promoters (Sanson 2018 Nat Comm).'''
    if strand == '+':
        return (tss_coord - 50, tss_coord + 300)
    return (tss_coord - 300, tss_coord + 50)

def crispra_window(tss_coord, strand='+'):
    '''Calabrese convention: -150 to -75 upstream of TSS.
    Reason: dCas9-VP64 (and SAM, SunTag) activate maximally when bound
    just upstream of Pol II loading. Horlbeck CRISPRa uses -400 to -50
    (broader, lower per-guide signal). For SAM, prefer Calabrese tightness;
    for SunTag, Horlbeck width is acceptable.'''
    if strand == '+':
        return (tss_coord - 150, tss_coord - 75)
    return (tss_coord + 75, tss_coord + 150)
```

**Critical nuance:** Cell-type-specific TSSs differ from the FANTOM5 consensus in ~15% of genes. For tissue-specific screens (e.g., neuron, hepatocyte), re-derive TSSs from a matched CAGE / GRO-seq / PRO-seq dataset before locking guide positions, or knockdown efficiency drops several-fold. The single most common cause of "weak" CRISPRi hits is mis-positioned guides against an alternative TSS.

## Genome-Wide Library Selection

| Library | Year | Modality | Size (genes x guides) | sgRNA rules | Notable |
|---------|------|----------|-----------------------|-------------|---------|
| GeCKOv2 | 2014 | Cas9 KO | ~19k x 6 (~123k) | Rule Set 1 + minimal off-target | Older; legacy datasets still use it |
| Avana | 2016 | Cas9 KO | ~18k x 4 (~73k) | Rule Set 2 + CFD | DepMap until 2021 |
| Brunello | 2016 | Cas9 KO | ~19k x 4 (~77k) | Rule Set 2 + CFD | Modern standard for new screens |
| TKOv3 | 2017 | Cas9 KO | ~18k x 4 (~71k) | Hart on/off-target | Bagel/BAGEL2-optimized |
| Humagne | 2021 | Cas9 KO | ~19k x ~3 (~58k) | Updated rules + paralog awareness | Minimal optimized |
| Horlbeck CRISPRi v2 | 2016 | dCas9-KRAB | ~18k x 5 (~104k) | Horlbeck CRISPRi rules | First-gen, still widely used |
| Dolcetto | 2018 | dCas9-KRAB | ~19k x 6 (~115k) | Horlbeck + Rule Set 2 | Modern CRISPRi standard |
| Horlbeck CRISPRa | 2016 | dCas9-VP64 | ~18k x 5 (~104k) | Horlbeck CRISPRa rules | Original CRISPRa |
| Calabrese | 2018 | dCas9-VP64 | ~19k x 6 (~117k) | Tight TSS window | Modern CRISPRa standard |
| Inzolia | 2024 | enAsCas12a | ~18k x 4 (~72k arrays) | enAsCas12a rules | Paralog-pair multiplex; ~30% smaller than Cas9 |
| in4mer | 2023 | Cas12a (4-guide) | Custom | enAsCas12a multiplex | Triple/quadruple KO per cassette |

**dAUC trajectory (essentiality benchmark):** GeCKOv2 < Avana < Brunello/TKOv3 (Doench 2016 + Hart 2017). Moving from 4 to 6 sgRNAs/gene gives diminishing returns; the larger gain is moving from Rule Set 1 to Rule Set 2.

**Cost-coverage tradeoff:** A 77k-guide Brunello at 500x cells/sgRNA needs 38.5M cells in pool, scalable. A 117k-guide Calabrese at 500x needs 59M cells -- often the deciding factor against CRISPRa for difficult-to-grow lines.

## PAM Variants and Alternative Cas Enzymes

| Enzyme | PAM | Spacer length | Best for |
|--------|-----|---------------|----------|
| SpCas9 (WT) | NGG | 20 nt | Standard pooled screens; broadest library support |
| eSpCas9, SpCas9-HF1 | NGG | 20 nt | Lower off-target rate; use for therapeutic-grade nomination |
| SpCas9-NG | NG | 20 nt | Expanded targeting (~4x coverage); accept lower activity per guide |
| SpRY | NRN / NYN | 20 nt | Near-PAMless; coverage at every position; ~50% lower per-guide activity |
| SaCas9 | NNGRRT | 21 nt | AAV-packageable (small ORF); rarely used in pooled screens |
| AsCas12a, LbCas12a | TTTV | 23 nt | AT-rich regions; staggered cut; lower expression noise |
| enAsCas12a (DeWeirdt 2021) | Expanded TTTV + several non-canonical | 23 nt | Combinatorial / paralog screens |

**Decision rule:** If the screen requires every possible TSS position (saturation tiling, dense regulatory dissection), use SpRY despite lower activity; otherwise, NGG is best because the on-target predictors were trained on it.

## Control Guides

A genome-wide library should include:

| Control type | Count | Purpose |
|--------------|-------|---------|
| Non-targeting (scrambled, no genomic match) | 500-1,000 (~1% of library) | Primary null distribution for CRISPRi/a; safe baseline for normalization |
| Safe-harbor (AAVS1, ROSA26-equivalent) | 50-100 | Cas9-only: absorbs cut-toxicity baseline (matters for amplicon-correction) |
| Olfactory receptors (presumed non-expressed) | 50-100 | Second null set for orthogonal normalization |
| Reference essentials (CEGv2 subset: e.g. RPS3, RPL11, EIF3A, POLR2A) | 50-100 | Internal positive control; QC dropout signal |
| Reference non-essentials (NEGv1 subset) | 50-100 | Internal negative control; BAGEL2 calibration |

**Critical pitfall:** Using only AAVS1 as the negative control in a Cas9 screen creates a normalization baseline biased toward "any cut is bad." Always add NTCs or non-essentials so that downstream median normalization and PR-AUC against CEGv2 work without baseline-shift artifacts.

## Library Composition for Specialized Screens

**Paralog buffering (Cas12a multiplex):** Build 4-guide arrays where positions 1-2 target gene A and positions 3-4 target paralog gene B. Inzolia covers ~4,000 paralog pairs at ~72k arrays. Singleton controls (gene A alone, gene B alone) must be included to score genetic interaction = double_KO_LFC - sum(single_KO_LFC).

**Base editor screens (Sanson GRACE-style):** Tile NGG-adjacent spacers across exons; ensure editing window (positions 4-8 from PAM-distal end) lands inside coding exons; flag bystander Cs/As in the window for downstream interpretation. Restrict to 50-90% editing efficiency a priori (filter out predicted low-efficacy guides) -- see [[base-editing-analysis]].

**Tiling / regulatory dissection:** Dense (every 5-10 bp) CRISPRi or CRISPRa guides across the candidate region; CRISPRi has broader signal width (good for enhancer discovery) but Cas9-indel tiling has sharper resolution (good for pinpointing critical bases). Pair with CRISPR-SURF deconvolution.

## Oligo Design for Pooled Synthesis

**Goal:** Generate the final oligo sequence (full length 60-200 nt) ready for chip-based synthesis (Twist 12K/92K/244K, Genscript, Agilent).

**Approach:** Add subpool PCR primers (so multiple sublibraries can share a synthesis array), the BsmBI/Esp3I overhang for golden-gate cloning into LentiGuide-Puro (Addgene 52963) or LentiCRISPRv2, and append the tracrRNA scaffold if the array length permits.

```python
def build_oligo(spacer, vector='lentiGuide-Puro', subpool_idx=None):
    '''Construct final oligo for pooled synthesis.

    LentiGuide-Puro / LentiCRISPRv2 use BsmBI (Esp3I) with these overhangs:
        forward: 5'-CACCG[spacer]-3'
        reverse: 5'-AAAC[revcomp(spacer)]C-3'
    For chip synthesis, the spacer is flanked by subpool-specific PCR primers.'''
    subpool_fwd = {
        1: 'GGAAAGGACGAAACACCG',   # subpool 1 forward primer + BsmBI overhang
        2: 'GAGGCACTGGGCAGGTACCG',
    }.get(subpool_idx, 'GGAAAGGACGAAACACCG')
    scaffold_short = 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG'  # truncated tracr for short oligos
    oligo = subpool_fwd + spacer + scaffold_short
    if len(oligo) > 200:
        raise ValueError(f'Oligo length {len(oligo)} exceeds Twist 244K limit (200 nt)')
    return oligo
```

**Subpool design:** Twist 92K can be partitioned into multiple sublibraries via subpool primers; each sub-PCR amplifies its subpool, allowing one synthesis batch to serve several screens. Typical subpool size: 10k-20k oligos.

## Library QC After Cloning

| Metric | Target | Failure mode if missed |
|--------|--------|------------------------|
| sgRNA detection (>25 reads/guide in plasmid pool) | ≥99% | Founder effect: missing guides cannot be screened; dropout impossible to distinguish from missing |
| Gini coefficient of plasmid pool | <0.1 | Synthesis defects or PCR bias; pool unfit for screening at standard 500x coverage |
| Skew ratio (top 10% / bottom 10%) | <2 (good), <5 (acceptable) | Skew >5 means underrepresented guides cannot generate statistical signal even at 1000x |
| % zero-count sgRNAs in plasmid pool | <0.5% | Plasmid bottleneck during cloning; re-amplify or re-clone |
| Replicate Pearson on plasmid pool (between sequencing technical replicates) | >0.99 | Sequencing artifact, not biology |

**Plasmid pool sequencing convention:** 200-500 reads per sgRNA before any biology (i.e. 15-40M reads for a 77k Brunello). This is the baseline against which all downstream depletion is computed; sequencing the plasmid is non-negotiable.

## Failure Modes

### Wrong TSS in CRISPRi/a library

**Trigger:** Using Ensembl/RefSeq TSS instead of empirical CAGE peak for genes with broad or non-canonical promoters.
**Mechanism:** dCas9-KRAB knockdown is maximal within ±100 bp of the actual Pol II loading site; canonical annotation can be off by 1-10 kb.
**Symptom:** "Easy" essentials (RPS, RPL, EIF) show normal dropout but newer genes don't; library validates poorly against CEGv2.
**Fix:** Re-derive TSS from FANTOM5 CAGE highest-rank peak; for tissue-specific lines, use matched CAGE or GRO-seq.

### Library skew from PCR bias during amplification

**Trigger:** Amplifying the cloned plasmid pool with too many PCR cycles (>20) or with high-GC-bias polymerase.
**Mechanism:** GC-extreme guides amplify nonlinearly; high-GC guides dominate, low-GC guides drop out.
**Symptom:** Gini >0.2 on plasmid pool; sgRNAs with GC <30% systematically depleted.
**Fix:** Cap PCR at 15 cycles; use Q5 or NEBNext Ultra II (low-bias); sequence at 500x post-amp to confirm Gini.

### Oligo-synthesis dropouts in low-complexity guides

**Trigger:** Chip-synthesis errors at homopolymer runs or guides starting with GGGG.
**Mechanism:** Synthesis chemistry has higher error rate at low-complexity regions; missing oligos cannot be cloned.
**Symptom:** Specific guides absent from plasmid pool despite no design-rule violation.
**Fix:** Re-design replacement guides; for production runs, request 2-3x synthesis depth so dropouts are buffered.

### Polyclonality from high MOI

**Trigger:** Infection at MOI >0.5 to "save cells."
**Mechanism:** Poisson math: at MOI 0.3, 26% of infected cells carry ≥1 sgRNA, 4% carry ≥2; at MOI 0.5, 39% of infected cells carry ≥1, 9% carry ≥2.
**Symptom:** Hits include neutral genes that co-infect with true essentials.
**Fix:** MOI 0.3 strict; titer Cas9-positive cells specifically; re-check by qPCR of integration.

### Wrong control proportion

**Trigger:** <100 non-targeting controls in a 70k library.
**Mechanism:** Null distribution for normalization and FDR rests on the NTC variance; too few NTCs yields unstable median and inflated FDR.
**Symptom:** Erratic gene-level p-values; MAGeCK FDR fluctuates wildly between runs.
**Fix:** ~1% of library (500-1,000) NTCs; supplement with non-essential-gene controls.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| GC content | 30-70% | Doench 2016 Nat Biotech: guides outside this range have low activity |
| Poly-T avoidance | ≤3 consecutive T | U6 Pol III terminator; ≥4 Ts terminates sgRNA transcription |
| Guides per gene (Cas9) | 4 (Brunello/TKOv3 standard); 6 (Avana, older) | Doench 2016: 4 vs 6 yields ~5% AUC gain; diminishing returns above 6 |
| CRISPRi window | -50 to +300 from FANTOM5 TSS | Sanson 2018 Nat Comm (Dolcetto); Horlbeck 2016 eLife |
| CRISPRa window | -150 to -75 from TSS | Sanson 2018 (Calabrese); narrower than Horlbeck (-400 to -50) |
| NTCs in library | ~1% (500-1,000 in a 70k library) | DepMap library design notes; rule-of-thumb for stable null |
| MOI | 0.3 | Poisson: P(>=2 sgRNAs/cell) = 4% at MOI 0.3 |
| Coverage at infection | 500 cells/sgRNA | DepMap convention; 200x minimum, 1000x for noisy / in-vivo |
| Off-target CFD aggregate | <0.5 (CRISPOR convention) | Doench 2016; lower = more specific |
| Library skew (top 10% / bottom 10%) | <2 ideal, <5 acceptable | Joung 2017 Nat Protoc |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| sgRNA fails to express | Poly-T in spacer terminates U6 | Filter `TTTT` in design; this is the #1 silent failure |
| CRISPRi guide gives no knockdown | Wrong TSS used | Re-derive TSS from FANTOM5 / matched CAGE |
| Library Gini >0.3 in plasmid pool | PCR over-amplification or synthesis defect | Cap at 15 cycles; re-sequence plasmid; consider re-synthesis |
| Hits include amplified loci (e.g. ERBB2 in HER2+) | Copy-number amplicon false-essentiality | See [[copy-number-correction]] |
| Paralog gene absent from hit list despite expression | Cas9 single-KO buffering | Switch to Cas12a multiplex; see [[combinatorial-screens]] |
| Cas12a oligo doesn't cut | Forgot Cas12a's TTTV PAM is 5' of spacer, not 3' | Re-orient: PAM-then-spacer for Cas12a, opposite of Cas9 |

## References

- Doench JG et al. 2014. *Nat Biotechnol* 32:1262. Rule Set 1.
- Doench JG et al. 2016. *Nat Biotechnol* 34:184. Rule Set 2, CFD, Brunello/Avana libraries.
- Sanjana NE et al. 2014. *Nat Methods* 11:783. GeCKOv2.
- Hart T et al. 2017. *G3* 7:2719. TKOv3 library; CEGv2/NEGv1 reference essentiality gene sets.
- Sanson KR et al. 2018. *Nat Commun* 9:5416. Dolcetto + Calabrese libraries; CRISPRi/a TSS rules.
- Horlbeck MA et al. 2016. *eLife* 5:e19760. CRISPRi/a design rules; Horlbeck v2 library.
- Kim HK et al. 2019. *Nat Biotechnol* 37:657. DeepSpCas9.
- Xiang X et al. 2021. *Nat Commun* 12:3238. CRISPRon.
- Tycko J et al. 2019. *Mol Cell* 73:1273. Bassik library design.
- DeWeirdt PC et al. 2021. *Nat Biotechnol* 39:94. enAsCas12a optimization.
- Bayle JH et al. 2024. *Nat Commun* 15:3375. Inzolia / in4mer paralog library.
- Dede M et al. 2020. *Genome Biol* 21:262. Paralog buffering invisible to Cas9 single-KO.
- Joung J et al. 2017. *Nat Protoc* 12:828. Genome-wide library screen protocol.
- Chen S et al. 2015. *Cell* 160:1246. Original GeCKO design conventions.

## Related Skills

- crispr-screens/screen-qc - Validate library skew, Gini, replicate correlation
- crispr-screens/mageck-analysis - Analyze screens run with the designed library
- crispr-screens/combinatorial-screens - Cas12a multiplex / paralog-pair library design
- crispr-screens/base-editing-analysis - GRACE-style base-editor library design
- crispr-screens/prime-editing-screens - PRIDICT2-optimized pegRNA libraries
- crispr-screens/copy-number-correction - Filter amplicon-driven artifacts in cancer-cell-line screens
