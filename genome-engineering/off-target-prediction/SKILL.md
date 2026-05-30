---
name: bio-genome-engineering-off-target-prediction
description: Nominates and assesses CRISPR off-target sites genome-wide. Enumerates candidate sites by mismatch and bulge tolerance with Cas-OFFinder/CRISPRitz, ranks them with the published CFD score (SpCas9-only, relative ranker) or MIT/CRISTA/energy models, runs variant-aware screening against gnomAD/individual genomes (CRISPRme), and frames the empirical genome-wide discovery assays (GUIDE-seq, CIRCLE-seq, CHANGE-seq, DISCOVER-seq, Digenome-seq) and high-fidelity nuclease choice (HiFi Cas9, Sniper-Cas9, eSpCas9, SpCas9-HF1). Use when assessing guide RNA specificity, choosing among candidate guides, screening a therapeutic guide against population variation, or planning empirical off-target validation. Distinguishes predicted vs detected vs validated. On-target activity scoring and deaminase (Cas-independent) base/prime-editor off-targets are separate skills.
tool_type: mixed
primary_tool: Cas-OFFinder
---

## Version Compatibility

Reference examples tested with: Cas-OFFinder 3.0+, pandas 2.2+, Python 3.10+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

Results depend on inputs far more than tool versions: the candidate list is bounded by the **reference genome build, the mismatch/bulge tolerance, and the PAM pattern** searched, not by the Cas-OFFinder version. The **CFD matrix is SpCas9/NGG-specific and a relative ranker, not a calibrated cutting probability**. Load the published CFD tables (Doench 2016 / CRISPOR distribution) rather than hand-typing values. Cas-OFFinder is the maintained `snugel/cas-offinder` repository (native DNA/RNA bulge support from v3.0.0).

# Off-Target Prediction

**"Check my guide for off-targets"** -> Enumerate candidate sites genome-wide by mismatch/bulge tolerance, rank them by a per-site score, decide whether in-silico is sufficient or empirical discovery is required, and report each claim at the right rung: predicted, detected, or validated.
- CLI: `cas-offinder input.txt G output.txt` enumerates sites (no ranking)
- Python: CFD scoring from the published mismatch/PAM tables; aggregate specificity
- Web/CLI: `CRISPRme` for variant-aware (gnomAD + individual) nomination; `CRISPOR` to aggregate

## The Single Most Important Modern Insight -- in-silico enumeration nominates *candidates*; it does not measure *which sites are cut*

The naive model -- "search the genome within N mismatches, score by CFD, the high scorers are my off-targets" -- is wrong in three structural ways no better scoring fixes:

1. **Mismatch count is not cleavage.** A 2-mismatch site in closed chromatin may never be cut; a 3-mismatch site in open chromatin near an active promoter is. Cellular cutting depends on chromatin, dose, and exposure time -- invisible to a sequence search.
2. **Bulges and non-canonical PAMs are routinely missed.** Real validated off-targets occur with 1-2 nt DNA/RNA bulges and at NAG/NGA PAMs; fixed-alignment mismatch-only search misses them. **The failure is silent** -- a clean report looks identical whether the guide is specific or the search just couldn't see the off-target.
3. **CFD is a narrow, SpCas9-only relative ranker.** A CFD of 0.08 is not "8% chance of cutting"; comparing two guides' aggregate scores is fine, reading an absolute CFD as a safety threshold is not.

The corollary, and the central professor-level point: **in-silico lists overlap only partially with empirically validated off-targets, and the empirical genome-wide assays disagree with *each other* too.** No single method is authoritative. Off-target evidence escalates: **predicted -> detected by an unbiased assay -> validated by targeted amplicon deep-seq.** Conflating these rungs is the field's most common error. Therapeutic-grade assessment is *triangulation* (variant-aware in-silico + >=2 orthogonal empirical assays + amplicon validation + a structural readout), never one tool's output.

## In-Silico Taxonomy -- enumerate, then score, then aggregate

| Layer | Tool | Citation | Role / caveat |
|-------|------|----------|---------------|
| Enumerate | **Cas-OFFinder** | Bae 2014 *Bioinformatics* 30:1473 | exhaustive, alignment-free, GPU; **DNA/RNA bulges** (native v3.0.0); returns sites, **no ranking** |
| Enumerate (variant) | CRISPRitz | Cancellieri 2020 *Bioinformatics* 36:2001 | enumerates against genome **+ a VCF** of variants, with bulges; backend of CRISPRme |
| Enumerate (scale) | GuideScan2 | Schmidt 2025 *Genome Biol* 26:41 | genome-wide specificity databases (NOT Nat Biotechnol) |
| Score (per-site) | **CFD** | Doench 2016 *Nat Biotechnol* 34:184 | position x mismatch-type matrix x PAM penalty; **SpCas9/NGG only**, poor on bulges; de facto standard |
| Score (legacy) | MIT/Hsu | Hsu 2013 *Nat Biotechnol* 31:827 | original; **deprecated/flawed** -- report, don't lead with it |
| Score (ML) | CRISTA; Elevation | Abadi 2017; Listgarten 2018 | Elevation folds in chromatin accessibility |
| Aggregate | **CRISPOR**; **CRISPRme** | Concordet 2018 *NAR* 46:W242; Cancellieri 2023 *Nat Genet* 55:34 | CRISPOR = research one-stop; CRISPRme = variant-aware therapeutic nominator |

The unifying caveat: every score is bounded by the enumerator's coverage -- **if the enumerator didn't propose a site (bulge, distal PAM, beyond the mismatch cutoff), no scorer will ever flag it.**

## Empirical Discovery Assays -- each has a characteristic bias; concordance is partial

| Assay | Citation | Class | Bias |
|-------|----------|-------|------|
| **CIRCLE-seq** | Tsai 2017 *Nat Methods* 14:607 | in-vitro (cell-free) | **over-calls** (no chromatin); most sensitive *candidate generator* |
| **CHANGE-seq** | Lazzarotto 2020 *Nat Biotechnol* 38:1317 | in-vitro | scalable CIRCLE-seq; same over-call caveat |
| Digenome-seq | Kim 2015 *Nat Methods* 12:237 | in-vitro (WGS) | unbiased but depth-limited, expensive |
| SITE-seq | Cameron 2017 *Nat Methods* 14:600 | in-vitro | concentration series ranks sensitivity |
| **GUIDE-seq** | Tsai 2015 *Nat Biotechnol* 33:187 | cell-based (dsODN tag) | physiological; **misses rare sites**, cell-type-specific, hard in primary/RNP |
| **DISCOVER-seq** | Wienert 2019 *Science* 364:286 | cell-based (MRE11 ChIP, in situ) | tag-free, works in vivo; depends on transient MRE11 occupancy |
| TTISS | Schmid-Burgk 2020 *Mol Cell* 78:794 | cell-based | high-throughput; benchmarks fidelity variants |

**The load-bearing reality:** in-vitro assays over-call (high sensitivity, low cellular specificity); cell-based assays under-call rare sites and are cell-type-dependent (K562 yields far more hits than HEK293 for the same guide). Cross-method discordance is *information*, not noise -- sites found by both are high-confidence; in-vitro-only sites are likely chromatin-protected. The defensible workflow is the **VIVO** logic (Akcakaya 2018): sensitive in-vitro generator -> cell-based assay in the relevant cell type -> amplicon validation.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Research knockout / screen (some off-target tolerable) | CRISPOR or GuideScan2 to pick the most specific guide; Cas-OFFinder (<=4 mm + bulges) to eyeball top sites | in-silico is sufficient when being wrong is cheap |
| Choosing among candidate guides | rank by **aggregate CFD specificity** (compare guides, not absolute safety) | specificity is a separate axis from on-target activity (-> grna-design) |
| Human therapeutic guide | **variant-aware** CRISPRme vs gnomAD (+ patient genome), bulges on | a common ancestry-enriched SNP can create a real off-target (rs114518452 / BCL11A) |
| Therapeutic, choosing the nuclease | high-fidelity variant **in the delivery format actually used** | RNP -> HiFi Cas9 (R691A) or Sniper-Cas9; plasmid-tuned variants can lose their edge as RNP |
| Therapeutic validation | >=2 orthogonal empirical assays -> amplicon deep-seq with stated LoD -> structural readout | predicted != detected != validated; amplicons miss large deletions/translocations |
| Base/prime-editor off-targets | this skill covers **Cas-dependent** only | deaminase (Cas-independent) DNA/RNA off-targets -> base-editing-design / prime-editing-design |

## High-Fidelity Nucleases -- often a bigger lever than guide reselection

| Variant | Citation | Note |
|---------|----------|------|
| eSpCas9(1.1) | Slaymaker 2016 *Science* 351:84 | neutralizes non-target-strand contacts; characterized mostly as plasmid |
| SpCas9-HF1 | Kleinstiver 2016 *Nature* 529:490 | weakens 4 Cas9-DNA H-bonds; plasmid-characterized |
| HypaCas9 | Chen 2017 *Nature* 550:407 | conformational proofreading gate |
| evoCas9 | Casini 2018 *Nat Biotechnol* 36:265 | ~79x fidelity; ~90% residual on-target |
| **Sniper-Cas9** | Lee 2018 *Nat Commun* 9:3048 | high specificity **and works as RNP** |
| **HiFi Cas9 (R691A)** | Vakulskas 2018 *Nat Med* 24:1216 | single mutation; **the RNP-favored therapeutic variant** |

Two tacit points: (1) **delivery format matters** -- eSpCas9/HF1 can lose their fidelity advantage delivered as a high transient RNP bolus; HiFi Cas9 and Sniper-Cas9 stay specific and active as RNP. (2) **Fidelity has a guide-dependent on-target tax** -- a variant clean and active on guide A can be nearly dead on guide B. Pick the variant, then **test it on the target guide in the intended delivery format**; transferability is not assumable.

## Variant-Aware Screening (reference-only is a clinical liability)

A patient is not GRCh38. A common SNP can restore a PAM or remove the protective mismatch at a near-target site, *creating* an off-target that exists only in some individuals -- and because variant frequencies differ by ancestry, reference-only screening systematically misses off-targets common in under-represented populations. For a human therapeutic guide, an off-target check must expand from "checked off-targets" to "checked off-targets **variant-aware, across ancestries**" -- run CRISPRme against gnomAD (and the treated individual's genome).

## Enumerate Candidate Sites with Cas-OFFinder

**Goal:** Generate the genome-wide candidate-site list for one or more guides, including bulges and relaxed PAMs.

**Approach:** Write the Cas-OFFinder input file -- genome path, an optional DNA/RNA bulge line (v3.0.0+), a pattern with N's at guide positions and the PAM (use `NRG` to also catch NAG/NGG), then one query line per guide (guide bases + N's for the PAM positions, same length as the pattern) with its mismatch tolerance. Run on GPU if available. The output is a flat site list with mismatch counts -- it is a hypothesis set to score downstream, not a verdict.

```bash
# input.txt
# /path/to/genome_dir            # directory of FASTA (Cas-OFFinder indexes it)
# 2 2                            # DNA bulge, RNA bulge (omit this line for no-bulge search)
# NNNNNNNNNNNNNNNNNNNNNRG        # 20 N (guide) + NRG PAM -> also catches NAG
# GGCCGACCTGTCGCTGACGCNNN 4      # query: 20 guide bases + NNN (PAM positions), <=4 mismatches
cas-offinder input.txt G output.txt   # G=GPU, C=CPU, A=auto
```

## Score Candidates with the Published CFD Tables

**Goal:** Rank candidate sites by relative cleavage propensity and compute an aggregate guide-specificity score for comparing guides.

**Approach:** Do NOT hand-type the CFD matrix. Load the published Doench 2016 tables (`mismatch_score.pkl`, `pam_scores.pkl` -- they ship with CRISPOR and the Doench code), take the product of per-position mismatch penalties x the PAM penalty for each site, and aggregate as `100/(1 + sum(CFD))` with per-site CFDs on a 0-1 scale (the CRISPOR specificity formulation; equivalently `10000/(100 + 100*sum)`). Compare aggregate scores *among candidate guides*; never read an absolute CFD as a safety guarantee. (See `examples/off_target_analysis.py`.)

```python
import pickle

def load_cfd_tables(mismatch_pkl, pam_pkl):
    '''Load the published Doench 2016 CFD tables (distributed with CRISPOR) -- do not fabricate.'''
    with open(mismatch_pkl, 'rb') as f:
        mismatch = pickle.load(f)   # keys like 'rA:dG,3' -> penalty
    with open(pam_pkl, 'rb') as f:
        pam = pickle.load(f)        # keys like 'AG' -> penalty
    return mismatch, pam
```

## Structural Consequences Amplicon Panels Miss

Validating only with a short amplicon at each predicted site systematically misses the large-scale outcomes that are often the real safety concern:
- **Large deletions / complex rearrangements at the on-target** (Kosicki 2018 *Nat Biotechnol* 36:765) -- kilobase deletions whose alleles often drop out of the PCR, so the amplicon reads back *more wild-type than it is*.
- **Chromosomal translocations** between on- and off-target (or multiplexed) cuts -- need junction-capture (PEM-seq, UDiTaS, HTGTS, CAST-seq), not amplicon panels. A "clean" amplicon panel does **not** certify the absence of these.

## The Evidence Ladder & Limit of Detection

| Rung | Meaning | Method | Floor |
|------|---------|--------|-------|
| Predicted | sequence-similar candidate | Cas-OFFinder/CRISPOR/CRISPRme | n/a |
| Detected | nuclease acts there (unbiased) | GUIDE-/CIRCLE-/DISCOVER-/CHANGE-seq | assay-dependent |
| Validated | confirmed editing + allele frequency | targeted amplicon deep-seq (rhAmpSeq) + CRISPResso2 | ~0.1-0.5% (~0.1% with UMI/duplex) |

**"Not detected" means "below the LoD," never "zero."** State the LoD: 0.05% editing is irrelevant for a research knockout but is ~50,000 mis-edited cells in a 10^8-cell therapy.

## Per-Method Failure Modes

### "I ran Cas-OFFinder, so I checked my off-targets"
**Trigger:** treating an in-silico mismatch list as a verdict. **Mechanism:** the search sees sequence homology, not cellular cutting; bulges/chromatin/sub-LoD editing are invisible. **Symptom:** clean report, real off-targets later. **Fix:** in-silico chooses *which guide to try*; validate empirically when being wrong matters.

### Clean amplicon panel read as "safe"
**Trigger:** amplicon-seq only at predicted sites. **Mechanism:** large deletions drop out of PCR (Kosicki 2018); the panel can't discover sites the in-silico search missed. **Symptom:** falsely clean. **Fix:** feed the panel from an unbiased discovery assay; add a structural/translocation readout; state the LoD.

### One assay treated as ground truth
**Trigger:** "CIRCLE-seq is the gold standard." **Mechanism:** in-vitro over-calls, cell-based under-calls rare/cell-type-specific sites; they disagree by design. **Symptom:** over- or under-stated risk. **Fix:** triangulate (in-vitro generator + cell-based in the relevant cell type + validation).

### High-fidelity nuclease recommended without delivery context
**Trigger:** "use eSpCas9 for specificity." **Mechanism:** plasmid-tuned variants can lose the advantage as RNP; the on-target tax is guide-dependent. **Symptom:** lost activity or lost specificity. **Fix:** RNP -> HiFi Cas9/Sniper-Cas9; test the variant on the target guide in the intended format.

### Reference-only screen for a therapeutic guide
**Trigger:** searching GRCh38 only. **Mechanism:** ancestry-enriched SNPs create/destroy off-targets. **Symptom:** a real, population-specific off-target missed. **Fix:** CRISPRme vs gnomAD + the individual's genome.

### Bulge / non-canonical-PAM off-target missed
**Trigger:** mismatch-only search at NGG. **Mechanism:** the mismatch-count abstraction can't represent a 1 nt bulge or an NAG site. **Symptom:** assay finds an off-target the search "missed." **Fix:** enable bulges (Cas-OFFinder v3) and search a relaxed PAM (NRG).

## Quantitative Thresholds

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Mismatch cutoff | <=4 typical (CRISPOR default); up to 6 for thorough | meaningful cutting rare beyond 4-5 mm, but bulges/variants rescue more-distant sites |
| Bulge size | up to ~2 (DNA + RNA) | real validated off-targets occur with 1-2 nt bulges |
| CFD per-site | relative ranker; attention >~0.1-0.2; high-risk near on-target | not a calibrated probability |
| Aggregate specificity (CRISPOR) | higher better; >~80 commonly "good" *for choosing guides* | research heuristic, NOT a clinical pass/fail |
| Amplicon LoD | ~0.1-0.5% (~0.1% with UMI/duplex) | below this, PCR/sequencer error dominates |
| High-fidelity on-target tax | guide- and format-dependent | always test the variant on the target guide |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Cas-OFFinder returns nothing | wrong genome path / pattern-query length mismatch | query length must equal pattern length; check the genome dir/FASTA |
| CFD scores look fabricated/wrong | hand-typed matrix | load the published `mismatch_score.pkl`/`pam_scores.pkl` |
| Assay finds an off-target the search missed | mismatch-only, NGG-only search | enable bulges; search NRG |
| "No detectable off-targets" claimed as zero | LoD not stated | report the limit of detection; absence is bounded, not absolute |

## References

- Bae S, Park J, Kim JS (2014). Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. *Bioinformatics* 30(10):1473-1475.
- Doench JG, Fusi N, Sullender M, et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nat Biotechnol* 34(2):184-191.
- Hsu PD, et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. *Nat Biotechnol* 31(9):827-832.
- Cancellieri S, et al. (2020). CRISPRitz: rapid, high-throughput and variant-aware in silico off-target site identification. *Bioinformatics* 36(7):2001-2008.
- Schmidt H, et al. (2025). Genome-wide CRISPR guide RNA design and specificity analysis with GuideScan2. *Genome Biol* 26:41.
- Abadi S, et al. (2017). A machine learning approach for predicting CRISPR-Cas9 cleavage efficiencies and patterns (CRISTA). *PLoS Comput Biol* 13(10):e1005807.
- Listgarten J, et al. (2018). Prediction of off-target activities for the end-to-end design of CRISPR guide RNAs (Elevation). *Nat Biomed Eng* 2(1):38-47.
- Concordet JP, Haeussler M (2018). CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids Res* 46(W1):W242-W245.
- Yan J, et al. (2020). Benchmarking and integrating genome-wide CRISPR off-target detection and prediction. *Nucleic Acids Res* 48(20):11370-11379.
- Tsai SQ, et al. (2015). GUIDE-seq enables genome-wide profiling of off-target cleavage by CRISPR-Cas nucleases. *Nat Biotechnol* 33(2):187-197.
- Kim D, et al. (2015). Digenome-seq: genome-wide profiling of CRISPR-Cas9 off-target effects in human cells. *Nat Methods* 12(3):237-243.
- Cameron P, et al. (2017). Mapping the genomic landscape of CRISPR-Cas9 cleavage (SITE-seq). *Nat Methods* 14(6):600-606.
- Tsai SQ, et al. (2017). CIRCLE-seq: a highly sensitive in vitro screen for genome-wide CRISPR-Cas9 nuclease off-targets. *Nat Methods* 14(6):607-614.
- Wienert B, et al. (2019). Unbiased detection of CRISPR off-targets in vivo using DISCOVER-seq. *Science* 364(6437):286-289.
- Lazzarotto CR, et al. (2020). CHANGE-seq reveals genetic and epigenetic effects on CRISPR-Cas9 genome-wide activity. *Nat Biotechnol* 38(11):1317-1327.
- Schmid-Burgk JL, et al. (2020). Highly Parallel Profiling of Cas9 Variant Specificity (TTISS). *Mol Cell* 78(4):794-800.e8.
- Akcakaya P, et al. (2018). In vivo CRISPR editing with no detectable genome-wide off-target mutations (VIVO). *Nature* 561:416-419.
- Scott DA, Zhang F (2017). Implications of human genetic variation in CRISPR-based therapeutic genome editing. *Nat Med* 23:1095-1101.
- Lessard S, et al. (2017). Human genetic variation alters CRISPR-Cas9 on- and off-targeting specificity at therapeutically implicated loci. *PNAS* 114(52):E11257-E11266.
- Cancellieri S, et al. (2023). Human genetic diversity alters off-target outcomes of therapeutic gene editing (CRISPRme). *Nat Genet* 55(1):34-43.
- Slaymaker IM, et al. (2016). Rationally engineered Cas9 nucleases with improved specificity (eSpCas9). *Science* 351(6268):84-88.
- Kleinstiver BP, et al. (2016). High-fidelity CRISPR-Cas9 nucleases with no detectable genome-wide off-target effects (SpCas9-HF1). *Nature* 529(7587):490-495.
- Chen JS, et al. (2017). Enhanced proofreading governs CRISPR-Cas9 targeting accuracy (HypaCas9). *Nature* 550(7676):407-410.
- Casini A, et al. (2018). A highly specific SpCas9 variant is identified by in vivo screening in yeast (evoCas9). *Nat Biotechnol* 36(3):265-271.
- Lee JK, et al. (2018). Directed evolution of CRISPR-Cas9 to increase its specificity (Sniper-Cas9). *Nat Commun* 9:3048.
- Vakulskas CA, et al. (2018). A high-fidelity Cas9 mutant delivered as a ribonucleoprotein complex enables efficient gene editing in human hematopoietic stem and progenitor cells (HiFi Cas9). *Nat Med* 24(8):1216-1224.
- Kosicki M, Tomberg K, Bradley A (2018). Repair of double-strand breaks induced by CRISPR-Cas9 leads to large deletions and complex rearrangements. *Nat Biotechnol* 36:765-771.
- Clement K, et al. (2019). CRISPResso2 provides accurate and rapid genome editing sequence analysis. *Nat Biotechnol* 37(3):224-226.

## Related Skills

- grna-design - Design and on-target-score guides before the specificity check (a separate axis)
- base-editing-design - Owns the deaminase (Cas-independent) DNA/RNA off-target classes
- prime-editing-design - pegRNA off-target considerations and PE3 nicking-guide specificity
- crispr-screens/crispresso-editing - Quantify and validate editing at candidate sites from amplicon reads
- variant-calling/variant-annotation - Annotate whether off-targets hit genes/pathogenic loci
- genome-intervals/bed-file-basics - Intersect off-target sites with exons/oncogenes for prioritization
- database-access/ncbi-datasets-cli - Download the reference genome for the search
