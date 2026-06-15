---
name: bio-genome-engineering-grna-design
description: Designs and ranks guide RNAs (sgRNAs) for CRISPR-Cas9/Cas12a gene knockout by scanning a target for PAM sites (NGG SpCas9, NNGRRT SaCas9, TTTV Cas12a, NG SpCas9-NG, near-PAMless SpRY), enumerating candidate spacers, applying hard filters (Pol-III TTTT terminator, 5' G, GC), ranking on-target activity with the context-appropriate model (Rule Set 2/Azimuth for U6/lentiviral, CRISPRscan for T7/embryo, DeepHF for high-fidelity variants, DeepCpf1 for Cas12a), and predicting the indel/frameshift outcome (Bae out-of-frame score, inDelphi, FORECasT, Lindel). Use when selecting sgRNAs to knock out a gene, choosing a nuclease/PAM for a constrained locus, picking which exon to target, or shortlisting guides before an off-target check. Off-target specificity, base/prime editing, and HDR donors are separate skills.
tool_type: python
primary_tool: CRISPOR
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, CRISPOR 5.0+ (web/CLI).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

Output depends on inputs more than tool versions: on-target scores are **model-specific and not interchangeable** (a 0.7 Azimuth score is not a 0.7 CRISPRscan score), and the **valid model is set by how the guide is delivered/transcribed**, not by preference. Record the nuclease, the delivery context (U6/lentiviral vs in-vitro T7/RNP), and the reference genome build used for any off-target step.

# Guide RNA Design

**"Design guide RNAs to knock out my gene"** -> Establish the delivery context, scan the target for the nuclease's PAM on both strands, drop guides that fail hard filters, rank survivors with the context-valid on-target model, choose the cut site by exon/transcript biology, and prefer guides whose predicted indel spectrum is frameshift-rich.
- Python: enumerate PAMs and apply hard filters with `Bio.Seq` + `re`; compute a Bae-style microhomology out-of-frame score
- CLI/web: `crispor.py <genome> in.fa out.tsv` aggregates the context-appropriate on-target score + off-target nomination per genome
- Web/code: inDelphi / FORECasT / Lindel for the full repair-outcome distribution

## The Single Most Important Modern Insight -- a guide produces a reproducible indel *distribution*, not "a cut", and knockout success is a property of that distribution

Two facts that naive design ignores and that pass review constantly:

1. **On-target efficiency scores are weak, context-locked predictors.** Rule Set 2, CRISPRscan, and DeepCas9 scores correlate with measured cutting at only **Spearman ~0.4 across realistic contexts** (~0.7 is the ceiling even within one matched context; the *same* guides re-tested in another cell line correlate ~0.37-0.48). Each was trained on one assay -- U6-Pol-III lentiviral vs in-vitro T7 vs RNP -- and **does not transfer** across nuclease, delivery, promoter, cell type, or temperature (Haeussler 2016). Using CRISPRscan (T7/zebrafish-trained) to rank guides for a U6 lentiviral screen is a category error. **Rank to shortlist, then design 3-6 guides and validate** -- never trust the rank as truth.

2. **Efficient editing is not knockout.** A cut yields a *characteristic, reproducible* set of indels (Shen 2018; Allen 2019; Chen 2019); roughly **1/3 of indels are in-frame**, so a 95%-efficient guide can still leave functional protein. Worse, even a confirmed frameshift may not eliminate protein -- translation reinitiation, exon skipping, NMD escape, and transcriptional adaptation rescue ~1/3 of verified knockouts (Smits 2019; Mou 2017; El-Brolosy 2019). So the modern question is **"which guide, at which site, produces a high out-of-frame fraction in an NMD-competent, constitutive transcript region?"** -- couple an *outcome model* to *exon biology*, not just an efficiency score. Verify the knockout at the **protein** level.

## On-Target Score Taxonomy -- each model is valid for ONE context

| Model | Citation | Trained on (valid for) | Notes |
|-------|----------|------------------------|-------|
| Rule Set 1 | Doench 2014 *Nat Biotechnol* 32:1262 | U6 mammalian | superseded; origin of GC/position rules |
| **Rule Set 2 / Azimuth** | Doench/Fusi 2016 *Nat Biotechnol* 34:184 | **U6/lentiviral mammalian KO -- the default for screens & cell lines** | gradient-boosted; best U6 predictor (Haeussler 2016) |
| **CRISPRscan** | Moreno-Mateos 2015 *Nat Methods* 12:982 | **in-vitro T7 / embryo injection -- NOT U6** | wrong tool for lentiviral screens |
| DeepSpCas9 | Kim 2019 *Sci Adv* 5:eaax9249 | SpCas9 mammalian; strong transfer | CNN |
| **DeepHF** | Wang 2019 *Nat Commun* 10:4284 | conditions on the **enzyme variant** (WT, eSpCas9, HF1) | use when using a high-fidelity Cas9 |
| **DeepCpf1 / Seq-deepCpf1** | Kim 2018 *Nat Biotechnol* 36:239 | **AsCas12a** (Deep adds chromatin) | use for Cas12a, not Cas9 |

Treat any score as a **rank-and-shortlist** signal (Spearman ~0.4 across context), never an oracle.

## Nuclease & PAM Taxonomy -- expanding PAM range trades away activity/specificity

| Nuclease | PAM | Guide | Cut | When |
|----------|-----|-------|-----|------|
| **SpCas9 (WT)** | 5'-NGG-3' | 20 nt | blunt, ~3 bp 5' of PAM | default workhorse; most data, most scores |
| SaCas9 | 5'-NNGRRT-3' | ~21 nt | blunt | ~1 kb smaller -> **fits a single AAV** (Ran 2015) |
| SpCas9-NG | 5'-NG-3' | 20 nt | blunt | relaxed PAM; lower activity at many sites (Nishimasu 2018) |
| xCas9 | NG, GAA, GAT | 20 nt | blunt | broad PAM, high specificity, site-variable/modest activity (Hu 2018) |
| SpRY | near-PAMless (NRN>NYN) | 20 nt | blunt | "target anywhere"; pays in activity + off-target breadth (Walton 2020) |
| AsCas12a / LbCas12a | 5'-TTTV-3' (5' PAM) | ~20-23 nt | staggered 5' overhang | AT-rich targets; self-processing crRNA array = easy multiplexing |
| enAsCas12a | expanded (TTTV + non-canonical) | ~20-23 nt | staggered | ~2x activity + broadened range (Kleinstiver 2019) |

Default to WT-SpCas9-NGG; escalate to NG/xCas9/SpRY only when no acceptable NGG sits in the required window, and expect to validate harder (the valid on-target score and the off-target burden both change).

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Single-gene KO, NGG in an early constitutive exon | SpCas9 + Rule Set 2/Azimuth shortlist -> outcome model -> off-target | frameshift in an NMD-competent exon kills all isoforms |
| In-vitro-transcribed / embryo / RNP injection | score with **CRISPRscan**, apply T7 (not U6) filters | Rule Set 2 is invalid here; TTTT/5'G Pol-III rules do not apply |
| AT-rich target, no good NGG; or multiplex KO | Cas12a (TTTV) + DeepCpf1 | PAM availability and crRNA-array multiplexing, not on-target score, are limiting |
| AAV in-vivo delivery | SaCas9 (NNGRRT) | packaging limit dictates the compact nuclease, which dictates the PAM set |
| Functional/negative-selection screen | tile sgRNAs across the **conserved functional domain** (Shi 2015) | domain indels are LoF even in-frame -> more true nulls than 5'-exon targeting |
| Have ranked candidates, need specificity | -> off-target-prediction | on-target score does not predict specificity |
| Scale to many genes | -> crispr-screens/library-design | pooled library construction |
| Single base change / no DSB tolerated | -> base-editing-design or prime-editing-design | scarless, DSB-free; KO-by-stop also avoids indels |

## Enumerate and Filter Candidate Guides

**Goal:** Return valid candidate spacers for a target, on both strands, dropping guides that cannot work in the chosen delivery context.

**Approach:** Scan both strands for the nuclease's PAM, extract the protospacer upstream (Cas9) or downstream (Cas12a) of each PAM, and apply hard filters -- reject `TTTT` (Pol-III terminator) for U6/H1 expression, flag a missing 5' G for U6 (prepend a G rather than replace the first base), and note GC outside ~40-70% as a soft penalty. Ranking comes from the context-valid model (route to CRISPOR), not from a hand-rolled score.

```python
from Bio.Seq import Seq
import re

GC_MIN, GC_MAX = 0.40, 0.70   # outside this band on-target activity falls off (Doench 2014); soft penalty

def find_guides(sequence, pam='NGG', guide_length=20):
    '''Enumerate SpCas9 (NGG) spacers on both strands; spacer is 5' of the PAM.'''
    seq = sequence.upper()
    guides = []
    for m in re.finditer(r'(?=([ACGT]GG))', seq):
        pos = m.start()
        if pos >= guide_length:
            guides.append({'spacer': seq[pos - guide_length:pos], 'pam': seq[pos:pos + 3],
                           'cut': pos - 3, 'strand': '+'})   # SpCas9 cuts ~3 bp 5' of the PAM
    rc = str(Seq(seq).reverse_complement())
    n = len(seq)
    for m in re.finditer(r'(?=([ACGT]GG))', rc):
        pos = m.start()
        if pos >= guide_length:
            guides.append({'spacer': rc[pos - guide_length:pos], 'pam': rc[pos:pos + 3],
                           'cut': n - (pos - 3), 'strand': '-'})
    return guides

def passes_u6_filters(spacer):
    '''Hard filters for U6/H1 Pol-III expression (NOT applicable to in-vitro T7/RNP).'''
    gc = sum(c in 'GC' for c in spacer) / len(spacer)
    return 'TTTT' not in spacer and GC_MIN <= gc <= GC_MAX   # TTTT terminates Pol III
```

## Rank On-Target Activity in the Valid Context

**Goal:** Shortlist guides by predicted cutting using the model that matches the delivery context.

**Approach:** Do NOT hand-roll a scoring matrix. Route to CRISPOR, which selects the context-appropriate score (Rule Set 2/Azimuth for U6/lentiviral, CRISPRscan for T7/embryo) per the Haeussler 2016 logic and also nominates off-targets against the chosen genome. Treat the returned score as a shortlist signal, then carry 3-6 candidates forward.

```bash
# CRISPOR: aggregates the context-valid on-target score + off-target nomination per genome
crispor.py hg38 target.fa guides.tsv --maxOcc 60000
# columns include the on-target score (context-selected) and off-target counts/specificity
```

## Choose the Cut Site by Exon Biology (the under-used lever)

KO success is mostly won here, and pure efficiency ranking fails:
- Target an **early, constitutive coding exon** (present in all protein-coding isoforms) -- but **not the start-ATG region** (downstream reinitiation can rescue an N-terminal truncation).
- **Avoid the last exon and the last ~50 nt of the penultimate exon** -- PTCs there **escape NMD**, leaving a stable, possibly-functional truncated protein.
- Keep the cut **away from splice donor/acceptor sites** unless splice disruption is the goal -- indels there cause **exon skipping** that can restore frame (Mou 2017).
- Confirm the exon is constitutive **in the cell type of interest** (an exon spliced out of the dominant isoform is a silent failure), and screen for **SNPs under the protospacer/PAM** in the actual background (mismatch/PAM loss -> allele dropout).
- For ruthless KO / screens: **tile the conserved functional domain** (Shi 2015), not the gene start.

## Predict the Editing Outcome (frameshift fraction decides KO)

**Goal:** Prefer guides whose predicted indel spectrum is frameshift-rich (and, for a single-genotype line, dominated by one outcome).

**Approach:** Cas9 repair outcomes are predictable from the ~30 bp of local sequence flanking the cut. The cheap, no-ML signal is the **Bae 2014 microhomology out-of-frame score**: enumerate microhomology pairs flanking the cut, weight each predicted MMEJ deletion, and report the fraction whose length is not a multiple of 3. For a full genotype distribution use **inDelphi** (Shen 2018), **FORECasT** (Allen 2019), or **Lindel** (Chen 2019). Rank by **(editing efficiency) x (out-of-frame fraction)** -- a 70%-efficient guide with frameshift fraction 0.9 beats a 90%-efficient guide at 0.5. (See `examples/grna_design.py` for a runnable Bae-style out-of-frame implementation.)

## Per-Method Failure Modes

### "We used the top-ranked guide" with no validation
**Trigger:** sorting by on-target score and taking #1. **Mechanism:** scores are Spearman ~0.4 across context. **Symptom:** confident ranking, poor empirical hit rate. **Fix:** design 3-6 guides per gene and validate; treat the score as triage.

### Score used out of its training context
**Trigger:** CRISPRscan for a lentiviral screen, or Rule Set 2 for embryo RNP. **Mechanism:** each model is an assay artifact (Haeussler 2016). **Symptom:** "principled" but wrong ranking. **Fix:** pick the score from the delivery context before reading any number.

### Efficient cut, no knockout phenotype
**Trigger:** ranking by editing efficiency. **Mechanism:** ~1/3 in-frame indels + reinitiation/exon-skipping/NMD-escape/compensation. **Symptom:** high indel %, residual protein, milder-than-knockdown phenotype. **Fix:** rank by frameshift fraction (Bae/inDelphi), target early constitutive NMD-competent exons, verify at protein level.

### Last-exon / splice-site guide
**Trigger:** "early exon" applied naively. **Mechanism:** late PTC escapes NMD; splice-site indel skips the exon. **Symptom:** stable truncated/reframed protein. **Fix:** retarget an early constitutive exon away from junctions.

### Poly-T or missing 5' G in a U6 construct
**Trigger:** spacer with `TTTT` or non-G 5' end expressed from U6/H1. **Mechanism:** Pol-III termination / poor initiation. **Symptom:** little or no sgRNA. **Fix:** reject TTTT; prepend (do not replace) a 5' G. (Irrelevant for in-vitro T7/RNP.)

### Allele dropout in a non-reference background
**Trigger:** designing against GRCh38 for a patient/hybrid/cancer line. **Mechanism:** a SNP in the seed or PAM blocks one allele. **Symptom:** heterozygous "knockout" with a retained functional allele. **Fix:** design against the actual genotype.

## Quantitative Thresholds

| Parameter | Value | Source / rationale |
|-----------|-------|--------------------|
| On-target score use | rank/shortlist only; ~0.4 Spearman across context | Haeussler 2016 |
| GC content | ~40-70% (soft penalty) | Doench 2014 |
| Pol-III terminator | reject `TTTT` (U6/H1 only) | Pol-III termination |
| 5' G (U6) | prepend a G if absent | Pol-III initiation preference |
| SpCas9 cut | ~3 bp 5' of NGG (blunt) | Jinek 2012 |
| Bae out-of-frame score | prefer **>66** | Bae 2014 frameshift-reliability recommendation |
| KO ranking | efficiency x out-of-frame fraction | frameshift fraction, not cutting, drives KO |
| Guides per gene | **3-6**, validate empirically | scores are weak; redundancy buys back error |
| Exon target | early, constitutive, NMD-competent (not last exon / last ~50 nt of penult.) | PTC must trigger NMD across all isoforms |
| Residual protein after frameshift | expect ~1/3 retain protein | Smits 2019 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No guides found | no PAM in window / wrong PAM for nuclease | try Cas12a (TTTV) for AT-rich; widen window; SpCas9-NG/SpRY as last resort |
| Guide cuts but no KO phenotype | last exon / 3'UTR / in-frame indels / compensation | retarget early constitutive exon; rank by frameshift; verify protein |
| Score looks low for a clearly good guide | score used outside its training context | use the context-valid model |
| Heterozygous result in a non-reference line | SNP under guide/PAM | design against the actual genotype |

## References

- Jinek M, et al. (2012). A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity. *Science* 337(6096):816-821.
- Doench JG, et al. (2014). Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation. *Nat Biotechnol* 32(12):1262-1267.
- Doench JG, Fusi N, Sullender M, et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nat Biotechnol* 34(2):184-191.
- Moreno-Mateos MA, et al. (2015). CRISPRscan: designing highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo. *Nat Methods* 12(10):982-988.
- Haeussler M, et al. (2016). Evaluation of off-target and on-target scoring algorithms and integration into the guide RNA selection tool CRISPOR. *Genome Biol* 17:148.
- Kim HK, et al. (2019). SpCas9 activity prediction by DeepSpCas9. *Sci Adv* 5(11):eaax9249.
- Wang D, et al. (2019). Optimized CRISPR guide RNA design for two high-fidelity Cas9 variants by deep learning (DeepHF). *Nat Commun* 10:4284.
- Kim HK, et al. (2018). Deep learning improves prediction of CRISPR-Cpf1 guide RNA activity (DeepCpf1). *Nat Biotechnol* 36(3):239-241.
- Bae S, Kweon J, Kim HS, Kim JS (2014). Microhomology-based choice of Cas9 nuclease target sites. *Nat Methods* 11(7):705-706.
- Shen MW, et al. (2018). Predictable and precise template-free CRISPR editing of pathogenic variants (inDelphi). *Nature* 563(7733):646-651.
- Allen F, et al. (2019). Predicting the mutations generated by repair of Cas9-induced double-strand breaks (FORECasT). *Nat Biotechnol* 37(1):64-72.
- Chen W, et al. (2019). Massively parallel profiling and predictive modeling of the outcomes of CRISPR-Cas9 double-strand break repair (Lindel). *Nucleic Acids Res* 47(15):7989-8003.
- Shi J, et al. (2015). Discovery of cancer drug targets by CRISPR-Cas9 screening of protein domains. *Nat Biotechnol* 33(6):661-667.
- Smits AH, et al. (2019). Biological plasticity rescues target activity in CRISPR knock outs. *Nat Methods* 16(11):1087-1093.
- Mou H, et al. (2017). CRISPR/Cas9-mediated genome editing induces exon skipping by alternative splicing or exon deletion. *Genome Biol* 18(1):108.
- El-Brolosy MA, et al. (2019). Genetic compensation triggered by mutant mRNA degradation. *Nature* 568(7751):193-197.
- Ran FA, et al. (2015). In vivo genome editing using Staphylococcus aureus Cas9. *Nature* 520(7546):186-191.
- Nishimasu H, et al. (2018). Engineered CRISPR-Cas9 nuclease with expanded targeting space (SpCas9-NG). *Science* 361(6408):1259-1262.
- Hu JH, et al. (2018). Evolved Cas9 variants with broad PAM compatibility and high DNA specificity (xCas9). *Nature* 556(7699):57-63.
- Walton RT, et al. (2020). Unconstrained genome targeting with near-PAMless engineered CRISPR-Cas9 variants (SpRY). *Science* 368(6488):290-296.
- Kleinstiver BP, et al. (2019). Engineered CRISPR-Cas12a variants with increased activities and improved targeting ranges (enAsCas12a). *Nat Biotechnol* 37(3):276-282.
- Concordet JP, Haeussler M (2018). CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids Res* 46(W1):W242-W245.

## Related Skills

- off-target-prediction - Check genome-wide specificity after on-target design (a separate axis from activity)
- base-editing-design - DSB-free knockout via premature stop / splice disruption when indels are unwanted
- prime-editing-design - Scarless small edits without a double-strand break
- hdr-template-design - Design the donor when the goal is a precise knock-in, not a knockout
- crispr-screens/library-design - Pool guides into a screening library (domain tiling, Rule Set 2 logic)
- crispr-screens/crispresso-editing - Quantify indel/editing outcomes from amplicon sequencing
- primer-design/primer-basics - Design validation/genotyping primers around the cut
- primer-design/primer-specificity - Confirm genotyping primers are unique near paralogs/off-targets
- genome-intervals/gtf-gff-handling - Get exon coordinates to restrict guide placement
