---
name: bio-workflows-crispr-editing-pipeline
description: Orchestrates an end-to-end CRISPR editing experiment design from target gene to delivery-ready, validatable constructs. Sequences guide design, off-target assessment, edit-modality selection (knockout, base editing, prime editing, HDR knock-in), and template/donor design, with a QC checkpoint at each handoff. Use when designing a complete CRISPR experiment for knockout, point correction, or tagging and the order of operations, the modality decision, and the cross-cutting traps are needed rather than a single step. Defers each step's mechanics to the genome-engineering skills.
tool_type: mixed
primary_tool: CRISPOR
workflow: true
depends_on:
  - genome-engineering/grna-design
  - genome-engineering/off-target-prediction
  - genome-engineering/base-editing-design
  - genome-engineering/prime-editing-design
  - genome-engineering/hdr-template-design
qc_checkpoints:
  - after_grna_design: "Context-valid on-target shortlist (CRISPOR: Rule Set 2 for U6/lentiviral, CRISPRscan for T7/embryo); reject TTTT and GC extremes; rank by predicted frameshift/out-of-frame fraction, not raw activity; carry 3-6 guides in an early constitutive NMD-competent exon"
  - after_offtarget: "Escalate predicted -> detected -> validated; reject guides with a low-mismatch high-CFD off-target in a gene; variant-aware (gnomAD) for therapeutic guides; high-fidelity nuclease in the delivery format used"
  - after_template: "Blocking (PAM/seed) mutation present AND codon-checked; edit within ~10 bp of the cut; donor format matches cell type (ssODN/lssDNA/dsDNA/AAV; HITI for post-mitotic); report edit:indel purity for base editing"
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, pandas 2.2+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

This workflow coordinates the five genome-engineering skills; it does not re-implement their scoring. Real on-target ranking comes from CRISPOR (context-valid model), off-target nomination from Cas-OFFinder/CRISPRme, base-editor outcomes from BE-Hive, and prime-editing ranking from PRIDICT/DeepPrime -- the embedded code is illustrative orchestration only.

# CRISPR Editing Pipeline

**"Design a complete CRISPR editing experiment for my target"** -> Run guide design -> off-target assessment -> edit-modality selection -> template/donor design -> validation, applying a QC checkpoint at each handoff and routing every mechanic to the relevant genome-engineering skill.
- Python: orchestrate the stages; enumerate/filter candidate guides with `Bio.Seq`
- CLI/web: CRISPOR (on-target + off-target), Cas-OFFinder/CRISPRme (off-target), BE-Hive, PRIDICT

## The Single Most Important Modern Insight -- the pipeline is a chain of handoffs, each with a checkpoint, and the pivotal decision is the edit modality

A CRISPR experiment fails most often not at one step but at a handoff where an unstated assumption carries through: a guide picked by on-target score that turns out non-specific, an "efficient" guide that never knocks out the protein, a base edit reported by efficiency that is a genotype soup, an HDR donor with no blocking mutation whose edit is silently re-cut. The workflow's job is to make each handoff explicit and gated. The pivotal branch is **which edit modality**: a transition (C->T/A->G) is usually a base-editing job; any other small precise edit is prime editing; a knockout is a plain nuclease; a large or non-transition insertion is HDR (or PE+integrase). Choosing the modality first reframes every downstream step. The cross-cutting traps the checkpoints exist to catch: **on-target activity != specificity** (two separate axes), **efficient editing != knockout** (frameshift fraction and NMD-competent exon biology decide it), **base-editor efficiency != purity** (bystanders), **a donor without a blocking mutation self-destructs** (re-cutting reads out as failed HDR), and **predicted != detected != validated** for off-targets.

## Edit-Modality Decision Tree (the pivotal branch)

| Goal / edit | Modality | Route to |
|-------------|----------|----------|
| Gene knockout (any frameshift) | nuclease + NHEJ | grna-design (rank by frameshift fraction) |
| Knockout without a DSB / non-dividing / multiplex | base-editor premature stop or splice disruption | base-editing-design |
| C*G->T*A or A*T->G*C transition | base editing (CBE/ABE) | base-editing-design |
| C->G transversion | CGBE | base-editing-design |
| Other transversion, small indel, combined edit | prime editing | prime-editing-design |
| Small precise edit, no DSB tolerated | prime editing (PE) | prime-editing-design |
| Tag / reporter / allele replacement (cycling cells) | HDR knock-in | hdr-template-design |
| Large insertion / post-mitotic cells | HDR (AAV/HITI) or PE+integrase (PASTE/twinPE) | hdr-template-design / prime-editing-design |

## Workflow Overview

```
Target gene / position
        |
        v
[1. Guide design] ----> CRISPOR (context-valid on-target) + outcome model (Bae microhomology / inDelphi)
        |                CHECKPOINT: shortlist 3-6, frameshift-rich, early constitutive exon
        v
[2. Off-target assessment] ----> Cas-OFFinder (+bulges) / CRISPRme (variant-aware) + CFD
        |                CHECKPOINT: no low-mm high-CFD off-target in a gene; predicted->detected->validated
        v
    DECISION: which edit modality?
        |
    +----------+-------------+--------------+-------------+
    v          v             v              v             v
[3a. KO]   [3b. Base edit] [3c. Prime edit] [3d. HDR knock-in]
 frameshift  window+purity   pegRNA panel     donor + codon-checked block
        |          |              |              |
        v          v              v              v
[4. Validation] ----> amplicon deep-seq (CRISPResso2); report purity/indels; state LoD
```

## Stage 1 -- Guide Design (-> grna-design)

**Goal:** A shortlist of 3-6 specificity-checkable guides whose predicted repair outcome is frameshift-rich, in an early constitutive NMD-competent exon.

**Approach:** Establish the delivery context (it sets the valid on-target model and the hard filters), enumerate PAMs on both strands, drop TTTT/GC-extreme guides, rank on-target with the context-valid model via CRISPOR (not a hand-rolled score), and rank knockout candidates by predicted frameshift/out-of-frame fraction (Bae microhomology / inDelphi). **Checkpoint:** carry 3-6 guides; do not commit on raw activity alone.

## Stage 2 -- Off-Target Assessment (-> off-target-prediction)

**Goal:** Reject promiscuous guides and, for therapeutics, establish an evidence-laddered specificity profile.

**Approach:** Enumerate candidates with Cas-OFFinder including bulges and a relaxed PAM; rank by CFD; for a research knockout this in-silico pass is sufficient. For a therapeutic, run variant-aware nomination (CRISPRme vs gnomAD + individual), choose a high-fidelity nuclease in the delivery format used, and plan empirical discovery + amplicon validation. **Checkpoint:** on-target score does not predict specificity; treat predicted/detected/validated distinctly.

## Stage 3 -- Modality-Specific Design

**Goal:** Produce the construct(s) for the chosen modality.

**Approach:** Branch by the decision tree. Knockout -> the frameshift-ranked guide. Base editing -> position the target base at the window peak, minimize bystanders, choose the editor variant, report the genotype spectrum (-> base-editing-design). Prime editing -> a PBS x RTT panel with PAM-disrupting/MMR-evading silent edits and a 3' motif, ranked by PRIDICT/DeepPrime (-> prime-editing-design). HDR -> the donor format for the cell type with a mandatory codon-checked blocking mutation and the cut within ~10 bp of the edit (-> hdr-template-design). **Checkpoint:** blocking mutation present and codon-checked; base-editing purity reported.

## Stage 4 -- Validation

**Goal:** Quantify the intended edit and its byproducts.

**Approach:** Design genotyping/amplicon primers around the edit (-> primer-design/primer-basics; keep both 3' ends off the cut site and any expected indel, and confirm the amplicon is unique near paralogs/pseudogenes -> primer-design/primer-specificity) and quantify outcomes by amplicon deep sequencing (CRISPResso2 / BE-Analyzer) -- intended-edit rate, indels, and (for base/prime editing) product purity -- stating the limit of detection (-> crispr-screens/crispresso-editing). The critical hand-off across the wet-lab gap: give CRISPResso2 the UNEDITED amplicon of the specific system as `--amplicon_seq` (the actual wild-type/pre-edit sequence -- matching the cell line's SNPs and primer product, NOT a mismatched canonical genome), the actual protospacer as `--guide_seq` so the quantification window centers on the cut, and for HDR/KI the intended edit as `--expected_hdr_amplicon_seq`. Reads are scored "unmodified" by matching `--amplicon_seq`, so supplying the EDITED sequence there makes real edits score as unmodified (~0%) with no error raised. **Checkpoint:** report purity and LoD, not a lone efficiency number.

## Common Errors (integration level)

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Top guide has a near-perfect off-target | picked by on-target score alone | re-rank by specificity; on-target and specificity are separate axes |
| Efficient editing, no knockout phenotype | in-frame indels / late-exon / compensation | rank by frameshift fraction; target an early constitutive exon; verify protein |
| Base edit "80% efficient" but messy genotypes | bystanders in the window | report the spectrum; reposition or use a narrowed-window editor |
| HDR gives only indels | donor lacks a blocking mutation | add a codon-checked PAM/seed block; the edit was re-cut |
| Validation shows ~0% editing on a working edit | EDITED (or wrong) sequence supplied as `--amplicon_seq`, so edited reads match the reference / wrong guide window | give CRISPResso2 the UNEDITED reference as `--amplicon_seq` (+ `--expected_hdr_amplicon_seq` for HDR) and the actual `--guide_seq` |
| "No off-targets" claimed | LoD not stated / reference-only | state the LoD; variant-aware for therapeutics |

## References

- Doench JG, Fusi N, Sullender M, et al. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nat Biotechnol* 34(2):184-191.
- Concordet JP, Haeussler M (2018). CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. *Nucleic Acids Res* 46(W1):W242-W245.
- Bae S, Park J, Kim JS (2014). Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. *Bioinformatics* 30(10):1473-1475. [off-target search]
- Bae S, Kweon J, Kim HS, Kim JS (2014). Microhomology-based choice of Cas9 nuclease target sites. *Nat Methods* 11(7):705-706. [microhomology/MMEJ frameshift-outcome predictor -- the exp(-deletionLength/20) length weight]
- Clement K, Rees H, Canver MC, et al. (2019). CRISPResso2 provides accurate and rapid genome editing sequence analysis. *Nat Biotechnol* 37(3):224-226.
- Paquet D, Kwart D, Chen A, et al. (2016). Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9. *Nature* 533(7601):125-129.

## Related Skills

- genome-engineering/grna-design - Guide design and outcome-aware knockout ranking
- genome-engineering/off-target-prediction - Specificity assessment and the evidence ladder
- genome-engineering/base-editing-design - CBE/ABE window, bystander purity, off-target classes
- genome-engineering/prime-editing-design - pegRNA panel design and PE system selection
- genome-engineering/hdr-template-design - Donor format and codon-checked blocking mutation
- crispr-screens/crispresso-editing - Quantify and validate editing outcomes from amplicon reads
- crispr-screens/library-design - Scale single-gene design to a pooled screen
- primer-design/primer-basics - Design the genotyping/amplicon primers around the edit
- primer-design/primer-specificity - Confirm the genotyping amplicon is unique near paralogs/pseudogenes
