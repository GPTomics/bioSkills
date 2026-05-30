---
name: bio-genome-engineering-base-editing-design
description: Designs cytosine (CBE, C-to-T) and adenine (ABE, A-to-G) base-editor guides by positioning the target base at the activity-peak of the editing window (protospacer positions ~5-7, PAM-distal numbering), minimizing bystander edits for product purity, reading dinucleotide context (APOBEC1 TC favored / GC disfavored), and selecting the editor variant (BE4max, ABEmax, ABE8e, YE1/SECURE, TadCBE, CGBE, SpG/SpRY-BE). Covers knockout by premature stop (CRISPR-STOP/iSTOP) and splice-site disruption, the three off-target classes (Cas-dependent, Cas-independent DNA, RNA), outcome prediction (BE-Hive/DeepBE), and the base-vs-prime-vs-HDR decision. Use when installing a transition mutation without a double-strand break, knocking out a gene without indels, or choosing CBE vs ABE. Generic guide scoring, prime editing, and HDR donors are separate skills.
tool_type: python
primary_tool: BioPython
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

BE-Hive is NOT a pip package -- it is a web tool (crisprbehive.design) plus clone-only repos (`maxwshen/be_predict_efficiency`, `be_predict_bystander`) pinning old dependencies (scikit-learn 0.20.3, BioPython 1.73). Use the web tool for a quick answer, the clone for batch use. Editing windows and outcome predictions are **editor-variant- and cell-type- specific** and are activity gradients, not hard boxes -- record which editor model was used.

# Base Editing Design

**"Install a transition mutation without a double-strand break"** -> Write the edit as a base-pair change to pick the editor family, find a PAM that lands the target base at the window peak while keeping bystanders out, read sequence context, choose the editor variant, and report the predicted genotype spectrum -- not a lone efficiency number.
- Python: window/bystander scan with `Bio.Seq` + `re`; per-base context reading
- Web/clone: BE-Hive / DeepBE for the per-base efficiency + bystander genotype spectrum
- Web: BE-Designer (guide enumeration), BE-Analyzer / CRISPResso2 (validate from NGS)

## The Single Most Important Modern Insight -- the job is product PURITY, not editing efficiency

A base editor produces a *distribution of genotypes at one site*, not a binary cut. "80% editing" can mean 80% of alleles carry the clean intended edit, or 80% carry *some* edit -- a soup where a bystander C/A two positions away is also converted half the time, so the exact desired genotype is only 30% of alleles. Both report as "80% efficient." The number that answers the biology is **precise editing: the fraction of alleles with the target edit AND no bystander edit.** Design is dominated by three coupled levers on one ~5-nt window: **(1) where the target base sits** (set by PAM/guide; drives efficiency), **(2) what else sits in the window** (bystander C's/A's; drives purity), **(3) the local dinucleotide context** of each base (drives which actually edit). The corollary trap: a *more active* editor (ABE8e over ABE7.10) raises headline efficiency while *lowering* purity (higher processivity sweeps more bystanders). Rank guides by the predicted **outcome spectrum**, never by efficiency alone.

The second field-defining insight: **base editors have THREE off-target classes, and the two that distinguish CBE from ABE are invisible to every Cas-off-target tool** (below). Most people only think of class 1 (guide-directed); the classes that matter for safety are Cas-independent and guide-invisible.

## Mechanism & the Position-Numbering Convention (get this right)

Both families are a **ssDNA deaminase fused to a Cas9 NICKASE (nCas9, D10A)**, guided by an ordinary sgRNA -- **no double-strand break**. When Cas9 binds, the non-target (PAM-containing) strand is displaced as ssDNA in the R-loop; the deaminase edits bases on that displaced strand within the window. **Canonical numbering (Komor 2016): position 1 = PAM-DISTAL (5' end of the spacer), position 20 = PAM-PROXIMAL (next to the PAM at 21-23). The editing window is ~positions 4-8, peak ~5-7.** (This is the opposite of what older tutorials sometimes say.)
- **CBE:** cytidine deaminase (rAPOBEC1) C->U, + **UGI** (blocks uracil excision -- the main CBE purity determinant), nCas9 nicks the unedited strand -> resolves to C->T (G->A other strand).
- **ABE:** lab-evolved TadA* deaminase A->inosine (read as G) -> resolves to A->G (T->C other strand). No UGI needed (inosine is not efficiently excised) -- a key reason early ABE was intrinsically cleaner than CBE.

## The Window Is a Gradient, Editor-Specific -- not a box

Activity inside the "box" is a steep gradient: a target at position 6 edits far better than one at 8; a bystander at 8 edits at a fraction of one at 5. The design move is "how close to the 5-6-7 peak is my target, and how far toward the cold edges can I push the bystanders?" The width is a property of the **editor variant**: ABE8e's high processivity widens it to ~3-11, so a guide that was bystander-clean with ABE7.10 becomes dirty when "upgraded" to ABE8e without re-checking. Carrying a window assumption across an editor switch is the most common silent failure.

## Editor Variant Selection

| Editor | Class | When |
|--------|-------|------|
| **BE4max / AncBE4max** | CBE | default modern CBE (Koblan 2018) |
| **ABEmax** | ABE | strong default ABE (Koblan 2018) |
| **ABE8e** | ABE | maximum activity (hard targets, screens) -- but WIDER window, MORE bystanders/off-target, LOSES context preference (Richter 2020; Lapinaite 2020) |
| YE1 / SECURE-BE3 | CBE | narrowed window / low RNA off-target (Grunewald 2019) |
| **TadCBE / TadDE** | CBE / dual | lowest Cas-independent DNA+RNA off-target CBE (Neugebauer 2023) |
| CGBE1 | C->G | the only transversion a base editor does well (Kurt 2021) |
| A&C-BEmax / SPACE | dual | simultaneous C->T and A->G (niche; adds bystander surface) (Grunewald 2020; Sakata 2020) |
| SpG / SpRY base editors | any | when no canonical NGG positions the base (Walton 2020) -- more class-1 off-target |

Default BE4max/ABEmax; reach for ABE8e for activity at a purity cost; a SECURE/TadCBE variant when off-target matters; a PAM-flexible editor when no NGG positions the base. The *deaminase* choice (rAPOBEC1 vs APOBEC3A vs eA3A vs evolved TadA) sets window width, context, and off-target far more than the BE3-vs-BE4 generation number.

## The Three Off-Target Classes -- ask "which of the three?" first

| Class | Mechanism | Detected by | Mitigation |
|-------|-----------|-------------|------------|
| **1. Cas-dependent DNA** | guide mismatch-tolerance, like ordinary Cas9 | GUIDE-/CIRCLE-seq; in-silico (-> off-target-prediction) | high-fidelity Cas; guide selection |
| **2. Cas-INDEPENDENT DNA** | deaminase edits transiently-exposed genomic ssDNA, **no guide** | GOTI (Zuo 2019); WGS (Jin 2019) | choose a low-activity-deaminase variant (YE1, **TadCBE**) |
| **3. RNA** | deaminase edits cellular mRNA, **no guide**, transient | RNA-seq (Grunewald/Rees/Zhou 2019) | **SECURE** variants (rAPOBEC1 R33A; ABE F148A); TadCBE |

Classes 2-3 have no protospacer, so they are **invisible to every guide-based predictor** -- the guide cannot be designed to avoid them; only a cleaner editor can. **CBE has historically been the dirtier family on classes 2-3** (Zuo/Jin: CBE >20x background Cas-independent SNVs, mostly C>T in transcribed DNA; early ABE did not) -- a genuine input to the CBE-vs-ABE choice. But this is *variant*-specific, not family-destiny: ABE8e raised it back up; TadCBE pulled CBE's down. Class-3 RNA edits are transient (RNA turns over) -- a real dose/exposure-dependent liability for chronic/AAV/therapeutic use, often tolerable for a transient cell-line transfection.

## Knockout by Base Editing (no DSB)

Base editing knocks out a gene without a double-strand break -- no indel lottery, no large deletions/translocations, no p53 response, works in **non-dividing cells**, and safe for **multiplex** (no translocations between simultaneous cut sites). Two routes:
- **Premature stop (CRISPR-STOP/iSTOP):** a CBE converts **CAA->TAA, CAG->TAG, CGA->TGA** (sense) or **TGG->stop** via the antisense strand. Only these four codons are reachable; target an **early** stop (before functional domains, NMD-competent). iSTOP precomputed sgRNAs cover 97-99% of genes (Billon 2017; Kuscu 2017).
- **Splice-site disruption:** edit the invariant splice-donor **GT** or acceptor **AG** -> mis- splicing / exon skipping (Kluesner 2021). Often the more robust KO (every gene has many junctions; only ~half have a well-placed early stop codon).

The quiet failure mode: a base-editor KO is only as complete as the editing -- an incompletely edited cell still makes wild-type protein, and a bystander can turn an intended silent KO into a missense allele. Prefer an early pmSTOP; fall back to splice disruption; verify protein.

## Decision Tree by Scenario

| Desired edit | Use | Why |
|--------------|-----|-----|
| Write edit as base-pair change first | -- | C:G->T:A = CBE; A:T->G:C = ABE; resolves strand confusion (a "G->A" edit is a CBE job) |
| C:G->T:A or A:T->G:C transition, base positionable | **CBE / ABE** | higher efficiency, cleaner, no DSB -- preferred over PE/HDR for transitions |
| C->G transversion | CGBE | the only transversion a base editor does |
| Any other transversion, small indel, multi-base | -> prime-editing-design | beyond base-editor chemistry |
| Large insertion / knock-in | -> hdr-template-design (or PE+integrase) | beyond base editing |
| Knockout, no DSB / non-dividing / multiplex | **BE pmSTOP or splice disruption** | no translocations, works in post-mitotic cells |
| Both CBE and ABE could make it (opposite strands) | break the tie on **purity + off-target** | ABE historically cleaner on classes 2-3 (per variant) |
| No NGG positions the base | SpG/SpRY base editor | accept more class-1 off-target |
| Off-target raised | ask "which of the three classes?" | class 1 -> off-target-prediction; 2-3 -> editor variant |
| Validate outcomes | -> crispr-screens/base-editing-analysis | amplicon NGS spectrum |

## Find Editable Guides and Read the Window

**Goal:** Find guides that place the target base near the window peak with the fewest in-window bystanders, and read the genotype spectrum the window implies.

**Approach:** Scan both strands for PAMs, compute where the target base lands in the spacer (PAM-distal = position 1), keep guides where it falls in the window, list bystander C's/A's, and read the 5' dinucleotide context of each (TC favored / GC disfavored for APOBEC1). The position-efficiency values below are **coarse illustrative gradients, not measurements** -- for a real outcome spectrum use BE-Hive/DeepBE, then validate by NGS. (See `examples/base_editing_design.py`.)

```python
from Bio.Seq import Seq
import re

CBE_WINDOW = (4, 8)   # activity gradient, peak ~5-7; PAM-distal numbering (position 1 = 5' end of spacer)
ABE_WINDOW = (4, 8)   # ABE7.10 tighter (~4-7); ABE8e WIDER (~3-11) -- editor-specific
```

## Per-Method Failure Modes

### Ranked guides by efficiency, shipped a genotype soup
**Trigger:** sorting by on-target editing %. **Mechanism:** efficiency hides the bystander spectrum. **Symptom:** "80% edited" but few alleles have the exact desired genotype. **Fix:** rank by predicted precise (bystander-free) editing; report the spectrum (BE-Hive).

### Carried a window across an editor switch
**Trigger:** "upgraded" to ABE8e, kept the ABE7.10 window. **Mechanism:** ABE8e's window is wider (~3-11). **Symptom:** new bystanders, dirtier product. **Fix:** use the variant-specific window; re-check bystanders after any editor change.

### Treated a base-editor off-target like a Cas9 off-target
**Trigger:** running GUIDE-seq/in-silico predictors and declaring it safe. **Mechanism:** those cover only class 1; classes 2-3 are guide-invisible. **Symptom:** clean class-1 report, genome/transcriptome-wide deaminase collateral. **Fix:** ask "which of the three?"; for 2-3 pick a SECURE/TadCBE variant.

### Picked the wrong editor (strand confusion)
**Trigger:** wanting a "G->A" change and reaching for ABE. **Mechanism:** G->A is C->T on the complement = a CBE job. **Symptom:** no editor can make it. **Fix:** write the edit as a base-pair change first.

### "I designed a stop, the gene is off"
**Trigger:** assuming pmSTOP = knockout. **Mechanism:** incomplete editing leaves WT protein; a late/NMD-escaping stop leaves functional product; a bystander makes a missense allele. **Fix:** target an early stop or splice site; verify at the protein level.

### Reached for SpRY by default
**Trigger:** using a PAM-flexible editor for convenience. **Mechanism:** relaxed PAM tolerates more mismatches -> more class-1 off-target. **Fix:** exhaust NGG first; use SpG/SpRY as a deliberate trade and check off-target harder.

## Quantitative Thresholds

| Parameter | Value | Source |
|-----------|-------|--------|
| CBE window | ~positions 4-8, peak 5-7 (PAM-distal numbering) | Komor 2016 |
| ABE7.10 window | ~4-7/4-8 | Gaudelli 2017 |
| ABE8e window | ~3-11 (wider, more bystanders) | Richter 2020; Lapinaite 2020 |
| APOBEC1 context | TC strongly preferred, GC strongly disfavored | Komor 2016 |
| iSTOP codons | CAA/CAG/CGA (sense) + TGG (antisense) | Billon 2017 |
| iSTOP gene coverage | 97-99% of genes (with some editor) | Billon 2017 |
| CBE Cas-independent DNA off-target | >20x background (early CBE; not early ABE) | Zuo 2019; Jin 2019 |
| Purity metric | report % target-edit-AND-no-bystander, or the full spectrum | BE-Hive (Arbab 2020) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No editable guide found | no PAM lands the base at ~5-7 | use SpG/SpRY base editor; check the other strand |
| Many in-window bystanders | multi-C/multi-A window | reposition to the peak; narrowed-window variant (YE1); exploit GC-context disfavor |
| "Fewest off-target edits" conflates bystanders with off-targets | mislabeling | bystanders are in-window on-locus; off-targets are elsewhere -- different fixes |
| Clean predicted spectrum but genomic collateral | model predicts on-target window only | classes 2-3 need editor choice + orthogonal assays |

## References

- Komor AC, Kim YB, Packer MS, Zuris JA, Liu DR (2016). Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. *Nature* 533(7603):420-424.
- Komor AC, Zhao KT, Packer MS, et al. (2017). Improved base excision repair inhibition and bacteriophage Mu Gam protein yields C:G-to-T:A base editors with higher efficiency and product purity (BE4/BE4-Gam). *Sci Adv* 3(8):eaao4774.
- Gaudelli NM, Komor AC, Rees HA, et al. (2017). Programmable base editing of A:T to G:C in genomic DNA without DNA cleavage (ABE7.10). *Nature* 551(7681):464-471.
- Koblan LW, Doman JL, Wilson C, et al. (2018). Improving cytidine and adenine base editors by expression optimization and ancestral reconstruction (BE4max/AncBE4max/ABEmax). *Nat Biotechnol* 36(9):843-846.
- Richter MF, Zhao KT, Eton E, et al. (2020). Phage-assisted evolution of an adenine base editor with improved Cas domain compatibility and activity (ABE8e). *Nat Biotechnol* 38(7):883-891.
- Lapinaite A, Knott GJ, Palumbo CM, et al. (2020). DNA capture by a CRISPR-Cas9-guided adenine base editor. *Science* 369(6503):566-571.
- Zuo E, Sun Y, Wei W, et al. (2019). Cytosine base editor generates substantial off-target single-nucleotide variants in mouse embryos (GOTI). *Science* 364(6437):289-292.
- Jin S, Zong Y, Gao Q, et al. (2019). Cytosine, but not adenine, base editors induce genome-wide off-target mutations in rice. *Science* 364(6437):292-295.
- Grunewald J, Zhou R, Garcia SP, et al. (2019). Transcriptome-wide off-target RNA editing induced by CRISPR-guided DNA base editors. *Nature* 569(7756):433-437.
- Zhou C, Sun Y, Yan R, et al. (2019). Off-target RNA mutation induced by DNA base editing and its elimination by mutagenesis. *Nature* 571(7764):275-278.
- Rees HA, Wilson C, Doman JL, Liu DR (2019). Analysis and minimization of cellular RNA editing by DNA adenine base editors (SECURE-ABE). *Sci Adv* 5(5):eaax5717.
- Grunewald J, Zhou R, Iyer S, et al. (2019). CRISPR DNA base editors with reduced RNA off-target and self-editing activities (SECURE-BE3). *Nat Biotechnol* 37(9):1041-1048.
- Billon P, Bryant EE, Joseph SA, et al. (2017). CRISPR-Mediated Base Editing Enables Efficient Disruption of Eukaryotic Genes through Induction of STOP Codons (iSTOP). *Mol Cell* 67(6):1068-1079.
- Kuscu C, Parlak M, Tufan T, et al. (2017). CRISPR-STOP: gene silencing through base-editing-induced nonsense mutations. *Nat Methods* 14(7):710-712.
- Kluesner MG, Lahr WS, Lonetree CL, et al. (2021). CRISPR-Cas9 cytidine and adenosine base editing of splice-sites mediates highly-efficient disruption of proteins in primary and immortalized cells. *Nat Commun* 12:2437.
- Arbab M, Shen MW, Mok BY, et al. (2020). Determinants of Base Editing Outcomes from Target Library Analysis and Machine Learning (BE-Hive). *Cell* 182(2):463-480.
- Song M, Kim HK, Lee S, et al. (2020). Sequence-specific prediction of the efficiencies of adenine and cytosine base editors (DeepBE). *Nat Biotechnol* 38(9):1037-1043.
- Hwang GH, Park J, Lim K, et al. (2018). Web-based design and analysis tools for CRISPR base editing (BE-Designer/BE-Analyzer). *BMC Bioinformatics* 19:542.
- Kurt IC, Zhou R, Iyer S, et al. (2021). CRISPR C-to-G base editors for inducing targeted DNA transversions in human cells (CGBE1). *Nat Biotechnol* 39(1):41-46.
- Neugebauer ME, Hsu A, Arbab M, et al. (2023). Evolution of an adenine base editor into a small, efficient cytosine base editor with low off-target activity (TadCBE/TadDE). *Nat Biotechnol* 41(5):673-685.
- Walton RT, Christie KA, Whittaker MN, Kleinstiver BP (2020). Unconstrained genome targeting with near-PAMless engineered CRISPR-Cas9 variants (SpG/SpRY). *Science* 368(6488):290-296.
- Grunewald J, Zhou R, Lareau CA, et al. (2020). A dual-deaminase CRISPR base editor enables concurrent adenine and cytosine editing (A&C-BEmax). *Nat Biotechnol* 38(7):861-864.
- Sakata RC, Ishiguro S, Mori H, et al. (2020). Base editors for simultaneous introduction of C-to-T and A-to-G mutations. *Nat Biotechnol* 38(7):865-869.

## Related Skills

- grna-design - Cas-domain on-target/class-1 guide scoring this skill builds the window/bystander layer on
- prime-editing-design - For transversions (except C->G), indels, and multi-base replacements
- off-target-prediction - Owns the class-1 Cas-dependent off-target assays/prediction
- crispr-screens/base-editing-analysis - Quantify editing outcomes from NGS after the experiment
- crispr-screens/crispresso-editing - Amplicon editing-spectrum quantification
- variant-calling/variant-annotation - Annotate the consequence of the installed edit
