---
name: bio-genome-engineering-hdr-template-design
description: Designs donor/repair templates for precise CRISPR knock-ins -- choosing the format (ssODN, long-ssDNA/Easi-CRISPR, dsDNA/plasmid, AAV6), sizing homology arms, placing the cut within ~10 bp of the edit, and adding a mandatory codon-checked blocking (PAM/seed) mutation so the edited allele is not re-cut. Frames the HDR-vs-NHEJ-vs-MMEJ pathway competition, the MMEJ (PITCh) and homology-independent (HITI/HMEJ) alternatives for post-mitotic cells, ssODN strand/asymmetry choice, phosphorothioate end-protection, and ranked HDR enhancers. Use when designing a donor for a point mutation, epitope/fluorophore tag, allele replacement, or knock-in, or when HDR efficiency is low. Guide design and base/prime editing are separate skills.
tool_type: python
primary_tool: primer3-py
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, primer3-py 2.0+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

primer3-py designs PRIMERS (use `primer3.bindings.design_primers(seq_args, global_args)`; the camelCase `designPrimers` is deprecated since 1.0.0), not homology arms -- arm extraction and codon-aware blocking are the skill's own BioPython code. Design arms and guide against the **actual cell line's sequence**, not GRCh38 -- a SNP in an arm reduces annealing and a SNP in the PAM/seed can mean the guide does not cut.

# HDR Template Design

**"Design a donor for my CRISPR knock-in"** -> Decide the format by edit size and whether the cell cycles, size the arms, confirm a guide cuts within ~10 bp of the edit, and add a codon-checked blocking mutation so the corrected allele cannot be re-cut.
- Python: arm extraction + codon-aware PAM/seed blocking with `Bio.Seq`; `primer3.bindings.design_primers()` for arm-amplification and junction-validation primers
- Decision: format/route by cell type (cycling vs post-mitotic) and insert size

## The Single Most Important Modern Insight -- HDR is the minority pathway, and a donor without a blocking mutation is a self-destructing one

A Cas9 double-strand break is repaired by whichever pathway wins a kinetic race, and in most cells the winner is **classical NHEJ** (fast, all cell-cycle phases). MMEJ (microhomology, S/G2) and **HDR (template-dependent, S/G2 only)** are minority players, so **unenhanced HDR knock-in is typically single-digit to low-double-digit percent -- that is normal, not a failure.** The donor is not "the sequence to insert"; it is the toolkit for tilting a race NHEJ is structurally favored to win. The corollary, and the field's most expensive misread: **a donor with perfect arms but no blocking mutation gets its successful edit erased** -- the corrected allele still has an intact protospacer + PAM, so Cas9 re-cuts it and NHEJ scars it, and the indel reads out as "HDR failed" (indistinguishable from low HDR). So when someone reports "low HDR, lots of indels," the first question is not "how long are the arms?" -- it is **"does the donor disrupt the PAM or seed?"** A blocking mutation is mandatory and must be codon-checked; without it the readout is re-cutting, not HDR.

## The Pathway Competition (everything follows from this)

| Pathway | Cell cycle | Template | Signature | Relevance |
|---------|-----------|----------|-----------|-----------|
| **c-NHEJ** | all phases (dominant) | none | indels | the competitor; the engine HITI exploits |
| **MMEJ / alt-EJ** (Pol theta) | S/G2 | 5-25 bp microhomology | microhomology-flanked deletions | the PITCh route |
| **HDR / HR** | **S/G2 only** | sister chromatid or **exogenous donor** | precise, scarless | the classic knock-in route; minority |
| SSA | S/G2 | repeats | deletion between repeats | nuisance |

End resection (cell-cycle-gated, licensed in S/G2; 53BP1-RIF1 protects ends/pro-NHEJ, BRCA1 antagonizes it/pro-HDR) decides the fork. Consequences: **post-mitotic cells barely do HDR** -> for neurons/muscle/in-vivo tissue, HITI/HMEJ (NHEJ-based) is the *correct first choice*, not a fallback. Timing RNP+donor delivery into S/G2 raises HDR (Lin 2014, up to ~38% in HEK293T) -- the donor is necessary but the cell-cycle state gates it.

## Donor Format Decision

| Format | Insert | Arms | Best for | Caveat |
|--------|--------|------|----------|--------|
| **ssODN** | <= ~50 bp edits | ~30-60 nt each (total ~120-200 nt) | point mutations, small tags, loxP | synthesis ceiling ~200 nt; strand choice contested |
| **long ssDNA (Easi-CRISPR)** | ~0.2-2 kb | ~50-100 nt | cassettes, floxed/conditional alleles, zygote KI | **less toxic & less random integration than dsDNA**; harder to make |
| dsDNA (PCR/linear) | ~0.1-1.5 kb | ~200-800 bp | medium cassettes | **dsDNA is toxic** (innate sensing) + random integration |
| plasmid / HMEJ | up to several kb | ~500-2000 bp | large insertions, conditional alleles | backbone integration risk; slowest |
| **AAV6** | <= ~4.5 kb (ITR-to-ITR, arms included) | ~400 bp-1 kb | hard-to-transfect primary cells (HSPC, T, iPSC), in vivo | manufacturing cost; cargo cap is a hard wall |

Heuristics: point/small edit -> ssODN; 0.2-2 kb -> **lssDNA over dsDNA** (cleaner) for animal/zygote work; large cassette -> plasmid/HMEJ (lines) or **AAV6+RNP** (primary cells); **post-mitotic -> HITI**.

## Homology Arms, Edit-to-Cut Distance, and the Blocking Mutation (the most-botched trio)

- **Arm length:** ssODN ~30-60 nt each (more does not help, costs synthesis); dsDNA/plasmid ~500-800 bp sweet spot. Match the **actual cell-line sequence**, not the reference.
- **Edit-to-cut distance drives guide choice:** HDR incorporation falls sharply with distance, so place the cut **within ~10 bp of the edit** (Paquet 2016). The guide and donor are a *joint* design -- a "great" guide cutting 25 bp away is worse than a mediocre one cutting 3 bp away. **If no guide cuts within ~10 bp, HDR is the wrong tool** -> reconsider base/prime editing.
- **Blocking mutation (mandatory, codon-checked):** disrupt the **PAM** synonymously (preferred -- change a G in the NGG at a wobble position); if the PAM has no synonymous option, introduce **silent seed-region** mutations (PAM-proximal ~10-12 nt; PAM-distal mismatches are tolerated and do not block). The blocking edit must sit within the ~10 bp incorporation window (which is also where it blocks best). Paquet 2016 (CORRECT) raises per-allele accuracy ~10-fold and allows zygosity control by distance.

## ssODN Strand & Asymmetry (an over-cited rule) + the reliable win

Richardson 2016 proposed an ssODN **complementary to the non-target strand**, **asymmetric with the longer arm PAM-proximal (~91 nt) and the shorter PAM-distal (~36 nt)**. Subsequent systematic work could **not** reproduce this as universal: the optimal strand flips by locus and the asymmetric advantage often vanishes once both arms are >=30 nt. Treat it as a **prior to test, not a law** -- generate both strands and symmetric+asymmetric variants and test them. By contrast, **phosphorothioate (PS) end-protection (2-3 terminal bases each end)** is a near-universal cheap win (exonuclease resistance) -- encode these at opposite confidence levels.

## MMEJ / Homology-Independent Routes (when HDR is the wrong tool)

- **PITCh / CRIS-PITCh (MMEJ, Nakade 2014):** ~5-25 bp microhomologies instead of long arms; Pol theta joins donor to genome. Appeal is purely donor-construction convenience (microhomologies are primer overhangs); cost is error-prone junctions.
- **HITI (Suzuki 2016):** **homology-INDEPENDENT, NHEJ-based -> works in non-dividing cells.** The donor carries the same Cas9 target site(s) in **reverse orientation** flanking the insert; wrong-orientation insertions reconstitute the site and get re-cut/ejected, right-orientation insertions destroy it and lock in. Junctions can carry small indels.
- **HMEJ (Yao 2017):** ~800 bp arms PLUS flanking gRNA sites that linearize the donor in vivo; higher KI than HR/NHEJ/MMEJ in some contexts but ties/loses in others (mESC, N2a) -- test at the target locus.

| Situation | Route |
|-----------|-------|
| Point/small edit, cycling cells | ssODN + HDR (with blocking mutation) |
| Medium/large cassette, cycling line | HDR (lssDNA/plasmid) or HMEJ |
| Large cassette, primary cells (HSPC/T/iPSC) | AAV6 donor + RNP + HDR |
| Clean zygote/animal KI, <=2 kb | lssDNA Easi-CRISPR + HDR |
| Trivial donor construction wanted | PITCh (MMEJ) |
| **Non-dividing / post-mitotic / in vivo** | **HITI** (or HMEJ) |
| Edit far from any cut / single base | -> base-editing-design or prime-editing-design (donor-free) |

## HDR Enhancers -- ranked experiments, not multipliers

Most enhancers are marginal, cell-type-specific, and frequently non-reproducible; the published fold-changes are line-specific maxima. A blocking mutation and a cut near the edit matter more than any small molecule.
- **First tier (try by default, low risk):** cell-cycle timing of RNP delivery (Lin 2014); cold shock (32 C, 24-48 h; Guo 2018); PS end-protection; RNP+ssODN co-delivery.
- **Second tier (test in the target cells, expect variability):** DNA-PKcs inhibition (M3814/nedisertib -- the most consistently potent small molecule); 53BP1 inhibition (i53 / Alt-R HDR Enhancer).
- **Bottom tier (mention with a reproducibility warning):** SCR7 (widely un-reproducible), RS-1.

## Generate the Donor, Block Re-cutting (codon-checked), and Design Validation Primers

**Goal:** Assemble a donor that incorporates the edit AND survives re-cutting, with primers to amplify the arms and genotype the junction.

**Approach:** Extract arms flanking the cut, insert the edit, then add a blocking mutation -- disrupt the PAM synonymously if a wobble option exists, else introduce silent seed mutations -- verifying the change does not alter the encoded amino acid. Use primer3-py for arm-amplification/junction primers. (See `examples/hdr_template_design.py` for codon-aware blocking and a primer3 call.)

```python
from Bio.Seq import Seq

def synonymous_pam_block(codon_table, pam_codon, alt_codon):
    '''Return True only if a PAM-disrupting codon swap keeps the same amino acid (silent).'''
    return codon_table.get(pam_codon) == codon_table.get(alt_codon)   # never mutate the PAM without this check
```

## Per-Method Failure Modes

### "I got an indel, so HDR failed"
**Trigger:** low edit, mostly indels, no blocking mutation. **Mechanism:** the corrected allele keeps an intact PAM -> Cas9 re-cuts -> NHEJ scar. **Symptom:** indels indistinguishable from no-HDR. **Fix:** add a codon-checked PAM/seed blocking mutation; the readout was re-cutting, not HDR.

### Edit far from the cut
**Trigger:** best-cutting guide is 25 bp from the edit. **Mechanism:** HDR incorporation falls with distance. **Symptom:** only the blocking mutation is incorporated (useless silent-only allele) or no edit. **Fix:** choose a guide cutting within ~10 bp; if none, switch to base/prime editing.

### Frame-unaware blocking mutation
**Trigger:** blindly changing the NGG's second G to A. **Mechanism:** the PAM may be in a coding frame. **Symptom:** an unintended missense/nonsense change. **Fix:** verify the swap is synonymous; else use silent seed mutations.

### dsDNA in sensitive cells
**Trigger:** a plasmid/PCR donor in iPSC/primary/zygotes. **Mechanism:** dsDNA toxicity + random integration. **Symptom:** low viability, random integrants. **Fix:** use lssDNA (Easi-CRISPR) or AAV6.

### HDR donor in post-mitotic cells
**Trigger:** ssODN/plasmid for neurons/in-vivo tissue. **Mechanism:** HDR runs only in S/G2. **Symptom:** essentially no knock-in. **Fix:** use HITI (NHEJ-based) or HMEJ.

### Arms designed against the reference
**Trigger:** GRCh38 arms for a passaged/cancer line. **Mechanism:** line-specific SNPs in the arm or PAM/seed. **Symptom:** poor annealing or no cut. **Fix:** design against the cell line's actual sequence; account for ploidy/zygosity.

## Quantitative Thresholds

| Parameter | Value | Source |
|-----------|-------|--------|
| ssODN total length | ~120-200 nt | synthesis ceiling |
| ssODN arm | ~30-60 nt each | below ~30 HDR drops; above ~60 diminishing returns |
| dsDNA/plasmid arm | ~200-800 bp (up to ~2 kb) | ~500-800 bp common sweet spot |
| lssDNA insert | ~0.2-2 kb | Easi-CRISPR range |
| AAV cargo | <= ~4.5 kb (arms included) | packaging limit (hard wall) |
| PITCh microhomology | ~5-25 bp | MMEJ working range |
| **Edit-to-cut distance** | **<= ~10 bp** | HDR incorporation falls with distance (Paquet 2016) |
| Phosphorothioate | 2-3 terminal bases each end | exonuclease resistance |
| Cold shock | 32 C, 24-48 h | G2/M accumulation (Guo 2018) |
| Typical raw HDR | single-digit to ~20% (up to ~38-60% optimized) | minority pathway |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low HDR, mostly indels | no blocking mutation (re-cutting) | add codon-checked PAM/seed block |
| Only the silent mutation incorporated | edit too far from cut | cut within ~10 bp or switch to base/prime editing |
| Toxicity / random integration | dsDNA in sensitive cells | lssDNA or AAV6 |
| No knock-in in neurons/in vivo | HDR donor in post-mitotic cells | HITI/HMEJ |
| AAV donor will not package | arms + insert exceed ~4.5 kb | shorten arms/insert; budget against the cap |

## References

- Richardson CD, Ray GJ, DeWitt MA, Curie GL, Corn JE (2016). Enhancing homology-directed genome editing by catalytically active and inactive CRISPR-Cas9 using asymmetric donor DNA. *Nat Biotechnol* 34(3):339-344.
- Lin S, Staahl BT, Alla RK, Doudna JA (2014). Enhanced homology-directed human genome engineering by controlled timing of CRISPR/Cas9 delivery. *eLife* 3:e04766.
- Paquet D, Kwart D, Chen A, et al. (2016). Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9 (CORRECT). *Nature* 533(7601):125-129.
- Quadros RM, Miura H, Harms DW, et al. (2017). Easi-CRISPR: a robust method for one-step generation of mice carrying conditional and insertion alleles using long ssDNA donors and CRISPR ribonucleoproteins. *Genome Biol* 18:92.
- Nakade S, Tsubota T, Sakane Y, et al. (2014). Microhomology-mediated end-joining-dependent integration of donor DNA in cells and animals using TALENs and CRISPR/Cas9 (PITCh). *Nat Commun* 5:5560.
- Suzuki K, Tsunekawa Y, Hernandez-Benitez R, et al. (2016). In vivo genome editing via CRISPR/Cas9 mediated homology-independent targeted integration (HITI). *Nature* 540(7631):144-149.
- Yao X, Wang X, Hu X, et al. (2017). Homology-mediated end joining-based targeted integration using CRISPR/Cas9 (HMEJ). *Cell Res* 27(6):801-814.
- Guo Q, Mintier G, Ma-Edmonds M, et al. (2018). 'Cold shock' increases the frequency of homology directed repair gene editing in induced pluripotent stem cells. *Sci Rep* 8:2080.
- Untergasser A, Cutcutache I, Koressaar T, et al. (2012). Primer3 -- new capabilities and interfaces. *Nucleic Acids Res* 40(15):e115.

## Related Skills

- grna-design - Choose the guide; HDR consumes its cut site (tightly coupled via edit-to-cut distance)
- base-editing-design - Donor-free alternative for single-base transitions far from any cut
- prime-editing-design - Donor-free precise small edits, and twinPE/PASTE for large insertions
- primer-design/primer-basics - PCR primers for arm amplification and junction genotyping
- primer-design/primer-validation - Check genotyping primers for dimers and hairpins
- primer-design/primer-specificity - Confirm the genotyping amplicon is unique (off-target/pseudogenes)
- sequence-io/read-sequences - Parse GenBank CDS/start/stop features for tag placement and codon-aware design
- variant-calling/variant-annotation - Confirm the installed edit and its consequence
