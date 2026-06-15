# HDR Template Design - Usage Guide

## Overview

Designs donor/repair templates for precise CRISPR knock-ins -- ssODN, long-ssDNA (Easi-CRISPR), dsDNA/plasmid, and AAV6 donors. The skill picks the format by edit size and whether the target cell cycles, sizes the homology arms, requires a guide that cuts within ~10 bp of the edit, and adds a mandatory codon-checked blocking (PAM/seed) mutation so the corrected allele is not re-cut. It frames the HDR-vs-NHEJ-vs-MMEJ pathway competition and the HITI/HMEJ/PITCh alternatives for cells where HDR is the wrong tool. Its central discipline: HDR is the minority pathway, and a donor without a blocking mutation is self-destructing.

## Prerequisites

```bash
pip install biopython primer3-py
# primer3-py designs primers (primer3.bindings.design_primers); arm extraction and codon-aware
#   blocking are the skill's own BioPython code. Design arms/guide against the ACTUAL cell-line
#   sequence, not GRCh38. lssDNA production (Easi-CRISPR) and AAV donors are wet-lab methods.
```

## Quick Start

Tell your AI agent what you want to do:
- "Design an ssODN to correct this point mutation, with a blocking mutation"
- "Add a GFP tag to MYC -- which donor format should I use?"
- "My HDR gives mostly indels and almost no edit -- what's wrong?"
- "I'm knocking in to neurons -- how should the donor differ?"
- "The nearest guide cuts 25 bp from my edit -- is HDR still the right tool?"

## Example Prompts

### Point mutations and tags

> "Design an ssODN to install the HBB sickle-cell correction. Confirm a guide cuts within 10 bp, add a synonymous PAM-blocking mutation so the allele isn't re-cut, give me both strands to test, and add phosphorothioate caps."

> "Add a C-terminal FLAG tag to my gene. Pick the donor format and explain why, and design junction-genotyping primers."

### Format and pathway

> "I'm knocking a 1.5 kb cassette into primary T cells -- ssODN, plasmid, lssDNA, or AAV?"

> "Why does dsDNA donor kill my iPSCs, and what should I use instead?"

> "I need a knock-in in post-mitotic neurons -- HDR isn't working. What now?"

### Troubleshooting

> "Low HDR, lots of indels -- walk me through the likely cause."

> "Only the silent mutation got incorporated, not my edit -- why?"

## What the Agent Will Do

1. Establish the cell type (cycling vs post-mitotic) and edit size -- this picks HDR vs HITI/HMEJ and the donor format.
2. Confirm a guide cuts within ~10 bp of the edit; if not, recommend a closer guide or base/prime editing.
3. Choose the format (ssODN / lssDNA / dsDNA / plasmid / AAV6) and size the arms.
4. Add a mandatory codon-checked blocking mutation (synonymous PAM disruption preferred; silent seed mutations otherwise).
5. For ssODN, emit both strands and symmetric/asymmetric variants to test, and recommend phosphorothioate end-protection.
6. Suggest HDR enhancers as ranked experiments (cell-cycle timing / cold shock / M3814 / 53BP1 first).
7. Design arm-amplification and junction-genotyping primers, and recommend screening enough clones.

## Tips

- A donor without a blocking mutation is self-destructing -- the corrected allele keeps an intact PAM, Cas9 re-cuts it, and the indel reads out as "HDR failed." Block the PAM (codon-checked) or the seed.
- HDR is the minority pathway: single-digit to low-double-digit percent unenhanced is normal; correct edits are recovered by screening clones, not by expecting bulk efficiency.
- Edit-to-cut distance couples donor to guide: place the cut within ~10 bp, or only the blocking mutation gets incorporated. If no guide cuts close, switch to base/prime editing.
- The Richardson asymmetric/non-target-strand ssODN rule is a prior to test, not a law -- generate both strands; phosphorothioate end-caps are the one reliable cheap win.
- dsDNA is toxic and integrates randomly -- prefer lssDNA (Easi-CRISPR) for 0.2-2 kb and AAV6 for primary cells.
- Post-mitotic / in-vivo cells barely do HDR -- use HITI (NHEJ-based, reverse-orientation target sites) or HMEJ, not an HDR donor.
- HDR enhancers are ranked experiments, not multipliers: cell-cycle timing / cold shock / PS caps / M3814 / 53BP1 first; SCR7 and RS-1 are unreliable.
- Design against the actual cell-line genome (and account for ploidy), not the reference.

## Related Skills

- grna-design - Choose the guide; HDR consumes its cut site (coupled via edit-to-cut distance)
- base-editing-design - Donor-free alternative for single-base transitions far from any cut
- prime-editing-design - Donor-free precise small edits, and twinPE/PASTE for large insertions
- primer-design/primer-basics - PCR primers for arm amplification and junction genotyping
- primer-design/primer-validation - Check genotyping primers for dimers and hairpins
- primer-design/primer-specificity - Confirm the genotyping amplicon is unique (off-target/pseudogenes)
- sequence-io/read-sequences - Parse GenBank CDS/start/stop features for tag placement
- variant-calling/variant-annotation - Confirm the installed edit and its consequence
