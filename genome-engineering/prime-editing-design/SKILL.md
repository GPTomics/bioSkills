---
name: bio-genome-engineering-prime-editing-design
description: Designs pegRNAs and nicking guides for prime editing (PE) -- choosing the nick/strand, tuning the primer-binding site (PBS) and reverse-transcription template (RTT) as a per-locus panel, selecting the PE system (PE2/PE3/PE3b/PE4/PE5/PEmax/PE7), adding MMR-evading and PAM-disrupting silent edits, appending epegRNA 3' motifs (tevopreQ1/mpknot), and ranking with PRIDICT/DeepPrime. Covers twinPE/PASTE for large insertions and the prime-vs-base-editing decision. Use when designing a scarless point mutation, small insertion/deletion, or any of the 12 base conversions without a double-strand break, when efficiency is low and MMR inhibition or pegRNA stabilization is needed, or when routing a large insertion to an integrase method. Generic guide scoring and base editing are separate skills.
tool_type: mixed
primary_tool: PrimeDesign
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, PrimeDesign 1.2+ (Docker), PRIDICT2.0 (web/code).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

PrimeDesign is Docker-only (no pip) and takes the edit inline in a single string with exact parenthesis notation (below) -- the most-hallucinated thing in PE tooling; verify it against the repo, never reconstruct from memory. Outcome-prediction models are trained mostly on HEK293T + small edits (<=3 bp); their scores are priors, not measurements, and degrade off-distribution. The PE *system* (MMR status, expressed vs synthetic pegRNA) drives efficiency more than any oligo tweak.

# Prime Editing Design

**"Install a precise small edit without a double-strand break"** -> Establish the edit, cell type, MMR status, and delivery; choose the nick position/strand; design a *panel* of PBS x RTT combinations; pick the PE system; add the free wins (PAM-disrupting + MMR-evading silent edits, a 3' motif); rank with a model; and test.
- CLI (Docker): `PrimeDesign` generates ranked pegRNA + nicking-guide components from a reference + edit string
- Python: assemble/sweep PBS x RTT panels with `Bio.Seq`; enforce the don't-end-on-C and 5'-G rules
- Web/code: PRIDICT2.0 / DeepPrime rank candidates by intended-edit and indel rate

## The Single Most Important Modern Insight -- there is no universal PBS/RTT optimum, and the *system* choice carries the order of magnitude

Two reframes:

1. **PBS and RTT length are parameters to optimize per locus, not constants to look up.** The PBS x RTT optimum is locus-specific -- it depends on local GC (which sets the PBS annealing Tm), the edit, the nick-to-edit distance, and chromatin. A high-GC target wants a *short* PBS; a low-GC target a *long* one; the "13/15" that is perfect at one locus is useless 200 bp away. A hard-coded default produces a sequence that *looks* valid, so nothing flags it until the data come back at 2%. The correct deliverable is **a ranked panel** (a few PBS x a few RTT x the viable nicks), tested or model-ranked -- emitting a single pegRNA is the tell of someone who has never run PE.

2. **Prime editing efficiency is a cellular-genetics problem, not just oligo design.** The cell's mismatch repair (MMR; MutSalpha/MutLalpha) detects the edit:original heteroduplex and excises the *edited* strand, reverting it and spawning indels. The biggest post-2019 jump was not a better PBS -- it was **inhibiting MMR (MLH1dn -> PE4/PE5, ~7.7x average)**. The second was **stopping the pegRNA 3' end from being degraded** (epegRNA motifs; PE7's La protein). Design now means choosing the *system* (PE2 vs PE3b vs PE5max+epegRNA vs PE7) as much as the sequence. First branch: what edit, what cell type, MMR-proficient or not, expressed or synthetic.

## Mechanism (the design rules fall out of it)

The prime editor (Anzalone 2019) is **Cas9 H840A nickase + engineered M-MLV reverse transcriptase**, programmed by a **pegRNA** = sgRNA (spacer + scaffold) with a 3' extension read 5'->3' as **[RTT][PBS]**. (1) The nickase cuts the protospacer (PAM) strand ~3 nt 5' of the PAM, exposing a free 3'-OH. (2) The **PBS** anneals to that nicked 3' end (the genomic strand becomes the primer). (3) The RT extends through the **RTT**, synthesizing a new 3' DNA flap that *encodes the edit*. (4) FEN1-type nucleases preferentially excise the unedited 5' flap, favoring incorporation of the edited 3' flap; ligation seals it. (5) The resulting heteroduplex is resolved by MMR -- which preferentially reverts the edit (hence the MMR section below). Consequences: PBS length is tuned by annealing Tm; **RTT length = nick-to-edit distance + edit + 3' homology tail (~10-16 nt)**; efficiency falls as the edit moves farther from the nick; the edit must lie within the RTT.

## The PE System Stack -- orthogonal axes, not a "bigger number is better" ranking

| System | Adds over previous | Acts on | Cite |
|--------|--------------------|---------|------|
| PE1 | Cas9 H840A + **wild-type** M-MLV RT | proof of concept | Anzalone 2019 |
| **PE2** | **engineered M-MLV RT** (pentamutant) | the workhorse enzyme | Anzalone 2019 |
| PE3 | + second **nicking sgRNA** on the non-edited strand (~1.5-4x) | flap resolution / MMR strand bias -- **but raises indels** (transient near-DSB) | Anzalone 2019 |
| **PE3b** | PE3 ngRNA matching only the **edited** sequence -> nick fires *after* the edit | near-eliminates PE3's indels; **only possible when the edit makes/breaks a protospacer** | Anzalone 2019 |
| PE4 | PE2 + **MLH1dn** (dominant-negative MMR) (~7.7x avg) | MMR globally | Chen 2021 |
| PE5 | PE3 + MLH1dn | second nick + MMR | Chen 2021 |
| **PEmax** | optimized protein (codon, NLS, R221K/N394K, linker); +MLH1dn = **PE4max/PE5max** | the protein | Chen 2021 |
| PE7 | PEmax-family + **La-protein** RBD capping the pegRNA 3' end | pegRNA stability | Yan 2024 |

The expert move is to reason about which axis the problem needs: low efficiency in an MMR-active cell -> add MLH1dn; too many indels -> drop to PE2 or design **PE3b** (not PE3); short pegRNA half-life -> epegRNA/PE7. Note: in **MMR-deficient lines** (HCT116, many tumor lines) PE2 already behaves like PE4, so MLH1dn adds nothing -- benchmark numbers from such lines overstate the gain in MMR-proficient primary cells. **PE5max + epegRNA is the modern default workhorse** for hard, MMR-active contexts.

## pegRNA Parameters & the Free Wins

- **PBS (~8-17 nt; start ~11-15):** tune to annealing Tm/GC, not a fixed length. pegFinder's starting heuristic is PBS ~= 24 - (GC%/5), clamped 8-17; test a small ladder (e.g. 10/13/15/17).
- **RTT:** = nick-to-edit + edit + ~10-16 nt 3' homology. Shorter RTT is usually more efficient -- use the shortest that spans the edit with adequate homology, then test a couple.
- **Don't end the synthesized flap on a C** (a C at the +1 templated position lowers efficiency; PrimeDesign exposes `--filter_c1_extension`).
- **5' G for U6:** prepend a G if the spacer lacks one -- **prepend, do not replace** the first base (replacing creates a spacer:target mismatch).
- **PAM-disrupting silent edit (free win):** if the edit (or an added silent change) destroys the protospacer/PAM, the editor cannot re-nick the edited strand -> fewer indels, and the change doubles as an MMR-evading mismatch. Always check whether the edit can be routed to disrupt the PAM.
- **MMR-evading bystander edits (free win):** add 1-2 *silent* substitutions next to the intended edit to make a >=3-bp edited "bubble" that MMR recognizes less efficiently -> higher correct-edit yield. Trivial in coding sequence (synonymous codons); the tactic of choice before reaching for MLH1dn.

## epegRNA 3' Motifs & pegRNA Stability

The pegRNA 3' extension (RTT+PBS) is single-stranded RNA that is **exonucleolytically degraded** before it can prime RT -- an invisible failure (the molecule is made, just chewed back). **epegRNAs** append a structured pseudoknot motif to the 3' end (Nelson 2022): use **tevopreQ1** by default (~3-4x gain, no added off-target); **mpknot** is larger and benefits most from a **pegLIT**-designed linker (tevopreQ1/evopreQ1 often work linker-free). **PE7** (La protein) attacks the same degradation from the protein side and is **partly redundant** with epegRNAs (PE7's gains are largest with plain pegRNAs) -- don't stack them as if independent. For synthetic (non-expressed) pegRNAs where a folded motif is awkward, PE7 / La-optimized 3' chemistry is the lever instead.

## Outcome Prediction (rank, but still test)

| Model | Predicts | Cite |
|-------|----------|------|
| **PRIDICT / PRIDICT2.0** | intended-edit + unintended (indel) rate; 2.0 is chromatin-aware across lines | Mathis 2023 *Nat Biotechnol* 41:1151; Mathis 2025 *Nat Biotechnol* 43:712 |
| **DeepPrime / DeepPrime-FT** | efficiency across 8 PE systems x 7 cell types, edits <=3 bp | Yu 2023 *Cell* 186:2256 |
| Easy-Prime | XGBoost pegRNA design with RNA-structure features | Li 2021 *Genome Biol* 22:235 |

Limits: trained mostly on HEK293T + small edits; scores degrade for large edits, untrained cell types, primary/iPS cells, and in vivo loci. A high score says "worth synthesizing," not "will work in the target cell." Report **edit:indel purity, not efficiency alone** (PE3's indel liability hides when only the intended-edit rate is reported).

## Large / Advanced Edits (single-pegRNA PE runs out of room)

| Strategy | Mechanism | Size | Cite |
|----------|-----------|------|------|
| **twinPE** | two pegRNAs template complementary flaps -> replacement/deletion/inversion | up to ~hundreds bp; +recombinase -> kb | Anzalone 2022 *Nat Biotechnol* 40:731 |
| GRAND editing | dual pegRNAs, RTTs complementary to each other (non-genomic) -> template-free insertion | up to a few hundred bp (drops sharply >~400 bp) | Wang 2022 |
| **PASTE** | PE writes a serine-integrase attB site, integrase drops in a donor | **~10-36 kb**, DSB-free | Yarnall 2023 *Nat Biotechnol* 41:500 |

Route "knock in a 2 kb reporter" to twinPE+integrase/PASTE (or HDR/HITI) -- a single giant-RTT pegRNA is a category error.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| C->T / G->A or A->G / T->C transition, base positionable in a window | -> base-editing-design | BE is higher-efficiency, cleaner, no flap/MMR competition for its transition |
| Any of the other small edits (other transversions, small indels, combined) | prime editing, panel of PBS x RTT | PE owns the precise-small-edit-without-a-DSB box |
| Low efficiency in an MMR-proficient cell | PE4/PE5 (MLH1dn) + MMR-evading silent edits | MMR is the dominant barrier |
| Indels unacceptable (therapeutic) | PE2 or **PE3b** (not PE3) | PE3's second nick raises indels |
| Expressed pegRNA | add a tevopreQ1 3' motif (PE5max+epegRNA default) | fixes invisible 3'-degradation |
| Large insertion (genes/tags, >~hundreds bp) | -> twinPE+integrase / PASTE / hdr-template-design | beyond single-pegRNA flap capacity |
| Knockout only (any frameshift) | -> grna-design (plain Cas9) | PE's precision is wasted; nuclease is simpler/more efficient |
| Validate edits | -> crispr-screens/crispresso-editing | quantify intended-edit and indel rates from amplicons |

## Generate Designs with PrimeDesign (verified notation)

**Goal:** Produce ranked pegRNA + nicking-guide candidates for a precise edit.

**Approach:** Encode the reference and edit in ONE inline string with PrimeDesign's exact parenthesis notation, then run the Docker CLI; it sweeps PBS/RTT, ranks pegRNAs (PAM-disrupted preferred), and nominates ngRNAs. Do not hand-roll the design as the only step.

```bash
# PrimeDesign edit-string notation (verify against the repo README; the most-hallucinated PE detail):
#   substitution:  ...AAACG(T/A)CTTCC...        # ref/edit, slash-separated
#   insertion:     ...AAACGT(+CTT)CTTCC...      # bare leading + (also (/CTT))
#   deletion:      ...AAAAC(-GTCT)TCCAAT...     # bare leading - (also (GTCT/))
#   combinatorial: GCCTGTGACTAACTGC(G/T)CCA(+ATCG)AAACGTC(-TTCC)AATCCCCTTATCCAATTTA
docker run -v ${PWD}/:/DATA -w /DATA pinellolab/primedesign primedesign_cli \
  -f edits.csv -pbs 10 12 14 -rtt 10 16 22 -nick_dist_min 0 -nick_dist_max 100 -out designs/
```

## Sweep a PBS x RTT Panel and Enforce the Hard Rules

**Goal:** Build a small, ordered panel of pegRNA extensions for one nick, applying the don't-end-on-C and 5'-G rules.

**Approach:** For each PBS length, take the reverse complement of the genomic sequence 5' of the nick; for each RTT length, build the edited 3' flap and reject extensions whose first templated base is C. Rank the panel by a model (PRIDICT/DeepPrime) for synthesis. (See `examples/prime_editing_design.py`.)

```python
from Bio.Seq import Seq

def prepend_u6_g(spacer):
    return spacer if spacer.startswith('G') else 'G' + spacer   # prepend, never replace
```

## Per-Method Failure Modes

### One pegRNA from a fixed PBS=13/RTT=15
**Trigger:** treating PBS/RTT as constants. **Mechanism:** the optimum is locus-specific (GC/Tm/nick distance/chromatin). **Symptom:** valid-looking pegRNA, ~2% editing. **Fix:** design and test a PBS x RTT panel; rank with PRIDICT2.0/DeepPrime.

### Designed for the edit, ignored the repair machinery
**Trigger:** installing only the literal intended base. **Mechanism:** MMR reverts the edit; an intact PAM lets the editor re-nick. **Symptom:** low yield + indels. **Fix:** add a PAM-disrupting silent edit and 1-2 MMR-evading silent edits; use PE4/PE5 (MLH1dn) in MMR-active cells.

### Reached for PE3 when PE3b was available
**Trigger:** reading the ladder as a scalar. **Mechanism:** PE3's second nick is a transient near-DSB. **Symptom:** good efficiency, unacceptable indels. **Fix:** if the edit makes/breaks a protospacer, design PE3b; otherwise drop to PE2/PE4.

### Reported % editing without % indels
**Trigger:** efficiency-only readout. **Mechanism:** PE yields a mix (edit/unedited/indel). **Symptom:** a "40%" pegRNA that throws 15% indels looks fine. **Fix:** report edit:indel purity (PRIDICT predicts both).

### Trusted a model score off-distribution / forgot the locus
**Trigger:** picking the top-scored pegRNA, skipping the panel, in a non-HEK293T context. **Mechanism:** models are trained on HEK293T + small edits; chromatin dominates and is invisible to sequence. **Symptom:** "designed perfectly, didn't work." **Fix:** weight the model less far from training; still test; a closed locus may sink any design.

### 5' G replaced, or flap ends on C, or large insert forced into one pegRNA
**Trigger:** `'G'+spacer[1:]`; RTT ending on C; 2 kb into one RTT. **Mechanism:** spacer:target mismatch; +1-C re-incorporation; flap can't template/resolve. **Fix:** prepend the G; shift RTT off a terminal C; route large inserts to twinPE/PASTE.

## Quantitative Thresholds

| Parameter | Value | Source |
|-----------|-------|--------|
| PBS length | ~8-17 nt, tuned to Tm/GC (start ~24-GC%/5) | Anzalone 2019; pegFinder heuristic |
| RTT | edit + ~10-16 nt 3' homology; shortest workable | Anzalone 2019 |
| Nick-to-edit | as small as possible; efficiency falls with distance | Anzalone 2019 |
| PE3 ngRNA distance | ~40-100 bp (sweet spot ~50-90), non-edited strand | Anzalone 2019 |
| Flap +1 base | not C | Anzalone 2019 / PrimeDesign `--filter_c1_extension` |
| MMR inhibition gain | ~7.7x avg (MMR-proficient cells only) | Chen 2021 |
| epegRNA 3' motif | tevopreQ1 default; ~3-4x | Nelson 2022 |
| Deliverable | a ranked panel, report edit:indel purity | field practice |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Editing ~2% despite a "perfect" pegRNA | fixed PBS/RTT, unfavorable locus | test a panel; consider MLH1dn/epegRNA; the locus may be closed |
| High indels with PE3 | second nick on non-edited strand | use PE3b (if the edit makes/breaks a protospacer) or PE2 |
| PrimeDesign mis-encodes the edit | wrong inline notation | use exact `(ref/edit)`/`(+ins)`/`(-del)`; verify against the repo |
| No benefit from MLH1dn | MMR-deficient cell line | PE2 already behaves like PE4 there |

## References

- Anzalone AV, Randolph PB, Davis JR, et al. (2019). Search-and-replace genome editing without double-strand breaks or donor DNA. *Nature* 576(7785):149-157.
- Chen PJ, Hussmann JA, Yan J, et al. (2021). Enhanced prime editing systems by manipulating cellular determinants of editing outcomes. *Cell* 184(22):5635-5652.
- Nelson JW, Randolph PB, Shen SP, et al. (2022). Engineered pegRNAs improve prime editing efficiency. *Nat Biotechnol* 40(3):402-410.
- Yan J, Oyler-Castrillo P, Ravisankar P, et al. (2024). Improving prime editing with an endogenous small RNA-binding protein. *Nature* 628(8008):639-647.
- Anzalone AV, Gao XD, Podracky CJ, et al. (2022). Programmable deletion, replacement, integration and inversion of large DNA sequences with twin prime editing. *Nat Biotechnol* 40(5):731-740.
- Yarnall MTN, Ioannidi EI, Schmitt-Ulms C, et al. (2023). Drag-and-drop genome insertion of large sequences without double-strand DNA cleavage using CRISPR-directed integrases (PASTE). *Nat Biotechnol* 41(4):500-512.
- Hsu JY, Grunewald J, Szalay R, et al. (2021). PrimeDesign software for rapid and simplified design of prime editing guide RNAs. *Nat Commun* 12:1034.
- Chow RD, Chen JS, Shen J, Chen S (2021). A web tool for the design of prime-editing guide RNAs (pegFinder). *Nat Biomed Eng* 5(2):190-194.
- Mathis N, Allam A, Kissling L, et al. (2023). Predicting prime editing efficiency and product purity by deep learning (PRIDICT). *Nat Biotechnol* 41(8):1151-1159.
- Mathis N, Allam A, Talas A, et al. (2025). Machine learning prediction of prime editing efficiency across diverse chromatin contexts (PRIDICT2.0). *Nat Biotechnol* 43(5):712-719.
- Yu G, Kim HK, Park J, et al. (2023). Prediction of efficiencies for diverse prime editing systems in multiple cell types (DeepPrime). *Cell* 186(10):2256-2272.
- Li Y, Chen J, Tsai SQ, Cheng Y (2021). Easy-Prime: a machine learning-based prime editor design tool. *Genome Biol* 22:235.

## Related Skills

- base-editing-design - Preferred for the single transition a base editor can make
- grna-design - Generic spacer scoring; plain-nuclease knockout when precision is unneeded
- off-target-prediction - pegRNA spacer and PE3 nicking-guide off-target considerations
- hdr-template-design - Large-insertion alternative (HDR/HITI) when PASTE/twinPE is not used
- crispr-screens/prime-editing-screens - Pooled prime-editing screen analysis
- crispr-screens/crispresso-editing - Quantify intended-edit vs indel rates from amplicons
- variant-calling/variant-annotation - Identify the pathogenic variant to correct or install
