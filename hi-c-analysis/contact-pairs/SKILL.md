---
name: bio-hi-c-analysis-contact-pairs
description: Turns Hi-C/Micro-C FASTQ into a deduplicated, filtered .pairs file with pairtools and decides whether the library worked. Covers the bwa mem -SP5M / bwa-mem2 / chromap --preset hic alignment idiom (mates mapped as independent single-end reads), pairtools parse vs parse2 and the walks-policy choice (5unique pairwise vs all for Pore-C/Micro-C concatemers), pair-type classification (keep UU and rescued UC), dedup (PCR vs optical/by-tile), select by pair_type/MAPQ/distance, restriction-fragment handling (restrict, Arima dual-enzyme, Micro-C/DNase fragment-free), and allele-specific phasing (pairtools phase to two coolers). The library-QC decision uses % long-range cis as the one-number quality metric, trans as the noise floor, orientation balance as fragment-map-free dangling-end/self-circle QC, and % duplicates as a complexity proxy. Use when processing Hi-C/Micro-C/Omni-C reads into pairs, judging library quality, handling multi-enzyme or restriction-agnostic protocols, or generating allele-specific contacts.
tool_type: cli
primary_tool: pairtools
---

## Version Compatibility

Reference examples tested with: pairtools 1.1+, bwa 0.7.17+ (or bwa-mem2 2.2+), chromap 0.2+, samtools 1.19+, cooler 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

pairtools defaults have shifted across releases (e.g. `parse --max-molecule-size` is 750 bp in 1.1.x; `dedup --backend` defaults to scipy). `parse` and `parse2` report DIFFERENT positions by default - confirm `--report-position` before mixing outputs. If a command errors, introspect with `pairtools <subcommand> --help` and adapt rather than retrying.

# Hi-C Contact Pairs

**"Turn my Hi-C reads into clean contacts and tell me if the library worked"** -> Align mates independently through ligation junctions, classify and deduplicate pairs to a 5'-canonical .pairs file, then read the cis/trans and orientation statistics to decide library quality before any matrix is built.
- CLI: `bwa mem -SP5M ref.fa R1.fq R2.fq | pairtools parse -c chrom.sizes | pairtools sort | pairtools dedup | pairtools stats`

## The Single Most Important Modern Insight -- The Read Count Is a Lie Until Pairs Are Classified; Library Quality Is Decided in pairtools, Not in the Aligner

A Hi-C library's *usable signal* is not "reads sequenced" - it is the **uniquely-mapped, deduplicated, long-range cis** contacts. Everything between FASTQ and the matrix exists to strip a *specific* artifact of proximity-ligation chemistry, and the diagnostic ratios from `pairtools stats` reveal whether the experiment succeeded **before** any compute is spent binning it. Three load-bearing consequences:

1. **% long-range cis is the one-number quality metric; trans is the noise floor.** True crosslink-ligation contacts are overwhelmingly cis and distance-decaying. Random ligation between two unrelated molecules in solution is as likely to be trans as cis-far, so **trans% is a direct readout of the spurious-ligation floor**. A good in-situ human library runs cis>=1kb ~50-65%, inter-chromosomal <10%. But these numbers are **genome-size dependent** - a yeast or bacterial genome legitimately has higher expected trans (more inter-chromosomal volume per cis distance). Never apply a human trans threshold to a microbe.

2. **The ligation junction lives INSIDE the read, so a local/split aligner aligning mates independently is required.** A single read often sequences *through* a ligation junction (locus A | locus B within one read). An end-to-end aligner soft-clips or mis-maps it and the contact is lost. `bwa mem -SP5M` aligns R1 and R2 as **independent single-end reads** (`-SP` skips mate rescue and pairing - proper-pair logic would destroy every long-range and trans contact) and marks the **5'-most chimeric segment primary** (`-5`, the anchor for pairtools' 5' convention). The chimera fraction rises with read length, so on 150bp PE and on Micro-C long reads this is a large, real chunk of contacts.

3. **Keep UU AND rescued UC; selecting only UU silently discards every rescued ligation.** A naive `select pair_type=="UU"` throws away the chimeric reads pairtools successfully reconstructed (UC = combined-unique) - a meaningful fraction on long reads. The 4DN standard keeps **UU and UC**.

## Aligner Taxonomy

| Aligner | Role | Hi-C invocation | When |
|---------|------|-----------------|------|
| bwa mem -SP5M | reference standard; local/split alignment reconstructs in-read junctions | `bwa mem -SP5M -t N ref.fa R1 R2` | default; best inter-contig accuracy in benchmarks |
| bwa-mem2 | drop-in faster reimplementation, identical output, same flags | `bwa-mem2 mem -SP5M -t N ref.fa R1 R2` | when speed matters and the larger index fits RAM |
| chromap --preset hic | ultrafast integrated align + dedup + pairs (4DN .pairs out) | `chromap --preset hic -x idx -r ref.fa -1 R1 -2 R2 -o out.pairs` | ~10x faster scans; trades fine walk-policy control for speed |

`-SP5M`, letter by letter: `-S` skip mate rescue; `-P` skip pairing (no proper-pair rescue); together `-SP` align mates as independent single-end reads. `-5` mark the 5'-most split segment primary (anchors the 5' convention). `-M` is legacy compatibility only (secondary flag 256 vs supplementary 2048); pairtools handles either - never agonize over `-M`, never drop `-SP5`.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Standard in-situ/Omni-C, FASTQ -> matrix | bwa mem -SP5M -> parse (5unique) -> sort -> dedup -> stats | the 4DN/distiller default; restriction-agnostic |
| Fastest scan, control not critical | chromap --preset hic | integrated align+dedup+pairs ~10x faster |
| Multi-way contacts (Pore-C, MC-3C, Micro-C walks) | `parse --walks-policy all` or `parse2 --expand` | 5unique COLLAPSES concatemers to pairwise silently |
| Micro-C / DNase Hi-C | NO fragment map; do NOT apply a 1kb min-distance cut | sub-1kb (nucleosome ladder) is signal, not artifact |
| Arima / Hi-C 3.0 dual-enzyme, fragment-level | fragment file must encode BOTH motifs (`restrict -f`) | a single-enzyme digest file is silently wrong |
| Repeat-heavy genome / stringent loop anchors | raise `parse --min-mapq 30` | reclassifies borderline reads M, drops repeat false anchors |
| Allele-specific / diploid folding | diploid ref + `-XA` suboptimal hits -> `pairtools phase` -> two coolers | needs the suboptimal-score gap to resolve haplotypes |
| Decide whether to sequence deeper | complexity curve from dup model / preseq lc_extrap | a bare dup% without depth is meaningless |
| Build the matrix from clean pairs | -> hic-data-io (`cooler cload pairs`) | binning happens after classification/dedup |
| Annotate boundary/anchor coordinates | -> genome-intervals/bed-file-basics | pairs are 1-based, half-open conventions differ |

## Align: Mates as Independent Single-End Reads

```bash
bwa index ref.fa                                            # or bwa-mem2 index ref.fa
bwa mem -SP5M -t 16 ref.fa R1.fq.gz R2.fq.gz | \            # -SP: independent SE; -5: 5' segment primary
    samtools view -b -@ 8 - > aligned.bam
# chromap fast path (integrated align + dedup + 4DN pairs, no separate pairtools needed):
# chromap -i -r ref.fa -o idx && \
# chromap --preset hic -x idx -r ref.fa -1 R1.fq.gz -2 R2.fq.gz -o sample.pairs
```

## Parse, Sort, Dedup, Select: the pairtools Core

```bash
# Parse alignments into a 5'-canonical .pairs. min-mapq 1 (default) is permissive: only MAPQ-0 is "multi".
pairtools parse -c chrom.sizes --walks-policy 5unique --min-mapq 1 \
    --add-columns mapq --drop-sam aligned.bam | \
pairtools sort --nproc 8 | \                                # block-sort; flips to upper-triangular (5'-canonical)
pairtools dedup --max-mismatch 3 --mark-dups \              # within 3bp on both sides = duplicate; tag DD
    --output-stats sample.dedup.stats | \
pairtools select '(pair_type=="UU") or (pair_type=="UC")' \ # keep both-unique AND rescued chimeric
    -o sample.valid.pairs.gz
```

Dedup MUST run on a flipped, 5'-canonical file (`sort` does the flip); on a non-canonical file dedup under-collapses and dup% reads falsely low. `--max-molecule-size` (750 bp in 1.1.x) governs single-ligation chimera rescue; `--max-inter-align-gap` (20 bp) sets when a coverage gap becomes a null alignment.

## Library QC: the Decision This Skill Owns

`pairtools stats` is the canonical readout. Read it as a funnel, not a single number.

```bash
pairtools stats --bytile-dups -o sample.stats.tsv sample.valid.pairs.gz
# Key fields: frac_dups; frac_cis; cis_1kb+/cis_20kb+; trans; pair_types; dist_freq orientation FF/FR/RF/RR.
```

- **% long-range cis (cis>=1kb, often cis>=20kb)** = signal quality. **trans = noise floor** (genome-size dependent).
- **Orientation vs distance = fragment-map-free dangling-end/self-circle QC.** Above ~1kb the four orientations FF/FR/RF/RR each converge to ~25% (random, the positive QC signal). A short-range **FR (inward) spike** = dangling ends / undigested / self-ligation; a short-range **RF (outward) spike** = self-circles / religation. The distance where orientations equalize is the **minimum reliable contact distance** - derive the min-distance cut from this plot, do not hardcode 1kb (Micro-C structure lives below 1kb).
- **% duplicates = complexity proxy**, but `--bytile-dups` separates OPTICAL dups (patterned NovaSeq flowcells, same tile, adjacent coordinates) from PCR dups. Only the PCR fraction reflects library complexity; reading total dup% as over-amplification wrongly condemns a good library. A bare dup% without the depth it was measured at is meaningless - use the complexity/yield curve (preseq lc_extrap) to decide whether deeper sequencing buys uniques or duplicates.

Apply the QC-derived distance cut without a fragment map:

```bash
pairtools select '(chrom1!=chrom2) or (abs(pos2-pos1) > 1000)' \   # keep trans + cis beyond the orientation-equalization distance
    -o sample.filtered.pairs.gz sample.valid.pairs.gz
```

## Restriction Fragments: Opt-In, Not Default

Modern pipelines SKIP fragment filtering on purpose - the distance cut + dedup + UU/UC filter + balancing absorb the residual, and a digest file is one genome-specific place to mis-specify the enzyme. `pairtools restrict -f frags.bed` is opt-in for: sub-kb / restriction-fragment-resolution maps, capture-Hi-C / 4C-style fragment analysis, and bench QC where the dangling-end *fraction* is the digestion-efficiency readout.

```bash
cooler digest --out frags.bed chrom.sizes ref.fa DpnII     # single-enzyme: DpnII ^GATC
pairtools restrict -f frags.bed -o restricted.pairs.gz parsed.pairs.gz
```

**Arima dual-enzyme has FOUR junction motifs** (GATCGATC, GANTGATC, GANTANTC, GATCANTC); a single-enzyme (DpnII-only) digest file silently mis-assigns fragments. **Micro-C / DNase Hi-C have NO fragment map** (MNase/DNase cut sequence-nonspecifically) - any tool requiring a restriction file cannot process them, which is exactly why the restriction-agnostic pairtools path became the field default.

## Allele-Specific Contacts: pairtools phase -> Two Coolers

**Goal:** Resolve each contact to maternal vs paternal haplotype for allele-specific 3D folding.

**Approach:** Align to a diploid reference reporting suboptimal hits, so each read carries its best alignment on *both* homologs; `pairtools phase` reads the gap between the two best alignment scores to tag each side resolved-hap1 / resolved-hap2 / non-resolved / multi-mapper - the score gap is what separates a genuinely uninformative read from an actual repeat.

```bash
bcftools consensus -H 1 -f ref.fa phased.vcf.gz > hap1.fa   # build a diploid (two-homolog) reference
bcftools consensus -H 2 -f ref.fa phased.vcf.gz > hap2.fa
bwa mem -SP5M ref_diploid.fa R1.fq R2.fq | \                # report suboptimal hits so both homologs are kept
pairtools parse -c chrom.sizes --add-columns XB,AS,XS | \
pairtools phase --phase-suffixes _hap1 _hap2 --tag-mode XB | \
pairtools sort | pairtools dedup -o phased.pairs.gz
# Then split into hap1/hap2/trans pairs and cload each into a SEPARATE cooler.
```

## Per-Method Failure Modes

### Aligned Hi-C as a normal PE library
**Trigger:** plain `bwa mem` without `-SP`, or proper-pair logic. **Mechanism:** mate rescue forces the shotgun insert model on mates from different loci. **Symptom:** trans% collapses, compartments vanish, sparse map. **Fix:** `bwa mem -SP5M`, mates aligned independently.

### Dropped `-5`
**Trigger:** copying a pre-2016 `-SP` command lacking `-5`. **Mechanism:** the 5'-most chimeric segment is not primary, so pairtools picks an inconsistent anchor. **Symptom:** degraded flip/dedup, smeared loops. **Fix:** always `-SP5M` (or `-SP5`).

### walks-policy 5unique on a multi-way protocol
**Trigger:** Pore-C/MC-3C/Micro-C walks parsed with the default. **Mechanism:** 5unique reports only the two 5'-most alignments, collapsing concatemers to pairwise. **Symptom:** "we found few multi-contacts." **Fix:** `--walks-policy all` or `parse2 --expand`.

### Mixing parse and parse2 outputs
**Trigger:** combining a parse2 (.pairs, default outer/junction-anchored) file with a parse (5'-anchored) file. **Mechanism:** the two report positions by different conventions. **Symptom:** coordinates shift by the alignment length; dedup under-collapses; loops smear. **Fix:** one parser/convention per project; prefer plain `parse` if walks are not needed.

### Selecting only UU
**Trigger:** `select pair_type=="UU"`. **Mechanism:** discards rescued chimeric (UC) pairs. **Symptom:** lower valid-pair yield, especially on long reads. **Fix:** keep `(pair_type=="UU") or (pair_type=="UC")`.

### Dedup on a non-canonical file
**Trigger:** dedup run before `sort`/flip, or on mixed parse/parse2 positions. **Mechanism:** duplicates are not in canonical coordinates. **Symptom:** dup% reads falsely low - library looks better than it is. **Fix:** `sort` (flips to 5'-canonical) before `dedup`.

### Optical dups read as PCR dups
**Trigger:** total dup% on a patterned flowcell taken as library complexity. **Mechanism:** optical (same-tile) dups inflate apparent PCR rate. **Symptom:** a good library condemned as over-amplified. **Fix:** `--bytile-dups` / `--output-bytile-stats`; judge complexity on the PCR fraction only.

### Fixed 1kb cut on Micro-C
**Trigger:** applying Hi-C's >1kb min-distance cut to Micro-C. **Mechanism:** Micro-C's nucleosome-ladder signal lives below 1kb. **Symptom:** the structure Micro-C exists to capture is erased. **Fix:** derive the cut from the orientation-vs-distance plot per library.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| cis>=1kb ~50-65% of nodup pairs (good in-situ human) | Dovetail/Arima QC guidance (~approx) | long-range cis is the signal; library- and genome-size dependent |
| inter-chromosomal (trans) <10% (clean), 20-30% acceptable | in-situ Hi-C practice | trans is the spurious-ligation floor; NEVER apply to small genomes |
| FF/FR/RF/RR -> ~25% each above ~1kb | random strand combination at true contacts | convergence is the fragment-map-free positive QC signal |
| `parse --min-mapq` 1 default, raise to 30 for stringency | pairtools default | 1 drops only MAPQ-0; 30 removes repeat-driven false anchors |
| `dedup --max-mismatch` 3 bp | pairtools default | tolerates mapping wobble; 0 over-splits, larger over-collapses complexity |
| `parse --max-molecule-size` 750 bp (1.1.x) | pairtools default | bound on single-ligation chimera rescue; revisit for unusual size selection |
| min-distance cut ~1kb (Hi-C), derive from orientation plot | orientation-equalization distance | the cut is protocol-specific; Micro-C signal is sub-1kb |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| trans% high, no compartments | aligned without `-SP` (proper-pair logic) | re-align `bwa mem -SP5M` |
| Few multi-way contacts on Pore-C/Micro-C | default `--walks-policy 5unique` collapsed walks | `--walks-policy all` / `parse2 --expand` |
| Loops smeared, dedup under-collapses | mixed parse/parse2 position conventions | one parser per project; `sort` before `dedup` |
| Valid-pair yield lower than expected | `select` kept only UU | keep UU and UC |
| dup% suspiciously low | dedup ran before `sort`/flip | sort to 5'-canonical first |
| dup% high on NovaSeq, complexity looks bad | optical dups counted as PCR | `--bytile-dups`; judge on PCR fraction |
| Micro-C structure disappears after filtering | fixed 1kb min-distance cut | derive cut from orientation-vs-distance |
| Fragment assignment wrong on Arima data | single-enzyme digest file | encode all four Arima junction motifs |
| `phase` resolves nothing | aligned to a haploid reference | diploid ref + suboptimal (`-XA`) alignments |

## References

- Open2C, Abdennur N, Fudenberg G, Flyamer IM, Galitsyna AA, Goloborodko A, Imakaev M, Venev SV. 2024. Pairtools: from sequencing data to chromosome contacts. *PLoS Comput Biol* 20(5):e1012164.
- Li H. 2013. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997.
- Zhang H, Song L, Wang X, et al. 2021. Fast alignment and preprocessing of chromatin profiles with Chromap. *Nat Commun* 12:6566.
- Durand NC, Shamim MS, Machol I, et al. 2016. Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. *Cell Syst* 3:95-98.
- Servant N, Varoquaux N, Lajoie BR, et al. 2015. HiC-Pro: an optimized and flexible pipeline for Hi-C data processing. *Genome Biol* 16:259.
- Akgol Oksuz B, Yang L, Abraham S, et al. 2021. Systematic evaluation of chromosome conformation capture assays. *Nat Methods* 18:1046-1055.
- Krietenstein N, Abraham S, Venev SV, et al. 2020. Ultrastructural details of mammalian chromosome architecture (Micro-C). *Mol Cell* 78:554-565.
- Ramani V, Cusanovich DA, Hause RJ, et al. 2016. Mapping 3D genome architecture through in situ DNase Hi-C. *Nat Protoc* 11:2104-2121.
- Abdennur N, Mirny LA. 2020. Cooler: scalable storage for Hi-C data and other genomically labeled arrays. *Bioinformatics* 36:311-316.
- Daley T, Smith AD. 2013. Predicting the molecular complexity of sequencing libraries (preseq). *Nat Methods* 10:325-327.

## Related Skills

- hic-data-io - Bins the deduplicated valid pairs into a .cool/.mcool matrix
- matrix-operations - Balancing and O/E that the binned pairs feed into
- hic-visualization - Render contact maps from the resulting cooler
- read-alignment/bwa-alignment - Aligner upstream; this skill adds the Hi-C `-SP5M` idiom
- alignment-files/duplicate-handling - General duplicate-marking context for the pairtools dedup step
- genome-intervals/bed-file-basics - Coordinate/digest BED handling for restriction fragments and anchors
- genome-assembly/scaffolding - Same Hi-C reads used to order contigs into chromosomes
