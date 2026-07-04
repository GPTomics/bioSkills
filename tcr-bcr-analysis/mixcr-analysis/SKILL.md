---
name: bio-tcr-bcr-analysis-mixcr-analysis
description: Align V(D)J reads and assemble TCR/BCR clonotypes with MiXCR, driven by a chemistry-matched preset. Use when choosing/auditing the preset for a library (5'RACE/template-switch vs multiplex-primer amplicon -> rigid vs floating boundaries; RNA vs gDNA -> --rna/--dna; bulk vs 10x single-cell; UMI vs no-UMI -> tag pattern and barcode collapse; kit presets Takara/NEBNext/QIAseq/BD/MiLaboratory); assembling clonotypes by CDR3 vs VDJRegion; setting the reads-vs-UMI-vs-cell quantitation denominator; exporting native MiXCR fields vs AIRR rearrangement TSV for downstream Immcantation/scirpy/VDJtools; and running alignment/chain-usage QC. Keywords: MiXCR, analyze, align, refineTagsAndSort, assemblePartial, assemble, assembleCells, exportClones, exportAirr, exportQc, CDR3, V(D)J, clonotype, UMI, cell barcode, 10x VDJ, license.
tool_type: cli
primary_tool: MiXCR
---

## Version Compatibility

Reference examples tested with: MiXCR 4.7+ (Java 17)

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mixcr --version` (also prints the JVM) then `mixcr <command> --help` to confirm flags

If a command errors with an unknown-flag, unknown-preset, or missing-license message, run
`mixcr <command> --help` and `mixcr exportPreset --preset-name <name>` and adapt rather than retrying.

Note: MiXCR 4.x is a rearchitecture of 3.x. The hand-built `mixcr analyze amplicon`/`analyze shotgun` pipelines are GONE, replaced by `mixcr analyze <preset>`; `correctAndSortTags` became `refineTagsAndSort`; AIRR export moved to a dedicated `mixcr exportAirr`. Current 4.x needs Java 17 and an activated license (see the licensing gate below). Any 3.x tutorial is stale.

# MiXCR Analysis

**"Extract TCR/BCR clonotypes from my sequencing data"** -> align raw reads to V/D/J/C germline, collapse molecules/cells by barcode, and assemble reads into clonotypes keyed on CDR3 + V + J.
- CLI: `mixcr analyze <preset>` runs the whole ordered pipeline; the underlying stages (`align` -> `refineTagsAndSort` -> `assemblePartial`/`extend` -> `assemble` -> `assembleCells` -> `exportClones`/`exportAirr` -> `qc`) can be run by hand for control.

## The governing principle: the preset IS the analysis, and the wrong one fails silently

In MiXCR 4.x there is no default-correct pipeline. Correctness is ~90% preset choice plus library-chemistry match. `mixcr analyze <preset>` expands the preset into an ordered stage list encoding material (RNA vs DNA), 5'/3' alignment-boundary behavior (rigid vs floating), the barcode tag pattern, the assembling feature, and species defaults. The wrong preset does NOT raise an error -- it emits plausible-but-wrong clonotypes: a mismatched RNA/DNA model or boundary model mis-places V/J boundaries and truncates CDR3; a UMI kit run without a tag pattern skips barcode collapse and inflates diversity with PCR/sequencing artifacts. Because the failure is silent, the load-bearing skill is choosing and AUDITING the preset. Dump exactly what a preset does with `mixcr exportPreset --preset-name <name>` (full resolved parameter YAML), and confirm chemistry with `mixcr exportQc align`/`chainUsage` after the run. A clonotype is an analyst choice, not a fact: it is CDR3 (+ V + J) at a chosen boundary (assembling feature), counted in a chosen denominator (reads vs UMIs vs cells) -- every downstream number depends on these.

## Licensing gate (do this first, or every run fails)

MiXCR 4.x refuses to run any analysis command until a license is activated -- the single most common reason a copied 3.x recipe fails today. Academic/non-profit use is free (obtain a key at platforma.bio/getlicense); for-profit use needs a business license. Activate by any one of:
- `mixcr activate-license` then paste the key (interactive).
- Place `mi.license` (or `~/.mi.license`) in `~/`, next to `mixcr.jar`, or next to the executable.
- Set `MI_LICENSE=<key content>` or `MI_LICENSE_FILE=/path/mi.license` (best for HPC/Docker/CI).

MiXCR also validates the key over the internet periodically. On air-gapped or firewalled compute nodes, whitelist IPv4 `75.2.96.100` and `99.83.215.63` (and the corresponding IPv6) or arrange an offline license, or a job silently stalls waiting on egress.

## Preset selection by library type

The preset must match the exact wet-lab chemistry. Inspect the built-in list with `mixcr exportPreset` and the docs; verify current names against `mixcr analyze --help` since MiLaboratories occasionally renames presets between minor releases.

| Library / chemistry | Preset (verified 4.7) | Material | 5' boundary | UMI/barcode | Biology consequence if mismatched |
|---|---|---|---|---|---|
| 5'RACE / template-switch bulk (e.g. SMARTer) | kit preset, or `generic-amplicon`/`-with-umi` with `--rigid-left-alignment-boundary` | `--rna` | RIGID (5' set by template-switch oligo) | kit UMI or `--tag-pattern` | Floating-left on RACE trims real 5' V sequence; missing tag pattern skips UMI collapse |
| Multiplex-primer amplicon (V/J or V/C primers) | kit preset, or `generic-amplicon` with `--floating-left-alignment-boundary` | `--rna` or `--dna` | FLOATING on the primer side | as designed | Rigid boundary counts primer bases as germline mismatch -> wrong V call, truncated CDR3 |
| gDNA multiplex (genomic template) | `--dna` variant preset | `--dna` (include introns) | floating on primer side | usually none | `--rna` on gDNA drops intron-containing alignments; gDNA count approximates cell count |
| Bulk RNA-seq mining (non-targeted) | `rna-seq` | `--rna` | n/a (fragmented) | none | Needs `assemblePartial` x2 + `extend`; judged by absolute yield, not % aligned |
| 10x single-cell V(D)J (TCR+BCR) | `10x-sc-xcr-vdj` | preset-set | preset-set | CELL+UMI (preset) | Missing cell/UMI pattern -> no pairing, fake diversity; count CELLS not reads |
| 10x 5' GEX repertoire mining | `10x-sc-5gex` | preset-set | preset-set | CELL+UMI | Shallow repertoire mined from GEX; not a substitute for enriched VDJ |
| Takara SMARTer human TCR/BCR | `takara-human-rna-tcr-umi-smarter-v2`, `takara-human-rna-bcr-umi-smarter`, `...-smartseq` | `--rna` | RIGID (template-switch) | 12nt UMI (preset) | Uses the correct RACE boundary + UMI pattern automatically |
| NEBNext immune-seq | `neb-human-rna-xcr-umi-nebnext` (`neb-mouse-...`) | `--rna` | preset-set | UMI (preset) | `xcr` = both TCR and BCR in one preset |
| QIAseq immune | `qiagen-human-rna-tcr-umi-qiaseq` (`...-mouse-...`) | `--rna` | preset-set | UMI (preset) | -- |
| BD Rhapsody single-cell | `bd-human-sc-xcr-rhapsody-cdr3`, `bd-sc-xcr-rhapsody-full-length` | preset-set | preset-set | CELL+UMI | full-length variant enables SHM/contig work |
| MiLaboratories kits | `milab-human-rna-tcr-umi-race`, `milab-human-rna-tcr-umi-multiplex`, `milab-human-dna-tcr-multiplex`, ... | per name | per name | per name | name decodes `<vendor>-<species>-<rna/dna>-<chain>-[umi]-<protocol>` |

For `generic-*` presets `--species <hsa|mmu|...>` is REQUIRED (forgetting it fails or misaligns). `xcr` presets cover TCR and BCR together; single-chain presets (`trb`, `ig`) cover one locus. For gamma-delta (TRG/TRD), use the same generic/kit presets and restrict chains at export with `-c TRG` / `-c TRD` (or a gd-specific kit preset if the wet-lab kit targets gd); note that a gd repertoire is invisible if the library only primed alpha-beta.

## Pipeline stages and where each one fails

`mixcr analyze <preset> R1.fastq.gz R2.fastq.gz out_prefix` runs the ordered stages below; the preset is embedded in the binary `.vdjca`/`.clns` files so hand-run commands only name the preset on `align`.

| Stage | Command | Purpose | Common failure |
|---|---|---|---|
| Align | `mixcr align -p <preset> --species hsa ...` | Reads -> V/D/J/C germline; extract barcodes if tag pattern set | Low alignment rate: wrong species/material/boundaries, untrimmed primers, reads too short to span CDR3 |
| Refine tags | `mixcr refineTagsAndSort` | UMI + cell-barcode error correction and sort | Skipped on a UMI library -> barcode errors become fake clonotypes; memory-heavy (~32 GB on large single-cell) |
| Assemble partial | `mixcr assemblePartial` (run x2) | Overlap fragmented mates that each cover part of CDR3 (RNA-seq/10x) | Needs `align --keep-non-CDR3-alignments` first; on amplicon reads that already span CDR3 it is wasted |
| Extend | `mixcr extend` | Impute unambiguous missing V/J germline ends | Safe for TCR; on BCR can fabricate germline over SHM-mutated ends |
| Assemble | `mixcr assemble` | Collapse alignments into clonotypes by the assembling feature; PCR/error correction, UMI consensus | Wrong assembling feature merges/splits clones; low-quality CDR3 filtered |
| Assemble cells | `mixcr assembleCells` | Single-cell: group per-chain clones by CELL barcode into paired cells | Needs cell tags; barcode contamination -> mispaired cells |
| Export | `mixcr exportClones` / `mixcr exportAirr` | Write clonotype TSV (native or AIRR) | Native field-name mistakes; forgetting `-c/--chains`; not filtering non-productive |
| QC | `mixcr qc`, `mixcr exportQc align`/`chainUsage` | Alignment rate, chain composition, tag coverage | Not run -> silent quality problems pass downstream |

From MiXCR 4.7, presets that do not intrinsically define an assembling feature REQUIRE `--assemble-clonotypes-by <feature>` (e.g. `CDR3`, `VDJRegion`); older tutorials that omit it now error. `CDR3` is the robust default on short reads; `VDJRegion` needs reads/contigs spanning V-through-J and keeps SHM variants separate (useful for BCR full-length).

## The quantitation denominator: reads vs UMIs vs cells

Clonotype abundance is only meaningful relative to the chemistry. Report the right unit or reintroduce the bias the chemistry was meant to remove:
- Non-UMI bulk: abundance = `readCount` (`cloneCount` is an alias). PCR-amplification biased -- not a molecule count.
- UMI bulk: after `refineTagsAndSort`, report `uniqueMoleculeCount` (generic form `uniqueTagCount Molecule`), NOT reads. Reporting reads on a UMI library re-adds the amplification bias the UMIs removed.
- Single-cell: the unit is the CELL (`uniqueTagCount Cell` / cellGroup), not reads or UMIs.

## Export: native MiXCR fields vs AIRR

MiXCR's native export headers are NOT AIRR or VDJtools names. Downstream renaming to a chosen schema is a user-side step; the field names to select from MiXCR are its own.

```bash
mixcr exportClones -c TRB \
    -cloneId -readCount -readFraction -uniqueMoleculeCount \
    -nSeqCDR3 -aaSeqCDR3 -bestVGene -bestJGene -allVHitsWithScore -isProductive VRegion \
    clones.clns clones_TRB.tsv
```

Key native fields: `cloneId`, `readCount`/`readFraction` (aliases `cloneCount`/`cloneFraction`), `uniqueMoleculeCount`, `nSeqCDR3`/`aaSeqCDR3` (CDR3 nt/aa -- the headline fields), `bestVGene`/`bestJGene`/`bestCGene` (gene-level), `bestVHit` (allele-level best), `allVHitsWithScore` (full hit list), `isProductive <feature>`. Filter flags: `-c/--chains TRB`, `-o` (drop out-of-frame), `-t` (drop stops), `--export-productive-clones-only`.

For AIRR-schema interchange (Immcantation, scirpy, any AIRR tool) use the dedicated command, which emits `sequence_id`, `v_call`, `d_call`, `j_call`, `junction`, `junction_aa`, `productive`, `duplicate_count`, `cell_id`:

```bash
mixcr exportAirr clones.clns clones.airr.tsv
```

Field-name traps (do NOT use as MiXCR selectors): `count`/`frequency`/`cdr3_aa`/`vGene` are VDJtools/AIRR conventions, not MiXCR headers. A downstream tool expecting AIRR names should be fed `exportAirr` output, not renamed native output.

## The D-gene caveat: never key or trust the D call in TRB/IGH

Dbeta and Dh segments are short (~12-16 nt) and heavily trimmed at both ends with N-additions between; the surviving germline-matchable D stretch is often 0-5 nt, statistically indistinguishable from random junctional nucleotides. A substantial fraction of TRB rearrangements have no detectable D at all (de Greef & de Boer 2021 *PNAS* 118:e2104367118), and any "longest germline D match" over-calls D by chance. Treat `bestDGene`/`allDHitsWithScore` as unreliable: never use the D call as a clonotype key, never stratify biology by D usage without heavy skepticism, and expect large tool-to-tool D disagreement. Clonotypes are keyed on CDR3 + V + J -- not D. This caveat is TRB- and IGH-specific: the TRD (delta) chain can incorporate one to two D segments in tandem, giving more germline D content than TRB's single heavily-trimmed D, so the D call is more informative for gamma-delta work (the junction is still highly diverse from N-additions).

## QC: match preset to chemistry, catch cross-contamination

```bash
mixcr qc clones.clns
mixcr exportQc align results/*.clns qc_align.pdf
mixcr exportQc chainUsage results/*.clns qc_chainUsage.pdf
```

Read `exportQc align`: targeted amplicon should align high (often >80-90%); a low rate signals wrong species/library/boundaries or untrimmed primers, and "absent CDR3" means reads too short or wrong boundaries. RNA-seq mining legitimately aligns a tiny fraction (only receptor-overlapping reads) -- judge it by absolute clonotype yield, not %. Read `chainUsage`: a TRB library showing appreciable IGH signals cross-contamination or index hopping on patterned flowcells. A huge reads-to-clonotypes drop (millions -> thousands, worse after UMI collapse) is normal; a tiny clone count with high alignment suggests over-aggressive filtering or a wrong assembling feature.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Every processing command (align/analyze/assemble) refuses to run | No activated license (4.x mandatory for the pipeline; `mixcr --version`/`exportPreset` still work) | `mixcr activate-license` or set `MI_LICENSE_FILE`; whitelist phone-home IPs on firewalled nodes |
| `mixcr analyze amplicon ...` unknown | 3.x command removed in 4.x | Use `mixcr analyze <preset>`; pick a chemistry-matched preset |
| Runs cleanly but clonotypes look wrong (truncated CDR3, odd V calls) | Wrong preset / boundary / material -- silent, no error | Match preset to chemistry; audit with `mixcr exportPreset`; check `exportQc align` |
| Diversity far too high, many near-identical clones | UMI kit run without tag pattern -> no barcode collapse | Use the UMI preset or add `--tag-pattern`; ensure `refineTagsAndSort` ran; report `uniqueMoleculeCount` |
| RNA-seq run yields almost no clonotypes | No `assemblePartial`/`extend`; partials filtered at align | `align --keep-non-CDR3-alignments`, `assemblePartial` twice, then `extend` (or use `rna-seq` preset) |
| `assemble` errors asking for an assembling feature | Preset lacks intrinsic feature (4.7+) | Add `--assemble-clonotypes-by CDR3` (or `VDJRegion`) |
| Downstream AIRR tool rejects the table | Fed native MiXCR headers, not AIRR | Export with `mixcr exportAirr`, not renamed `exportClones` |
| D-gene usage plot looks meaningless / irreproducible | Trusting the near-unassignable D call in TRB/IGH | Drop D from keys and usage; report V/J only |
| `--species` missing on a generic preset | Generic presets require species | Add `--species hsa` (or `mmu`, taxon id) |
| refineTagsAndSort out-of-memory on single-cell | Barcode-heavy step needs large heap | `mixcr -Xmx32g refineTagsAndSort ...` |

## Related Skills

- vdjtools-analysis - Downstream diversity and overlap on bulk clonotypes
- immcantation-analysis - BCR clonal clustering, SHM and lineage from AIRR output
- scirpy-analysis - Single-cell VDJ integration with gene expression
- repertoire-visualization - Plot V/J usage and clonal structure
- specificity-annotation - Antigen-specificity clustering and database lookup
- read-qc/adapter-trimming - Upstream read QC and adapter handling
- workflows/tcr-pipeline - End-to-end orchestration

## References

- Bolotin DA, et al. MiXCR: software for comprehensive adaptive immunity profiling. *Nat Methods* 12:380-381 (2015).
- Bolotin DA, et al. Antigen receptor repertoire profiling from RNA-seq data. *Nat Biotechnol* 35:908-911 (2017).
- de Greef PC, de Boer RJ. TCRbeta rearrangements without a D segment are common, abundant, and public. *PNAS* 118:e2104367118 (2021).
- Vander Heiden JA, et al. AIRR Community standardized representations for annotated immune repertoires. *Front Immunol* 9:2206 (2018).
- MiXCR documentation. https://mixcr.com/mixcr/ (presets, mixins, exportClones/exportAirr, licensing, QC).
