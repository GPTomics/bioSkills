---
name: bio-genome-assembly-hifi-assembly
description: Assembles haplotype-resolved diploid and telomere-to-telomere (T2T) genomes from PacBio HiFi reads with hifiasm (HiFi-only, Hi-C, or trio phasing) and verkko (HiFi + ultralong ONT for T2T), extracting contigs from GFA and routing phasing QC to k-mer/trio metrics. Covers why a primary assembly is a haplotype mosaic that exists in no cell, partial-vs-full phasing (the .bp. vs .dip. filename convention), the purge-default trap on inbred samples, the --hom-cov coverage-estimate alarm, and verkko-vs-hifiasm for T2T. Use when assembling a diploid eukaryote from HiFi, phasing haplotypes with parents (trio) or Hi-C, deciding whether to chase T2T, or diagnosing switch errors invisible to N50/BUSCO/QV.
tool_type: cli
primary_tool: hifiasm
---

## Version Compatibility

Reference examples tested with: hifiasm 0.25.0+, yak 0.1+, verkko 2.3+, gfatools 0.5+, meryl 1.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

hifiasm output filenames AND the default purge level are version-dependent: confirm the `.bp.*`/`.dip.*` prefixes and `-l` defaults against `hifiasm --help` and the man page for the installed build before scripting around them. The k-mer-size convention also matters more than the binary version: hifiasm/yak trio uses k=31; verkko/Merqury hap-mers commonly use k=30 - mixing them silently produces garbage hap-mer matching. If a command errors, introspect the installed tool and adapt rather than retrying.

# HiFi Assembly

**"Assemble a diploid genome from HiFi reads"** -> Build two cleanly phased haplotypes (not one mosaic primary) from accurate long reads, choosing the phasing mechanism from the sample type and available data, and validate phasing with k-mer QC the headline metrics cannot see.
- CLI: `hifiasm -o prefix -t 32 reads.hifi.fq.gz` (HiFi-only), `--h1/--h2` (Hi-C), `-1/-2 *.yak` (trio); `verkko --hifi ... --nano ...` (T2T)

## The Single Most Important Modern Insight -- A Primary Assembly Is a Mosaic Chimera, Not a Haplotype

A "primary assembly" is not a genome that exists in any cell. At every heterozygous block the assembler picks one allele, and which one it picks switches arbitrarily from block to block - so the primary is a stitched chimera: maternal here, paternal there, matching no gamete, no parent, no individual. It is a fine *haploid representation* for "roughly where is gene X" and a terrible substrate for anything allele-aware (phased variant calling, allele-specific expression, HLA/KIR typing, compound-heterozygote analysis). HiFi's combination of length AND ~Q30 accuracy made *phased diploid* assembly routine, so the field's deliverable shifted from one collapsed reference to **two phased haplotypes** (and onward to the **pangenome**: each HPRC node is a per-sample phased diploid assembly, Liao 2023 *Nature* 617:312).

**The killer corollary:** the thing that is wrong with a mosaic - that it is a haplotype mix - is exactly the thing none of the metrics people check can see. A switch error does not break a contig (N50 unchanged), delete a gene (BUSCO unchanged), or introduce a wrong base (QV unchanged - both alleles are real sequence, just assigned to the wrong haplotype). **Switch errors and haplotype mosaicism are structurally invisible to N50, BUSCO, and even base-level QV.** Only k-mer/trio QC - Merqury hap-mer blob plots and switch/hamming error against parental k-mers - can see them. "Reported N50 and BUSCO but no switch/hamming or hap-mer plot" means the phasing was never validated. Subtle trap: the HiFi-only `.bp.hap1/hap2` are *partially* phased (locally phased, switch errors between blocks) - "hap1" in a filename does NOT certify global phasing; the phasing *data* supplied (trio or Hi-C) does.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| hifiasm | Cheng 2021 *Nat Methods* | phased string-graph HiFi assembler; built-in purging | the default for diploid HiFi; HiFi-only / Hi-C / trio / `--ul` |
| hifiasm (Hi-C) | Cheng 2022 *Nat Biotechnol* | global phasing from proximity-ligation, no parents | diploid, Hi-C available, parents unavailable (the broad default) |
| yak | (Li, hifiasm suite) | parental haplotype-specific k-mer DBs for trio | feeds hifiasm `-1/-2`; trio gold-standard phasing |
| verkko | Rautiainen 2023 *Nat Biotechnol* | HiFi + ultralong-ONT graph assembler (MBG -> GraphAligner -> rukki) | T2T-grade; trio or Hi-C phasing; reference-quality |
| verkko2 | Antipov 2025 *Genome Res* | adds proximity-ligation phasing into the De Bruijn graph | T2T with Hi-C; ~doubled T2T-scaffold yield |
| HiCanu | Nurk 2020 *Genome Res* | HiFi Canu fork; segmental dups/satellites/allelic variants | legacy/SD-focused; superseded by hifiasm for routine diploid |
| meryl/Merqury | Rhie 2020 *Genome Biol* | hap-mer DBs; k-mer QV, completeness, switch/hamming, blob plots | the only QC that sees phasing (-> assembly-qc) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Sample type/heterozygosity unknown | profile first: k-mer spectrum (GenomeScope) for genome size + heterozygosity | purge setting, n-hap, and whether to phase are all downstream of outbred-vs-inbred |
| Outbred diploid, HiFi only, quick draft | hifiasm default | primary + partially-phased `.bp.hap1/hap2`; partial phasing is the cost of no linkage data |
| Diploid, HiFi + Hi-C, no parents | hifiasm `--h1/--h2` (Cheng 2022) | global phasing from the sample itself; the pragmatic default for the non-human bestiary |
| Diploid, HiFi + both parents | hifiasm trio `-1 pat.yak -2 mat.yak` | gold standard - read-level phasing, no switch ambiguity where parents are informative |
| Inbred / doubled-haploid / mole | hifiasm `-l0` (purging OFF) | nothing to phase; default purging would DELETE real segmental duplications |
| Unbalanced hap1/hap2 size after a run | re-run with `--hom-cov` set to the k-mer peak | the size imbalance is a mis-estimated coverage alarm, not a phasing result |
| T2T-grade reference, have HiFi + ultralong ONT + trio/Hi-C | verkko (or hifiasm `--ul` for speed) | only graph + UL spanning reaches gapless centromeres; T2T is a project, not a flag |
| No ultralong ONT but want a reference | hifiasm-HiFi-only is the sensible stop | neither tool reaches T2T without UL reads to span repeats |
| One of 50 individuals in a diversity panel | phased diploid (primary may suffice) | do NOT chase T2T per-sample; match grade to question |
| Reads not yet QC'd | -> long-read-sequencing/long-read-qc | HiFi length/QV/contamination cap assembly quality |

## hifiasm Invocations (the fragile commands - run as written)

```bash
# HiFi-only (default): primary mosaic + PARTIALLY-phased hap1/hap2
hifiasm -o prefix -t 32 reads.hifi.fq.gz

# Hi-C phased (no parents) - Cheng 2022; global phasing
hifiasm -o prefix -t 32 --h1 hic_R1.fq.gz --h2 hic_R2.fq.gz reads.hifi.fq.gz

# Trio phased (gold standard) - build parental k-mer DBs first (k=31)
yak count -k31 -b37 -t16 -o pat.yak paternal_R1.fq.gz paternal_R2.fq.gz
yak count -k31 -b37 -t16 -o mat.yak maternal_R1.fq.gz maternal_R2.fq.gz
hifiasm -o prefix -t 32 -1 pat.yak -2 mat.yak reads.hifi.fq.gz

# Inbred / homozygous / mole: DISABLE purging or real duplications are deleted
hifiasm -o prefix -t 32 -l0 reads.hifi.fq.gz

# Unbalanced haplotype sizes: pin the homozygous-coverage peak read off the k-mer histogram
hifiasm -o prefix -t 32 --hom-cov 38 reads.hifi.fq.gz

# Ultralong ONT toward T2T (combine with Hi-C or trio for phasing)
hifiasm -o prefix -t 32 --ul ul_ont.fq.gz reads.hifi.fq.gz
```

**Key flags:** `-l` purge level (`0` none, `1` light, `2`/`3` aggressive; **default 3 in HiFi-only, 0 in trio**); `--h1/--h2` Hi-C R1/R2; `-1/-2` paternal/maternal yak DBs; `--ul` ultralong ONT; `--hom-cov INT` force the homozygous-coverage peak; `--n-hap INT` ploidy (default 2); `--primary` emit primary + alternate (`a_ctg`) instead of dual hap1/hap2.

## Output Filenames: the .bp. vs .dip. Convention (read the prefix, it encodes the mode)

The prefix is not cosmetic - it tells how the assembly was phased:
- **`.bp.` ("balanced phasing") = HiFi-only OR Hi-C mode.** `prefix.bp.p_ctg.gfa` (the mosaic primary), `prefix.bp.hap1.p_ctg.gfa`/`prefix.bp.hap2.p_ctg.gfa` (the two haplotypes - *partially* phased in HiFi-only, *fully* phased with Hi-C). Also `prefix.bp.r_utg.gfa`/`p_utg.gfa` (unitig graphs).
- **`.dip.` ("diploid") = trio mode.** `prefix.dip.hap1.p_ctg.gfa` (paternal/hap1), `prefix.dip.hap2.p_ctg.gfa` (maternal/hap2).
- **`--primary` mode:** `prefix.p_ctg.gfa` (primary) + `prefix.a_ctg.gfa` (alternate). The alternate is incomplete by construction (only heterozygous loci produce alt contigs).
- **Reusable binaries** `prefix.ec.bin`/`ovlp.*.bin` let a re-run with different `-l`/phasing skip error-correction.

People grep `hap1.p_ctg` and are confused when a trio run made `dip.hap1.p_ctg` and a Hi-C run made `bp.hap1.p_ctg`. Verify the actual emitted names against the installed version.

## GFA Is Not FASTA (extract S lines before anything downstream)

hifiasm and verkko emit assembly **graphs** (GFA), and downstream tools want FASTA. Contig sequences live in GFA `S` (segment) lines:

```bash
gfatools gfa2fa prefix.bp.hap1.p_ctg.gfa > hap1.fa          # preferred
awk '/^S/{print ">"$2"\n"$3}' prefix.bp.hap1.p_ctg.gfa > hap1.fa   # dependency-free fallback
```

The graph also carries bubbles/alternate paths the FASTA throws away - keep the GFA. Verkko works internally in homopolymer-compressed coordinates (its `.gfa` is HPC); its final `assembly.fasta` is in normal space.

## verkko (T2T, when the goal genuinely needs it)

T2T is a project, not a flag. Reach for verkko only when T2T completeness is the actual goal AND the data exist (deep HiFi PLUS good ultralong ONT N50 PLUS trio or Hi-C). It is slower, heavier (hundreds of CPU-hours, high RAM, Snakemake-orchestrated), and more fragile than hifiasm. CHM13 - the first T2T human - was a hydatidiform mole precisely because a mole is effectively homozygous, decoupling repeat-resolution from phasing (Nurk 2022 *Science* 376:44).

```bash
# Trio-phased T2T: build parental + CHILD meryl DBs (k=30), derive hap-mers, then run.
# --hap-kmers needs HAP-MER DBs (haplotype-specific k-mers), not raw parental count DBs.
meryl count k=30 paternal.fq.gz output pat.meryl
meryl count k=30 maternal.fq.gz output mat.meryl
meryl count k=30 child.fq.gz    output child.meryl
$MERQURY/trio/hapmers.sh mat.meryl pat.meryl child.meryl   # hapmers.sh takes maternal first; emits mat/pat.hapmer.meryl
verkko -d asm_out --hifi hifi.fq.gz --nano ul_ont.fq.gz \
  --hap-kmers pat.hapmer.meryl mat.hapmer.meryl trio       # paternal first -> haplotype1=paternal (matches hifiasm -1 pat)

# Hi-C-phased (no parents)
verkko -d asm_out --hifi hifi.fq.gz --nano ul_ont.fq.gz --hic1 hic_R1.fq.gz --hic2 hic_R2.fq.gz
```

`--hifi` accurate reads; `--nano` ultralong ONT (the *spanning* data, strongly recommended); `--hap-kmers <hap1> <hap2> trio` (argument ORDER sets which becomes haplotype1/haplotype2 - it is not inferred from biology) or `--hic1/--hic2` for phasing. Outputs `assembly.fasta`, `assembly.haplotype1.fasta`, `assembly.haplotype2.fasta`, `assembly.homopolymer-compressed.gfa`. verkko-vs-hifiasm(`--ul`): verkko closes more chromosomes automatically at higher cost; hifiasm `--ul` brings much of the spanning benefit far cheaper but needs more manual finishing. Neither reaches T2T without ultralong reads.

## Per-Method Failure Modes

### Primary assembly used as if it were a haplotype
**Trigger:** running the primary `.bp.p_ctg`/`p_ctg` into an allele-aware analysis. **Mechanism:** the primary is a per-block haplotype mosaic. **Symptom:** garbage phased variant / ASE / HLA results; "best" metrics. **Fix:** use trio/Hi-C-phased hap1/hap2; never the primary for allele-aware work.

### Partial phasing mistaken for full phasing
**Trigger:** treating HiFi-only `.bp.hap1/hap2` as trio-grade haplotypes. **Mechanism:** no linkage data -> locally phased with inter-block switch errors. **Symptom:** phasing breaks at block boundaries; high switch rate vs trio. **Fix:** add Hi-C (`--h1/--h2`) or trio (`-1/-2`) for global phasing; validate with hap-mers.

### Default purging on an inbred/homozygous sample
**Trigger:** running hifiasm at default `-l3` on an inbred line, doubled-haploid, or mole. **Mechanism:** real segmental duplications/paralogs look like duplicate haplotigs in a homozygous genome and get purged. **Symptom:** assembly shrinks below true genome size; real duplications collapsed. **Fix:** `-l0` (purging off) for low-heterozygosity samples.

### Unbalanced hap1/hap2 size read as biology
**Trigger:** one haplotype much larger than the other. **Mechanism:** hifiasm mis-estimated the homozygous-coverage peak (bimodal coverage, odd ploidy, contamination), so it over-/under-purged. **Symptom:** size ratio far from 1. **Fix:** read the peak off the k-mer/coverage histogram (GenomeScope) and set `--hom-cov`.

### Phasing reported without phasing QC
**Trigger:** N50 + BUSCO + QV reported, no switch/hamming or hap-mer blob plot. **Mechanism:** those metrics are structurally blind to mosaicism. **Symptom:** a "great" assembly that is a phasing disaster. **Fix:** Merqury hap-mer blob plot + switch/hamming vs trio/Hi-C hap-mers (-> assembly-qc). No hap-mers (no trio/Hi-C) = phasing unvalidated by construction.

### verkko recommended reflexively / without ultralong reads
**Trigger:** reaching for verkko (or `--ul`) without ultralong ONT or without a T2T goal. **Mechanism:** UL N50 - not HiFi coverage - is the binding T2T constraint; 60x HiFi cannot span a 3 Mb satellite array. **Symptom:** months of compute, no T2T gain. **Fix:** confirm UL N50 (a meaningful fraction >100 kb) and a real T2T need first; else hifiasm-HiFi-only.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| HiFi coverage >=13x per haplotype; ~30-40x diploid sweet spot | hifiasm FAQ + community | below ~13x/hap het bubbles fragment; past ~40x diminishing returns and worse false-dup |
| HiFi >40x is not free quality | community convention | extra coverage can confuse the `--hom-cov` estimate and inflate false duplications |
| Ultralong ONT N50: the longer the better, meaningful fraction >100 kb (T2T targets >100 kb-1 Mb) | T2T practice | UL value is entirely about *spanning* centromeres/satellites; short "long" reads add little |
| `-l` purge: default 3 (HiFi-only), 0 (trio), set 0 for inbred | hifiasm man page | aggressive purging cleans outbred het-dups but deletes real sequence in homozygous samples |
| k=31 (yak/hifiasm trio) vs k=30 (verkko/Merqury hap-mers) | tool conventions | mixing k silently corrupts hap-mer matching - match each tool |
| Assembly QV: HiFi diploid routinely Q40-Q50; T2T reference ~Q60-Q70+ | Merqury practice | base QV is a SEPARATE axis from phasing accuracy (switch/hamming) |
| Do NOT reflexively polish HiFi | over-polishing risk | HiFi is already ~Q30+ at read level; polishing an accurate assembly often lowers QV (-> assembly-polishing) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Downstream tool rejects hifiasm output | GFA, not FASTA | extract `S` lines (`gfatools gfa2fa` or awk) |
| Assembly below expected genome size, real duplications gone | default purging on an inbred/mole sample | re-run with `-l0` |
| hap1 and hap2 wildly different sizes | mis-estimated homozygous-coverage peak | set `--hom-cov` from the k-mer histogram |
| "hap1" file but phasing breaks across blocks | HiFi-only partial phasing | add Hi-C/trio for global phasing |
| Allele-aware analysis (ASE/phased variants) looks wrong | used the primary mosaic | use phased hap1/hap2, not primary |
| QV drops after polishing the HiFi assembly | over-polishing an already-accurate assembly | stop; HiFi rarely needs short-read polish |
| verkko run never reaches T2T | no/short ultralong ONT | UL N50 is the binding constraint; add long UL or stop at hifiasm |
| Trio run produced `dip.*`, scripts expected `bp.*` | mode encodes the prefix | match filename prefix to phasing mode |

## References

- Cheng H, Concepcion GT, Feng X, Zhang H, Li H. 2021. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods* 18:170-175.
- Cheng H, Jarvis ED, Fedrigo O, et al. 2022. Haplotype-resolved assembly of diploid genomes without parental data. *Nat Biotechnol* 40:1332-1335.
- Nurk S, Walenz BP, Rhie A, et al. 2020. HiCanu: accurate assembly of segmental duplications, satellites, and allelic variants from high-fidelity long reads. *Genome Res* 30:1291-1305.
- Koren S, Rhie A, Walenz BP, et al. 2018. De novo assembly of haplotype-resolved genomes with trio binning. *Nat Biotechnol* 36:1174-1182.
- Rautiainen M, Nurk S, Walenz BP, et al. 2023. Telomere-to-telomere assembly of diploid chromosomes with Verkko. *Nat Biotechnol* 41:1474-1482.
- Antipov D, Rautiainen M, Nurk S, et al. 2025. Verkko2 integrates proximity-ligation data with long-read de Bruijn graphs. *Genome Res* 35:1583-1594.
- Nurk S, Koren S, Rhie A, et al. 2022. The complete sequence of a human genome (T2T-CHM13). *Science* 376:44-53.
- Liao WW, Asri M, Ebler J, et al. 2023. A draft human pangenome reference. *Nature* 617:312-324.
- Rhie A, Walenz BP, Koren S, Phillippy AM. 2020. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.

## Related Skills

- genome-profiling - Decide outbred/inbred/ploidy and the homozygous-coverage peak before setting purge level and `--hom-cov`
- assembly-qc - Merqury hap-mer blob plots and switch/hamming are the only QC that sees phasing
- assembly-polishing - HiFi is already accurate; deciding whether to polish at all
- scaffolding - Hi-C used here for phasing; chromosome-scale scaffolding (YaHS/SALSA2) is the adjacent step
- contamination-detection - Screen contigs for foreign sequence after assembly
- long-read-sequencing/long-read-qc - HiFi read length/QV/contamination before assembly
- workflows/genome-assembly-pipeline - End-to-end profile -> assemble -> phase -> QC -> scaffold
