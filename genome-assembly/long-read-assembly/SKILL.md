---
name: bio-genome-assembly-long-read-assembly
description: Assembles genomes de novo from noisy long reads (Oxford Nanopore R9/R10/Dorado, PacBio CLR) with Flye (repeat graph), Canu (correct-trim-assemble OLC), NextDenovo, Shasta, Raven, wtdbg2, or miniasm, and reconciles bacterial assemblies into a consensus with Trycycler/Autocycler. Covers matching the input flag to the basecaller era (--nano-hq vs --nano-raw), why a raw long-read assembly is contiguous but low-QV and not finished until polished, haplotig false-duplication and purge_dups, coverage and read-N50 as non-substitutable inputs, and mid-read adapter de-chimerization. Use when assembling a bacterial or eukaryotic genome from ONT or PacBio noisy reads, choosing a long-read assembler, or diagnosing an over-collapsed or duplicated assembly. For PacBio HiFi use hifi-assembly instead.
tool_type: cli
primary_tool: Flye
---

## Version Compatibility

Reference examples tested with: Flye 2.9+, Canu 2.2+, NextDenovo 2.5+, Shasta 0.11+, Raven 1.8+, wtdbg2 2.5+, miniasm 0.3+, minimap2 2.26+, purge_dups 1.2+, Porechop_ABI 0.5+, Trycycler 0.5+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Two behaviors are version-dependent and load-bearing: Flye `--genome-size` became optional (auto-estimated) and `--scaffold` flipped OFF by default in recent releases - confirm against the installed Flye. Shasta `--config` names are dated and chemistry-specific (e.g. `Nanopore-R10-Fast-Nov2022`) - list with `shasta --command listConfigurations` rather than hard-coding one. Canu read-type flag spellings and `correctedErrorRate` defaults changed across major versions. If a command errors, introspect the installed tool and adapt rather than retrying.

# Long-Read Assembly (Noisy Reads)

**"Assemble a genome from Nanopore/PacBio-CLR long reads"** -> Build a contiguous de novo assembly whose input mode matches the basecaller error regime, then hand off to polishing because the raw consensus is contiguous but not yet accurate.
- CLI: `flye --nano-hq reads.fq.gz --out-dir out -t 16` (modern ONT R10/Dorado), `canu -p asm -d out genomeSize=4.6m -nanopore reads.fq.gz` (thorough), `wtdbg2 -x ont -g 4.6m -i reads.fq.gz -fo asm && wtpoa-cns -i asm.ctg.lay.gz -fo asm.fa` (fast draft)

## The Single Most Important Modern Insight -- The Basecaller Era Picks the Flag, and the Assembly Is Only Half-Done Until Polished

Two load-bearing facts no assembler README states:

1. **The basecaller era dictates the input flag, and a mismatch silently wrecks the assembly while the report looks finished.** Noisy long reads are not one error regime: ONT R9.4.1+Guppy is ~5-10% error, ONT R10.4.1+Dorado SUP is Q20+ (~1-2%), PacBio CLR is ~10-15%. Each assembler has separate modes tuned to each (Flye `--nano-raw` / `--nano-hq` / `--pacbio-raw`). The dangerous direction is **telling the assembler the reads are noisier than they are** - feeding R10/Dorado data to `--nano-raw` makes Flye treat real repeat-copy and allele differences as noise and **over-collapse** them. The result is *fewer contigs and a HIGHER N50* as the assembly gets *worse* - the most dangerous failure mode because the headline metric improved and nothing crashes. The chemistry/kit/basecaller model is an assembly *parameter*, not optional metadata; if it is unknown, the flag cannot be chosen and the assembly is uninterpretable - find out, do not guess.

2. **The bottleneck flipped from contiguity to consensus accuracy.** A raw noisy-long-read assembly emits a FASTA with a spectacular N50 that *looks* finished, but per-base accuracy is often Q20-Q30 (one error every ~100-1000 bp). The dominant ONT error is the **indel**, especially in homopolymers, and a single indel frameshifts a protein - so a contiguous assembly can have every gene model broken. Contiguity and correctness are orthogonal axes; N50 is blind to QV. The assembler is the **halfway point**: the deliverable is a *polished* assembly with a *measured* QV (Merqury), not a high-N50 FASTA. Flye does one internal polishing round and stops - that is not "polished." Hand off to assembly-polishing and assembly-qc.

## Tool Taxonomy

| Tool | Citation | Paradigm | When |
|------|----------|----------|------|
| Flye | Kolmogorov 2019 *Nat Biotechnol* | repeat graph (disjointigs -> explicit repeat structure) | the fast general-purpose default; bacteria -> eukaryote |
| Canu | Koren 2017 *Genome Res* | correct -> trim -> assemble (OLC, MHAP) | maximum single-assembler quality; slow, grid-oriented |
| NextDenovo | Hu 2024 *Genome Biol* | correct (NextCorrect) -> string-graph (NextGraph) | large/repetitive plant and animal genomes; contiguity-first |
| Shasta | Shafin 2020 *Nat Biotechnol* | run-length-encoded marker graph | ONT human-scale, speed/cost critical |
| Raven | Vaser & Sikic 2021 *Nat Comput Sci* | OLC string graph, near parameter-free | fast simple draft; common Trycycler input |
| wtdbg2 | Ruan & Li 2020 *Nat Methods* | fuzzy de Bruijn graph (+ mandatory wtpoa-cns) | fastest, lowest RAM; lowest accuracy -> polish hard |
| miniasm | Li 2016 *Bioinformatics* | string-graph layout, NO consensus | layout demo / pipeline component only; output = raw read error |
| Trycycler / Autocycler | Wick 2021 *Genome Biol* / Wick 2025 | consensus of multiple independent assemblies | bacterial reliability; catches single-assembler structural errors |
| purge_dups | Guan 2020 *Bioinformatics* | read-depth + self-alignment | remove haplotig false-duplication from a diploid primary |

PacBio CLR (`--pacbio-raw`) is legacy - superseded by HiFi for all new PacBio work; treat CLR support as maintenance for archival data, do not recommend generating new CLR.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| ONT R10.4.1 / Dorado HAC-SUP, any genome | Flye `--nano-hq` | modern default mode; matches the Q20+ error regime |
| ONT R9.4.1 SUP (Guppy5+/Dorado) | Flye `--nano-hq --read-error 0.05` | hq mode with the error floor raised to R9-SUP |
| ONT R9.4.1 legacy fast/HAC | Flye `--nano-raw` | the genuinely-noisy mode; do not use on R10 |
| PacBio CLR (archival) | Flye `--pacbio-raw` or Canu `-pacbio` | legacy noisy mode; polish hard (arrow/GCpp) |
| Bacterial isolate, want a *correct* finished genome | Trycycler (interactive) / Autocycler (automated) | consensus across assemblers; fixes structural errors polishing can't |
| Large repetitive plant/animal, contiguity-first | NextDenovo | memory-efficient, top contiguity on big genomes |
| ONT human-scale, speed-critical | Shasta `--config <era-matched>` | RLE marker graph; human genomes in days |
| Quick draft / compute is the bottleneck | wtdbg2 or Raven | fastest; accept lower accuracy then polish |
| PacBio HiFi (Q30+, CCS) | -> hifi-assembly | hifiasm phased haplotypes; wrong tool here |
| After assembling (always) | -> assembly-polishing then assembly-qc | raw consensus is low-QV; not finished until polished + QV-measured |
| Reads not yet QC'd / unknown chemistry | -> long-read-sequencing/long-read-qc, long-read-sequencing/basecalling | garbage-in caps the assembly; basecaller model sets the flag |
| Genome size / coverage unknown | GenomeScope2 on accurate short reads (not raw ONT) | k-mer histograms from noisy reads inflate unique k-mers |

## Flye (the default)

```bash
flye --nano-hq reads.fq.gz --out-dir out -t 16                              # ONT R10 / Dorado SUP
flye --nano-hq reads.fq.gz --read-error 0.05 --out-dir out -t 16            # ONT R9 SUP
flye --nano-raw reads.fq.gz --out-dir out -t 16                             # legacy ONT R9 fast/HAC
flye --pacbio-raw reads.fq.gz --out-dir out -t 16                           # PacBio CLR (legacy)
flye --nano-hq reads.fq.gz --genome-size 3g --asm-coverage 40 -o out -t 32  # large genome: use longest 40x for initial assembly
```

`--genome-size` is optional in recent Flye (auto-estimated) but **required when paired with `--asm-coverage`**, which downsamples to the longest N-coverage of reads for the initial disjointig step (cuts runtime/RAM on deep large-genome data; the rest are still used). `--iterations` defaults to 1 polishing round (`0` to skip); `--keep-haplotypes` retains alt bubble paths for diploid awareness; `--meta` is metaFlye for uneven-coverage communities (-> metagenome-assembly). Output: `assembly.fasta`, `assembly_info.txt` (per-contig length/coverage/circularity), `assembly_graph.gfa`.

## Canu (thorough), wtdbg2 / Raven (fast), miniasm (layout only)

```bash
canu -p asm -d out genomeSize=4.6m -nanopore reads.fq.gz useGrid=false maxThreads=16   # correct->trim->assemble
wtdbg2 -x ont -g 4.6m -t 16 -i reads.fq.gz -fo asm && wtpoa-cns -t 16 -i asm.ctg.lay.gz -fo asm.ctg.fa  # consensus step is MANDATORY
raven -t 16 reads.fq.gz > asm.fasta                                                     # near parameter-free
```

Canu's master meta-parameter is `correctedErrorRate` (max expected difference between two corrected reads): defaults ~0.144 (Nanopore) / ~0.045 (PacBio); raise for heterozygosity/divergence, lower for clean high-coverage data. Read-type flags `-nanopore` / `-pacbio` / `-pacbio-hifi` (there is no `-nanopore-hifi`; high-accuracy ONT still uses `-nanopore`); `useGrid=false` forces a single machine. wtdbg2 needs the separate `wtpoa-cns` consensus call - the `-x` preset (`ont`/`sq`/`rs`/`ccs`) is set FIRST, and `-L` discards short reads (default 5000 for `ont`/`sq`). miniasm does **no consensus at all** - its output carries the full raw read error rate and is unusable until polished, so use it only as a fast layout inside a polished pipeline.

## Haplotig False-Duplication and purge_dups

A diploid assembly from noisy reads either collapses heterozygous loci or emits both haplotypes as separate primary contigs (**haplotig duplication**). The diagnostic trifecta: (1) **assembly size 1.5-2x the expected genome size**; (2) a **bimodal read-depth histogram with a half-coverage peak** (haplotigs split reads between two copies); (3) **inflated BUSCO-Duplicated**. Fix with purge_dups (read-depth + self-alignment):

```bash
minimap2 -xmap-ont asm.fa reads.fq.gz | gzip > aln.paf.gz
pbcstat aln.paf.gz && calcuts PB.stat > cutoffs          # INSPECT the coverage histogram before trusting auto-cuts
split_fa asm.fa > asm.split && minimap2 -xasm5 -DP asm.split asm.split | gzip > self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov self.paf.gz > dups.bed
get_seqs -e dups.bed asm.fa                              # -> purged.fa (primary) + hap.fa (haplotigs)
```

purge_dups **cannot tell a haplotig from a real recent segmental duplication/paralog** - both look like similar sequence at fractional depth. On a genome with known recent WGD or high SD content (many plants), over-purging deletes real genes; eyeball the histogram and validate against a related assembly. True haplotype-*resolved* assembly is a HiFi capability - do not promise phasing from noisy reads (-> hifi-assembly).

## Pre-Assembly: Mid-Read Adapter De-Chimerization

```bash
porechop_abi -abi -i reads.fq.gz -o trimmed.fq.gz -t 16   # ab initio adapter detection; SPLITS internal-adapter chimeras
```

ONT occasionally sequences two molecules as one read with an **internal adapter** - an untrimmed chimera becomes a structural mis-join (a layout error polishing cannot fix). Porechop_ABI detects adapters ab initio (no fixed DB, which matters because kit adapter sequences change) and splits chimeric reads. Use it, not the original Porechop (unmaintained since 2018, frozen adapter DB). PacBio CLR handles adapter/scrap removal upstream on the instrument, so this is an ONT-specific concern.

## Per-Method Failure Modes

### Flag noisier than reads actually are
**Trigger:** `--nano-raw` (or low-error-rate omission) on R10/Dorado-SUP data. **Mechanism:** assembler treats real repeat-copy and allele differences as noise and over-collapses. **Symptom:** fewer contigs, HIGHER N50, lost repeats/SVs - looks better. **Fix:** match the flag to the basecaller era (`--nano-hq` for R10); record pore+kit+model before assembling.

### "It assembled, so it's done"
**Trigger:** shipping the raw assembler FASTA on its high N50. **Mechanism:** raw consensus is Q20-Q30; indels frameshift genes. **Symptom:** broken gene models, failed variant calling, BLAST misses present genes. **Fix:** polish (long-read Medaka/Racon, model-matched), then *measure* QV with Merqury; an assembly without a stated QV is a draft.

### Trusting polishing to fix a structural error
**Trigger:** expecting polish to repair a mis-resolved repeat, inversion, collapsed tandem array, or dropped plasmid. **Mechanism:** polishers correct per-base consensus (QV) only; they never change contig layout. **Symptom:** a Q50 assembly that is still structurally wrong. **Fix:** catch layout errors with read-back coverage uniformity, multi-assembler consensus (Trycycler/Autocycler), or Hi-C - not QV.

### Haplotig duplication read as completeness
**Trigger:** "my genome is bigger than expected - more complete!" **Mechanism:** both haplotypes kept as primary contigs. **Symptom:** size 1.5-2x expected, half-coverage depth peak, inflated BUSCO-Duplicated. **Fix:** purge_dups (inspect cutoffs); do not over-purge real SDs.

### High coverage from a short-N50 library
**Trigger:** 200x of reads with a 5 kb N50, expecting contiguity. **Mechanism:** reads physically cannot span long repeats; depth buys consensus, not spanning. **Symptom:** fragmented assembly that more depth never fixes. **Fix:** coverage and read-N50 are non-substitutable; spend effort on read length (extraction, size selection) or ultra-long ONT.

### Over-filtering by length
**Trigger:** aggressive Filtlong/chopper length filter to "keep the best reads." **Mechanism:** the longest reads are often not the highest quality; the filter discards the long-but-lower-Q reads that span hard repeats. **Symptom:** lost contiguity. **Fix:** filter conservatively; spanning reads are precious.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Coverage ~30-60x | field convention (approx) | below ~20-30x consensus too thin (per-base accuracy collapses); beyond ~60x more of the *same* reads adds little contiguity and slows overlap |
| `--asm-coverage 40` with `--genome-size` | Flye usage | downsample to longest 40x for the initial assembly on deep large genomes |
| Canu `correctedErrorRate` ~0.144 ONT / ~0.045 PacBio | Canu 2.2 reference | master knob; raise for heterozygosity, lower for clean high-coverage |
| Assembly size 1.5-2x expected = red flag | diploid norm | haplotig false-duplication; confirm with half-coverage peak + BUSCO-Duplicated -> purge_dups |
| Merqury QV40 (~1 error/10 kb), Q50 reference-grade | Rhie 2020 *Genome Biol* | the QV stop signal; report it, never an N50 alone |
| R9 raw ~Q15-17, R10 simplex Q20+, duplex Q30+ | Wick-blog / ONT (approx) | sets the input flag and whether short-read polishing is even needed |
| wtdbg2 `-x ont` `-L` default 5000 | wtdbg2 preset | silently discards reads shorter than 5 kb |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Fewer contigs, higher N50, but lost variation | `--nano-raw` on R10/Dorado data (over-collapse) | use `--nano-hq`; match the basecaller era |
| Gene prediction frameshifts everywhere | unpolished low-QV assembly | polish then measure QV (Merqury); -> assembly-polishing |
| Assembly ~2x expected size, BUSCO-Duplicated high | uncollapsed haplotigs | purge_dups; inspect coverage cutoffs first |
| wtdbg2 output is empty/short | forgot the `wtpoa-cns` consensus step | run `wtpoa-cns` on `.ctg.lay.gz` |
| miniasm assembly full of errors | miniasm does no consensus | polish (Racon/Medaka) or use a consensus assembler |
| Chimeric contigs / structural mis-joins | internal-adapter chimeric reads | de-chimerize with Porechop_ABI before assembly |
| Canu runs for days, huge RAM | normal for the correction stage on large genomes | `useGrid=true` on a cluster, or use Flye |

## References

- Kolmogorov M, Yuan J, Lin Y, Pevzner PA. 2019. Assembly of long, error-prone reads using repeat graphs (Flye). *Nat Biotechnol* 37:540-546.
- Koren S, Walenz BP, Berlin K, Miller JR, Bergman NH, Phillippy AM. 2017. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. *Genome Res* 27:722-736.
- Hu J, Wang Z, Sun Z, et al. 2024. NextDenovo: an efficient error correction and accurate assembly tool for noisy long reads. *Genome Biol* 25:107.
- Shafin K, Pesout T, Lorig-Roach R, et al. 2020. Nanopore sequencing and the Shasta toolkit enable efficient de novo assembly of eleven human genomes. *Nat Biotechnol* 38:1044-1053.
- Vaser R, Sikic M. 2021. Time- and memory-efficient genome assembly with Raven. *Nat Comput Sci* 1:332-336.
- Ruan J, Li H. 2020. Fast and accurate long-read assembly with wtdbg2. *Nat Methods* 17:155-158.
- Li H. 2016. Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences. *Bioinformatics* 32:2103-2110.
- Guan D, McCarthy SA, Wood J, Howe K, Wang Y, Durbin R. 2020. Identifying and removing haplotypic duplication in primary genome assemblies (purge_dups). *Bioinformatics* 36:2896-2898.
- Wick RR, Judd LM, Cerdeira LT, et al. 2021. Trycycler: consensus long-read assemblies for bacterial genomes. *Genome Biol* 22:266.
- Rhie A, Walenz BP, Koren S, Phillippy AM. 2020. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.

## Related Skills

- genome-profiling - k-mer genome-size estimate feeds Flye `--genome-size` and the haplotig sanity check
- hifi-assembly - PacBio HiFi (Q30+) phased haplotype-resolved assembly with hifiasm; the right tool for accurate reads
- assembly-polishing - Polishes the contiguous-but-low-QV contigs this skill produces; assembly is not finished without it
- assembly-qc - QUAST/BUSCO/Merqury; QV is the deliverable, N50 alone is a vanity metric
- long-read-sequencing/long-read-qc - Read length/quality QC and conservative filtering before assembly
- long-read-sequencing/basecalling - Basecaller era and model that determine the assembler input flag
- read-qc/contamination-screening - Screen reads for host/vector contamination before assembling
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
