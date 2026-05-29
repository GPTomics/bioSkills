---
name: bio-genome-assembly-short-read-assembly
description: Assembles a genome de novo from Illumina short reads with SPAdes (isolate/careful/sc/meta/plasmid/rna modes), MEGAHIT (low-memory, huge datasets), Unicycler (bacterial finishing/hybrid), MaSuRCA (large hybrid), ABySS (Bloom-filter), and Platanus (heterozygous diploids), using multi-k de Bruijn graphs. Covers the repeat-resolution limit, why N50 plateaus at the genome not the depth, GenomeScope2 k-mer profiling first, the heterozygosity/haplotig trap, error-correction erasing rare alleles, GC dropout, and NG50/auN/BUSCO reporting. Use when assembling a bacterial isolate, fungal, small-eukaryotic, single-cell, or metagenome genome from Illumina reads, or when deciding whether short reads can even produce the assembly being asked for.
tool_type: cli
primary_tool: SPAdes
---

## Version Compatibility

Reference examples tested with: SPAdes 4.0+, MEGAHIT 1.2+, Unicycler 0.5+, ABySS 2.3+, GenomeScope2 2.0+, KMC 3.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

SPAdes 4.0 (June 2024) requires Python >=3.8, defaults to GFA v1.2 output with circular-path tags (`TP:Z:circular`), and is stated by the maintainers to be the last major feature release - pin the version in methods, because k auto-selection and the `--isolate`/`--careful` recommendation flipped across 3.x. MEGAHIT preset k-lists (`meta-sensitive`, `meta-large`) have changed across versions - confirm against `megahit -h` before hard-coding. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Short-Read Assembly

**"Assemble a genome from Illumina reads"** -> Build a multi-k de Bruijn graph from short reads and walk it into contigs, knowing the contiguity ceiling is set by the genome's repeat structure, not the depth.
- CLI: `spades.py --isolate -1 R1.fq.gz -2 R2.fq.gz -o out -t 16` (bacterial isolate, the genuine sweet spot)
- CLI: `megahit -1 R1.fq.gz -2 R2.fq.gz -o out -t 16` (huge/low-memory, metagenome)
- CLI: profile first (see below): `genomescope2 -i kmer.hist -o gs -k 21 -p 2` (size/het/repeat before parameters)

## The Single Most Important Modern Insight -- Short Reads Cannot Span Repeats, So Contiguity Is Capped by the Genome, Not the Depth

A short read (~150 bp) or paired insert (~300-600 bp) cannot resolve any repeat longer than itself: the de Bruijn graph enters a repeat from multiple unique flanks, traverses an identical internal path, and the assembler's only safe move is to break the contig at the repeat boundary or collapse the copies. Genomic repeats - rRNA operons (~5 kb), transposons, segmental duplications, satellites - are routinely kilobases to megabases. So a short-read assembly is **structurally fragmented at every long repeat**, and that is biology, not a software defect. Three load-bearing corollaries:

1. **N50 plateaus at the repeat structure, not the sequencing depth.** Past ~50-100x for an isolate (~50-60x for a eukaryote), more Illumina does NOT lengthen contigs - it only adds error-correction signal, and above ~150x it actively hurts by amplifying GC bias and duplicate-driven coverage distortion. The reflex "add more reads to get a better assembly" is wrong: the limit is read length.

2. **The assembler is almost never the bottleneck - the input DNA and the genome's repeat/heterozygosity structure are.** A clean, high-molecular-weight, PCR-free library at the right depth assembles beautifully with defaults; a degraded/contaminated library assembles badly with every setting and every k. The expert's first move on a bad assembly is to look at the k-mer spectrum, GC-vs-coverage, and duplication rate (GenomeScope2 profiling, see the profiling section below), not to swap assemblers or sweep k.

3. **Short reads are no longer the frontier instrument for a finished genome.** For any genome that must be complete or finished, the answer is long reads that span the repeats (`-> long-read-assembly`, `-> hifi-assembly`), with short reads relegated to polishing, k-mer-spectrum QC, and the cheap bacterial-isolate/surveillance workhorse. A skill that lets an agent attempt a de novo short-read assembly of a large, repeat-rich, heterozygous eukaryote and then blame parameters is doing harm.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| SPAdes | Bankevich 2012 *J Comput Biol* | multi-k de Bruijn assembler with mode dispatch | bacterial/fungal/small-euk isolate; the field default |
| MEGAHIT | Li 2015 *Bioinformatics* | succinct (compressed) de Bruijn, ultra-low memory | huge short-read datasets, metagenomes, memory-constrained |
| Unicycler | Wick 2017 *PLoS Comput Biol* | SPAdes wrapper + graph bridging + circularization | bacterial finishing; hybrid Illumina+long-read |
| MaSuRCA | Zimin 2013 *Bioinformatics* | super-reads + mega-reads (hybrid OLC/DBG) | large eukaryotic genomes with mixed data |
| ABySS 2.0 | Jackman 2017 *Genome Res* | Bloom-filter de Bruijn, MPI-parallel | large genomes on memory-constrained HPC |
| Platanus / Platanus-allee | Kajitani 2014 *Genome Res* | bubble-aware DBG for high heterozygosity | highly heterozygous diploids (a real niche) |
| GenomeScope2 / Smudgeplot | Ranallo-Benavidez 2020 *Nat Commun* | k-mer-histogram model: size/het/repeat/ploidy | run FIRST; sets the ceiling reference-free |
| Velvet / SOAPdenovo2 | Zerbino 2008 *Genome Res* / Luo 2012 *GigaScience* | single-k DBG, superseded | reproduce old papers only; do not start new projects |

Every mainstream short-read assembler is a de Bruijn assembler because all-vs-all overlap (OLC) is infeasible for hundreds of millions of short reads. The single most consequential parameter is k: small k over-connects (collapses repeats, chimeras), large k is more unique (resolves repeats up to length k) but demands higher coverage and is error-fragile. No single k is optimal everywhere in a genome, which is exactly why SPAdes and MEGAHIT iterate over a k-series - do NOT hand-pick a single k.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Reads not yet profiled | GenomeScope2 on a k-mer histogram (profiling section below) | the spectrum sets size/het/repeat ceiling before any parameter |
| Reads not yet QC'd | -> read-qc/quality-reports, read-qc/adapter-trimming | over-trimming creates coverage holes that fragment the graph |
| Bacterial/archaeal/viral isolate, Illumina only | SPAdes `--isolate` | multi-k DBG + paired-end repeat resolution; the genuine sweet spot |
| Small bacterial genome, want mismatch/indel cleanup | SPAdes `--careful` (small genomes ONLY) | runs MismatchCorrector; NOT for large eukaryotes; incompatible with --meta/--rna |
| Bacterial isolate, want a finished/circular assembly | Unicycler (short-only or hybrid) | graph bridging + circularization + rotation to dnaA |
| Single-cell / MDA-amplified | SPAdes `--sc` | designed for the wildly uneven MDA coverage |
| Highly heterozygous diploid, short reads | Platanus-allee, then purge | bubble-aware; but consider -> hifi-assembly instead |
| Metagenome (recover MAGs) | -> metagenome-assembly (metaSPAdes / MEGAHIT) | community co-assembly; MAG recovery, not single-genome N50 |
| Huge dataset, RAM-constrained | MEGAHIT | succinct DBG, tiny memory footprint |
| Large/repeat-rich/heterozygous eukaryote, finished genome wanted | -> hifi-assembly or -> long-read-assembly | short reads mathematically cannot span the repeats |
| Assembly built, judging quality | -> assembly-qc (NG50 + auN + BUSCO + Merqury QV) | never report N50 alone; contiguity is not correctness |

## SPAdes -- the dominant short-read assembler (the mode is a graph choice, not a cosmetic flag)

```bash
spades.py --isolate -1 R1.fq.gz -2 R2.fq.gz -o out -t 16 -m 64   # recommended default for isolates; does NOT run MismatchCorrector
spades.py --careful -1 R1.fq.gz -2 R2.fq.gz -o out               # SMALL genomes only; runs MismatchCorrector; NOT eukaryotes; not with --meta/--rna
spades.py --sc      -1 R1.fq.gz -2 R2.fq.gz -o out               # single-cell/MDA (default k 21,33,55)
spades.py --plasmid -1 R1.fq.gz -2 R2.fq.gz -o out              # plasmidSPAdes (coverage-based extraction)
spades.py -1 R1.fq.gz -2 R2.fq.gz -o out -k 21,33,55,77         # explicit multi-k LIST (odd, <= read length) only if overriding
spades.py --only-assembler -1 R1.fq.gz -2 R2.fq.gz -o out       # skip BayesHammer (reads pre-corrected, or het/pooled data)
```

Picking the wrong mode is the most common SPAdes error. `--isolate` and `--careful` solve different problems: modern SPAdes pushes `--isolate` as the isolate default (tuned for high, even coverage), while `--careful` runs BWA-based MismatchCorrector to cut mismatches/short indels on **small** genomes only - it is explicitly slow, memory-heavy, and low-benefit on medium/large eukaryotes, and is incompatible with `--meta` and `--rna`. `--isolate` and `--careful` are mutually exclusive - SPAdes aborts (`cannot specify --careful in isolate mode`) if both are given; pick one. `--meta`, `--rna`, `--bio` dispatch to genuinely different pipelines (metaSPAdes, rnaSPAdes, biosyntheticSPAdes) with different graph models - they are NOT genome assembly. BayesHammer (the built-in Illumina error-corrector) runs by default; `--only-assembler` skips it. Outputs: `contigs.fasta`, `scaffolds.fasta`, `assembly_graph.gfa`. Always report contig metrics, not just scaffold metrics (scaffold gaps are estimated N-runs, not sequence).

## MEGAHIT, Unicycler, and the others

```bash
megahit -1 R1.fq.gz -2 R2.fq.gz -o out -t 16                    # default k 21,29,39,59,79,99,119,141; min-count 2 (drops singleton error k-mers)
megahit -1 R1.fq.gz -2 R2.fq.gz --presets meta-sensitive -o out # min-count 1; verify k-list against megahit -h
unicycler -1 R1.fq.gz -2 R2.fq.gz -o out -t 16                  # short-read-only bacterial finishing
unicycler -1 R1.fq.gz -2 R2.fq.gz -l long.fq.gz -o out          # hybrid (short for accuracy + long to bridge repeats)
abyss-pe name=asm k=96 B=2G in='R1.fq.gz R2.fq.gz'              # B = Bloom-filter size (the 2.0 low-memory mode); single k
```

MEGAHIT's succinct DBG gives a tiny memory footprint at some cost in contiguity/accuracy versus metaSPAdes - the default for huge datasets and a metagenome alternative. Unicycler wraps SPAdes, cleans the graph, bridges repeats, and rotates circular replicons to `dnaA`; hybrid mode (short + ONT/PacBio) is the bacterial finishing standard but a short-only run still cannot beat the repeat limit. ABySS 2.0's Bloom-filter mode enables large-genome assembly on modest hardware but is single-k and niche today.

## Pre-Assembly Genome Profiling (run this FIRST)

```bash
kmc -k21 -t16 -m64 -ci1 -cs10000 @reads.lst kmcdb tmp/         # k=21: long enough to be mostly unique, short enough for k-mer coverage
kmc_tools transform kmcdb histogram kmer.hist -cx10000
genomescope2 -i kmer.hist -o gs_out -k 21 -p 2                 # -p ploidy; estimates size, heterozygosity, repeat content
```

GenomeScope2 fits a model to the k-mer frequency histogram and returns, reference-free and before assembly, the genome size, heterozygosity rate, repeat content, and ploidy. A single histogram peak ~ haploid/homozygous; two peaks (a het peak at ~half the homozygous-peak coverage) ~ a heterozygous diploid that will fragment and inflate a short-read assembly. This is the ceiling check: it predicts whether short reads can even produce the assembly being asked for, and it gives the expected size to report assembly total against (assembly >> estimate = un-purged haplotigs, not a big genome).

## Per-Method Failure Modes

### Adding more coverage to fix contiguity
**Trigger:** re-sequencing deeper because contigs are short. **Mechanism:** contiguity is capped by read length vs repeat length, not depth. **Symptom:** N50 unchanged past ~50-100x; cost wasted. **Fix:** accept the repeat-determined ceiling, or switch to long reads (`-> hifi-assembly`).

### Hand-picking a single large k for highest N50
**Trigger:** "we used k=127 because it gave the best N50". **Mechanism:** large k demands high effective k-mer coverage and is error-fragile; climbing k until N50 peaks is an N50-gaming search that rewards chimeric mega-contigs. **Symptom:** gappy assembly, or a suspiciously contiguous mis-joined one. **Fix:** let multi-k auto-select; if overriding, give a small->large list, never a point.

### `--careful` on a eukaryote
**Trigger:** cargo-culting `--careful` from a bacterial tutorial. **Mechanism:** MismatchCorrector is small-genome-only; on a large genome it is slow, memory-hungry, low-benefit. **Symptom:** the run stalls/OOMs for little quality gain. **Fix:** `--isolate` (no `--careful`) for isolates; defer base accuracy to dedicated polishing (`-> assembly-polishing`).

### Error-correcting heterozygous/pooled/metagenomic data
**Trigger:** running BayesHammer/BFC/Lighter on a diploid, pooled, or community sample. **Mechanism:** k-mer-spectrum correctors assume rare k-mers are errors and edit them toward the consensus - but the minor allele / rare strain IS those rare k-mers. **Symptom:** real low-frequency variation silently flattened before assembly. **Fix:** `--only-assembler` for het/pooled/meta data; correct only clonal isolates.

### Reporting the inflated heterozygous assembly size as the genome size
**Trigger:** assembly total ~1.5-2x the GenomeScope estimate, high BUSCO-Duplicated. **Mechanism:** both haplotypes assembled separately (un-purged haplotigs), not a real duplication. **Symptom:** "genome size" wrong by up to 2x. **Fix:** purge_dups / redundans, or -> hifi-assembly; report size against the profiling estimate.

### Treating a short-read bacterial assembly's contig count as a failure
**Trigger:** "why isn't my isolate one contig?" **Mechanism:** ~7 rRNA operons and IS elements each exceed the insert and force a break - the contig count is roughly a count of long repeats. **Symptom:** 30-100 contigs from a perfect run. **Fix:** accept it for typing/AMR/phylogenetics; for a finished single contig, use Unicycler hybrid or long reads.

### Chasing GC-dropout gaps with more depth
**Trigger:** re-assembling to close a gap whose flanks have normal coverage. **Mechanism:** PCR/chemistry under-represents GC-extreme regions; the bias, not the depth, is the problem. **Symptom:** an interior coverage hole at extreme GC that stays missing at 30x and 300x. **Fix:** upstream (PCR-free prep, alternative polymerase) or different chemistry, not a parameter.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Bacterial isolate coverage ~50-100x | field convention | below ~30x error/true k-mers blur and gaps proliferate; 50-100x = clean correction + strong k-mer peak |
| Above ~100-150x: diminishing/negative | field convention | does not fix repeats; amplifies GC bias; duplicates distort the coverage model and GenomeScope fit |
| Eukaryote short-read coverage ~50-60x | field convention | higher wastes money on unresolvable repeats; lower starves the graph |
| GenomeScope/profiling k = 21 | field convention | long enough to be mostly unique, short enough for adequate k-mer coverage |
| k constraint: odd, <= read length | DBG property | even k can self-loop on palindromes; k > read is impossible |
| Heterozygous assembly inflation ~1.5-2x haploid | het diploid norm | both haplotypes assembled separately; do not report as genome size |
| MEGAHIT min-count default 2 | MEGAHIT default | filters singleton (error) k-mers; presets drop to 1 for low-abundance metagenome members |
| Report NG50 + auN + BUSCO + a reference-free correctness check | reporting convention | N50 alone is gamed by dropping sequence and inflated by mis-joins (auN = sum(L_i^2)/G) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Out of memory | large/complex dataset, high coverage | lower `-m`; use MEGAHIT (succinct DBG) or ABySS Bloom-filter mode |
| Assembly ~2x expected size, high BUSCO-Duplicated | un-purged haplotigs (heterozygosity) | purge_dups/redundans; Platanus-allee; or -> hifi-assembly |
| Many contigs from a clean bacterial isolate | rRNA operons / IS elements exceed the insert (biology) | accept for draft uses; Unicycler hybrid or long reads for a finished genome |
| Extra contigs at odd coverage | contamination/index hopping, not a novel replicon | coverage-vs-GC screen (`-> contamination-detection`) |
| Poor assembly after heavy trimming | over-trimming created coverage holes | trim adapters + bad tails only; keep reads long |
| N50 high but downstream results wrong | mis-join / N50-gaming; high N50 != correct | report NG50 + auN + BUSCO; reference-free correctness (`-> assembly-qc`) |
| `--careful` rejected / errors | used with `--meta` or `--rna` | drop `--careful`; it is small-genome isolate-only |

## References

- Bankevich A, et al. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *J Comput Biol* 19:455-477.
- Compeau PEC, Pevzner PA, Tesler G. 2011. How to apply de Bruijn graphs to genome assembly. *Nat Biotechnol* 29:987-991.
- Li D, et al. 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 31:1674-1676.
- Wick RR, et al. 2017. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. *PLoS Comput Biol* 13:e1005595.
- Jackman SD, et al. 2017. ABySS 2.0: resource-efficient assembly of large genomes using a Bloom filter. *Genome Res* 27:768-777.
- Zerbino DR, Birney E. 2008. Velvet: algorithms for de novo short read assembly using de Bruijn graphs. *Genome Res* 18:821-829.
- Luo R, et al. 2012. SOAPdenovo2: an empirically improved memory-efficient short-read de novo assembler. *GigaScience* 1:18.
- Zimin AV, et al. 2013. The MaSuRCA genome assembler. *Bioinformatics* 29:2669-2677.
- Kajitani R, et al. 2014. Efficient de novo assembly of highly heterozygous genomes from whole-genome shotgun short reads (Platanus). *Genome Res* 24:1384-1395.
- Ranallo-Benavidez TR, Jaron KS, Schatz MC. 2020. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nat Commun* 11:1432.
- Lander ES, Waterman MS. 1988. Genomic mapping by fingerprinting random clones: a mathematical analysis. *Genomics* 2:231-239.

## Related Skills

- genome-profiling - Estimate genome size, heterozygosity, and ploidy from a k-mer spectrum before assembling
- assembly-qc - NG50 + auN + BUSCO + Merqury QV; never report N50 alone
- assembly-polishing - Base-accuracy polishing belongs here, not in `--careful` reflexes
- metagenome-assembly - Community co-assembly and MAG recovery (metaSPAdes/MEGAHIT)
- long-read-assembly - The real fix for repeats short reads cannot span
- hifi-assembly - Phased, repeat-spanning assembly for heterozygous/complex eukaryotes
- read-qc/quality-reports - Garbage-in caps assembly quality; QC reads before assembling
- read-qc/adapter-trimming - Light adapter/quality trimming; over-trimming fragments the graph
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
