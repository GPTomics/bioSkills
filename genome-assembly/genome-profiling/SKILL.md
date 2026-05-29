---
name: bio-genome-assembly-genome-profiling
description: Profiles a genome from raw reads BEFORE assembly with a k-mer spectrum (KMC or Jellyfish histogram), then models it with GenomeScope2 to estimate genome size, heterozygosity, repeat content, and ploidy, and Smudgeplot to infer ploidy from heterozygous k-mer pairs (diploid AB vs triploid AAB vs tetraploid AABB). Covers choosing k via Merqury best_k.sh, the k-mer-coverage vs sequencing-coverage confusion, reading het/repeat/contamination/organelle peaks, why noisy ONT must not be used for counting, and how the estimate becomes the NG50 denominator, the Flye -g value, the hifiasm --hom-cov/purge setting, and the 1.5-2x-too-big haplotig sanity check. Use when starting any de novo assembly, deciding whether short reads can work, estimating genome size for an unknown organism, diagnosing ploidy, or sanity-checking an assembly's size against expectation.
tool_type: cli
primary_tool: GenomeScope2
---

## Version Compatibility

Reference examples tested with: GenomeScope2 2.0+, KMC 3.2+, Jellyfish 2.3+, meryl 1.4+ (Merqury 1.3+), Smudgeplot 0.2.5+, KAT 2.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Smudgeplot changed its backend: classic releases (0.2.x) run `smudgeplot.py hetkmers`/`smudgeplot.py plot` on a KMC dump; newer releases (0.3+) run `smudgeplot hetmers`/`smudgeplot all` on a FastK database. Confirm which interface is installed (`smudgeplot.py --version` or `smudgeplot --version`) before scripting. GenomeScope2 ships as `genomescope2` and as `genomescope.R`; both take the same flags. meryl `best_k.sh` lives in the Merqury install. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Genome Profiling

**"What am I about to assemble, and what should I expect?"** -> Build a k-mer spectrum from raw accurate reads and model it to estimate genome size, heterozygosity, repeat content, and ploidy, which set every downstream assembly expectation and parameter.
- CLI: `kmc -k21 ... reads kmc_db tmp/ && kmc_tools transform kmc_db histogram reads.histo` (count) then `genomescope2 -i reads.histo -o gs_out -k 21 -p 2` (model); `smudgeplot.py hetkmers` / `smudgeplot.py plot` (ploidy)

## The Single Most Important Modern Insight -- Profile the Genome Before Assembling, or the Assembler Guesses For Itself

A k-mer spectrum built from the raw reads -- reference-free, before a single contig exists -- estimates genome size, heterozygosity, repeat content, and ploidy, and those four numbers set the rest of the project: the NG50 denominator (assembly-qc), Flye's `-g`/`--genome-size`, hifiasm's `--hom-cov`/purge level, and whether short reads can produce the assembly being asked for at all. Skipping it leaves the assembler to infer the homozygous-coverage peak itself; when it mis-estimates (heterozygosity, odd ploidy, contamination, a bimodal coverage spectrum) it over- or under-purges, and that is exactly why people publish genomes inflated 1.5-2x by uncollapsed haplotigs. Three load-bearing moves:

1. **The estimate is the denominator AND the sanity check.** NG50 is N50 against the *expected* genome size, not the assembly size, so without the estimate NG50 cannot be reported honestly. After assembly, the same number is the haplotig test: an assembly 1.5-2x the GenomeScope size with high BUSCO-Duplicated is uncollapsed haplotypes, not a big genome -- purge before believing the size (see hifi-assembly, assembly-qc).
2. **The spectrum reads as a diagnostic, not just a size estimate.** A single homozygous peak ~= haploid/inbred; a distinct half-coverage (AB) peak *left* of the homozygous (AA) peak is heterozygous diploid, and the het peak's area gives the heterozygosity rate. A heavy high-multiplicity tail is repeat content; a spike at very high multiplicity is organelle or high-copy repeat; a left-shoulder near multiplicity 1 is sequencing error or a low-coverage contaminant. Read it before trusting any number from it.
3. **Count with accurate reads only.** GenomeScope's negative-binomial mixture model assumes errors are rare and Poisson-like. Noisy ONT (raw/HAC) injects so many unique error k-mers that the error shoulder swamps the real peaks and the fit fails. Count from Illumina or PacBio HiFi; use the ONT reads to *assemble*, never to *profile*.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| KMC | Kokot 2017 *Bioinformatics* | disk-based k-mer counter -> histogram | default counter; frugal on RAM, fast on large genomes |
| Jellyfish | Marcais & Kingsford 2011 *Bioinformatics* | in-memory k-mer counter -> histogram | alternative counter; classic GenomeScope input |
| meryl | Rhie 2020 *Genome Biol* | k-mer counter + `best_k.sh` | derives k from genome size; feeds Merqury QV downstream |
| GenomeScope2 | Ranallo-Benavidez 2020 *Nat Commun* | model: size, het, repeat, ploidy from the histogram | the profiling model for diploids and polyploids |
| Smudgeplot | Ranallo-Benavidez 2020 *Nat Commun* | ploidy from het k-mer-pair coverage ratios | unknown ploidy; cross-check GenomeScope's `-p` |
| KAT | Mapleson 2017 *Bioinformatics* | spectra plots, reads-vs-assembly k-mer comparison | contamination triage; post-assembly completeness/spectra-cn |

## Decision Tree by Scenario

| Scenario | What the profile indicates | Path |
|----------|---------------------------|------|
| Single sharp peak, size ~ expected, low het | haploid/inbred or clonal; clean | proceed -> short-read-assembly or hifi-assembly with default purge |
| Two peaks (AB at ~half AA) | heterozygous diploid; het rate from AB area | HiFi+phasing best; if short reads only, expect haplotigs -> short-read-assembly (Platanus) |
| High het + Illumina only | short-read DBG will fragment and inflate 1.5-2x | get long reads, or plan purge_dups; do not report inflated size |
| GenomeScope `-p 2` fits poorly; Smudgeplot shows AAB/AABB | triploid/tetraploid; ploidy not 2 | re-run GenomeScope with correct `-p`; -> hifi-assembly haplotype expectations |
| Coverage peak < ~15-20x | too shallow for a stable model fit | sequence more, or treat size/het as lower-confidence; -> read-qc/quality-reports |
| Extra peak at odd multiplicity, or bimodal spectrum | contamination / organelle / mixed sample | KAT spectra triage; screen reads -> read-qc/quality-reports before assembling |
| Only noisy ONT available | cannot profile reliably from error-dominated spectrum | assemble first, then estimate size from the assembly + Merqury -> assembly-qc |
| Need the genome-size denominator for QC | GenomeScope haploid length | feed as NG50 expected size and Flye `-g` -> assembly-qc, long-read-assembly |

## Choosing k

**Goal:** Pick a k that is large enough that most k-mers are genomically unique but small enough to keep per-k-mer depth high.

**Approach:** Derive k from the expected genome size with Merqury's `best_k.sh` (formula `k = log4(G(1-p)/p)`, default tolerable collision rate `p=0.001`); a vertebrate-scale ~3 Gb genome returns k=21 (the de facto GenomeScope default), ~1 Gb returns k=20, and a ~12 Mb yeast returns k=17.

```bash
sh $MERQURY/best_k.sh 3100000000          # ~3.1 Gb vertebrate -> k=21 (the common GenomeScope default)
sh $MERQURY/best_k.sh 1000000000          # ~1 Gb -> k=20
sh $MERQURY/best_k.sh 12000000            # ~12 Mb yeast -> k=17
```

Too small a k saturates: nearly every k-mer recurs across the genome by chance, the unique/repeat peaks merge, and size is overestimated. Too large a k loses depth (per-k-mer coverage is `c*(L-k+1)/L`, so it drops as k rises) and the peaks blur into the error shoulder. k=21 is the long-standing default at vertebrate (~3 Gb) scale because it sits in this window; smaller genomes want smaller k (best_k.sh returns ~20 at 1 Gb, ~17 at 12 Mb). Use the SAME k for counting and for the `-k` passed to GenomeScope2.

## Counting K-mers and Running GenomeScope2

```bash
# KMC: -ci1 keeps singletons (the error shoulder GenomeScope models), -cs10000 caps the histogram tail
kmc -k21 -t16 -m64 -ci1 -cs10000 @fastq_list.txt kmc_db tmp/
kmc_tools transform kmc_db histogram reads.histo -cx10000

# GenomeScope2: -p 2 = diploid; raise for known/Smudgeplot-suggested polyploidy
genomescope2 -i reads.histo -o gs_out -k 21 -p 2

# Jellyfish alternative (-C canonical k-mers, mandatory for unstranded WGS)
jellyfish count -C -m 21 -s 4G -t 16 reads_*.fastq -o reads.jf
jellyfish histo -t 16 reads.jf > reads_jf.histo
```

GenomeScope2 writes a model fit (`model.txt`), the linear/log spectrum plots, and a summary with haploid genome length, heterozygosity, and a "% unique" (inverse of repeat content). `-cs10000`/`-cx10000` cap the histogram so a single organelle/repeat spike at multiplicity 100000+ does not dominate the file; raise the cap only for very high-coverage data.

## K-mer Coverage vs Sequencing Coverage (the confusion that breaks the fit)

GenomeScope reports `kmercov` (lambda) -- the mean coverage *per k-mer*, the x-position of the homozygous peak. This is NOT per-base sequencing coverage. They differ by the `(L-k+1)/L` factor: at read length L=150 and k=21, a k-mer is covered `130/150 ~= 0.87x` as often as a base, so a 50x-sequenced genome shows a homozygous peak near multiplicity 43, not 50. Passing the sequencing coverage where GenomeScope expects the k-mer-coverage peak (e.g. as an `-l` initial guess), or reading lambda back as sequencing depth, mis-scales the model and the size estimate. Read the peak off the plot; let GenomeScope estimate lambda unless the fit fails, then seed `-l` with the observed peak position.

## Ploidy with Smudgeplot

```bash
# Classic (0.2.x, KMC backend): pick L/U coverage cutoffs, dump het k-mer pairs, plot
L=$(smudgeplot.py cutoff kmc_db.histo L)     # ~0.5x the haploid peak (errors below)
U=$(smudgeplot.py cutoff kmc_db.histo U)     # ~8.5x the haploid peak (repeats above)
kmc_tools transform kmc_db -ci"$L" -cx"$U" dump -s kmc_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmer_pairs < kmc_L"$L"_U"$U".dump
smudgeplot.py plot kmer_pairs_coverages.tsv -o sample
```

Smudgeplot infers ploidy from the coverage RATIO of heterozygous k-mer pairs, independent of GenomeScope's model: a diploid AB pair sits at `CovB/(CovA+CovB) ~= 0.5`, a triploid AAB near 0.33, a tetraploid AABB shows smudges at 0.25 and 0.5. When Smudgeplot's inferred ploidy disagrees with the `-p` that fit GenomeScope best, that disagreement is itself the finding -- re-examine for polyploidy, aneuploidy, or contamination rather than forcing one answer.

## Per-Method Failure Modes

### Profiling from noisy ONT
**Trigger:** counting k-mers from raw/HAC Nanopore reads. **Mechanism:** ~5-10% per-base error makes nearly every error k-mer unique, burying the real peaks under the multiplicity-1 shoulder. **Symptom:** GenomeScope fit fails or returns absurd size/het. **Fix:** count from Illumina or HiFi; profile from the assembly + Merqury if only ONT exists.

### Too-low coverage
**Trigger:** homozygous peak below ~15-20x. **Mechanism:** the error shoulder and the real peak overlap; the negative-binomial mixture cannot separate them. **Symptom:** wide confidence intervals, unstable size, no clean peaks. **Fix:** sequence more depth, or report size/het as low-confidence.

### Reading lambda as sequencing depth
**Trigger:** treating GenomeScope `kmercov` as per-base coverage, or seeding `-l` with sequencing depth. **Mechanism:** k-mer coverage = depth * `(L-k+1)/L` < depth. **Symptom:** mis-scaled model, size off by the `(L-k+1)/L` factor. **Fix:** read lambda off the peak; do not substitute sequencing coverage.

### k too small (saturation)
**Trigger:** k=15-17 on a large/repeat-rich genome. **Mechanism:** most k-mers recur by chance; unique and repeat components merge. **Symptom:** inflated size, blurred peaks. **Fix:** use `best_k.sh` for the expected size (k=21 at vertebrate ~3 Gb scale; smaller for smaller genomes).

### Contamination/organelle distorting the spectrum
**Trigger:** an unexplained extra peak or a spike at very high multiplicity. **Mechanism:** a second organism's k-mers add their own peak; organelle/high-copy DNA spikes the tail. **Symptom:** GenomeScope size too large or a multi-modal spectrum. **Fix:** KAT spectra triage; screen reads (-> read-qc/quality-reports) before assembling.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| k from best_k.sh (21 at ~3 Gb, 20 at ~1 Gb, 17 at ~12 Mb) | best_k.sh `k=log4(G(1-p)/p)`, p=0.001 | balances uniqueness vs per-k-mer depth; k=21 is the common GenomeScope default at vertebrate scale |
| Homozygous peak >= ~15-20x | GenomeScope model stability | below this the error shoulder and real peak cannot be separated |
| AB peak at ~0.5x the AA peak | diploid k-mer theory | heterozygous k-mers are half-covered (one haplotype); het rate from AB area |
| Smudgeplot CovB/(CovA+CovB): 0.5 / 0.33 / 0.25 | k-mer-pair coverage ratio | AB diploid / AAB triploid / AABB tetraploid |
| Assembly size 1.5-2x GenomeScope estimate | haplotig norm | uncollapsed heterozygous haplotypes; purge before reporting size |
| KMC `-cs`/GenomeScope `-cx` cap ~10000 | histogram tail control | prevents an organelle/repeat spike from dominating the file |
| Smudgeplot L ~0.5x, U ~8.5x haploid peak | Smudgeplot guidance | excludes errors (below) and high-copy repeats (above) from pairing |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| GenomeScope size far larger than expected | k too small, contamination, or organelle spike | raise k via best_k.sh; cap the tail; screen reads |
| Model fit fails / "cannot fit" | too-low coverage or ONT error k-mers | more depth; count from Illumina/HiFi, not noisy ONT |
| Histogram peak at multiplicity ~43 for "50x" data | k-mer coverage = depth * (L-k+1)/L, not depth | expected; do not pass sequencing depth as `-l` |
| GenomeScope `-p 2` fits poorly | sample is not diploid | run Smudgeplot; re-run with the correct `-p` |
| Jellyfish histogram looks halved | counted without `-C` (canonical) | recount with `-C`; WGS is unstranded |
| Assembly 1.8x the profiled size | uncollapsed haplotigs | purge_dups / hifiasm purge; report the haploid size |

## References

- Ranallo-Benavidez TR, Jaron KS, Schatz MC. 2020. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nat Commun* 11:1432.
- Vurture GW, et al. 2017. GenomeScope: fast reference-free genome profiling from short reads. *Bioinformatics* 33:2202-2204.
- Kokot M, Dlugosz M, Deorowicz S. 2017. KMC 3: counting and manipulating k-mer statistics. *Bioinformatics* 33:2759-2761.
- Marcais G, Kingsford C. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers (Jellyfish). *Bioinformatics* 27:764-770.
- Rhie A, Walenz BP, Koren S, Phillippy AM. 2020. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.
- Mapleson D, et al. 2017. KAT: a K-mer analysis toolkit to quality control NGS datasets and genome assemblies. *Bioinformatics* 33:574-576.

## Related Skills

- short-read-assembly - High het from the profile predicts short-read fragmentation and haplotig inflation
- hifi-assembly - The profiled size/het sets hifiasm --hom-cov and the purge level
- long-read-assembly - The profiled genome size feeds Flye -g and ONT/PacBio expectations
- assembly-qc - The estimate is the NG50 denominator and the 1.5-2x haplotig sanity check
- read-qc/quality-reports - QC and contamination-screen reads before profiling and assembling
- workflows/genome-assembly-pipeline - Profiling is the first step of the end-to-end assembly workflow
