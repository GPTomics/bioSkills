---
name: bio-genome-assembly-assembly-qc
description: Evaluates genome assembly quality across the three orthogonal axes - contiguity (QUAST auN/NG50/NGx, not bare N50), completeness (BUSCO/compleasm gene-space plus Merqury k-mer completeness), and correctness (reference-free Merqury QV, Inspector/CRAQ structural errors, asmgene false-duplication/collapse). Covers why N50 is the most-gamed metric, why QV measured on the polishing reads is circular, distinguishing uncollapsed haplotigs from real WGD, and the EBP/VGP 6.C.Q40 standard. Use when judging whether an assembly is good enough to annotate or publish, comparing assemblers, diagnosing a fragmented or duplicated assembly, or assessing a phased diploid assembly.
tool_type: cli
primary_tool: QUAST
---

## Version Compatibility

Reference examples tested with: QUAST 5.2+, BUSCO 5.5+ (and 6.x for odb12 lineages), compleasm 0.2.6+, Merqury 1.3+, meryl 1.4+, minimap2 2.26+, Inspector 1.2+, CRAQ 1.0+, merfin 1.0+, GenomeScope2 2.0+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

Results depend on inputs that outlive the binary version - record them:
- BUSCO/compleasm depend on the **lineage dataset** and OrthoDB generation. `_odb10` (BUSCO 5) and `_odb12` (BUSCO 6 default) gene sets are not comparable across the version boundary; a 99% on the shallow `eukaryota_odb10` (~255 genes) is a different claim from 99% on a deep clade set (~5,500+).
- Merqury QV/completeness depend on the **k-mer size** (from `best_k.sh <genome_size>`, not hardcoded) and the **read set** used for the k-mer DB (use accurate reads; see the circularity warning below).
- NG50/NGx/auNG depend on the **expected genome-size estimate** (GenomeScope2 / flow cytometry / a congener).

If code throws an error, introspect the installed tool and adapt rather than retrying.

# Assembly QC

**"Is my genome assembly any good?"** -> Measure all three orthogonal axes - contiguity, completeness, correctness - with reference-free methods, because no single number (least of all N50) is quality.
- CLI: `quast.py asm.fa --large --eukaryote -o out` (contiguity + reference-based structure), `busco -i asm.fa -m genome -l <lineage>` or `compleasm run -a asm.fa -l <lineage>` (gene completeness), `merqury.sh reads.meryl asm.fa out` (reference-free QV + k-mer completeness), `inspector.py -c asm.fa -r reads.fq` (reference-free structural errors)

## The Single Most Important Modern Insight -- Quality Is Three Orthogonal Axes; N50 Is the Most-Gamed One

Assembly quality is **three genuinely orthogonal axes - contiguity, completeness, correctness - and a single number on any one is not quality.** The axes do not predict each other, and the diagnostic failure modes prove it:

- **Contiguous + wrong:** a single-contig "chromosome" that is three chromosomes misjoined. Perfect N50, catastrophic correctness. Only Hi-C / a same-species reference / read-discordance catches it.
- **Complete + shredded:** BUSCO 99%, but repeats collapsed, segmental duplications merged, intergenic space wrong. BUSCO is gene-space-only and cannot see it.
- **Accurate + incomplete:** QV60 over the 92% that assembled, with the hard 8% (centromeres, rDNA, satellites) simply absent. QV is silent about what is not there.

The field's historical sin is reporting **contiguity alone** because it is cheapest to compute and easiest to game. **N50 is the most-gamed metric in genomics:** it rises when sequence is thrown away (N50 is computed on what survives), when misjoins are *not* broken (a misjoined contig is a long contig), and when haplotigs are retained. A bigger N50 is louder, not better. Three load-bearing moves:

1. **Report auN/NGx, not bare N50.** auN = the area under the Nx curve = length-weighted mean contig length; it integrates the whole curve and is continuous where N50 jumps discontinuously (the small-L50 / T2T regime). NGx/auNG normalize to the **expected genome size**, coupling contiguity to completeness (an assembly that drops half the genome gets a great N50 but a terrible NG50). Always report contig AND scaffold N50 - if scaffold >> contig, the contiguity is glue (Ns), not sequence.
2. **Default to reference-free.** For a *novel* genome there is no trusted reference; QUAST against a divergent relative reports real inversions/SVs as "misassemblies" and real SNPs as "mismatches". Use Merqury QV (accuracy) + Inspector/CRAQ (structure) + asmgene (false dup/collapse). QUAST is the special case "I have a same-organism reference," not the default.
3. **Report a Merqury QV - and never compute it on the polishing reads.** QV is the reference-free accuracy standard reviewers now demand; an assembly paper with no QV is a red flag. But QV from the same reads used for polishing is **circular** - the polisher already made the assembly agree with those reads, so the QV measures convergence, not correctness. Build the k-mer DB from accurate, ideally independent reads (HiFi/Illumina, not noisy ONT).

## Tool Taxonomy

| Tool | Citation | Axis / Role | When |
|------|----------|-------------|------|
| QUAST / calN50 | Gurevich 2013 *Bioinformatics*; auN = Li blog (no journal) | contiguity (auN/NGx, N50/L50) + reference-based structure (NA50, misassemblies) | always for contiguity; structure only vs a same-organism reference |
| BUSCO | Manni 2021 *Mol Biol Evol*; Simão 2015 *Bioinformatics* | gene-space completeness (C/S/D/F/M) | universal; conservative on good genomes |
| compleasm | Huang & Li 2023 *Bioinformatics* | gene-space completeness, miniprot-based | faster + more sensitive; default on HiFi/T2T-era genomes |
| Merqury | Rhie 2020 *Genome Biol* | reference-free QV + k-mer completeness + spectra-cn + phasing | always; the accuracy gold standard |
| merfin | Formenti 2022 *Nat Methods* | multiplicity-corrected QV / polishing | refine QV biased by k-mer multiplicity |
| Inspector | Chen 2021 *Genome Biol* | reference-free structural + base errors (long reads) | novel genomes; can also correct |
| CRAQ | Li 2023 *Nat Commun* | reference-free structural/regional errors (clipped alignments) | novel genomes; flags misjoins to split |
| asmgene | Li (minimap2, no separate journal) | gene collapse / false duplication | high-quality genomes where BUSCO saturates |
| GenomeScope2 | Ranallo-Benavidez 2020 *Nat Commun* | genome size / het / repeat % from k-mers | the size estimate NG50/auNG/spectra-cn need (-> genome-profiling) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Novel genome, no trusted reference | Merqury QV + k-mer completeness + BUSCO/compleasm + auN/NGx + Inspector/CRAQ | reference-free triad; QUAST structure is uninterpretable here |
| Same-species (near-isogenic) reference available | add QUAST `--large --eukaryote -r ref.fa` (NA50, misassemblies) | reference-based structure is trustworthy only vs the same organism |
| High-quality HiFi/T2T-era assembly, BUSCO looks low | compleasm | BUSCO under-reports good genomes (its predictor, not the assembly, misses genes) |
| Contiguity claim must resist gaming | auN/auNG via `calN50.js -L <size>` | N50 is a single unstable order-statistic and is gameable |
| High BUSCO-Duplicated, size > expected | spectra-cn + asmgene + GenomeScope2 size -> purge_dups | distinguish uncollapsed haplotigs (purge) from real WGD (keep) |
| Phased diploid / trio assembly | Merqury hap-mers: switch/hamming error + blob plot | phasing accuracy is the extra axis |
| Need genome size / het before NG50 | -> genome-profiling (GenomeScope2) | NG/auN/spectra-cn all need a size estimate |
| Reads not yet QC'd | -> read-qc/quality-reports | garbage-in caps assembly quality |
| Bacterial isolate / MAG completeness+contamination | -> contamination-detection (CheckM2/GUNC/MIMAG) | marker-gene completeness/contamination is a different problem |

## Contiguity -- auN/NGx (not bare N50)

```bash
k8 calN50.js -L <genome_size> asm.fa       # N50/L50 + NG50/NGx + auN/auNG; -L sets genome size for NG/auNG (ships with minimap2)
quast.py asm.fa --large --eukaryote -t 16 -o quast_out   # N50/L50, GC, # contigs; NG50 only with -r or --est-ref-size; structure only if -r given
```

`--large` implies `--eukaryote --min-contig 3000 --min-alignment 500 --extensive-mis-size 7000`. Report contig AND scaffold N50; a scaffold N50 far above the contig N50 means the contiguity is scaffolding Ns, and every gap is a join hypothesis that could be a misassembly. NA50 (QUAST, contigs broken at misassemblies) far below N50 means the contiguity is partly fictional.

## Completeness -- BUSCO / compleasm + Merqury k-mer completeness

```bash
busco -i asm.fa -m genome -l vertebrata_odb10 -o busco_out -c 16   # metaeuk predictor (default)
busco -i asm.fa -m genome --auto-lineage -o busco_out -c 16        # auto-pick if clade unknown
compleasm run -a asm.fa -l vertebrata -o compleasm_out -t 16        # miniprot-based; faster + more sensitive
```

Reported as `C:[S,D],F,M,n`. **Read C, F, and M together, never C alone:** high Fragmented with high Complete signals a contiguity/base-quality problem hidden behind the headline. Use the **deepest applicable clade dataset** (a 99% on the shallow `eukaryota_odb10` ~255-gene set is trivially easy and not comparable to a deep clade set), and record the lineage + OrthoDB generation. On a high-quality assembly, BUSCO reported ~95.7% complete where compleasm reported ~99.6% on the same human genome (Huang & Li 2023) - the missing ~4% was missing from BUSCO's *predictor*, not the genome - so prefer compleasm on good genomes. Both share gene-space blindness: they say nothing about intergenic/repeat/regulatory sequence. **Merqury k-mer completeness** scores the whole genome (reliable read k-mers found in the assembly / reliable read k-mers in the reads), catching missing sequence BUSCO cannot see; it is blind to structure (a scrambled-but-present genome scores 100%).

## Correctness -- Merqury QV (reference-free) and structural validation

**Goal:** Get a reference-free per-base accuracy (QV) plus a copy-number/false-duplication picture, then structural errors without a reference.

**Approach:** Build a meryl k-mer DB at the `best_k.sh`-derived k from accurate reads, run Merqury for QV + completeness + spectra-cn, and map raw long reads back with Inspector/CRAQ for structural errors. Refine QV with merfin where multiplicity bias matters.

```bash
best_k.sh <genome_size>                          # prints recommended k (NOT hardcoded); ~18-21 for Gbp genomes
meryl count k=21 reads.fastq output reads.meryl  # k from best_k.sh; use ACCURATE reads (HiFi/Illumina)
merqury.sh reads.meryl asm.fa out                # -> out.qv (per-scaffold + overall), out.completeness.stats, spectra-cn

inspector.py -c asm.fa -r reads.fq -o insp_out --datatype hifi -t 16   # reference-free structural + base errors
craq -g asm.fa -sms long_reads.bam -ngs short_reads.bam -o craq_out    # R-AQI/S-AQI; CRE (regional)/CSE (structural)
```

QV: with `E = K_asm-only / K_total`, per-base error `P = 1 - (1 - E)^(1/k)` and `QV = -10*log10(P)` (the `^(1/k)` converts a k-mer error rate to per-base, since one wrong base breaks k overlapping k-mers). QV40 = 1 error/10 kb (the EBP/VGP floor), QV50 strong, ~QV60 = T2T-grade (1/Mb). The **spectra-cn plot** reads completeness and false duplication in one figure: a black "missing" peak at homozygous depth = real content absent; 2-copy k-mers under the 1-copy peak = uncollapsed haplotigs; error k-mers sit far left.

## EBP/VGP standards and phased QC

The EBP minimum is **6.C.Q40**: `x.y.z` where x = log10 of contig NG50 (6 = 1 Mb), y = scaffold level (C = chromosome-scale), z = QV (40 = <1 error/10 kb). It is literally the triad turned into a label, and the bar moves - T2T pushed the achievable frontier to ~Q60/gapless, so bragging about QV40 in 2026 is hitting the floor. Match the bar to the organism (the relaxed "5" tier, >100 kb contig NG50, exists for low-input species). For phased diploid/trio assemblies, Merqury hap-mers (parental k-mers, or Hi-C) give the **switch error** (local haplotype flips within a block) and **hamming error** (global mis-assignment fraction) - report both, and read the hap-mer blob plot (cleanly phased contigs sit on one axis).

## Per-Method Failure Modes

### Leading with N50 (and stopping)
**Trigger:** reporting a single N50 as the quality verdict. **Mechanism:** N50 rises on thrown-away sequence, unbroken misjoins, and retained haplotigs; it is also a single unstable order-statistic. **Symptom:** big N50, unstated QV/completeness/NG50. **Fix:** report auN/NGx + BUSCO/compleasm + Merqury QV; treat N50-only claims as untrustworthy.

### QUAST against a divergent reference
**Trigger:** `quast.py -r congener.fa` on a novel genome. **Mechanism:** real inversions/SVs and SNPs between organism and reference are scored as "misassemblies"/"mismatches"; the count scales with divergence, not error. **Symptom:** "hundreds of misassemblies" on a correct assembly. **Fix:** use reference-free correctness (Inspector/CRAQ/Merqury); reserve QUAST structure for a same-organism reference.

### QV computed on the polishing reads
**Trigger:** QV from the exact reads used to polish. **Mechanism:** the polisher made the assembly agree with those reads by construction. **Symptom:** impressively high QV that rose after polishing with the QV reads. **Fix:** build the k-mer DB from accurate, ideally independent reads; consider merfin for multiplicity-corrected QV.

### High BUSCO-Duplicated read as success
**Trigger:** treating high BUSCO-D as "extra coverage / more complete". **Mechanism:** uncollapsed haplotigs (both alleles kept as separate primary contigs) vs real WGD vs split models. **Symptom:** D in high single digits to tens, assembly size >> GenomeScope2 estimate. **Fix:** triangulate size + spectra-cn 2-copy peak + asmgene false-dup; if no WGD -> purge_dups, then re-QC (watch for over-purge: size dropping below the estimate deletes real segmental duplications).

### QV/BUSCO accepted on an incomplete genome
**Trigger:** QV60 + BUSCO 99% taken as "done". **Mechanism:** QV is measured on what assembled; BUSCO scores the conserved easy core only. **Symptom:** high accuracy and gene-completeness while 8-15% of sequence (repeats/centromeres) is absent. **Fix:** add Merqury k-mer completeness (whole-genome) and inspect read mapping-rate/coverage uniformity.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Merqury QV >= 40 | EBP/VGP minimum (Rhie 2021) | 1 error/10 kb; QV50 strong, ~Q60 T2T-grade; report the actual value |
| QV from polishing reads | circularity trap | always biased high; use independent/accurate reads, prefer merfin |
| BUSCO Complete >= 95%, Fragmented < 5% | field convention | read F+M with C; high F = contiguity/base-quality problem behind a good C% |
| BUSCO Duplicated ~1-3% (clean haploid); >5-8% no WGD | assembly norm | uncollapsed haplotigs -> purge_dups; cross-check size + spectra-cn + asmgene |
| compleasm preferred on good genomes | Huang & Li 2023 | BUSCO under-reports (~95.7% vs ~99.6% on human) due to its predictor |
| k-mer completeness (Merqury) >= 95% | field convention | lower = sequence absent that the BUSCO gene set cannot see |
| Contig NG50 >= 1 Mb (EBP "6") | Rhie 2021 / EBP standards | the 6.C.Q40 contig bar; relaxed "5" (>100 kb) for low-input species |
| Report auN/NGx, not bare N50 | Li (auN blog) | N50 is gameable and a single unstable order-statistic |
| NG/auN need a genome-size estimate | by definition | GenomeScope2/flow cytometry; couples contiguity to completeness |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Big N50, no QV reported | leading with the most-gamed metric | add Merqury QV, k-mer completeness, auN/NGx |
| Hundreds of QUAST "misassemblies" on a novel genome | divergent reference; biology scored as error | reference-free (Inspector/CRAQ); QUAST only vs same organism |
| QV suspiciously high, rose after polishing | QV computed on the polishing reads (circular) | independent/accurate-read k-mer DB; merfin |
| Assembly ~1.5-2x expected size, high BUSCO-D | uncollapsed haplotigs (false duplication) | purge_dups; verify with spectra-cn + asmgene |
| Size drops below GenomeScope2 estimate after purging | over-purged real segmental duplications | back off purge stringency; check asmgene collapse direction |
| BUSCO low-90s on a HiFi/T2T assembly | BUSCO predictor misses present genes | re-run compleasm before concluding incompleteness |
| Scaffold N50 >> contig N50 reported as contiguity | contiguity is gap-Ns, not sequence | report contig N50 too; each gap is a join hypothesis |

## References

- Gurevich A, Saveliev V, Vyahhi N, Tesler G. 2013. QUAST: quality assessment tool for genome assemblies. *Bioinformatics* 29:1072-1075.
- Manni M, et al. 2021. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Mol Biol Evol* 38:4647-4654.
- Simão FA, et al. 2015. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. *Bioinformatics* 31:3210-3212.
- Huang N, Li H. 2023. compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics* 39:btad595.
- Rhie A, Walenz BP, Koren S, Phillippy AM. 2020. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.
- Formenti G, et al. 2022. Merfin: improved variant filtering, assembly evaluation and polishing via k-mer validation. *Nat Methods* 19:696-704.
- Chen Y, et al. 2021. Accurate long-read de novo assembly evaluation with Inspector. *Genome Biol* 22:312.
- Li K, et al. 2023. CRAQ: identification of errors in draft genome assemblies at single-nucleotide resolution for quality assessment and improvement. *Nat Commun* 14:6556.
- Ranallo-Benavidez TR, Jaron KS, Schatz MC. 2020. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nat Commun* 11:1432.
- Rhie A, et al. 2021. Towards complete and error-free genome assemblies of all vertebrate species (VGP). *Nature* 592:737-746.
- Li H. 2020. auN: a new metric to measure assembly contiguity. Blog post (lh3.github.io); auN/asmgene tools ship in minimap2/calN50 (Li 2018 *Bioinformatics* 34:3094-3100).
- Guan D, et al. 2020. Identifying and removing haplotypic duplication in primary genome assemblies (purge_dups). *Bioinformatics* 36:2896-2898.

## Related Skills

- short-read-assembly - Short-read assemblies plateau at the repeat structure; QC shows it
- long-read-assembly - Produces the contiguous-but-error-prone contigs this QC evaluates
- hifi-assembly - Phased diploid output whose false duplication and switch/hamming error this QC checks
- assembly-polishing - Merqury QV plateau is the honest stop signal; never QV on the polishing reads
- scaffolding - Validate chromosome-scale joins (Hi-C/contact map) before trusting scaffold NG50
- contamination-detection - MAG completeness/contamination (CheckM2/GUNC/MIMAG) is a separate problem
- genome-profiling - GenomeScope2 genome-size estimate that NG50/auNG/spectra-cn require
- genome-annotation/annotation-qc - Assembly-side completeness; purge haplotigs before annotating
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
