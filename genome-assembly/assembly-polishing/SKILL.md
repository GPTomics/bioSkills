---
name: bio-genome-assembly-assembly-polishing
description: Decides whether and how to polish a draft genome assembly to raise consensus accuracy (QV) with read-type-matched tools - Racon and medaka (ONT consensus), dorado polish, Polypolish and pypolca (Illumina, repeat-aware), Pilon (legacy short-read), NextPolish/NextPolish2, Hapo-G (haplotype-aware), ntEdit, and DeepPolisher/PEPPER-Margin-DeepVariant for human. Covers the do-not-polish-HiFi rule, the medaka basecaller-model footgun, held-out Merqury QV as the only honest stop signal, and the haplotype-collapse trap. Use when correcting homopolymer indels or residual SNPs in a long-read assembly, deciding if a HiFi assembly needs polishing, or choosing an ONT vs hybrid vs short-read polishing chain.
tool_type: cli
primary_tool: Pilon
---

## Version Compatibility

Reference examples tested with: Racon 1.5+, medaka 2.0+, minimap2 2.26+, bwa 0.7.17+, samtools 1.19+, Polypolish 0.6+, pypolca 0.3+, Pilon 1.24+, Merqury 1.3+, meryl 1.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

The **medaka consensus model** matters more than the medaka binary version: the model must match the basecaller + pore chemistry + caller mode + version (`-m r1041_e82_400bps_sup_v5.0.0`); a mismatched model silently degrades the consensus. List options with `medaka tools list_models`. medaka's status is a moving target - ONT now steers toward `dorado polish` and has deprecated medaka's diploid-variant workflow in favour of Clair3; verify current guidance and the `medaka_consensus` wrapper vs newer subcommands. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Assembly Polishing

**"Polish my genome assembly"** -> Decide whether base-level consensus correction is warranted, then apply a read-type-matched polisher and measure the gain with held-out k-mers - or, for an already-accurate assembly, decline.
- CLI: `racon reads.fq aln.sam draft.fa` then `medaka_consensus -i reads.fq -d draft.fa -o out -m <model>` (ONT); `polypolish polish draft.fa f1.sam f2.sam` + `pypolca run` (hybrid/short); `merqury.sh reads.meryl asm.fa out` (measure QV before/after)

## The Single Most Important Modern Insight -- Polishing Is Being Engineered Out, and Reflexive Polishing Does Net Harm

Polishing exists to fix the characteristic error of *noisy* long reads: indels in homopolymers and low-complexity tracts (the pore/RT stutters on `AAAAAA`) plus residual substitutions. As read accuracy rose (PacBio HiFi ~Q40-50 reads, ONT R10.4.1 + Dorado `sup`/duplex), the assembly consensus is already near-perfect, and a mapping-based polisher's read-pileup step injects more errors than it removes. Two load-bearing rules dominate:

1. **Do NOT polish HiFi assemblies by default.** A hifiasm/HiCanu HiFi assembly starts at ~Q40+; mapping-based polishers (Racon, Pilon, GCpp) routinely *lower* QV and introduce haplotype-switch errors. Only polish HiFi with a tool explicitly built not to overcorrect (NextPolish2, DeepPolisher) and only if Merqury QV says there is a real deficit. For HiFi the burden of proof is on polishing, not against it.
2. **Measure the gain with reference-free, HELD-OUT Merqury k-mers - never the reads polished with, never BUSCO.** A polisher's literal objective is to maximize agreement with its input reads, so scoring it with those same reads is circular and always looks like an improvement. The honest stop signal is a Merqury QV *plateau*, not a fixed iteration count - and polishing can drive QV *down* silently. Polish with one platform, evaluate with another (polish with ONT, measure with Illumina k-mers).

Polishing is transitional infrastructure: indispensable in the CLR/early-ONT error era, shrinking toward irrelevance as raw accuracy climbs upstream (better basecalling, hifiasm's own consensus). The skill's job is to say *when not to polish* as firmly as how to.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| Racon | Vaser 2017 *Genome Res* | POA read-to-assembly consensus; fast workhorse | first-pass ONT/CLR self-polish (1-4 rounds); needs minimap2 SAM/PAF |
| medaka | nanoporetech/medaka (software) | ONT neural-network consensus | one pass *after* Racon; model MUST match basecaller |
| dorado polish | ONT/dorado (software) | modern ONT consensus, owns basecalling context | the emerging medaka replacement; verify current status |
| Polypolish | Wick & Holt 2022 *PLoS Comput Biol* | Illumina polish using ALL alignments (`bwa mem -a`) | repeat-rich long-read assemblies; best-in-class short-read polish |
| pypolca | Bouras 2024 *Microb Genom* (POLCA: Zimin 2020) | fast Illumina SNP/indel correction | pair with Polypolish (modern bacterial practice) |
| Pilon | Walker 2014 *PLoS ONE* | all-in-one SNP/indel/local short-read fix | legacy; small genomes only; mismaps in repeats, ~1 GB heap/Mb |
| NextPolish / NextPolish2 | Hu 2020 *Bioinformatics* / Hu 2024 *GPB* | iterative LR+SR polish / HiFi repeat+phase-aware | NextPolish2 is the HiFi-safe successor |
| Hapo-G | Aury & Istace 2021 *NARGAB* | haplotype-aware short-read polish | heterozygous diploids; preserves both alleles |
| ntEdit | Warren 2019 *Bioinformatics* | Bloom-filter k-mer polish, no mapping | very large (>3 Gb) genomes; scales where alignment can't |
| DeepPolisher | Mastoras 2025 *Genome Res* | transformer correction (PHARAOH ONT phasing) | HiFi/T2T human; halves errors without overcorrection |
| PEPPER-Margin-DeepVariant | Shafin 2021 *Nat Methods* | haplotype-aware ONT consensus via variant calling | human ONT-only polishing |
| DeepConsensus | Baid 2023 *Nat Biotechnol* | improves CCS *reads* UPSTREAM of assembly | NOT a polisher - operates before assembly (common category error) |

## Decision Tree by Scenario

| Scenario (data type) | Recommended | Why |
|----------------------|-------------|-----|
| PacBio HiFi only | **Do NOT polish by default** -> check Merqury QV first | already ~Q40-50; Racon/Pilon lower QV + collapse haplotypes |
| HiFi with a real Merqury QV deficit | NextPolish2 or DeepPolisher | repeat/phase-aware; won't overcorrect het sites |
| ONT-only, noisy (R9 / R10 hac) | Racon (1-4 rounds) -> medaka (1 pass), or dorado polish | canonical ONT chain; neural consensus on the ONT error model |
| ONT-only, human/diploid | PEPPER-Margin-DeepVariant or DeepPolisher | haplotype-aware; preserves both alleles |
| Long-read draft + Illumina (hybrid) | long-read pass, then Polypolish + pypolca | SR fixes residual homopolymer indels Racon/medaka missed |
| Short-read-only assembly (SPAdes) | usually no polishing needed | Illumina is already ~Q40+; structural sins survive anyway |
| HiFi + Illumina | use Illumina k-mers for *evaluation*, not correction | Merqury hybrid QV, not a polishing pass |
| Heterozygous / repeat-rich genome | Polypolish (`-a`), Hapo-G, NextPolish2 - or none | non-haplotype-aware polishers erase het / homogenize paralogs |
| Reads not yet QC'd / wrong basecaller | -> read-qc/quality-reports, -> long-read-sequencing/long-read-alignment | garbage-in or model-mismatch makes polishing worse than skipping |
| Complaint is "fragmented" not "low QV" | -> not polishing (scaffolding/gap-filling) | polishing fixes bases, not contiguity |
| Map reads before polishing | short -> read-alignment/bwa-alignment; long -> long-read-sequencing/long-read-alignment | polishers consume a BAM/SAM/PAF, not raw reads |

## The Canonical ONT Chain: Racon -> medaka

**Racon** does the bulk cheap consensus from the long reads themselves; it is a standalone consensus module and does NOT map - the SAM/PAF is supplied. Re-map every round (1-4 rounds; gains decay fast).

```bash
minimap2 -t 16 -ax map-ont draft.fasta reads.fq.gz > aln.sam   # map-pb for PacBio CLR
racon -t 16 reads.fq.gz aln.sam draft.fasta > racon1.fasta
# re-map reads to racon1.fasta and repeat for round 2... measure QV each round, stop at plateau
```

**medaka** applies an ONT-trained neural model for the final consensus - exactly ONE pass after the Racon rounds (it is not an iterate-many-times tool; running medaka twice is a tell). For medaka mechanics and model tables see long-read-sequencing/medaka-polishing; this skill owns the *strategy*.

```bash
medaka_consensus -i reads.fq.gz -d racon_final.fasta -o medaka_out -t 16 \
  -m r1041_e82_400bps_sup_v5.0.0          # model MUST match basecaller+chemistry+caller+version
```

Recent basecallers embed the model in the FASTQ so medaka auto-selects; if the data was basecalled with an old/unknown caller, pick manually from `medaka tools list_models`, and treat a stale/deprecated model name as a reason to rebasecall rather than proceed.

## Hybrid / Short-Read Polishing: Polypolish + pypolca

**Polypolish** aligns short reads to *all* locations (`bwa mem -a`) so it can disambiguate which repeat copy a read belongs to instead of forcing one placement - this is *how* it beats Pilon in repeats and why it almost never introduces errors.

```bash
bwa index draft.fasta
bwa mem -t 16 -a draft.fasta reads_1.fq.gz > aln_1.sam     # -a = ALL alignments (the whole point)
bwa mem -t 16 -a draft.fasta reads_2.fq.gz > aln_2.sam
polypolish filter --in1 aln_1.sam --in2 aln_2.sam --out1 filt_1.sam --out2 filt_2.sam
polypolish polish draft.fasta filt_1.sam filt_2.sam > polypolish.fasta   # --careful (v0.6+) for low depth
pypolca run -a polypolish.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz -o pypolca_out -t 16 --careful
```

Bouras 2024 depth-tiered bacterial recommendation: depth <5x -> Polypolish `--careful` alone; 5-25x -> Polypolish `--careful` + pypolca `--careful`; >25x -> Polypolish (default) + pypolca `--careful`. Modern bacterial best practice is Polypolish + pypolca, NOT Pilon.

## Pilon (Legacy Short-Read)

```bash
bwa mem -t 16 draft.fasta r1.fq r2.fq | samtools sort -o frags.bam
samtools index frags.bam
java -Xmx16G -jar pilon.jar --genome draft.fasta --frags frags.bam --output pilon --fix all --changes
```

`--fix` modes: `snps`, `indels`, `bases` (=snps+indels), `gaps`, `local`, `all` (**default**), `none`. Pilon needs a sorted+indexed BAM (route mapping to read-alignment/bwa-alignment). It is legacy: it uses best-placement alignments so it *mismaps in repeats* (miscorrects toward paralogs), and its ~1 GB heap per Mb of genome OOMs on large eukaryotes - both reasons Polypolish/pypolca displaced it.

## Measuring the Gain: Held-Out Merqury QV

**Goal:** Decide whether a polish actually helped, without fooling yourself.

**Approach:** Build a meryl k-mer DB from an *independent / different-platform* read set (k from Merqury's `best_k.sh`, not hardcoded), then run Merqury on the pre- and post-polish assemblies and compare QV. A polish that does not raise QV did not help; one that lowers it must be reverted.

```bash
K=$(sh $MERQURY/best_k.sh 5000000 | tail -n1 | awk '{print int($1+0.5)}')   # K from genome size; round float->int
meryl count k=$K output reads.meryl illumina_reads.fq.gz               # eval reads != polishing reads
merqury.sh reads.meryl draft.fasta    qv_before                        # QV of the input
merqury.sh reads.meryl polished.fasta qv_after                         # QV after polishing
```

Use a held-out set or a *different platform* than was polished with (polish with ONT, evaluate with Illumina k-mers); the CHM13 effort triangulated with both HiFi and Illumina k-mers (Mc Cartney 2022). Do NOT use BUSCO as the polishing metric - it measures gene-space completeness and barely moves with the homopolymer-indel QV that polishing changes, so a flat BUSCO masks a silent QV drop. A per-gene internal-stop / frameshift count is a useful secondary readout. See assembly-qc for the full QV/spectra-cn workflow.

## Per-Method Failure Modes

### medaka run with the wrong basecaller model
**Trigger:** specifying or defaulting to a model that does not match the actual basecaller+chemistry+version. **Mechanism:** the network applies corrections calibrated for errors that aren't there and misses the ones that are. **Symptom:** medaka completes with no warning, but Merqury QV is *lower* than the input. **Fix:** let medaka auto-detect from the FASTQ; if manual, derive the model from the real caller and confirm in `medaka tools list_models`; treat a stale model as a reason to rebasecall.

### Polishing a HiFi assembly with a mapping-based polisher
**Trigger:** running Racon/Pilon/GCpp on a hifiasm assembly. **Mechanism:** at het sites the pileup mixes both true alleles; the polisher overwrites toward the majority and homogenizes near-identical paralogs. **Symptom:** lost heterozygosity, haplotype-switch errors, QV unchanged or down. **Fix:** don't; if Merqury shows a real deficit, use NextPolish2/DeepPolisher only.

### Validating with the polishing reads (circularity)
**Trigger:** measuring QV with the same read set used to polish. **Mechanism:** the polisher already maximized agreement with those reads. **Symptom:** every polish "improves QV," monotonically. **Fix:** evaluate with held-out / different-platform k-mers.

### Over-polishing past the plateau
**Trigger:** "a few more rounds to be safe." **Mechanism:** once resolvable errors are fixed, extra rounds flip correct bases at het/repeat sites. **Symptom:** "N changes made" keeps reporting while QV oscillates or falls. **Fix:** track Merqury QV per round; stop at plateau. Treat a large change count on an accurate assembly as a risk signal, not success.

### Pilon mismapping in repeats
**Trigger:** Pilon on a repeat-rich genome. **Mechanism:** best-placement alignment forces a read onto one repeat copy; Pilon corrects that copy toward the wrong one. **Symptom:** repeat copies homogenized, paralog differences erased. **Fix:** Polypolish (`bwa mem -a`, uses all alignments).

### Treating DeepConsensus as a polisher
**Trigger:** running DeepConsensus on the assembly. **Mechanism:** DeepConsensus improves CCS *reads* upstream, before assembly. **Symptom:** tool expects subreads/CCS, not a FASTA. **Fix:** run it in the read-prep stage; for assembly polishing use a post-assembly tool.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| HiFi assembly ~Q40-50 already | HiFi read accuracy | starting point too high for mapping-based polishing to help |
| Racon 1-4 rounds, stop at QV plateau | Vaser 2017 + practice | gains decay after ~1-2; extra rounds flip correct bases |
| medaka exactly 1 pass after Racon | medaka design | trained-model pass, not an iterative tool |
| Polypolish `--careful` <5x; +pypolca 5-25x; default >25x | Bouras 2024 *Microb Genom* | depth-tiered to avoid false-positive repeat edits |
| QV40 (~1 err/10 kb) "reference-grade"; Q50 (~1/100 kb) modern aspiration; CHM13 ~Q73 | field convention; Mc Cartney 2022 | QV = -10*log10(error rate); report measured QV + method |
| Merqury k from `best_k.sh` (k=21 human-scale) | Rhie 2020 *Genome Biol* | wrong k silently degrades QV/completeness |
| Pilon ~1 GB heap per Mb of genome | Walker 2014 | OOM risk on large eukaryotes; set `-Xmx` |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| QV drops after polishing a HiFi assembly | over-polishing an already-accurate assembly | stop; HiFi rarely needs short-read polish |
| medaka output worse than input, no warning | wrong/stale basecaller model | auto-detect or match model exactly; else rebasecall |
| Every polishing round "improves" QV | measured with the polishing reads (circular) | evaluate with held-out / different-platform k-mers |
| Lost heterozygosity / switch errors | non-haplotype-aware polisher on a diploid | Hapo-G / NextPolish2 / Polypolish `-a`, or skip |
| Repeat copies homogenized | Pilon best-placement mismapping | Polypolish with `bwa mem -a` |
| Pilon OOM / crash on large genome | ~1 GB/Mb heap blowup | raise `-Xmx`; prefer Polypolish/ntEdit |
| Polishing didn't fix fragmentation | wrong operation | fragmentation is scaffolding/gap-filling, not polishing |

## References

- Vaser R, Sović I, Nagarajan N, Šikić M. 2017. Fast and accurate de novo genome assembly from long uncorrected reads (Racon). *Genome Res* 27:737-746.
- Walker BJ, et al. 2014. Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement. *PLoS ONE* 9:e112963.
- Wick RR, Holt KE. 2022. Polypolish: short-read polishing of long-read bacterial genome assemblies. *PLoS Comput Biol* 18:e1009802.
- Zimin AV, Salzberg SL. 2020. The genome polishing tool POLCA makes fast and accurate corrections in genome assemblies. *PLoS Comput Biol* 16:e1007981.
- Bouras G, et al. 2024. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. *Microb Genom* 10:001254.
- Hu J, et al. 2020. NextPolish: a fast and efficient genome polishing tool for long-read assembly. *Bioinformatics* 36:2253-2255.
- Hu J, et al. 2024. NextPolish2: a repeat-aware polishing tool for genome assemblies with HiFi long reads. *Genomics Proteomics Bioinformatics* 22:qzad009.
- Aury JM, Istace B. 2021. Hapo-G, haplotype-aware polishing of genome assemblies with accurate reads. *NAR Genom Bioinform* 3:lqab034.
- Warren RL, et al. 2019. ntEdit: scalable genome sequence polishing. *Bioinformatics* 35:4430-4432.
- Baid G, et al. 2023. DeepConsensus improves the accuracy of sequences with a gap-aware sequence transformer. *Nat Biotechnol* 41:232-238.
- Shafin K, et al. 2021. Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads. *Nat Methods* 18:1322-1332.
- Mastoras M, et al. 2025. Highly accurate assembly polishing with DeepPolisher. *Genome Res* 35:1595-1608.
- Rhie A, et al. 2020. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.
- Mc Cartney AM, et al. 2022. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. *Nat Methods* 19:687-695.

## Related Skills

- long-read-assembly - Produces the contiguous-but-error-prone contigs this skill polishes
- short-read-assembly - Source of Illumina reads for hybrid/short-read polishing
- hifi-assembly - HiFi assemblies that usually should NOT be polished
- assembly-qc - Merqury QV before/after is the polishing stop signal
- read-alignment/bwa-alignment - Map short reads to the draft for Polypolish/Pilon
- long-read-sequencing/long-read-alignment - minimap2 mapping of long reads for Racon
- long-read-sequencing/medaka-polishing - medaka mechanics and model tables; this skill owns the strategy
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
