---
name: bio-genome-annotation-repeat-annotation
description: Discovers, classifies, and masks repetitive elements and transposable elements with RepeatModeler2 (de novo family library), RepeatMasker (masking against a library), EDTA (plant/structural TEs), or EarlGrey (auto-curating wrapper), and quantifies TE expression from RNA-seq with TEtranscripts/SQuIRE. Covers de-novo-library-as-curation-project, soft-vs-hard masking, the domesticated-gene over-masking massacre, Dfam-vs-RepBase, TE classification (Class I/II, family-vs-copy), Kimura repeat landscapes, LAI, and the RNA-seq multimapping problem. Use when masking repeats before gene prediction, building a TE library for a non-model genome, or analyzing transposable-element content or expression.
tool_type: cli
primary_tool: RepeatMasker
---

## Version Compatibility

Reference examples tested with: RepeatModeler 2.0.5+, RepeatMasker 4.1.5+, EDTA 2.1+, EarlGrey 4.0+, TEtranscripts 2.2+, matplotlib 3.8+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

The **library database version matters as much as the binary**: RepeatMasker now ships with Dfam (open); RepBase has been paywalled since May 2019, so any pipeline that "requires RepBase" is a reproducibility/access hazard - record the Dfam release and library provenance. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Repeat and Transposable Element Annotation

**"Mask repeats in my genome assembly"** -> Build a de novo repeat-family library, annotate copies genome-wide, and soft-mask them as a prerequisite for gene prediction.
- CLI: `RepeatModeler -database mydb -LTRStruct` (library), `RepeatMasker -lib lib.fa -xsmall assembly.fa` (soft-mask), or `EarlGrey`/`EDTA.pl` (wrappers)

## The Single Most Important Modern Insight -- The Library Is the Experiment, and the Assembly Caps It

Two load-bearing truths the masker hides:

1. **De novo library construction is a curation research project, not a button.** A RepeatModeler2 run emits `mydb-families.fa` overnight - a *draft of a draft*: consensi are routinely 5'-truncated (L1 looks 1.5 kb when the active element is 6 kb), boundary-bled into flanking unique sequence, chimeric (two families merged), and 30-60% "Unknown" on a non-model genome. The dominant error in published TE annotations is the *library*, not the masker engine. Crucially, **masking percentage is robust to a bad library** (a chimeric consensus still masks roughly the right real estate), so the headline number survives while everything downstream rots: inflated family counts, wrong classification, distorted age landscapes, and - the killer - host-gene-contaminated consensi that silently mask real genes. Masking + gross % can use an automated library; any *per-family biological claim* (this family is young/active/novel) needs curation (Goubert 2022 *Mob DNA* 13:7; TE-Aid; MCHelper).

2. **Annotation quality is capped by assembly quality.** Short-read de Bruijn assemblers *collapse* near-identical TE copies and *drop* the youngest (most identical, most biologically active) ones, so short-read assemblies systematically under-count TEs and bias the age distribution toward "old" - which masquerades as the real signal "this lineage has no recent activity." Software cannot recover what the assembler threw away. Always ask what assembly a "% repeat" came from; HiFi/T2T raised the ceiling (LAI measures it) but T2T satellite/centromere repeats still exceed what the standard TE toolchain can annotate.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| RepeatModeler2 | Flynn 2020 *PNAS* | de novo family discovery -> consensus library | discover genome-specific families (run `-LTRStruct`) |
| RepeatMasker | Smit/Hubley/Green (software) | annotate/mask a genome **against** a library | the masking step; does not discover families |
| EarlGrey | Baril 2024 *MBE* | wraps RepeatModeler2 + auto consensus-elongation + RepeatMasker + plots | **non-model default**; minimal hand-work |
| EDTA | Ou 2019 *Genome Biol* | structural LTR/TIR/Helitron discovery + filtering | **plant / structurally-rich** genomes |
| LTR_retriever | Ou & Jiang 2018 *Plant Physiol* | isolate intact LTR-RTs; feeds LAI | LTR focus / assembly-quality (LAI) |
| TRF | Benson 1999 *NAR* | tandem/satellite repeats | a different algorithm class from TE maskers |
| RepeatClassifier / DeepTE / TERL | Flynn 2020; Yan 2020 | classify unknown consensi | attack the "Unknown" fraction (validate; can mislabel) |

Dfam (open) vs RepBase (paywalled since 2019) is the database schism - modern open pipelines build on Dfam + de novo. Engine `-e`: `rmblast` (default, consensus FASTA) vs `nhmmer` (Dfam profile HMMs, more sensitive for ancient repeats, slower) - the same genome reads a higher % with HMM detection.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Non-model eukaryote, defensible answer, minimal hand-work | EarlGrey | RepeatModeler2 + auto-curation + clean outputs |
| Plant / structurally-rich TE genome | EDTA | best-in-class LTR/TIR/Helitron structural annotation + host-gene filtering |
| Well-covered vertebrate, just need masking | RepeatMasker `-species` against Dfam | curated families already exist |
| Mask before gene prediction | RepeatMasker `-xsmall` (soft) + decontaminated library | predictors need soft-masking |
| Publication-grade TE biology claim | de novo -> manual curation (Goubert protocol, TE-Aid) | automated library is the start, not the end |
| Tandem/satellite/centromeric repeats | TRF + satellite tools (not RepeatMasker) | library-based TE tools don't see tandem arrays |
| TE expression from RNA-seq | -> TEtranscripts/SQuIRE (EM multimapper handling) | see expression section |
| TE insertion polymorphisms from reads | -> variant-calling (MELT/TEPID) | out of scope here |

## RepeatModeler2 -> RepeatMasker (the canonical pair)

```bash
# 1. De novo family discovery -> mydb-families.fa
BuildDatabase -name mydb assembly.fa
RepeatModeler -database mydb -threads 16 -LTRStruct   # -LTRStruct enables the LTR structural pipeline

# 2. (Recommended) decontaminate the library against host proteins, then UNION with Dfam clade
#    -> pull any consensus whose best hit is a host gene with no transposase/RT/integrase domain

# 3. Soft-mask against (custom library) for gene prediction
RepeatMasker -lib mydb-families.fa -xsmall -gff -e rmblast -pa 16 -dir rm_out assembly.fa
```

`-xsmall` = **soft-mask (lowercase)** - the key flag, the one people get wrong. Default `.masked` output hard-masks with N; `-x` masks with X. `-nolow` skips low-complexity/simple repeats (often wanted before gene prediction - see below). Outputs: `.masked`, `.out`, `.tbl` (summary %), `.align` (needed for the landscape).

## Soft vs Hard Masking (Critical, and the Over-Masking Massacre)

- **Gene prediction needs SOFT-masking.** Modern predictors (AUGUSTUS/GeneMark/BRAKER) avoid *nucleating* models in lowercase but let exons extend into repeats - real genes have TE-derived exons and TE-filled introns. **Hard-masking (N) destroys sequence**: any gene overlapping a repeat is truncated or never called, and the predictor reports nothing - no error, no log. The gene is simply absent.
- **Over-aggressive soft-masking is also a failure.** A gene whose promoter sits in an LTR, or a young gene inside a recent TE burst, gets suppressed because the predictor won't start a model in a heavily-lowercased locus. Before gene prediction, soft-mask *interspersed repeats only* and skip low-complexity (`-nolow`) - simple repeats overlap real coding microsatellites and low-complexity protein domains.
- **The domesticated-gene massacre.** A de novo library contains fragments of real multicopy gene families because they *look* repetitive - and masking them deletes the genome's most interesting genes. Named casualties: **RAG1/RAG2** (domesticated Transib transposase; V(D)J recombination), **syncytins** (captured retroviral env; placentation), **SETMAR/Metnase** (Hsmar1 mariner + SET domain; DNA repair), **CENP-B** (pogo transposase; centromere), and the **KRAB-ZNF arrays** (~350+ primate zinc-finger genes that exist *to repress TEs* - masking them deletes the genome's anti-TE machinery). Defense: **decontaminate the library against a protein DB before masking** (EDTA does a version; RepeatModeler2 does not by default). Treat any gene model falling entirely inside a masked "TE" as a flag to investigate, not a finished call.

## TE Classification (Wicker-compatible)

- **Class I (retrotransposons, copy-and-paste via RNA + RT):** LTR-RTs (Ty3/Gypsy, Ty1/Copia - dominate large plant genomes), LINEs (autonomous, often 5'-truncated), SINEs (non-autonomous, e.g. Alu).
- **Class II (DNA transposons, mostly cut-and-paste):** TIR superfamilies (hAT, Tc1/Mariner, CACTA, PIF/Harbinger, Mutator), Helitrons (rolling-circle, can capture host genes), MITEs (non-autonomous TIR derivatives, EDTA reclassifies ≤600 bp).
- **Family vs copy:** a *family* is one consensus in the library (~thousands); a *copy/insertion* is one genomic locus matching it (millions). "% genome masked" counts copies; classification is family-level; "10,000 TEs" is ambiguous between the two.
- **Full-length vs decayed:** most copies are dead, truncated, point-mutated relics; only a tiny fraction are intact. **Solo-LTR : full-length ratio** is real biology (recombinational LTR-RT removal rate), measurable only if the assembly resolved full-length elements. Wicker 2007 *Nat Rev Genet* is the reference scheme; present Gypsy/Copia and the ICTV Metaviridae/Pseudoviridae names both.

## Repeat Statistics and Age Landscape with Python

**Goal:** Summarize masked content by class and plot the Kimura-divergence landscape (a relative within-genome age readout).

**Approach:** Parse the RepeatMasker `.out` file, group by class for bp and genome fraction, then histogram percent divergence stratified by major TE class (x = divergence-from-consensus ~ relative age).

```python
import pandas as pd

def parse_repeatmasker_out(out_file):
    records = []
    with open(out_file) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            records.append({'perc_div': float(parts[1]), 'seqid': parts[4],
                            'repeat_class': parts[10], 'length': int(parts[6]) - int(parts[5]) + 1})
    return pd.DataFrame(records)

def repeat_summary(rm_df, genome_size):
    by_class = rm_df.groupby('repeat_class')['length'].sum().sort_values(ascending=False)
    total = rm_df['length'].sum()
    print(f'Total masked: {total/genome_size:.1%} of genome (a LOWER bound; ancient copies decay past detection)')
    return by_class / genome_size * 100
```

The landscape is **right-censored** - the most ancient TEs decayed past alignment detection, so "no old activity" can mean "old activity is invisible." A sharp left (low-divergence) peak is a recent/ongoing burst; treat *presence* of a recent peak as informative and *absence* of an old hump cautiously. A truncated/chimeric consensus distorts the whole x-axis (another reason curation matters); never compare landscapes across genomes annotated with different libraries.

## TE Expression from RNA-seq (the Multimapping Minefield)

A read from a young high-copy family maps equally to hundreds of near-identical loci. **Unique-only mapping** (standard RNA-seq QC) discards most TE signal and biases toward old, uniquely-mappable copies - measuring the *least* active elements. Use EM/probabilistic reassignment: **TEtranscripts/TElocal** (Jin 2015), **SQuIRE** (Yang 2019), **Telescope** (Bendall 2019). Subfamily-level (TEtranscripts: "L1 went up", high power, no locus) vs locus-level (SQuIRE/TElocal/Telescope: "this HERV-K on chr7 is on", noisy, mappability-sensitive) changes the conclusion, not just the resolution. The dominant false positive: **a TE in an intron or downstream of an expressed gene is not "expressed"** - read-through/intron-retention piles reads on it; distinguish autonomous transcription from passenger signal by strand and continuity (TEspeX filters embedded-TE reads). Be skeptical of any "TEs reactivated in disease/aging" headline that used unique-only mapping.

## Per-Method Failure Modes

### Hard-masking before gene prediction
**Trigger:** running RepeatMasker without `-xsmall` (default hard-masks with N). **Mechanism:** masked sequence is destroyed. **Symptom:** genes overlapping repeats silently absent from the GFF. **Fix:** `-xsmall`; hand a soft-masked genome to the predictor.

### Host-gene-contaminated library
**Trigger:** masking with an uncurated de novo library. **Mechanism:** multicopy gene families look repetitive and enter the library. **Symptom:** suspiciously few NLR/ZNF/OR genes; domesticated genes (RAG1, CENP-B) missing. **Fix:** BLAST the library against a protein DB; drop consensi hitting host genes with no TE domain.

### Trusting % repeat from a short-read assembly
**Trigger:** comparing TE content across studies/assemblies. **Mechanism:** short reads collapse/drop young copies; % depends on library+engine+assembly. **Symptom:** "low TE, all ancient" or non-comparable cross-study tables. **Fix:** check LAI/assembly type; report method + assembly with every number; never compare published % across papers.

### Discarding multimappers in TE RNA-seq
**Trigger:** unique-only TE quantification. **Mechanism:** young high-copy families are not uniquely mappable. **Symptom:** most TE signal lost, bias to old elements. **Fix:** EM tools (TEtranscripts/SQuIRE/Telescope); separate read-through from autonomous transcription.

### "Unknown" passed off as a result
**Trigger:** shipping a 40%-Unknown library without inspection. **Mechanism:** classification is the hardest, last, most-skipped step. **Symptom:** weak biological annotation; possible gene-family contamination hiding in Unknown. **Fix:** RepeatClassifier/DeepTE to triage; curate; note DB-coverage limits.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `-xsmall` soft-mask before gene prediction | predictor requirement | hard-mask truncates repeat-overlapping genes |
| TE content scales with genome size (human ~50%, maize ~85%, Arabidopsis ~20-25%, fungi ~1-20%) | clade norms (approx) | main driver of the C-value enigma; sanity-check vs genome size |
| "Unknown" ~<15% (mammal) vs 30-50% (non-model) | DB coverage | high Unknown bounds biological claims; very low on non-model = over-assignment |
| LAI <10 draft / 10-20 reference / >20 gold | Ou 2018 *NAR* | LTR-RT-resolution metric; only valid for LTR-rich genomes |
| Report library + engine + assembly with any % | reproducibility | % masked is non-comparable across methods |
| 80-80-80 (≥80% id over ≥80% length over ≥80 bp) | Wicker lineage | dereplication threshold, NOT a quality check |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Gene prediction finds too few genes | hard-masked, or over-masked low-complexity | `-xsmall`; `-nolow` before gene prediction |
| Suspiciously few NLR/ZNF/OR genes | uncurated library masked gene families | decontaminate library against a protein DB |
| Low masking percentage | novel repeats absent from DB | run RepeatModeler2 first; union de novo + Dfam |
| RepeatModeler very slow | normal for large genomes | `-threads`; consider EDTA (plants) or EarlGrey |
| "TE re-activated" result looks too clean | unique-only mapping / read-through | EM tools; check strand + continuity from neighbor |
| Cross-study % repeat disagree | different library/engine/assembly | re-annotate uniformly; report method |

## References

- Flynn JM, et al. 2020. RepeatModeler2 for automated genomic discovery of transposable element families. *PNAS* 117:9451-9457.
- Storer J, et al. 2021. The Dfam community resource of transposable element families, sequence models, and genome annotations. *Mob DNA* 12:2.
- Ou S, et al. 2019. Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline (EDTA). *Genome Biol* 20:275.
- Ou S, Jiang N. 2018. LTR_retriever: a highly accurate and sensitive program for identification of long terminal repeat retrotransposons. *Plant Physiol* 176:1410-1422.
- Ou S, Chen J, Jiang N. 2018. Assessing genome assembly quality using the LTR Assembly Index (LAI). *Nucleic Acids Res* 46:e126.
- Wicker T, et al. 2007. A unified classification system for eukaryotic transposable elements. *Nat Rev Genet* 8:973-982.
- Baril T, Galbraith J, Hayward A. 2024. Earl Grey: a fully automated user-friendly transposable element annotation and analysis pipeline. *Mol Biol Evol* 41:msae068.
- Goubert C, et al. 2022. A beginner's guide to manual curation of transposable elements. *Mob DNA* 13:7.
- Jin Y, et al. 2015. TEtranscripts: a package for including transposable elements in differential expression analysis of RNA-seq datasets. *Bioinformatics* 31:3593-3599.
- Yang WR, et al. 2019. SQuIRE reveals locus-specific regulation of interspersed repeat expression. *Nucleic Acids Res* 47:e27.
- Bendall ML, et al. 2019. Telescope: characterization of the retrotranscriptome by accurate estimation of transposable element expression. *PLoS Comput Biol* 15:e1006453.
- Benson G. 1999. Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Res* 27:573-580.

## Related Skills

- eukaryotic-gene-prediction - Receives the soft-masked genome; over-masking is a shared failure
- annotation-qc - LAI and repeat-content sanity in the assembly-to-annotation handoff
- genome-assembly/assembly-qc - Assembly type sets the TE-annotation ceiling (LAI)
- differential-expression/deseq2-basics - Differential TE expression from TEtranscripts/SQuIRE counts
- copy-number/recurrent-cnv - Segmental duplications, a distinct phenomenon from interspersed repeats
