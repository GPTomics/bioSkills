---
name: bio-workflows-genome-assembly-pipeline
description: Orchestrates an end-to-end de novo genome assembly project, routing each step to the right genome-assembly skill rather than restating it. Profiles the genome first (k-mer spectrum -> size, heterozygosity, ploidy), QCs reads, chooses an assembly path by data type (SPAdes for Illumina, Flye for noisy long reads, hifiasm for HiFi, metaFlye for communities), polishes only when needed, decontaminates, scaffolds with Hi-C, and finishes with three-axis QC (contiguity + completeness + correctness). Use when assembling a genome from raw reads and deciding which assembler, whether to polish, and how to prove the result is good.
tool_type: cli
primary_tool: Flye
workflow: true
depends_on:
  - genome-assembly/genome-profiling
  - read-qc/fastp-workflow
  - long-read-sequencing/long-read-qc
  - genome-assembly/short-read-assembly
  - genome-assembly/long-read-assembly
  - genome-assembly/hifi-assembly
  - genome-assembly/metagenome-assembly
  - genome-assembly/assembly-polishing
  - genome-assembly/contamination-detection
  - genome-assembly/scaffolding
  - genome-assembly/assembly-qc
qc_checkpoints:
  - after_profiling: "Genome-size and heterozygosity estimate obtained; sets NG50 denominator, purge level, assembler choice"
  - after_assembly: "Total length within ~10-20% of profiled size; contig count plausible for read type"
  - after_polishing: "Merqury QV improved or plateaued (do not over-polish HiFi); k-mers from accurate reads, not the polishing reads"
  - after_decontamination: "Single-organism: FCS-GX/BlobToolKit clean; MAG: CheckM2 >90% complete, <5% contam, GUNC pass"
  - after_scaffolding: "Contact map shows clean diagonal; off-diagonal blocks inspected/broken before calling chromosome-scale"
  - final_three_axis_qc: "Contiguity (auN/NG50 vs profiled size) + completeness (BUSCO/compleasm) + correctness (Merqury QV) all reported; never N50 alone"
---

## Version Compatibility

Reference examples tested with: GenomeScope2 2.0+, meryl 1.4+, Merqury 1.3+, fastp 0.23+, SPAdes 4.0+, Flye 2.9+, hifiasm 0.25+, metaFlye 2.9+, Racon 1.5+, medaka 2.0+, minimap2 2.26+, FCS-GX 0.5+, CheckM2 1.0+, GUNC 1.0+, YaHS 1.2+, QUAST 5.2+, BUSCO 5.5+, samtools 1.19+. Each owning genome-assembly skill is the source of truth for its tool's pinned version.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Tool outputs are driven by more than the binary version: medaka consensus quality depends on the basecaller MODEL string (must match the basecaller, e.g. `-m r1041_e82_400bps_sup_v5.0.0`); BUSCO/compleasm results depend on the lineage dataset and OrthoDB generation (record them); CheckM2/GTDB-Tk results track the reference DATABASE release; hifiasm output filenames and default purge behaviour change across versions (verify against the installed build). If a command errors, introspect the installed tool and adapt rather than retrying.

# Genome Assembly Pipeline

**"Assemble a genome from my sequencing reads and prove it is good"** -> Profile the genome, QC reads, route to the right assembler by data type, polish only if needed, decontaminate, scaffold if Hi-C exists, and finish with three-axis QC. This skill ORCHESTRATES the genome-assembly category; it routes each step to the owning skill and encodes the cross-cutting decisions, not each tool's full option set.

## The Single Most Important Modern Insight -- Assembly Is Three Orthogonal Questions, and Each Step Answers One

A genome project fails when one number stands in for the whole. Profiling sets expectations (how big, how heterozygous, how many haplotypes) BEFORE assembling, the assembler answers contiguity, polishing answers per-base accuracy, decontamination answers provenance, scaffolding answers arrangement, and QC must independently address all three of contiguity, completeness, and correctness. The orchestration job is to keep these separate and route each to its skill: a high N50 says nothing about whether the bases are right (Merqury QV) or whether the sequence is the organism's (contamination), and skipping profiling means the assembler guesses the parameters that profiling would have set.

## Decision Flow (Step 0 -> 6)

```
Raw reads (+ optional Hi-C, trio, short reads)
    |
    v
[0. Profile the genome] --> genome-assembly/genome-profiling
    |   k-mer spectrum (GenomeScope2) -> genome size, heterozygosity, ploidy.
    |   Sets NG50 denominator, expected haplotype count, hifiasm purge level,
    |   and which assembly path is even sensible. Do this BEFORE assembling.
    v
[1. QC reads] -----------> short: read-qc/fastp-workflow
    |                       long:  long-read-sequencing/long-read-qc
    |   Garbage-in caps assembly quality; record platform + basecaller era
    |   (it is an assembly PARAMETER, see step 2), trim internal adapters.
    v
[2. Choose path BY DATA TYPE]
    |  Illumina-only small/isolate -> genome-assembly/short-read-assembly (SPAdes)
    |  noisy ONT/CLR              -> genome-assembly/long-read-assembly (Flye --nano-hq for R10)
    |  PacBio HiFi                -> genome-assembly/hifi-assembly (hifiasm, phased)
    |  community sample           -> genome-assembly/metagenome-assembly (metaFlye/metaSPAdes + binning)
    |  large/heterozygous euk     -> long-read or HiFi, NOT short reads
    v
[3. Polish IF needed] ---> genome-assembly/assembly-polishing
    |   noisy long-read assemblies: Racon -> medaka (model MUST match basecaller).
    |   Do NOT polish HiFi reflexively (often net-harmful). Measure with Merqury QV,
    |   not the reads polished with. Skip entirely for SPAdes/HiFi when QV is already high.
    v
[4. Decontaminate] ------> genome-assembly/contamination-detection
    |   single organism: FCS-GX (GenBank-mandatory) + BlobToolKit blob plot.
    |   MAG:              CheckM2 + GUNC (chimerism). Two disjoint problems (see below).
    v
[5. Scaffold IF Hi-C] ---> genome-assembly/scaffolding
    |   automated YaHS produces a DRAFT; manual contact-map curation is the standard.
    |   Scaffold N50 != contig N50 (gaps are Ns). Skip if no Hi-C.
    v
[6. Three-axis QC] ------> genome-assembly/assembly-qc
        contiguity (auN/NG50 vs profiled size) + completeness (BUSCO/compleasm)
        + correctness (Merqury QV). Report the triad; NEVER N50 alone.
```

## Routing Table by Scenario

| Scenario | Path | Routes to |
|----------|------|-----------|
| Bacterial isolate, ONT R10 only | profile -> QC -> Flye `--nano-hq` -> medaka -> FCS-GX -> QC | long-read-assembly, assembly-polishing, contamination-detection |
| Bacterial isolate, Illumina only | profile -> fastp -> SPAdes `--isolate` -> FCS-GX -> QC | short-read-assembly |
| Small genome, ONT, max quality | profile -> QC -> multi-assembler consensus (Trycycler/Autocycler) -> medaka -> QC | long-read-assembly |
| Diploid eukaryote, HiFi (+Hi-C/trio) | profile -> QC -> hifiasm (hap1/hap2) -> purge check -> decontam -> scaffold -> QC | hifi-assembly, scaffolding, contamination-detection |
| Large heterozygous eukaryote, ONT | profile -> QC -> Flye -> purge_dups -> medaka -> decontam -> scaffold -> QC | long-read-assembly, scaffolding |
| Community / microbiome sample | QC -> metaFlye/metaSPAdes -> binning -> CheckM2 + GUNC | metagenome-assembly, contamination-detection |
| Hi-C reads available | after contigs+polish: scaffold, curate contact map | scaffolding |
| Reads not yet QC'd | start at step 1 | read-qc/fastp-workflow, long-read-sequencing/long-read-qc |

## Cross-Cutting Gotchas (surface these at every project)

- **Basecaller era must match the assembler flag.** `--nano-raw` on R10/Dorado-SUP reads silently collapses real repeats while RAISING N50; `--nano-hq` is the R10 default. The platform + basecaller model is an assembly parameter, not metadata.
- **A primary assembly is not a haplotype.** The hifiasm primary is a maternal/paternal mosaic that exists in no cell; for any allele-aware downstream use hap1/hap2 phased with trio or Hi-C, and treat HiFi-only hap1/hap2 as only partially phased.
- **N50 is gamed.** It rises when an assembly gets WORSE (misjoins, collapsed repeats, retained haplotigs). Report the triad (auN/NG50 + BUSCO + Merqury QV), never N50 alone.
- **A MAG is a population consensus, not a genome.** The unit of success is a binned, MIMAG-gated MAG, and "% contamination" conflates foreign-organism mixing, strain mixing, and assembly artifacts.
- **"Contamination" is two disjoint problems.** Single-organism cross-kingdom foreign sequence (FCS-GX, blob plot) is a different question from intra-domain MAG contamination/chimerism (CheckM2 + GUNC); do not apply one tool's question to the other's input.
- **Scaffold N50 >> contig N50 because gaps are Ns.** Scaffold contiguity is glue, not sequence; every join is a hypothesis a contact map must confirm. Report contig N50 alongside scaffold N50.

## Step 0: Profile the Genome (do this first)

Route the full treatment to genome-assembly/genome-profiling. The minimal orchestration step:

```bash
# k-mer count from ACCURATE reads (Illumina/HiFi, NEVER noisy ONT), then GenomeScope2 for size / heterozygosity / ploidy
meryl count k=21 output reads.meryl accurate_reads.fq.gz
meryl histogram reads.meryl > reads.hist
genomescope2 -i reads.hist -o gscope_out -k 21
# read off: estimated haploid genome size, heterozygosity %, and (with -p) ploidy.
# These set the NG50 denominator, the expected number of haplotypes, and the purge decision.
```

## Step 1: QC Reads

Short reads route to read-qc/fastp-workflow; long reads to long-read-sequencing/long-read-qc.

```bash
fastp -i R1.fq.gz -I R2.fq.gz -o t_R1.fq.gz -O t_R2.fq.gz \
    --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 --html qc.html
```

## Step 2: Assemble (route by data type)

Give the assembler the exact preset for the chemistry; the wrong preset is silent. Detailed options live in the owning skills.

```bash
# Illumina-only small/isolate genome -> short-read-assembly
spades.py --isolate -1 t_R1.fq.gz -2 t_R2.fq.gz -o spades_out -t 16
# NOTE: --careful is small-genome-only; do NOT use it on large eukaryote genomes.

# Noisy ONT (R10/Dorado-SUP) -> long-read-assembly. --nano-hq is the modern default.
flye --nano-hq ont.fq.gz --out-dir flye_out --threads 16     # --genome-size optional in recent Flye

# PacBio HiFi -> hifi-assembly (phased by default; verify output filenames per version)
hifiasm -o asm -t 16 hifi.fq.gz                              # add --h1/--h2 (Hi-C) or -1/-2 (trio) to phase

# Community sample -> metagenome-assembly
flye --nano-hq ont.fq.gz --meta --out-dir metaflye_out --threads 16    # --meta is a modifier; still need a read-type selector. Then bin + CheckM2/GUNC
```

## Step 3: Polish IF Needed

Polishing is read-type-matched and conditional. Route to genome-assembly/assembly-polishing.

```bash
# Noisy long-read assembly: medaka with the MATCHING model. medaka_consensus does its own
# read-to-assembly alignment from -i/-d (no separate minimap2/BAM step needed); add a Racon
# round upstream only if the assembler did not already polish - see assembly-polishing.
medaka_consensus -i ont.fq.gz -d flye_out/assembly.fasta -o medaka_out -t 16 \
    -m r1041_e82_400bps_sup_v5.0.0   # MUST match the basecaller model used to call the reads
# For a BACTERIAL isolate, prefer the methylation-aware bacterial model (medaka 2.0+):
#   medaka_consensus -i ont.fq.gz -d assembly.fasta -o medaka_out --bacteria
```

Do NOT reflexively polish a HiFi assembly (already ~Q30+; over-polishing lowers QV). SPAdes output needs no separate long-read polish. The stop signal is a Merqury QV plateau, not a fixed iteration count, and the QV must be measured against reads independent of those used to polish.

## Step 4: Decontaminate (route by sample type)

```bash
# Single-organism assembly (GenBank-mandatory foreign screen + blob plot)
python3 ./fcs.py screen genome --fasta assembly.fa --out-dir gx_out/ --gx-db "$GXDB/gxdb" --tax-id <taxid>
# acts on EXCLUDE/TRIM/FIX cross-kingdom contigs; keep host-integrated foreign sequence (see contamination-detection)

# MAG (intra-domain contamination + chimerism)
checkm2 predict --input bins/ --output-directory checkm2_out --threads 16
gunc run --input_dir bins/ --out_dir gunc_out                            # chimerism, orthogonal to CheckM2
```

## Step 5: Scaffold IF Hi-C Is Available

YaHS produces a draft; the contact map is the QC, not decoration. Route to genome-assembly/scaffolding.

```bash
# Map Hi-C to contigs, then YaHS; inspect the contact map (PretextMap/Juicer) and break misjoins.
yahs assembly.fasta hic_to_contigs.bam -o yahs_out          # output scaffolds + AGP; curate before publishing
```

## Step 6: Three-Axis QC (contiguity + completeness + correctness)

Route the full treatment to genome-assembly/assembly-qc. Report all three axes; lead with the QV.

```bash
# Contiguity vs the PROFILED genome size (NG50/auN, not bare N50)
quast.py final.fasta -o quast_out -t 16 --est-ref-size <profiled_size>

# Completeness on the DEEPEST applicable clade (compleasm on good genomes; BUSCO otherwise)
busco -i final.fasta -l <clade>_odb10 -o busco_out -m genome -c 16

# Correctness: Merqury QV from ACCURATE reads (k from best_k.sh, not hardcoded)
K=$(sh $MERQURY/best_k.sh <genome_size_bp> | tail -n1 | awk '{print int($1+0.5)}')   # round float->int
meryl count k=$K output reads.meryl accurate_reads.fq.gz
merqury.sh reads.meryl final.fasta merqury_out          # QV + k-mer completeness + spectra-cn
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Assembly ~1.5-2x profiled size, high BUSCO-Duplicated | uncollapsed haplotigs (false duplication) | purge_dups; check half-coverage depth peak; do not over-purge real segmental duplications |
| Contiguous but gene models frameshift | noisy long-read assembly not polished | Racon -> medaka (matched model); measure QV |
| QV drops after polishing | over-polishing an already-accurate (HiFi) assembly | stop polishing; HiFi rarely needs short-read polish |
| medaka consensus worse than input | wrong basecaller model string | set `-m` to the model the reads were basecalled with |
| Fewer contigs than expected but repeats collapsed | `--nano-raw` used on R10 reads | re-run Flye with `--nano-hq` |
| CheckM2 says clean but bin looks mixed | chimera with disjoint markers | run GUNC; CheckM2 marker redundancy cannot see chimerism |
| Scaffold N50 huge, contig N50 small | scaffolding glue, not sequence | inspect contact map, break off-diagonal misjoins |

## Related Skills

- genome-assembly/genome-profiling - Step 0: k-mer spectrum for size, heterozygosity, ploidy; sets expectations before assembling
- genome-assembly/short-read-assembly - SPAdes path for Illumina-only small/isolate genomes
- genome-assembly/long-read-assembly - Flye/Canu path for noisy ONT/CLR reads
- genome-assembly/hifi-assembly - hifiasm phased path for PacBio HiFi
- genome-assembly/metagenome-assembly - metaFlye/metaSPAdes + binning for community samples
- genome-assembly/assembly-polishing - Racon/medaka/Pilon, applied only when needed
- genome-assembly/contamination-detection - FCS-GX/BlobToolKit (single organism) vs CheckM2/GUNC (MAG)
- genome-assembly/scaffolding - YaHS Hi-C scaffolding and contact-map curation
- genome-assembly/assembly-qc - Three-axis QC: auN/NG50 + BUSCO + Merqury QV
- read-qc/fastp-workflow - Short-read QC before assembly
- long-read-sequencing/long-read-qc - Long-read length/quality QC and basecaller-era awareness
- workflows/genome-annotation-pipeline - Downstream: only a decontaminated, QC-passed FASTA should hand off to annotation

## References

- Rhie A, Walenz BP, Koren S, Phillippy AM (2020) Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biology* 21:245. DOI 10.1186/s13059-020-02134-9. (QV/completeness from k-mers.)
- Manni M, Berkeley MR, Seppey M, Simao FA, Zdobnov EM (2021) BUSCO update: novel and streamlined workflows. *Molecular Biology and Evolution* 38:4647-4654. DOI 10.1093/molbev/msab199.
- Rhie A, McCarthy SA, Fedrigo O, et al (2021) Towards complete and error-free genome assemblies of all vertebrate species. *Nature* 592:737-746. DOI 10.1038/s41586-021-03451-0. (three-axis / VGP standard.)
- Ranallo-Benavidez TR, Jaron KS, Schatz MC (2020) GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. *Nature Communications* 11:1432. DOI 10.1038/s41467-020-14998-3.
