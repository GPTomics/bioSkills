---
name: bio-genome-annotation-ncrna-annotation
description: Identifies non-coding RNAs (tRNA, rRNA, snoRNA, snRNA, riboswitches, sRNAs) using Infernal covariance-model search against Rfam, tRNAscan-SE 2.0 for tRNA, barrnap for rRNA, and ARAGORN for tmRNA, plus the small-RNA-seq boundary for miRNA and the transcript-assembly boundary for lncRNA. Covers the structure-conserved-not-sequence-conserved principle (why BLAST fails), GA-threshold and clan-competition correctness, tRNAscan-SE domain modes and pseudogene flags, rDNA copy-number collapse, and why homology annotation is a recall floor. Use when performing genome-wide ncRNA annotation, choosing the right tool for an RNA class, or interpreting ncRNA counts.
tool_type: cli
primary_tool: Infernal
---

## Version Compatibility

Reference examples tested with: Infernal 1.1.4+, Rfam 15+ (CM library), tRNAscan-SE 2.0.12+, barrnap 0.9+, ARAGORN 1.2+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

The **Rfam release version** drives results (Rfam 15 has ~4,200+ families); record it. The downloaded `Rfam.cm` ships pre-calibrated, so only `cmpress` is needed before `cmscan`. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Non-Coding RNA Annotation

**"Find non-coding RNAs in my genome"** -> Scan an assembly for structured ncRNA families using covariance models (sequence + secondary structure jointly), with specialist detectors for tRNA and the expression boundary for miRNA/lncRNA.
- CLI: `cmscan --cut_ga --rfam --nohmmonly --fmt 2 --clanin Rfam.clanin Rfam.cm genome.fa` (Infernal), `tRNAscan-SE -B genome.fa`

## The Single Most Important Modern Insight -- ncRNA Homology Is Structure-Conserved, Not Sequence-Conserved

Protein annotation rides on sequence/ORF signal; structured-ncRNA annotation rides on **base-pairing covariation**. A G-C pair can become A-U across evolution while the *structure* is preserved - both positions mutate together (compensatory substitution). To a sequence-only tool these look like two mismatches; to a covariance model (a profile SCFG, the Rfam/Infernal engine) the *correlated* change is the strongest possible evidence of homology. Three consequences:

1. **BLAST is categorically the wrong tool for structured ncRNA.** Two RNase P RNAs may share <60% identity (BLAST sees noise) yet have unmistakable shared structure. Use Infernal covariance models, not BLAST. R-scape (Rivas 2017 *Nat Methods* 14:45) is the statistical test for whether a structure is *actually* conserved - it found no significant covariation support for the proposed structures of HOTAIR, SRA, and Xist-RepA.
2. **A homology-based annotation is a recall floor, never a count.** A CM can only exist for a family whose structure is conserved across enough divergent sequences to seed a model - so fast-evolving and lineage-specific ncRNAs are invisible *by design*, and the floor drops further for organisms far from the curation spotlight. Report "at least N conserved-family loci," never "the genome has N ncRNAs."
3. **Whole classes need expression evidence, not genomic search.** miRNAs are hypotheses until small-RNA-seq confirms the Dicer processing signature; lncRNAs are not annotatable by homology at all (no conserved structure) - they are transcript catalogs. Use specialist tools where the biology has a sharper signal (tRNAscan-SE, miRDeep2-with-reads).

## Tool Taxonomy

| RNA class | Tool | Citation | Method |
|-----------|------|----------|--------|
| tRNA | tRNAscan-SE 2.0 | Chan 2021 *NAR* | isotype-specific Infernal CMs + pseudogene/high-confidence logic |
| rRNA (fast) | barrnap | Seemann (software) | nhmmer HMM profiles; kingdom flag; prokaryotic-pipeline default |
| rRNA (structure-aware) | Infernal + Rfam SSU/LSU | Nawrocki 2013 | CM; better boundaries / unusual taxa |
| tmRNA (+ bacterial tRNA) | ARAGORN | Laslett 2004 | heuristic cloverleaf + tmRNA models |
| miRNA | miRDeep2 (+ small-RNA-seq) | Friedländer 2012 | Dicer-processing model on read pileups |
| C/D, H/ACA snoRNA | snoscan / snoReport | Lowe 1999 | guide-target complementarity / SVM |
| Everything else structured | Infernal `cmscan` vs Rfam | Nawrocki 2013 | covariance models; the general engine |
| lncRNA | StringTie + CPC2/CPAT/FEELnc | - | transcript assembly + coding-potential, NOT CM search |

RNAmmer is the legacy rRNA tool (license-encumbered, HMMER2) - use barrnap instead unless reproducing old annotations.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Prokaryote, fast complete annotation | barrnap + tRNAscan-SE `-B`/ARAGORN + Infernal/Rfam | what Bakta/Prokka/PGAP wrap |
| tRNA is the question | always tRNAscan-SE 2.0 (not Rfam's generic tRNA model) | isotype, pseudogene, intron, high-confidence logic |
| rRNA, speed matters | barrnap | seconds per genome |
| Broad ncRNA sweep of a new genome | Infernal `cmscan` vs full Rfam.cm (GA + clan competition) | structure-aware, family-typed |
| miRNA | demand small-RNA-seq; miRDeep2 | genomic hairpin prediction is unreliable |
| lncRNA | transcript assembly + coding-potential | not structurally conserved; no CM |
| Claim a conserved structure | R-scape covariation test (report power) | thermodynamic fold != selected structure |
| Bacterial AMR/CRISPR arrays | -> prokaryotic-annotation / CRISPRCasFinder | array detection is a separate tool class |

## Infernal / cmscan (the General ncRNA Engine)

```bash
# Rfam ships pre-calibrated; press it once, then scan (cmscan = many models vs one genome)
cmpress Rfam.cm
cmscan -Z <dbsize_Mb> --cut_ga --rfam --nohmmonly \
       --tblout out.tblout --fmt 2 --clanin Rfam.clanin \
       Rfam.cm genome.fa > out.cmscan
grep -v " = " out.tblout > out.deoverlapped.tblout   # drop within-clan overlaps
```

- `--cut_ga` (gathering threshold): the single most important correctness flag. Each family has a curator-set, per-family bit-score threshold; a fixed E-value would treat a 70-nt tRNA model and a 2,900-nt LSU model identically, which is wrong. **Overriding GA to "find more" imports the false positives the curator deliberately excluded.**
- `--nohmmonly` forces full CM (structure-aware) scoring so scores are GA-comparable.
- `-Z <dbsize_Mb>` = total_residues x 2 / 1e6 (both strands), making E-values run-comparable.
- `--fmt 2 --clanin Rfam.clanin` + the `grep -v " = "` deoverlap step is **mandatory, not a nicety**: clans group related families (the tRNA models, SSU/LSU rRNA), so one locus hits several models and the raw table double-counts (a 16S locus becomes "several rRNA genes").

## tRNAscan-SE 2.0

```bash
tRNAscan-SE -B -o trna.out -f trna.ss -m trna.stats --gff trna.gff3 genome.fa   # bacterial
```

Modes: `-E` eukaryotic (default), `-B` bacterial, `-A` archaeal, `-G` general (mixed/metagenome), `-M mammal`/`-M vert` mitochondrial, `-O` organellar (disables pseudogene checking). **Domain choice is not cosmetic** - the wrong mode mis-scores and miscalls isotypes; there is no auto-detect. Report the **high-confidence set**, not raw hits (raw counts include pseudogenes/SINEs and can be 2-10x inflated in eukaryotes). The pseudogene flag is reliable in eukaryotic nuclear genomes but false-positive-prone in organelles/odd mito-tRNAs (truncated arms read as "decayed") - hence `-O`/`-D`.

## barrnap (rRNA)

```bash
barrnap --kingdom bac genome.fa > rrna.gff3   # bac | arc | euk | mito
```

Reports partial rRNA at contig edges as `(partial)`. **rDNA copy number from an assembly is essentially always wrong** - near-identical rRNA arrays collapse in short-read assemblies, so the annotated count is a floor (off by orders of magnitude in eukaryotes); use long reads or read depth for true copy number.

## Parsing and Combining ncRNA Calls with Python

**Goal:** Merge Infernal and tRNAscan-SE into one ncRNA annotation, preferring the specialist for tRNA.

**Approach:** Parse the deoverlapped Infernal table, drop its tRNA rows (tRNAscan-SE is the authority for tRNA), and combine with the tRNAscan-SE high-confidence set; keep evidence provenance per class.

```python
import pandas as pd

def parse_infernal_tbl(tbl_file):
    rows = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 18:
                continue
            rows.append({'rfam_name': p[1], 'seqid': p[3], 'strand': p[11],
                         'score': float(p[16]), 'evalue': float(p[17])})
    df = pd.DataFrame(rows)
    return df[~df['rfam_name'].str.contains('tRNA', case=False, na=False)]   # tRNAscan-SE owns tRNA
```

## The miRNA Disaster and lncRNA Non-Annotatability

- **Most computationally predicted miRNAs are false positives, and much of miRBase is contaminated.** Any genome folds into astronomically many hairpins; foldability is not evidence of a miRNA. The discriminating signal - a precise homogeneous 5' end, a detectable star strand, a ~22 nt mode - is visible only in small-RNA-seq read pileups (miRDeep2 scores this geometry). Fromm 2015 found <1/3 of human and ~16% of metazoan miRBase entries are bona fide; closely related species differing by >1,000 miRNAs in miRBase is annotation noise, not biology. Animal-trained tools mis-call **plant** miRNAs (DCL1 biogenesis, 21/24 nt) - use plant-specific tools. Report homology-only hits as "miRNA candidates," never genes.
- **lncRNAs are not annotatable by homology in principle** - no conserved structure (R-scape), poor sequence conservation. "lncRNA annotation" is transcript cataloguing (assemble RNA-seq, subtract coding potential), inheriting every transcript-assembly pathology, so catalog size tracks sequencing depth and filter policy. The same human genome "has" ~16,000 to >100,000 lncRNAs depending only on the annotation source (RefSeq vs GENCODE vs NONCODE) - these encode different answers to the transcription-vs-selection function debate (Graur 2013), not different biology. Never quote a bare lncRNA count without source + version.

## Per-Method Failure Modes

### BLAST for structured ncRNA
**Trigger:** using BLAST to find ncRNA homologs. **Mechanism:** BLAST scores sequence identity, blind to covariation. **Symptom:** misses structure-conserved homologs; ranks pseudogenes above real homologs. **Fix:** Infernal/Rfam covariance models.

### Missing GA / clan competition
**Trigger:** `cmscan` without `--cut_ga`/`--nohmmonly`, or without `--clanin`/`--fmt 2` deoverlap. **Mechanism:** wrong thresholds; HMM-only scores not GA-comparable; within-clan overlaps uncollapsed. **Symptom:** inflated/redundant family calls (one locus counted several times). **Fix:** the full canonical command + `grep -v " = "`.

### Wrong tRNAscan-SE domain mode
**Trigger:** `-E` on a bacterium, or default on a metagenome. **Mechanism:** domain CMs differ; no auto-detect. **Symptom:** mis-scored/miscalled isotypes. **Fix:** `-B`/`-A`/`-G` per source; report the high-confidence set.

### Counting raw tRNA/rRNA hits
**Trigger:** reporting raw tRNAscan-SE hits or genomic rRNA as copy number. **Mechanism:** SINEs/pseudogenes/NUMTs inflate tRNA; rDNA arrays collapse. **Symptom:** 150+ bacterial "tRNAs"; "5 copies of 18S". **Fix:** high-confidence set; treat copy numbers as floors; flag NUMT-type organellar tRNAs in the nucleus.

### miRNA/lncRNA from genome alone
**Trigger:** calling miRNAs from hairpins or lncRNAs from CM search. **Mechanism:** no genomic signal suffices. **Symptom:** false-positive "genes". **Fix:** small-RNA-seq (miRNA); transcript assembly + coding-potential (lncRNA); label homology-only as candidates.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `--cut_ga` per-family GA threshold | Rfam curation | family-specific noise floor; never override genome-wide |
| `-Z` = residues x 2 / 1e6 | Infernal | run-comparable E-values (both strands) |
| Bacterial tRNA ~28-90 (E. coli ~89; symbionts ~28-35) | GtRNAdb | 150+ implies fragments/pseudogenes |
| Eukaryotic tRNA ~170-570 (amplified) | tRNA literature | raw hits far higher (SINEs/pseudogenes); use high-confidence set |
| Bacterial rRNA operons 1-15 (E. coli ~7) | copy-number norm | annotated count is a floor (rDNA collapse) |
| miRNA needs small-RNA-seq Dicer signature | Fromm 2015 | foldability is not a miRNA; >1/3 miRBase is artifact |
| R-scape covariation + power before claiming structure | Rivas 2017/2020 | thermodynamic fold != selected structure |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| cmscan slow | full Rfam scan | `--rfam` preset; split genome; parallelize |
| Redundant overlapping calls | no clan competition | `--fmt 2 --clanin`; `grep -v " = "` |
| Missing expected ncRNAs | stale Rfam / cmpress not run | check Rfam version; verify `.i1{f,i,m,p}` files |
| Too many tRNA "pseudogenes" | normal in eukaryotes; false-positive in organelles | high-confidence set; `-O`/`-D` for organelles |
| Implausible rRNA copy number | rDNA array collapse | report as floor; use read depth/long reads |
| Two species differ by 1000s of miRNAs | miRBase contamination | use MirGeneDB; demand processing evidence |

## References

- Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics* 29:2933-2935.
- Ontiveros-Palacios N, et al. 2025. Rfam 15: RNA families database in 2025. *Nucleic Acids Res* 53:D258-D267.
- Chan PP, et al. 2021. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Res* 49:9077-9096.
- Lowe TM, Eddy SR. 1997. tRNAscan-SE: a program for improved detection of transfer RNA genes in genomic sequence. *Nucleic Acids Res* 25:955-964.
- Laslett D, Canback B. 2004. ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. *Nucleic Acids Res* 32:11-16.
- Lagesen K, et al. 2007. RNAmmer: consistent and rapid annotation of ribosomal RNA genes. *Nucleic Acids Res* 35:3100-3108.
- Friedländer MR, et al. 2012. miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. *Nucleic Acids Res* 40:37-52.
- Lowe TM, Eddy SR. 1999. A computational screen for methylation guide snoRNAs in yeast (snoscan). *Science* 283:1168-1171.
- Fromm B, et al. 2015. A uniform system for the annotation of vertebrate microRNA genes and the evolution of the human microRNAome (MirGeneDB). *Annu Rev Genet* 49:213-242.
- Rivas E, Clements J, Eddy SR. 2017. A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs (R-scape). *Nat Methods* 14:45-48.
- Graur D, et al. 2013. On the immortality of television sets: "function" in the human genome according to the evolution-free gospel of ENCODE. *Genome Biol Evol* 5:578-590.

## Related Skills

- prokaryotic-annotation - Bakta/Prokka wrap barrnap + tRNAscan-SE + Infernal for prokaryotic ncRNA
- eukaryotic-gene-prediction - Protein-coding prediction does not find ncRNAs
- annotation-qc - tRNA/rRNA count sanity in the annotation QC panel
- rna-structure/ncrna-search - Targeted covariance-model homology searches
- rna-structure/secondary-structure-prediction - Fold and visualize an annotated ncRNA
