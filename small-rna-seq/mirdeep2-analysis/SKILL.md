---
name: bio-small-rna-seq-mirdeep2-analysis
description: Discovers novel miRNAs and quantifies known miRNAs with miRDeep2 by scoring genome-mapped read stacks against the Dicer/Drosha biogenesis signature. Use when deciding whether a study needs de novo discovery at all versus known-miRNA quantification; choosing the species and related-species miRBase references; reading the miRDeep2 score as a signal-to-noise hypothesis rather than a fixed cutoff; or filtering novel candidates against tRNA/rRNA loci to reject the classic false positives.
tool_type: cli
primary_tool: miRDeep2
---

## Version Compatibility

Reference examples tested with: miRDeep2 2.0.1.3+, bowtie 1.3+ (NOT bowtie2), ViennaRNA 2.5+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# miRDeep2 Analysis

**"Discover novel miRNAs from my small RNA-seq data"** -> Map collapsed reads to the genome, excise candidate hairpins, fold them, and score how well the observed read stacks match the Dicer/Drosha processing signature.
- CLI: `mapper.pl` (map to genome, emit ARF) -> `miRDeep2.pl` (discover + quantify) -> `quantifier.pl` (known-only quantification)

## The governing principle: a miRDeep2 score is a biogenesis hypothesis, not a validated miRNA

miRDeep2 does not detect miRNAs by sequence; it asks whether the reads piled on a genomic hairpin look like the product of Dicer/Drosha processing: a sharp, abundant MATURE arm, a lower-abundance STAR (passenger) arm with the correct ~2-nt 3' overhang geometry, a depleted loop, and a thermodynamically stable fold whose minimum free energy is lower than shuffled controls (the randfold p-value). A log-odds model converts that fit into a score (Friedländer 2012). The decisive consequence is that any locus producing a stacked, hairpin-foldable read pile can mimic the signature, so novel discovery is intrinsically high false-positive. The textbook failure is contaminating tRNA and rRNA fragments: tRNAs fold into stable cloverleaf arms and throw sharp, abundant read stacks that score as "novel miRNAs." A high score is a structural and expression hypothesis that demands orthogonal validation, never a finding.

There is no universal score cutoff. `survey.pl` sweeps cutoffs and reports, at each, the estimated true positives, false positives, signal-to-noise ratio, and an estimated FDR derived from permuted controls; Friedländer 2012 chose, per analysis, the lowest cutoff giving signal-to-noise >= 5. Asserting "score > 10 = high confidence" as a fixed rule is folklore: read the survey output, pick a cutoff for an acceptable estimated FDR, and report it.

## Decision: is miRDeep2 the right tool?

| Goal | Use | Why |
|------|-----|-----|
| Discover NOVEL miRNAs in an animal genome | miRDeep2 (full discovery) | The dedicated probabilistic biogenesis model; genome-anchored |
| Quantify KNOWN miRNAs + isomiRs + tRFs on a supported species | mirge3-analysis | Faster, isomiR-aware; discovery machinery is expensive and high-FP |
| Quantify KNOWN miRNAs only, no discovery | `quantifier.pl` (miRDeep2) or mirge3 | Skip the discovery engine when discovery is not needed |
| Profile tRFs / piRNAs (not miRNAs) | trf-pirna-profiling | tRF/rRF stacks are miRDeep2 false positives, not the target |
| Plant small RNAs | ShortStack (see trf-pirna-profiling) | Plant hairpins and 24-nt siRNA biology break the animal model |
| Animal with NO genome assembly (non-model, single-cell) | Mirnovo (genome-free ML) | miRDeep2 is genome-anchored and cannot run without an assembly |

miRDeep2 requires a reference GENOME and bowtie 1 (not bowtie2). The species and related-species miRBase references are load-bearing: the same-species mature/hairpin define "known," and the other-species mature provides conservation evidence that raises confidence in novel calls.

## Workflow overview

```
collapsed reads (FASTA, _xN counts)
    |
    v   mapper.pl  --> bowtie align to genome, emit ARF
    v
miRDeep2.pl  --> excise hairpins, fold (RNAfold), randfold, score read stacks
    |
    v   quantifier.pl  --> known-miRNA counts (run alone if no discovery needed)
```

## Step 1: Build the genome index (bowtie 1)

```bash
# miRDeep2 uses bowtie 1, NOT bowtie2
bowtie-build genome.fa genome_index
```

## Step 2: Map reads with mapper.pl

```bash
mapper.pl reads.fastq \
    -e -h -i -j \
    -k TGGAATTCTCGGGTGCCAAGG \
    -l 18 -m \
    -p genome_index \
    -s reads_collapsed.fa \
    -t reads_vs_genome.arf \
    -v

# -e: input is FASTQ   -h: parse to FASTA   -i: convert RNA to DNA
# -j: remove reads with non-ACGTN   -k: clip 3' adapter   -l 18: discard < 18 nt
# -m: collapse identical reads   -p: bowtie index   -s/-t: collapsed FASTA + ARF
```

## Step 3: Prepare miRBase references

```bash
# miRBase distributes RNA (U) sequences; miRDeep2 needs DNA and no whitespace.
# Pin the miRBase version - accessions and sequences change between releases.
wget https://www.mirbase.org/download/mature.fa
wget https://www.mirbase.org/download/hairpin.fa

# Same-species mature + hairpin (here human, hsa) and a related species for conservation
grep -A1 '>hsa-' mature.fa | grep -v '^--$' > mature_hsa.fa
grep -A1 '>hsa-' hairpin.fa | grep -v '^--$' > hairpin_hsa.fa
grep -A1 '>mmu-' mature.fa | grep -v '^--$' > mature_mmu.fa
# Convert U->T and strip spaces if the tool's extract_miRNAs.pl is not used:
# sed '/^>/!s/U/T/g; /^>/!s/u/t/g' in.fa
```

## Step 4: Run discovery with miRDeep2.pl

```bash
miRDeep2.pl \
    reads_collapsed.fa \
    genome.fa \
    reads_vs_genome.arf \
    mature_hsa.fa \
    mature_mmu.fa \
    hairpin_hsa.fa \
    -t Human \
    2> report.log

# Positional args (ORDER is fixed): collapsed reads, genome, ARF,
#   same-species mature, other-species mature (or 'none'), same-species hairpin
# -t: species for miRBase labelling
```

## Step 5: Known-miRNA quantification only (skip discovery)

```bash
quantifier.pl \
    -p hairpin_hsa.fa \
    -m mature_hsa.fa \
    -r reads_collapsed.fa \
    -t hsa
# Output: miRNAs_expressed_all_samples_*.csv
# Note: quantifier.pl and miRDeep2.pl counts can differ (different mapping logic)
```

## Output files

| File | Description |
|------|-------------|
| result_*.csv | Ranked candidates: miRDeep2 score, randfold p, mature/star, miRBase match, estimated probability TP |
| result_*.html | Interactive report with read-stack and structure plots |
| miRNAs_expressed_all_samples_*.csv | Known-miRNA expression matrix |
| mirdeep_runs/, expression_analyses/, pdfs_*/ | Intermediate read-stack alignments (.mrd) and structures |

## Reading and filtering results

```python
import pandas as pd

def parse_mirdeep2_results(csv_path, score_cutoff):
    # score_cutoff is NOT universal: choose it from survey.pl signal-to-noise / FDR,
    # then report the value. There is no fixed 'score > 10' rule.
    df = pd.read_csv(csv_path, sep='\t', skiprows=1)
    return df[df['miRDeep2 score'] >= score_cutoff]

def reject_structured_rna_false_positives(candidates, trna_rrna_bed):
    # The classic miRDeep2 false positive is a tRNA/rRNA fragment hairpin.
    # Require: (a) no overlap with tRNA/rRNA/snoRNA loci, (b) some star-arm read
    # support, (c) reproducibility across replicates, before trusting a novel call.
    return candidates  # intersect coordinates against trna_rrna_bed with bedtools upstream
```

## Calling a novel miRNA real: the community criteria

A miRDeep2 score is a prefilter, not a verdict. A genuine novel miRNA must satisfy the community annotation criteria (Ambros 2003; MirGeneDB), and the deliverable should be a per-candidate criteria table, not a score-ranked list:
- CONSISTENT 5' processing of BOTH the mature and star arms across reads - this 5'-end homogeneity is the single most discriminating signal (a precise 5' end is what defines the seed; degradation gives smeared ends).
- A mature/star duplex with the ~2-nt 3' overhang geometry of Dicer cleavage.
- Star-arm read support (real miRNAs usually show some passenger reads).
- A ~22-nt mature length and a hairpin without large internal loops/bulges.
- Conservation or Dicer/Drosha-dependence (loss of signal on knockdown), and reproducibility across replicates.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "novel miRNAs" cluster at tRNA/rRNA loci | Structured-RNA fragments fold into scoring hairpins | Intersect candidates against GtRNAdb/rRNA annotations and discard overlaps |
| mapper.pl fails or maps almost nothing | bowtie2 index supplied, or genome not indexed with bowtie 1 | Rebuild with `bowtie-build` (bowtie 1); confirm reads were adapter-trimmed |
| miRDeep2.pl errors on the reference FASTA | miRBase U-containing or whitespace-laden sequences | Convert U->T and strip header whitespace, or use the bundled extraction script |
| Treating score > 10 as truth | No universal cutoff exists | Use survey.pl signal-to-noise/FDR to set and report a cutoff |
| Very few known miRNAs detected | Wrong species `-t`, or reads not collapsed (`_xN`) | Set the correct species code; collapse reads in mapper.pl (`-m`) |
| Novel call has no star-arm reads | Real miRNAs usually show some passenger reads | Down-weight single-arm candidates; require duplex evidence |

## Related Skills

- smrna-preprocessing - Adapter trimming and read collapsing before mapping
- mirge3-analysis - Faster known-miRNA + isomiR quantification when discovery is not needed
- differential-mirna - Differential expression of the resulting count matrix
- trf-pirna-profiling - For tRF/piRNA biology, which would otherwise appear as miRDeep2 false positives
- genome-annotation/ncrna-annotation - Annotating tRNA/rRNA/snoRNA loci to filter false positives

## References

- Friedländer MR, Mackowiak SD, Li N, Chen W, Rajewsky N. 2012. miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. *Nucleic Acids Res* 40:37-52. doi:10.1093/nar/gkr688
- Friedländer MR, Chen W, Adamidi C, et al. 2008. Discovering microRNAs from deep sequencing data using miRDeep. *Nat Biotechnol* 26:407-415. doi:10.1038/nbt1394
- Bonnet E, Wuyts J, Rouzé P, Van de Peer Y. 2004. Evidence that microRNA precursors, unlike other non-coding RNAs, have lower folding free energies than random sequences. *Bioinformatics* 20:2911-2917. doi:10.1093/bioinformatics/bth374
- Kozomara A, Birgaoanu M, Griffiths-Jones S. 2019. miRBase: from microRNA sequences to function. *Nucleic Acids Res* 47:D155-D162. doi:10.1093/nar/gky1141
- Fromm B, Domanska D, Høye E, et al. 2020. MirGeneDB 2.0: the metazoan microRNA complement. *Nucleic Acids Res* 48:D1172-D1180. doi:10.1093/nar/gkz885
- Ambros V, Bartel B, Bartel DP, et al. 2003. A uniform system for microRNA annotation. *RNA* 9:277-279. doi:10.1261/rna.2183803
