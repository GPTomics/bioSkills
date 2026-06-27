---
name: bio-small-rna-seq-trf-pirna-profiling
description: Profiles non-miRNA small RNAs - tRNA-derived fragments (tRFs/tsRNAs), piRNAs, and rRNA/snoRNA-derived species - with MINTmap, unitas, SPORTS, and proTRAC. Use when annotating all small-RNA classes in a library; quantifying tRFs at locus resolution where tRNA loci are redundant (exclusive vs ambiguous); testing the piRNA ping-pong signature; deciding whether a species is a processed functional RNA or a degradation fragment; or judging whether the prep could even capture 5'-OH/cyclic-phosphate classes.
tool_type: mixed
primary_tool: MINTmap
---

## Version Compatibility

Reference examples tested with: MINTmap 2.0+, unitas 1.7+, SPORTS1.0, proTRAC 2.4+, Python 3.10+ (numpy 1.26+, pandas 2.2+)

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# tRF and piRNA Profiling

**"Profile the tRFs and piRNAs in my small RNA-seq"** -> Annotate every small-RNA class, quantify tRNA-derived fragments at locus resolution, and test whether a piRNA population is real and active.
- CLI: `MINTmap` (tRFs), `unitas` / `SPORTS1.0` (all classes), `proTRAC` (piRNA clusters)

## The governing principle: small RNA-seq is a size cut over 8+ classes, and the kit decides what is captured

A small-RNA library is a ~18-40 nt size selection that pools miRNAs, piRNAs, endo-siRNAs, tRNA-derived fragments, rRNA-derived fragments, snoRNA-derived RNAs, and Y-RNA fragments - treating the output as "a miRNA dataset" is the field's most common error. Two facts dominate the analysis. First, end chemistry decides capture: standard TruSeq-style ligation requires a 5'-monophosphate and a 3'-OH, so 5'-OH and 2',3'-cyclic-phosphate species (angiogenin-cleaved tRNA halves, many tRFs and rRFs) are SILENTLY ABSENT, not lowly expressed - their absence in a TruSeq library is an assay artifact until proven otherwise, and capturing them needs T4 PNK pre-treatment (cP-RNA-seq / PANDORA-seq). Second, detection is not function: any abundant structured RNA sheds breakdown products into the 18-40 nt window, so high read count proves nothing. Functionality must be EARNED by precise reproducible ends, strand bias, phasing, the piRNA ping-pong signature, or AGO/PIWI loading - and database membership (piRBase, MINTbase) is annotation, not proof.

Multimapping is the other defining hazard: tRNA and piRNA loci are highly redundant (many genomic copies), so a read often cannot be assigned to one locus. This is why tRF tools report EXCLUSIVE versus AMBIGUOUS counts, and why piRNAs are quantified at the cluster/family level rather than per sequence.

Two end-chemistry and biogenesis facts change piRNA conclusions specifically. piRNAs (and plant miRNAs) carry a 3' 2'-O-methyl (HENMT1) that suppresses standard 3'-adapter ligation, so they are systematically UNDER-counted - low piRNA yield can be a 3'-end-chemistry artifact, not low abundance. And the ping-pong signature evidences the SECONDARY (slicer-driven, transposon) pathway only: primary piRNAs are PHASED (1U, Zucchini-dependent trail biogenesis), not ping-pong, and adult mammalian testis is >95% pachytene piRNAs that are repeat-depleted and largely non-transposon. A flat ping-pong z-score therefore does NOT mean "no piRNAs" - test phasing as well.

## Decision: which tool for which class

| Goal | Tool | Why |
|------|------|-----|
| tRFs/tsRNAs at locus resolution | MINTmap | deterministic, mapping-free; separates exclusive vs ambiguous tRF reads; MINTplate license-plate IDs |
| All small-RNA classes annotated hierarchically | unitas or sRNAbench | universal annotation (miRNA/piRNA/tRF/rRF/snoRNA) across ~800 species |
| tRF + rRF-centric biology (sperm/stress/aging) | SPORTS1.0 | finer tRF/rRF classification than miRNA-centric tools |
| piRNA clusters and ping-pong | proTRAC (+ a ping-pong test) | probabilistic cluster detection from mapped reads |
| Plant small RNAs (24-nt siRNA, phasiRNA) | ShortStack | DicerCall, phasing/PHAS-locus detection; animal tools misperform on plants |
| Known miRNAs only | mirge3-analysis | wrong tool for tRFs/piRNAs; miRNA-specific |

## tRF quantification with MINTmap

```bash
# MINTmap maps trimmed reads against a tRNA-space lookup and emits two tables:
# EXCLUSIVE tRFs (reads that map only within tRNA space) and AMBIGUOUS tRFs.
# Trust exclusive counts; ambiguous reads are shared with non-tRNA loci.
MINTmap -f trimmed.fastq -p sample_out
# Outputs: sample_out-MINTmap_v2-exclusive-tRFs.expression.txt
#          sample_out-MINTmap_v2-ambiguous-tRFs.expression.txt
# tRF type (tRF-5/tRF-3/tRF-1/i-tRF/tRNA-half) and the source tRNA are reported per row.
```

The tRF subtype carries a biogenesis tell: tRF-1 comes from the pre-tRNA 3' trailer (RNase Z/ELAC2, ending at the Pol III terminator), tRF-3 includes the post-transcriptional CCA (a marker of mature-tRNA origin), and tRNA halves are angiogenin-cleaved and stress-induced. Mitochondrially-encoded tRFs (mse-tRFs) are lost or misassigned if reads are mapped only to the nuclear genome.

## All-class annotation with unitas

```bash
# Hierarchical annotation: each read assigned to the first matching class.
# Reading the class composition is the first interpretation step.
unitas -input trimmed.fastq -species human
# Output: a UNITAS folder with per-class read fractions (miRNA / piRNA / tRF / rRF / snoRNA / ...)
```

## piRNA cluster detection with proTRAC

```bash
# Map reads (e.g. with sRNAmapper/bowtie), then call clusters probabilistically.
proTRAC_2.4.4.pl -genome genome.fa -map reads.map -format SAM
# A real primary-piRNA cluster shows strand asymmetry, 1U bias, and phased 3' ends.
```

## Test the ping-pong signature (is this an active piRNA pathway?)

**Goal:** Decide whether a putative piRNA population shows the slicer-driven ping-pong amplification signature.

**Approach:** For sense/antisense read pairs, count 5'-5' overlaps; an active pathway shows a sharp excess at exactly 10 nt (with 1U on primary and 10A on secondary piRNAs).

```python
import numpy as np
from collections import defaultdict

def ping_pong_zscore(plus_5p, minus_5p, max_overlap=30):
    # plus_5p / minus_5p: dict mapping genomic 5' coordinate -> read count, per strand.
    # A sense read at position i and an antisense read whose 5' end sits at i+overlap-1
    # overlap by 'overlap' nt at their 5' ends. Score the overlap histogram; a 10-nt
    # spike (z >> 0) is the ping-pong signature, evidence of an active piRNA pathway.
    hist = np.zeros(max_overlap + 1)
    for pos, n in plus_5p.items():
        for overlap in range(1, max_overlap + 1):
            partner = pos + overlap - 1
            if partner in minus_5p:
                hist[overlap] += n * minus_5p[partner]
    others = np.concatenate([hist[1:10], hist[11:]])
    z10 = (hist[10] - others.mean()) / (others.std() + 1e-9)
    return hist, z10


def phasing_zscore(same_strand_5p, period=27, max_dist=60):
    # Primary piRNAs are produced head-to-tail, so adjacent SAME-strand 5' ends are
    # spaced ~one piRNA length apart. Score the 5'-to-5' distance histogram: a peak at
    # the modal piRNA length (~26-28 nt) is the phasing signal (the primary-pathway
    # complement to ping-pong; proTRAC reports it natively). Test BOTH, not just ping-pong.
    pos = sorted(same_strand_5p)
    hist = np.zeros(max_dist + 1)
    for a in pos:
        for d in range(1, max_dist + 1):
            if (a + d) in same_strand_5p:
                hist[d] += same_strand_5p[a] * same_strand_5p[a + d]
    others = np.delete(hist[1:], period - 1)
    zp = (hist[period] - others.mean()) / (others.std() + 1e-9)
    return hist, zp
```

## Separate functional species from degradation

**Goal:** Avoid reporting random tRNA/rRNA breakdown as regulatory small RNAs.

**Approach:** Require end precision (a sharp, reproducible 5' terminus across replicates), strand bias, and class-appropriate length modality before trusting a non-miRNA species; rRFs are the hardest case because rRNA is so abundant that even tiny decay yields huge counts.

```python
def end_precision(read_5p_positions):
    # read_5p_positions: list of 5' coordinates for reads at a candidate locus.
    # A processed species has a dominant 5' end; random decay gives a smeared
    # distribution. Fraction of reads at the modal 5' end is a cheap discriminator.
    from collections import Counter
    c = Counter(read_5p_positions)
    return max(c.values()) / sum(c.values())   # near 1.0 = precise; low = decay-like
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "No tRNA halves / no tRFs" from a TruSeq library | 5'-OH / 2',3'-cyclic-phosphate ends are not ligatable | Absence is an assay artifact; use T4 PNK prep (cP-RNA-seq/PANDORA-seq) to capture them |
| tRF counts unstable across samples | Counting ambiguous (multimapped) tRF reads | Use MINTmap EXCLUSIVE counts; report ambiguous separately |
| Abundant "piRNAs" in a somatic/plasma sample | piRBase match by chance (often tRFs or Y-RNA fragments) | piRNAs are scarce in soma; require ping-pong/phasing, not database membership |
| Huge rsRNA counts called a discovery | rRNA is so abundant that minor decay dominates | Demand end precision and reproducibility before treating an rRF as a species |
| Plant data gives few "miRNAs" | Animal tools misread 24-nt siRNA / phasiRNA biology | Use ShortStack with DicerCall and phasing |
| Ping-pong test is flat | No active SECONDARY pathway, or primary/pachytene piRNAs (which are phased, not ping-pong) | Test phasing too; flat ping-pong does not mean no piRNAs (testis is >95% pachytene/phased) |
| Low piRNA yield despite a capable prep | 3' 2'-O-methyl blocks standard adapter ligation | Treat low piRNA counts as a possible end-chemistry artifact; use periodate/2'-OMe-tolerant chemistry |

## Related Skills

- smrna-preprocessing - Wider size windows and end chemistry that determine class capture
- mirdeep2-analysis - tRF/rRF stacks are miRDeep2 false positives; this skill targets them instead
- mirge3-analysis - Known miRNAs (and a basic tRF module)
- differential-mirna - The same count-based DE framework applies to tRF/piRNA matrices
- genome-annotation/ncrna-annotation - tRNA/rRNA/snoRNA locus annotation underlying these tools

## References

- Loher P, Telonis AG, Rigoutsos I. 2017. MINTmap: fast and exhaustive profiling of nuclear and mitochondrial tRNA fragments from short RNA-seq data. *Sci Rep* 7:41184. doi:10.1038/srep41184
- Pliatsika V, Loher P, Magee R, et al. 2018. MINTbase v2.0: a comprehensive database for tRNA-derived fragments. *Nucleic Acids Res* 46:D152-D159. doi:10.1093/nar/gkx1075
- Gebert D, Hewel C, Rosenkranz D. 2017. unitas: the universal tool for annotation of small RNAs. *BMC Genomics* 18:644. doi:10.1186/s12864-017-4031-9
- Shi J, Ko EA, Sanders KM, Chen Q, Zhou T. 2018. SPORTS1.0: a tool for annotating and profiling non-coding RNAs optimized for rRNA- and tRNA-derived small RNAs. *Genomics Proteomics Bioinformatics* 16:144-151. doi:10.1016/j.gpb.2018.04.004
- Rosenkranz D, Zischler H. 2012. proTRAC - a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics* 13:5. doi:10.1186/1471-2105-13-5
- Brennecke J, Aravin AA, Stark A, et al. 2007. Discrete small RNA-generating loci as master regulators of transposon activity in Drosophila. *Cell* 128:1089-1103. doi:10.1016/j.cell.2007.01.043
- Shi J, Zhang Y, Tan D, et al. 2021. PANDORA-seq expands the repertoire of regulatory small RNAs by overcoming RNA modifications. *Nat Cell Biol* 23:424-436. doi:10.1038/s41556-021-00652-7
