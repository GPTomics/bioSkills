---
name: bio-small-rna-seq-smrna-preprocessing
description: Trims kit-specific 3' adapters, strips UMIs or 4N degenerate ends, size-selects, and collapses small RNA-seq reads (miRNA, piRNA, tRF) with cutadapt or fastp. Use when choosing the kit's 3' adapter; setting the size window (18-26 nt miRNA vs 24-32 nt piRNA); deciding whether a library carries a true UMI (QIAseq) versus a 4N debiasing spacer (NEXTflex); reading the read-length histogram to judge library quality; or deciding whether to collapse identical reads before mapping.
tool_type: cli
primary_tool: cutadapt
---

## Version Compatibility

Reference examples tested with: cutadapt 4.4+, fastp 0.23+, seqkit 2.6+, umi_tools 1.1+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Small RNA Preprocessing

**"Preprocess my small RNA-seq reads"** -> Remove the 3' adapter, remove any UMI or 4N degenerate bases, size-select to the target class window, and collapse identical reads to a counted FASTA before quantification or discovery.
- CLI: `cutadapt -a ADAPTER -m 18 -M 30 --discard-untrimmed` then class-specific UMI/4N handling and `seqkit rmdup -s`

## The governing principle: the insert IS its ends, so adapter handling decides everything

In small RNA-seq the molecule is shorter than the read (insert ~18-32 nt, read 50-75 nt), so the 3' adapter is sequenced through on EVERY real insert. A read with no adapter is therefore not a complete small RNA (the insert was too long, or it is an adapter dimer or junk), which inverts the genomic-DNA intuition: here `--discard-untrimmed` is the correct default, not an aggressive one. Because the insert is defined by its exact 5' and 3' ends, ligation bias never averages out the way fragmentation does in mRNA-seq: T4 RNA ligase captures some miRNA ends 10-100x more efficiently than others, so absolute, cross-miRNA abundance WITHIN a sample is not trustworthy (only the same miRNA compared ACROSS samples, where the per-sequence bias cancels, is reliable; Giraldez 2018). Preprocessing cannot fix ligation bias, but mishandling adapters, UMIs, or size windows manufactures artifacts on top of it.

The read-length histogram after trimming is the primary QC readout, not an afterthought: a sharp peak at 21-23 nt is a healthy miRNA library; a 26-32 nt peak is piRNA (expected in germline, suspicious in soma/plasma); a broad 30+ nt smear with no 22 nt peak is degradation or tRNA/rRNA-fragment contamination or failed size selection; a spike near insert length 0 is adapter dimer eating flowcell capacity. Read the histogram before trusting any downstream count. This degradation heuristic assumes a standard ligation library: for T4-PNK / PANDORA-seq / phospho-RNA-seq preps (which deliberately capture 5'-OH and cyclic-phosphate tRFs and rRFs) the heuristic INVERTS - broad ~18-35 nt tRF/rRF peaks are the expected signal, not contamination.

## Decision: how to handle the 3' end depends on the kit

| Kit | 3' adapter | Degenerate / UMI design | Preprocessing consequence |
|-----|-----------|--------------------------|----------------------------|
| Illumina TruSeq | TGGAATTCTCGGGTGCCAAGG | invariant ends (high ligation bias) | trim adapter only; do NOT PCR-dedup |
| NEBNext | AGATCGGAAGAGCACACGTCT | invariant ends | trim adapter only; do NOT PCR-dedup |
| NEXTflex (Bioo/PerkinElmer) | TGGAATTCTCGGGTGCCAAGG | 4 random nt on each adapter end (debiasing spacer) | trim adapter, then STRIP 4 nt from each insert end (`-u 4 -u -4`); the 4N is NOT a UMI, discard it |
| QIAseq miRNA | AACTGTAGGCACCATCAAT | true 12-nt UMI 3' of the adapter | EXTRACT the UMI (keep it), align, then UMI-dedup; never position-dedup |
| SMARTer / CATS (template-switching) | no ligation adapter | adds a 3' poly-A/tail, not a 4N or UMI | trim the 3' poly-tail, NOT a ligation adapter; different bias profile; ultra-low input |
| RealSeq (circularization) | single adapter, one ligation | sidesteps the two-junction ligation bias | low-input; one-ligation chemistry, not two |

The single most damaging error in this table is conflating the NEXTflex 4N debiasing spacer (only 4^4=256 combinations, must be DISCARDED) with the QIAseq 12-nt UMI (must be KEPT and used to separate PCR duplicates from biological duplicates). Using the 4N as a pseudo-UMI saturates instantly and undercounts abundant miRNAs.

## Adapter trimming with cutadapt

```bash
# Standard ligation-based small-RNA library (TruSeq adapter shown)
cutadapt \
    -a TGGAATTCTCGGGTGCCAAGG \
    -m 18 \
    -M 30 \
    -q 20 \
    --discard-untrimmed \
    -j 8 \
    -o trimmed.fastq.gz \
    input.fastq.gz

# -a: 3' adapter (cutadapt finds it even when only a prefix is sequenced)
# -m 18 / -M 30: keep the small-RNA window; -m drops adapter dimers (trim to ~0)
# -q 20: light 3' quality trim, applied BEFORE adapter removal (cutadapt orders it internally)
# --discard-untrimmed: a read with no adapter is not a complete small RNA
```

## Class-specific size windows

```bash
# miRNA-focused window (mature miRNAs cluster at 21-23 nt)
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 26 --discard-untrimmed -o mirna.fastq.gz input.fastq.gz

# piRNA / tRNA-half window (widen -M; do not clip the very class of interest)
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 24 -M 35 --discard-untrimmed -o pirna.fastq.gz input.fastq.gz
```

## Removing 4N degenerate bases (NEXTflex / high-definition adapters)

```bash
# ORDER MATTERS: trim the adapter FIRST, then strip the 4 random nt from each insert end.
# Stripping a fixed 4 nt before adapter removal would corrupt the adapter search.
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 30 --discard-untrimmed -o adapter_trimmed.fastq.gz input.fastq.gz
cutadapt -u 4 -u -4 -o final.fastq.gz adapter_trimmed.fastq.gz

# -u 4: remove 4 nt from the 5' end; -u -4: remove 4 nt from the 3' end (negative = 3')
```

## Extracting and using a true UMI (QIAseq)

```bash
# The 12-nt UMI sits immediately 3' of the QIAGEN adapter. Capture it into the read name,
# align, then collapse reads sharing sequence+position+UMI (PCR duplicates) but keep reads
# that differ in UMI (distinct biological molecules). Position-only dedup is WRONG for small RNA.
umi_tools extract --extract-method=regex \
    --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT)(?P<umi_1>.{12}).*' \
    -I input.fastq.gz -S umi_extracted.fastq.gz
# ... adapter-trim, map ...
umi_tools dedup --method=directional -I aligned.bam -S deduped.bam
# directional models 1-edit UMI sequencing errors; raw unique-UMI counting overcounts
```

## Using fastp as an alternative

```bash
fastp \
    --in1 input.fastq.gz \
    --out1 trimmed.fastq.gz \
    --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
    --length_required 18 \
    --length_limit 30 \
    --json report.json --html report.html

# --length_limit caps the small-RNA window; do NOT use fastp --dedup on small RNA
# (it is sequence-based and deletes real biological duplicates)
```

## Collapse identical reads to a counted FASTA

```bash
# seqkit rmdup -s only DEDUPLICATES identical sequences; it does NOT append the _xN
# count that miRDeep2 needs. Use it to shrink the file, but generate the counted FASTA
# with the awk/Python helper below. For a UMI library, collapse on sequence+UMI (or skip
# collapsing) so the UMI survives dedup.
seqkit rmdup -s trimmed.fastq.gz -o dedup.fasta
```

```python
import gzip
from collections import Counter

def collapse_reads(fastq_path, lo=18, hi=30):
    counts = Counter()
    with gzip.open(fastq_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            if lo <= len(seq) <= hi:
                counts[seq] += 1
    return counts

def write_collapsed_fasta(counts, output_path):
    # miRDeep2 reads the _xN suffix as the read count; preserve it
    with open(output_path, 'w') as f:
        for i, (seq, count) in enumerate(counts.most_common()):
            f.write(f'>seq_{i}_x{count}\n{seq}\n')
```

## QC and contamination gate with miRTrace

```bash
# Run miRTrace BEFORE quantifying. Beyond length/complexity, it reports the RNA-class
# composition (miRNA vs rRNA/tRNA/artifact) AND fingerprints clade-specific miRNAs to
# detect cross-species / reagent / sample-swap contamination (found in >7% of public
# datasets) that a good genome mapping rate hides. It has kit presets via --protocol.
mirtrace qc --species hsa --protocol illumina -o mirtrace_out *.fastq.gz
# Read: a miRNA-dominant composition is healthy; rRNA/tRNA-dominant means poor size
# selection, degraded input, or low real miRNA; a foreign-clade signal flags contamination.
```

For plasma/serum specifically, hemolysis is the dominant QC: red blood cells are loaded with miR-451a, so even slight hemolysis floods the sample with erythroid miRNAs and corrupts the circulating profile. Flag it with the miR-451a (RBC-enriched, rises with hemolysis) vs miR-23a-3p (hemolysis-insensitive) relationship - an elevated miR-451a fraction (or delta-Cq(miR-23a-3p - miR-451a) > ~7 by qPCR) marks a hemolyzed sample. Exclude or model hemolyzed samples before differential analysis (see differential-mirna). Use exogenous spike-ins (cel-miR-39) for low-biomass technical normalization.

## Read the length distribution as QC

```python
import gzip
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_length_distribution(fastq_path, out_png):
    lengths = Counter()
    with gzip.open(fastq_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                lengths[len(line.strip())] += 1
    xs = sorted(lengths)
    plt.bar(xs, [lengths[x] for x in xs])
    plt.axvspan(21, 23, color='green', alpha=0.15)  # healthy miRNA peak
    plt.xlabel('read length (nt)')
    plt.ylabel('count')
    plt.savefig(out_png)
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Almost nothing maps; reads ~50-75 nt | 3' adapter never removed (wrong sequence or step skipped) | Set the kit's exact adapter; reads must shrink to ~18-30 nt after trimming |
| Length histogram peaks at ~8 random nt over the real peak | 4N spacer not stripped after adapter removal | Add `cutadapt -u 4 -u -4` as a second pass |
| Abundant miRNAs look flat / undercounted | Position-based PCR dedup on non-UMI data, or 4N used as a UMI | Do not dedup without a real UMI; for QIAseq use `umi_tools dedup` |
| Broad 30+ nt smear, no 22 nt peak | Degraded input / tRNA-rRNA fragments / failed size selection | Inspect RNA quality (DV200, not RIN); rerun size selection; expect mostly non-miRNA classes |
| Huge spike at insert length ~0 | Adapter dimers (no-insert ligation), common at low input | `-m 18` discards them; report the dimer fraction as a library-quality flag |
| Cross-sample counts incomparable | Libraries built with different kits/protocols (bias is protocol-specific) | Never merge or compare counts across kits; rebuild with one protocol |
| High mapping rate but odd composition / foreign reads | Cross-species or reagent contamination that mapping rate hides | Run miRTrace clade fingerprinting; exclude or investigate contaminated samples |
| Good phospho/PANDORA library flagged as "degraded" | Standard length-histogram heuristic applied to a 5'-OH/cP-capture prep | Expect broad ~18-35 nt tRF/rRF peaks for these preps; the heuristic inverts |
| Plasma profile dominated by a few miRNAs across all samples | Hemolysis: red-cell miR-451a contamination | Flag with miR-451a:miR-23a-3p; exclude/model hemolyzed samples; spike-in normalize |

## Related Skills

- mirdeep2-analysis - Novel miRNA discovery; consumes collapsed reads
- mirge3-analysis - Fast known-miRNA + isomiR quantification; has its own trimming
- trf-pirna-profiling - tRF and piRNA profiling, where wider size windows and 5'-OH/cP end chemistry matter
- read-qc/adapter-trimming - General adapter trimming concepts and tool behavior
- read-qc/umi-processing - UMI extraction and deduplication mechanics

## References

- Martin M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 17:10-12. doi:10.14806/ej.17.1.200
- Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34:i884-i890. doi:10.1093/bioinformatics/bty560
- Giraldez MD, Spengler RM, Etheridge A, et al. 2018. Comprehensive multi-center assessment of small RNA-seq methods for quantitative miRNA profiling. *Nat Biotechnol* 36:746-757. doi:10.1038/nbt.4183
- Smith T, Heger A, Sudbery I. 2017. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Res* 27:491-499. doi:10.1101/gr.209601.116
- Sorefan K, Pais H, Hall AE, et al. 2012. Reducing ligation bias of small RNAs in libraries for next generation sequencing. *Silence* 3:4. doi:10.1186/1758-907X-3-4
- Kang W, Eldfjell Y, Fromm B, et al. 2018. miRTrace reveals the organismal origins of microRNA sequencing data. *Genome Biol* 19:213. doi:10.1186/s13059-018-1588-9
- Shi J, Zhang Y, Tan D, et al. 2021. PANDORA-seq expands the repertoire of regulatory small RNAs by overcoming RNA modifications. *Nat Cell Biol* 23:424-436. doi:10.1038/s41556-021-00652-7
