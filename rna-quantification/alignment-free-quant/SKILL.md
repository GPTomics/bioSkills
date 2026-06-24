---
name: bio-rna-quantification-alignment-free-quant
description: Quantify transcript expression from FASTQ with Salmon (selective alignment) or kallisto (pseudoalignment), bypassing genome mapping. Use when quantifying RNA-seq without alignment, deciding whether a decoy-aware index is required, detecting and verifying library strandedness, enabling GC and sequence bias correction, or choosing whether to generate inferential replicates (bootstraps/Gibbs) for transcript-level downstream testing.
tool_type: cli
primary_tool: salmon
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: Salmon 1.10+, kallisto 0.50+, fastp 0.23+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment-Free Quantification

**"Quantify gene expression without aligning to the genome"** -> Estimate transcript abundances directly from FASTQ reads by mapping to a transcriptome, then resolving multi-mapping reads (paralogs, shared isoform sequence) with an EM/variational model.
- CLI: `salmon quant -i index -l A -1 R1.fq.gz -2 R2.fq.gz -o quant/`, `kallisto quant -i index -o out R1.fq.gz R2.fq.gz`

These tools do not just count reads; they run probabilistic inference. Most fragments are compatible with several transcripts, so the abundance of each transcript is a latent quantity estimated by EM (Salmon offline phase, RSEM) or variational Bayes (Salmon default). The two load-bearing decisions that determine whether the numbers are trustworthy are the index (decoy-aware or not) and the library type. Get either wrong and the run completes silently with biased output.

## Decision 1: the decoy-aware index is not optional

A transcriptome-only index has no target for reads that originate from introns, unannotated transcription, or intergenic DNA. Those reads do not vanish; they are force-fit onto whatever transcript shares enough sequence, inflating its count. Adding the genome as a set of decoy sequences fixes this: if a fragment aligns better to a decoy than to any transcript, all of its mappings are discarded rather than misassigned. Always build the decoy-aware index when the genome is available (Srivastava et al. 2020).

```bash
# Decoy names = every genome sequence header
grep "^>" genome.fa | cut -d " " -f 1 | sed 's/>//g' > decoys.txt

# gentrome = transcripts FIRST, then genome
cat transcripts.fa genome.fa > gentrome.fa

# Build (k=31 is right for reads >=75 bp; lower to ~23-25 for ~50 bp reads)
salmon index -t gentrome.fa -d decoys.txt -i salmon_index -k 31 -p 8
```

`--validateMappings` is deprecated and has no effect: selective alignment has been the default since Salmon 1.0.0. Do not pass it. This is the accuracy difference from pure pseudoalignment (kallisto): pseudoalignment commits a read to its compatible transcript set without scoring base-level mismatches, so it over-assigns intron-, pseudogene-, and error-derived reads, whereas selective alignment computes an actual alignment score around each candidate and drops low-scoring spurious mappings. kallisto has no equivalent decoy index; its closest analog is `--d-list genome.fa` (a distinguishing-k-mer filter, different mechanism) to partially compensate.

## Decision 2: library type drives strandedness

```bash
# Salmon: auto-detect, then VERIFY
salmon quant -i salmon_index -l A \
    -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    -o sample_quant --gcBias --seqBias -p 8

# Single-end
salmon quant -i salmon_index -l A -r sample.fastq.gz -o sample_quant -p 8
```

`-l A` auto-detects the type and writes the inferred format to `lib_format_counts.json`. Inspect it: a library with weak strand signal can be miscalled, and the wrong type makes correctly oriented fragments incompatible, collapsing or randomizing abundances. The dominant modern chemistry (dUTP / Illumina TruSeq Stranded / NEBNext Directional) is `ISR` for Salmon and maps to `featureCounts -s 2` (reverse). Map: unstranded `IU` <-> `-s 0`; forward `ISF`/`SF` <-> `-s 1`; reverse `ISR`/`SR` <-> `-s 2`. 3'-tag protocols carry their own strandedness (Lexogen QuantSeq FWD is forward, `SF`/`-s 1`; QuantSeq REV is reverse), so rely on `-l A` and `lib_format_counts.json` rather than assuming reverse.

## Bias correction

`--gcBias` and `--seqBias` learn sample-specific fragment-GC and random-hexamer-priming biases and reweight the read-to-transcript probabilities. They cost little and protect against the case that produces false positives: when library-prep batch is confounded with the biological condition, an uncorrected GC bias becomes a condition effect. Enable both as a near-default; reserve `--posBias` for degraded or visibly 3'-biased libraries.

## Inferential replicates: only when transcript-level uncertainty matters

Transcripts that share sequence are not individually identifiable, so their point estimates carry inferential (quantification) uncertainty on top of biological variance. This uncertainty cancels when isoforms are summed to the gene, so gene-level DESeq2/edgeR via tximport needs no replicates. It does not cancel at the transcript level: differential transcript expression (DTE) and usage (DTU) require propagating it.

```bash
# Generate inferential replicates ONLY for transcript-level downstream testing
salmon quant -i salmon_index -l A --gcBias --seqBias \
    --numGibbsSamples 20 \
    -1 R1.fq.gz -2 R2.fq.gz -o sample_quant -p 8

# kallisto bootstraps (for sleuth)
kallisto quant -i kallisto_index -o sample_quant -b 100 R1.fq.gz R2.fq.gz
```

Use `--numGibbsSamples 20` (or `--numBootstraps 30`) for Salmon; ~100 bootstraps for kallisto. Downstream consumers: swish (alternative-splicing/isoform-switching), sleuth (expression-matrix/counts-ingest), edgeR catchSalmon (differential-expression/edger-basics). Do not pay this cost for gene-level work.

## kallisto Workflow

```bash
kallisto index -i kallisto_index transcripts.fa

# Paired-end learns the fragment-length distribution from mate distances
kallisto quant -i kallisto_index -o sample_quant R1.fastq.gz R2.fastq.gz

# Single-end CANNOT observe fragment length -> must supply mean (-l) and sd (-s),
# which set effective lengths and therefore TPM; wrong values bias every TPM
kallisto quant -i kallisto_index -o sample_quant --single -l 200 -s 20 sample.fastq.gz
```

## Output

`quant.sf` (Salmon) columns: `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads`. kallisto `abundance.tsv`: `target_id`, `length`, `eff_length`, `est_counts`, `tpm`; `abundance.h5` holds bootstraps. EffectiveLength is the transcript length convolved with the fragment-length distribution; a transcript shorter than the mean fragment length has a tiny, unstable effective length, so its TPM is hypersensitive to small count changes. Treat short-transcript TPMs with suspicion and filter low-count features before testing.

Import the estimated counts (`NumReads`/`est_counts`), not TPM, into DESeq2/edgeR via tximport, which adds the length offset (rna-quantification/tximport-workflow). TPM is a within-sample proportion and is invalid for cross-sample differential expression.

## Salmon vs kallisto vs RSEM

| Tool | Speed | Accuracy | Best when |
|------|-------|----------|-----------|
| Salmon (selective alignment + decoy) | Fast | Highest among lightweight | Default for bulk RNA-seq; decoy absorbs intron/pseudogene reads |
| kallisto | Fastest | Excellent | Speed-critical or sleuth-based DTE; add `--d-list` to mitigate intron over-assignment |
| RSEM | Slowest (needs a separate aligner) | Reference standard | Defensible benchmark accuracy; runs on a transcriptome BAM (STAR `--quantMode TranscriptomeSAM`) |

Methodology evolves; confirm current defaults against the Salmon and kallisto docs before relying on a flag.

## Combine TPM / counts for inspection

```python
import pandas as pd
from pathlib import Path

samples = ['sample1', 'sample2', 'sample3']
tpm, counts = {}, {}
for s in samples:
    df = pd.read_csv(Path(f'{s}_quant/quant.sf'), sep='\t', index_col=0)  # kallisto: abundance.tsv
    tpm[s], counts[s] = df['TPM'], df['NumReads']                          # kallisto: tpm, est_counts
pd.DataFrame(tpm).to_csv('tpm_matrix.csv')
pd.DataFrame(counts).to_csv('counts_matrix.csv')
```

## Quality Checks

```bash
grep "Mapping rate" sample_quant/logs/salmon_quant.log   # expect > ~70% for a matched reference
cat sample_quant/lib_format_counts.json                  # confirm one consistent library type
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Low mapping rate (<50%) | Wrong/old transcriptome version, contamination, or no decoy | First confirm a decoy-aware index and a transcriptome matching the GTF release; then run FastQ Screen for contamination |
| One gene/transcript implausibly high | Intron or pseudogene reads forced onto it (transcriptome-only index) | Rebuild with genome decoys |
| Counts halve or look random for stranded data | Library type miscalled by `-l A` | Read `lib_format_counts.json`; for dUTP/TruSeq it should be `ISR` |
| Inconsistent library types across samples | Mixed library preps or a sample swap | Verify metadata; quantify suspect samples separately and compare |
| `swish: no inferential replicates found` downstream | Salmon/kallisto run without Gibbs/bootstraps | Re-run with `--numGibbsSamples 20` (or kallisto `-b 100`) |

## Related Skills

- rna-quantification/tximport-workflow - Import counts with the length offset for DESeq2/edgeR
- rna-quantification/featurecounts-counting - Alignment-based counting alternative
- read-qc/fastp-workflow - Upstream adapter/quality trimming
- alternative-splicing/isoform-switching - swish DTE/DTU using Salmon Gibbs samples
- expression-matrix/counts-ingest - sleuth on kallisto bootstraps
- differential-expression/edger-basics - catchSalmon transcript-level DTE
- differential-expression/deseq2-basics - Gene-level downstream analysis

## References

- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods 14(4):417-419. doi:10.1038/nmeth.4197
- Bray NL, Pimentel H, Melsted P, Pachter L. 2016. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34(5):525-527. doi:10.1038/nbt.3519
- Srivastava A, Malik L, Sarkar H, et al. 2020. Alignment and mapping methodology influence transcript abundance estimation. Genome Biol 21:239. doi:10.1186/s13059-020-02151-8
- Love MI, Hogenesch JB, Irizarry RA. 2016. Modeling of RNA-seq fragment sequence bias reduces systematic errors in transcript abundance estimation. Nat Biotechnol 34(12):1287-1291. doi:10.1038/nbt.3682
