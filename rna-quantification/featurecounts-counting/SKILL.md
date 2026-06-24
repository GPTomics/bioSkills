---
name: bio-rna-quantification-featurecounts-counting
description: Count reads per gene from aligned BAM files using Subread featureCounts. Use when turning STAR/HISAT2 BAMs into a gene-level count matrix for DESeq2/edgeR, deciding library strandedness, handling paired-end fragment counting, choosing how to treat multi-mapping and multi-overlapping reads, or diagnosing a low assignment rate from the summary file.
tool_type: cli
primary_tool: featureCounts
---

## Version Compatibility

Reference examples tested with: Subread 2.0+, STAR 2.7.11+, HISAT2 2.2.1+, DESeq2 1.42+, edgeR 4.0+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# featureCounts Counting

**"Count reads per gene from my BAM files"** -> Assign each aligned read to at most one gene by overlap with a GTF, discarding ambiguous reads, to produce an integer gene-by-sample matrix for differential expression.
- CLI: `featureCounts -a genes.gtf -o counts.txt sample1.bam sample2.bam`

featureCounts is bookkeeping, not inference: it tallies reads to genes and discards anything ambiguous. That is correct for gene-level DE, where almost every read's gene of origin is unambiguous even when its isoform is not. The two settings that silently corrupt the matrix if wrong are strandedness (`-s`) and, for paired-end data, fragment counting (`--countReadPairs`).

## Basic Usage

```bash
# Multiple samples in one run -> a single aligned matrix (recommended)
featureCounts -a annotation.gtf -o counts.txt sample1.bam sample2.bam sample3.bam

# Defaults: -t exon -g gene_id (count reads over exons, aggregate by gene)
```

## Decision 1: strandedness (`-s`) is the load-bearing setting

| `-s` | Meaning | Read 1 |
|------|---------|--------|
| 0 | Unstranded | strand ignored |
| 1 | Forward stranded | read 1 is sense |
| 2 | Reverse stranded | read 1 is antisense; read 2 is sense |

The dominant chemistry, dUTP / Illumina TruSeq Stranded / NEBNext Directional, is reverse, `-s 2`. Setting the wrong strand does not error; it silently destroys the matrix. For a truly stranded library, the correct `-s` assigns ~80-90% of reads while the opposite setting collapses to ~5-20% (counting only antisense background). Do not trust the kit name; determine it empirically:

```bash
# Method A: RSeQC reports the strand pattern fractions
infer_experiment.py -r genes.bed -i sample.bam

# Method B: run all three and pick the one that maximizes Assigned in the .summary
for s in 0 1 2; do featureCounts -s $s -a annotation.gtf -o counts_s$s.txt sample.bam; done
```

If `-s 1` and `-s 2` give wildly different Assigned fractions, the data are stranded (use the higher); if both are roughly equal and about half of `-s 0`, the data are unstranded. STAR `--quantMode GeneCounts` provides a free cross-check (see below).

## Decision 2: paired-end fragment counting

```bash
# Subread >= 2.0.2: -p only declares paired input; --countReadPairs is REQUIRED to count fragments
featureCounts -p --countReadPairs -a annotation.gtf -o counts.txt *.bam

# Stricter: require both ends mapped, exclude chimeric/discordant pairs
featureCounts -p --countReadPairs -B -C -a annotation.gtf -o counts.txt *.bam
```

Omitting `--countReadPairs` on paired-end data counts each mate separately, roughly doubling counts and breaking the count model. `-B` requires both ends aligned; `-C` excludes pairs mapping across chromosomes or in the wrong orientation.

## Decision 3: multi-mapping and multi-overlap reads

```bash
# Default (recommended for gene-level DE): discard both -> uniquely, unambiguously assigned reads only
featureCounts -a annotation.gtf -o counts.txt *.bam

# Count multimappers fractionally (1/N) or fully (1 each) -- NOT recommended for DE
featureCounts -M --fraction -a annotation.gtf -o counts.txt *.bam
featureCounts -M -a annotation.gtf -o counts.txt *.bam

# Count reads overlapping >1 gene in all of them
featureCounts -O -a annotation.gtf -o counts.txt *.bam
```

Discarding multimappers is the right default for gene-level DE. `-M --fraction` looks principled but biases exactly the genes where resolution matters: a read truly from gene A that also maps to paralog A' is split 0.5/0.5, diluting both. This is the regime where alignment-free EM quantifiers (rna-quantification/alignment-free-quant) outperform featureCounts, because they reassign by full likelihood rather than a flat split.

## Quality and feature options

```bash
featureCounts -Q 10 -a annotation.gtf -o counts.txt *.bam      # min MAPQ (aligner-specific scale)
featureCounts --primary -a annotation.gtf -o counts.txt *.bam  # primary alignments only
featureCounts -t CDS -g gene_id -a annotation.gtf -o counts.txt *.bam  # count CDS instead of exon
```

`-Q` thresholds mapping quality, but MAPQ conventions are aligner-specific (STAR assigns 255 to unique reads, low values to multimappers), so confirm the scheme before choosing a cutoff. Do NOT add `--ignoreDup` for standard RNA-seq: high duplication is expected from highly expressed genes, and position-based deduplication discards real signal. Deduplicate only with UMIs. For exon-level usage testing (DEXSeq), use a flattened annotation rather than gene-level counting (alternative-splicing/isoform-switching).

## Output

```
counts.txt:           Geneid Chr Start End Strand Length sample1.bam sample2.bam ...
counts.txt.summary:   Status              sample1.bam  sample2.bam
                      Assigned            1523456      1678234
                      Unassigned_NoFeatures 234567     245678
```

Reading the `.summary` is the primary QC step. A good poly-A library assigns ~70-90% of mapped reads (rRNA-depletion libraries run lower).

| Dominant unassigned category | Likely cause | Action |
|------------------------------|--------------|--------|
| Unassigned_NoFeatures high | GTF/genome mismatch (chr naming `1` vs `chr1`, wrong release), DNA contamination | Match GTF release and chromosome naming to the BAM |
| Unassigned_MultiMapping high | rRNA carryover or repetitive content | Check rRNA depletion; inspect with FastQ Screen |
| Unassigned_Ambiguity high | Overlapping/nested gene models or wrong feature level | Expected in gene-dense regions; reconsider `-O`/feature type |
| Assigned low, others spread thin | Wrong strandedness | Re-test `-s` (see Decision 1) |

## STAR cross-check

If aligned with STAR `--quantMode GeneCounts`, `ReadsPerGene.out.tab` gives a free independent count: column 2 = unstranded (≈ `-s 0`), column 3 = forward (≈ `-s 1`), column 4 = reverse (≈ `-s 2`). The larger of columns 3 vs 4 reveals the strand directly, and the per-gene counts should track featureCounts at the matching `-s`.

## Extract the matrix

```bash
cut -f1,7- counts.txt | tail -n +2 > count_matrix.txt   # drop the 6 annotation columns
```

```python
import pandas as pd
counts = pd.read_csv('counts.txt', sep='\t', comment='#')
mat = counts.set_index('Geneid').iloc[:, 5:]
mat.columns = [c.replace('.bam', '') for c in mat.columns]
mat.to_csv('count_matrix.csv')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Assigned ~half of expected, no error | Wrong `-s`, or paired-end without `--countReadPairs` (double-counting) | Determine strand empirically; add `--countReadPairs` for paired-end |
| Near-zero counts for known genes | `gene_id` attribute or feature type mismatch with the GTF | Confirm `-t`/`-g` match the annotation; check the GTF attribute names |
| Counts much higher than read count | Paired-end mates counted separately | Add `-p --countReadPairs` |
| Inflated correlated paralog counts | `-M`/`-O` fractional counting enabled | Drop `-M`/`-O` for DE; use alignment-free EM for paralog-heavy genes |
| Low Assigned across all `-s` values | GTF does not match the aligned genome | Use the GTF release and contig names matching the alignment reference |

## Related Skills

- alignment-files/sam-bam-basics - Input BAM handling and filtering
- read-alignment/star-alignment - Producing BAMs and ReadsPerGene.out.tab strand cross-check
- genome-intervals/gtf-gff-handling - GTF/GFF annotation files
- rna-quantification/alignment-free-quant - EM-based alternative; better for multimappers/paralogs
- rna-quantification/count-matrix-qc - QC the resulting matrix before DE
- differential-expression/deseq2-basics - Gene-level DE from these counts

## References

- Liao Y, Smyth GK, Shi W. 2014. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30(7):923-930. doi:10.1093/bioinformatics/btt656
- Liao Y, Smyth GK, Shi W. 2013. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Res 41(10):e108. doi:10.1093/nar/gkt214
