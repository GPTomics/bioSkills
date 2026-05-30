# Interval Arithmetic - Usage Guide

## Overview
Interval arithmetic is the set-algebra layer of genomic analysis: intersect (overlaps), subtract (set difference), merge/cluster (collapse), complement (gaps), map/groupby (value transfer and aggregation), and multiinter/unionbedg (multi-sample stacks), plus jaccard/fisher as quick similarity/screen mechanics. The operations are exact and deterministic - the failures come from silent preconditions: forgetting to `sort` before `merge`, trusting `-sorted` without `-g`, omitting `-split` on spliced features, and choosing the wrong overlap fraction. This skill owns the mechanics; whether an overlap count is *more than chance* belongs to overlap-significance.

## Prerequisites
```bash
# bedtools (the reference implementation; pybedtools needs this binary on PATH)
conda install -c bioconda bedtools

# Python engines (choose what fits the workflow)
pip install pybedtools pyranges bioframe
```

## Quick Start
Tell your AI agent what you want to do:
- "Find which of my peaks overlap promoter regions"
- "Remove ENCODE blacklist regions from my peak calls"
- "Merge overlapping peaks from my three replicates into a consensus"
- "Transfer mean signal from a bedGraph onto each gene"
- "Which of these two SV call sets are the same event by 50% reciprocal overlap?"

## Example Prompts

### Finding and counting overlaps
> "Report each peak in peaks.bed once if it overlaps any gene in genes.bed, and separately list the peaks that overlap nothing."
> "Annotate every peak with the genes it overlaps using a -wa -wb join, then count B hits per peak."
> "Intersect my spliced transcript models against the peak BED at exon resolution (use -split) so intronic positions don't count."

### Removing and combining regions
> "Subtract exons from gene bodies to get introns, keeping fragments."
> "Concatenate my three replicate peak files, sort them, and merge peaks within 100 bp into a consensus peakset."
> "Drop any peak entirely if any part of it overlaps the blacklist."

### Value transfer and aggregation
> "Map the mean of column 4 from my bedGraph onto each gene interval (sort both first)."
> "Intersect peaks against genes with -wo, then groupby the gene columns and sum the overlap bp."

### Multi-sample and similarity
> "Build a presence/absence matrix of which of my three samples cover each sub-interval with multiinter."
> "Compute the Jaccard similarity between two peak sets for all-vs-all clustering (this is a similarity scalar, not a p-value)."
> "Run bedtools fisher as a quick screen, but tell me to route a real enrichment test to overlap-significance."

## What the Agent Will Do
1. Harmonize chromosome naming across inputs and the genome file (`chr1` vs `1`).
2. Sort every input with the same chromosome order before `merge`/`cluster`/`map`/`groupby`, and pass `-g genome.txt` whenever `-sorted`, `complement`, or `shuffle` is involved.
3. Add `-split` when an operand is BED12 or a spliced BAM and exon-level truth matters.
4. Select the overlap-output flag (`-u`/`-v`/`-c`/`-wa -wb`/`-loj`/`-wo`/`-wao`) and overlap fraction (`-f`/`-F`/`-r`) that match the question.
5. Execute the bedtools / pybedtools / pyranges / bioframe operation and save results.
6. For overlap *significance*, hand off to overlap-significance rather than reporting a raw count.

## Choosing the Right Operation

| Goal | Operation | Notes |
|------|-----------|-------|
| Which A intervals overlap B? | `intersect -u` | whole A once, no per-hit duplication |
| Which A intervals don't overlap B? | `intersect -v` | feature-level set difference |
| Annotate A with its B hits | `intersect -wa -wb` / `-loj` | join; `-loj` keeps A with no hit |
| Per-A count of B hits | `intersect -c` | 0 for no overlap |
| Remove B portions from A | `subtract` (`-A` for all-or-nothing) | A can fragment without `-A` |
| Collapse overlapping intervals | `sort \| merge` (`-d N`) | sort is mandatory |
| Gaps in coverage | `complement -g genome.txt` | genome file required |
| Transfer/aggregate B values onto A | `map -c -o` | both inputs sorted |
| Same event (SV/CNV)? | `intersect -f 0.5 -r` | 50% reciprocal convention |
| Is the overlap more than chance? | -> overlap-significance | needs a matched permutation null |

## Tips
- Sort each file once and reuse it; pipe operations (`bedtools sort -i in.bed | bedtools merge`).
- Whenever you use `-sorted`, also pass `-g genome.txt` - it turns a silent wrong answer (chromosome-order mismatch) into a loud crash.
- A surprisingly *small* overlap result usually means a sort/order problem, not biology - check sorting before interpreting.
- The output flags (`-wa`, `-wb`, `-wo`, `-loj`...) change only what is printed, not what overlaps; pick by the columns you need downstream.
- `-f` thresholds the fraction of A, `-F` the fraction of B - threshold the smaller set; the default is 1 bp (`-f 1e-9`), rarely the right biology.
- pybedtools writes many temp files; call `pybedtools.cleanup()` at the end of a script.
- pyranges 0.x and 1.0 differ - check `pyranges.__version__` before pasting idioms; pyranges/bioframe need no bedtools binary.
- jaccard is a similarity scalar and fisher is a weak screen; for a real enrichment p-value use overlap-significance.

## Resources
- [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
- [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html)
- [pybedtools documentation](https://daler.github.io/pybedtools/)
- [pyranges documentation](https://pyranges.readthedocs.io/)
- [bioframe documentation](https://bioframe.readthedocs.io/)

## Related Skills

- bed-file-basics - BED format, coordinate systems, and the conversions this skill depends on
- overlap-significance - Whether an overlap count exceeds a matched null (permutation, GAT/regioneR/LOLA/GREAT)
- proximity-operations - closest, window, flank, slop for adjacency rather than membership
- coverage-analysis - per-base depth and bedGraph signal feeding map/unionbedg
- gtf-gff-handling - exon/feature models whose `-split` behavior this skill depends on
- chip-seq/peak-calling - source of the peak BED files these operations consume
- atac-seq/consensus-peakset - replicate merge via merge/multiinter
