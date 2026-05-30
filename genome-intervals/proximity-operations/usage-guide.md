# Proximity Operations - Usage Guide

## Overview
Proximity operations compute the geometric relationships between genomic intervals: the nearest feature and its signed, strand-aware distance (`bedtools closest`), all features within a search radius (`window`), strand-aware interval extension (`slop`), and the regions beside a feature (`flank`). The defining judgment in this skill is interpretive, not mechanical: `closest` answers "what is the nearest annotated TSS?" (a coordinate fact), which is routinely misread as "which gene does this element regulate?" (a biology claim). For distal enhancers those disagree the majority of the time, so the agent does the bedtools arithmetic honestly and routes real enhancer-to-gene linking to activity/contact/QTL methods - while recognising that for GWAS locus-to-gene the nearest protein-coding gene is a defensible prior, the opposite advice from the same command.

## Prerequisites
```bash
# bedtools (required for all operations)
conda install -c bioconda bedtools

# pybedtools (optional Python wrapper; needs the bedtools binary on PATH)
pip install pybedtools

# samtools, only to build a genome (chrom-sizes) file from a FASTA index
conda install -c bioconda samtools
```

A chrom-sizes file (`genome.txt`, two columns chrom<TAB>length) is REQUIRED for `slop` and `flank`. Build it with `cut -f1,2 reference.fa.fai > genome.txt` or download UCSC `*.chrom.sizes`. `closest` requires both inputs coordinate-sorted.

## Quick Start
Tell your AI agent what you want to do:
- "Find the nearest gene to each ChIP-seq peak, signed by the gene's strand"
- "Build strand-aware promoters at TSS minus 2 kb / plus 200 bp from my gene model"
- "List all genes within 50 kb of each ATAC peak as candidates"
- "Which gene does this distal enhancer regulate?" (the agent will flag this as a linking question, not a proximity one)
- "Get 1 kb flanking regions on each side of every exon"

## Example Prompts

### Nearest Feature
> "Sort my peaks and the gene model, then run bedtools closest with signed distance by the gene's strand, ignoring overlaps and resolving ties to the first feature, and drop the no-feature -1 rows before reporting."
> "Report the distance from each variant to the nearest exon and give me the distribution, not just a promoter/distal flag."

### Enhancer-to-Gene (the interpretive case)
> "I have distal H3K27ac peaks - assign each to its target gene." (Expect the agent to explain that nearest-gene is wrong the majority of the time for distal elements, produce a candidate set with closest/window, and route the real linking to atac-seq/enhancer-gene-linking via ABC/PCHi-C/eQTL.)
> "For this fine-mapped GWAS credible-set SNP, what gene does the locus implicate?" (Here nearest protein-coding gene is a fair ~50-65% prior; the agent should say so and point to causal-genomics for the rigorous version.)

### Promoters and Windows
> "Collapse my BED6 genes to TSS, then build strand-aware promoter windows of 2 kb upstream and 200 bp downstream, and warn me about any windows clipped at chromosome ends."
> "Count how many genes fall within a 50 kb window of each peak."

### Flanking Regions
> "Give me the 1 kb regions on each side of every exon for splice context, strand-aware."

## What the Agent Will Do
1. Identify the regime - promoter-proximal mark, distal enhancer, or GWAS locus - because it changes whether nearest-gene is an answer or only a candidate.
2. Prepare inputs: build/verify the chrom-sizes file for slop/flank; coordinate-sort both inputs for closest; harmonize chromosome naming across all files.
3. For promoters, collapse genes to TSS first (strand-aware), then slop strand-aware - never slop a gene body.
4. Run the operation with the correct strand and sign flags (`-D b` for upstream/downstream biology, `-s` for slop/flank), resolving ties explicitly.
5. Filter the `-1` no-feature sentinel and verify slop/flank widths near chromosome ends.
6. For distal-enhancer linking, return candidates plus a clear hand-off to ABC/PCHi-C/eQTL methods rather than a single "nearest gene" call.

## Tips
- Use `-D b` (sign by the gene's strand), never `-D ref`, for any upstream/downstream claim - `-D ref` silently mis-signs every minus-strand gene.
- A promoter is a definition you impose on a TSS, not an annotated feature; report the window (e.g. -2000/+200) and treat it as a tunable parameter.
- `slop -b 2000` on a gene body is the wrong promoter; collapse to TSS, then `slop -s -l 2000 -r 200`.
- `slop` and `flank` clip silently at chromosome ends - verify `end-start` equals the requested width; treat contig-end features as edge cases.
- Default `closest -t all` reports every tie; counting rows as peaks double-counts equidistant peaks (concentrated at bidirectional promoters). Use `-t first` or aggregate by distinct peak, or use GREAT/rGREAT for peak-set enrichment.
- `flank` and `slop` need a genome file; `closest` and `window` do not.
- pybedtools shells out to the bedtools binary - ensure it is on PATH, and call `.sort()` before `.closest()`.

## Resources
- [bedtools closest](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html)
- [bedtools window](https://bedtools.readthedocs.io/en/latest/content/tools/window.html)
- [bedtools slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html)
- [bedtools flank](https://bedtools.readthedocs.io/en/latest/content/tools/flank.html)

## Related Skills
- bed-file-basics - BED coordinate systems and the sort/conversion this skill depends on
- gtf-gff-handling - Extract TSS and gene models from GTF/GFF for promoter construction
- interval-arithmetic - intersect/merge/subtract; window -w 0 is approximately intersect
- chip-seq/peak-annotation - Assigns peaks to genes via the same closest-TSS logic and caveats
- atac-seq/enhancer-gene-linking - The real enhancer-to-gene science this skill routes distal calls to
- data-visualization/genome-tracks - Render the promoter/proximity intervals built here
