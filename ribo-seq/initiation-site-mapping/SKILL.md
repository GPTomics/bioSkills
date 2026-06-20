---
name: bio-ribo-seq-initiation-site-mapping
description: Map translation initiation sites, including non-AUG and alternative starts, from initiation-drug ribosome profiling (TI-seq). Use when locating start codons, detecting near-cognate or upstream initiation, or analyzing harringtonine, lactimidomycin (GTI-seq/QTI-seq), or retapamulin (Ribo-RET) data.
tool_type: mixed
primary_tool: Ribo-TISH
---

## Version Compatibility

Reference examples tested with: Ribo-TISH 0.2.7+, PRICE/GEDI 1.0.5+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Translation Initiation Site Mapping

**"Map where translation starts in my Ribo-seq data"** -> Locate translation initiation sites (TIS) at single-nucleotide resolution, including non-AUG and upstream starts, from initiation-drug profiling experiments.
- CLI: `Ribo-TISH` for TIS detection from harringtonine/LTM data; `PRICE` for EM-based cryptic-start detection

This is a distinct analysis from elongation ORF detection: it asks WHERE initiation occurs (which start codon), not which ORF bodies are translated. It typically requires a dedicated initiation-drug library paired with a standard elongation library.

## Initiation-drug data types (which experiment produced the data)

| Method | Drug(s) | Signal | Citation |
|--------|---------|--------|----------|
| Harringtonine TIS | harringtonine | binds free 60S, blocks the first peptide bond; broad start peak | Ingolia 2011 |
| GTI-seq | lactimidomycin (LTM) + CHX in parallel | LTM blocks translocation at the assembled 80S; sharp start peak | Lee 2012 |
| QTI-seq | LTM then puromycin (sequential) | puromycin strips elongating ribosomes; quantitative, low background | Gao 2015 |
| Ribo-RET (bacteria) | retapamulin | arrests initiating 70S at start codons | Meydan 2019 |

All three drugs CREATE the initiation signal by halting or removing elongation; the data is a deliberate artifact read out at the start codon. LTM gives sharper peaks than harringtonine because it cannot act on elongating ribosomes whose E-site is occupied. QTI-seq (LTM then puromycin) is analyzed on the same LTM path below; the puromycin step only strips elongating ribosomes to lower the background, so the TIS library is still passed as the LTM-type `-t` input. Without an initiation-drug library, start codons can only be inferred indirectly from elongation periodicity (see orf-detection).

## Near-cognate and alternative starts

Initiation occurs at AUG and near-cognate codons differing by one base; the biologically used set is CUG, GUG, ACG, UUG, AUU, AUC, AUA (AAG/AGG also differ by one base but initiate negligibly). CUG is the dominant near-cognate start (~16% of mapped sites in GTI-seq; AUG remains >50%). uORFs especially use near-cognate starts, so initiation mapping must enable alternative start codons to recover them; an AUG-only search misses most upstream initiation.

## Tool selection

| Situation | Tool | Why |
|-----------|------|-----|
| TIS from harringtonine/LTM data, with QC | Ribo-TISH | quality + predict modes; near-cognate via --alt; differential TIS |
| Cryptic/near-cognate starts, EM model | PRICE | per-codon EM; handles near-cognate; designed for cryptic events |
| Bacterial initiation (Ribo-RET) | dedicated retapamulin analysis | prokaryote initiation; eukaryote periodicity tools fit poorly |

## QC the initiation library

**Goal:** Confirm the drug enriched start-codon signal and pick P-site offsets before predicting.

**Approach:** Run Ribo-TISH quality, which reports the metagene profile and writes a per-length offset parameter file.

```bash
# Writes a <ribo.bam>.para.py offset file and a QC figure
ribotish quality -b ribo_elongation.bam -g annotation.gtf -o ribo_quality.txt -f ribo_qc.pdf
ribotish quality -b ribo_tis.bam -g annotation.gtf -o tis_quality.txt -f tis_qc.pdf
```

## Predict initiation sites with Ribo-TISH

**Goal:** Call TIS, including non-AUG starts, using the initiation-drug library.

**Approach:** Run `ribotish predict` with the elongation BAM (-b) and the TIS/harringtonine/LTM BAM (-t), enabling alternative start codons.

```bash
# --harr marks the TIS library as harringtonine-type; --alt enables near-cognate starts
ribotish predict \
    -b ribo_elongation.bam \
    -t ribo_tis.bam \
    -g annotation.gtf \
    -f genome.fa \
    --harr --harrwidth 15 --alt \
    -o tis_predictions.txt
```

The output lists initiation sites with the start codon, ORF type, and significance. For differential initiation across conditions, `ribotish tisdiff` compares two TIS libraries.

## Alternative: cryptic starts with PRICE

**Goal:** Detect cryptic and near-cognate initiation with an EM model.

**Approach:** Prepare the genome and run the Price tool in GEDI on the Ribo-seq reads.

```bash
gedi -e Price -reads ribo_elongation.bam -genomic prepared_genome -prefix price_out
```

PRICE reports a per-ORF p-value from a generalized binomial model (not multiple-testing corrected); codon-level activity is written to `price_out.codons.cit`.

## Interpreting initiation sites

A called TIS is strongest when it shows a sharp drug-induced start peak, a downstream in-frame elongation signal in the standard library, and (for novel sites) conservation or peptide support. Alternative N-terminal starts and uORF starts frequently use near-cognate codons; report the start codon identity, not just the position. Initiation at a uORF does not guarantee a stable protein product (see orf-detection validation).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Only AUG starts found | Alternative starts not enabled | Add `--alt` (Ribo-TISH) or use PRICE for near-cognate |
| Broad, smeared start peaks | Harringtonine data treated as sharp LTM data | Use `--harr`; expect broader peaks than LTM |
| `predict` gives weak calls | Missing the paired elongation BAM (-b) | Provide both -b (elongation) and -t (TIS) libraries |
| TIS analysis on elongation-only data | No initiation-drug library present | Initiation mapping needs harringtonine/LTM/RET data; otherwise infer from periodicity |
| Bacterial data mis-called | Eukaryote TIS tool on Ribo-RET data | Use a retapamulin/prokaryote initiation workflow |

## Related Skills

- orf-detection - Call and validate the ORF bodies downstream of mapped starts
- ribosome-periodicity - Calibrate P-site offsets for both libraries
- riboseq-preprocessing - Align the elongation and initiation-drug libraries
- ribosome-stalling - Initiation drugs are not for elongation pausing

## References

- Ingolia NT, Lareau LF, Weissman JS. 2011. Ribosome profiling of mouse embryonic stem cells reveals the complexity and dynamics of mammalian proteomes. Cell 147(4):789-802. doi:10.1016/j.cell.2011.10.002
- Lee S, Liu B, Lee S, Huang SX, Shen B, Qian SB. 2012. Global mapping of translation initiation sites in mammalian cells at single-nucleotide resolution. Proc Natl Acad Sci USA 109(37):E2424-E2432. doi:10.1073/pnas.1207846109
- Gao X, Wan J, Liu B, Ma M, Shen B, Qian SB. 2015. Quantitative profiling of initiating ribosomes in vivo. Nat Methods 12(2):147-153. doi:10.1038/nmeth.3208
- Zhang P, He D, Xu Y, et al. 2017. Genome-wide identification and differential analysis of translational initiation. Nat Commun 8:1749. doi:10.1038/s41467-017-01981-8
- Erhard F, Halenius A, Zimmermann C, et al. 2018. Improved Ribo-seq enables identification of cryptic translation events. Nat Methods 15(5):363-366. doi:10.1038/nmeth.4631
- Meydan S, Marks J, Klepacki D, et al. 2019. Retapamulin-assisted ribosome profiling reveals the alternative bacterial proteome. Mol Cell 74(3):481-493. doi:10.1016/j.molcel.2019.02.017
