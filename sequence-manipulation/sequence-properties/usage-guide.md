# Sequence Properties - Usage Guide

## Overview

This skill enables AI agents to help you calculate physical and chemical properties of nucleotide and protein sequences using Biopython. It covers GC content and GC skew, molecular weight, melting temperature, and protein properties like isoelectric point, instability, and hydropathy, with emphasis on the unit and default traps that make these functions silently wrong.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate the GC content of every sequence in this FASTA file"
- "What is the double-stranded molecular weight of this DNA?"
- "Compute the nearest-neighbor Tm for this PCR primer at 50 mM Na and 250 nM primer"
- "Find the replication origin from cumulative GC skew on this bacterial genome"
- "Give me a full biophysical report for this protein: MW, pI, instability, GRAVY"

## Example Prompts

### GC Content
> "What is the GC content of each sequence in my FASTA file, as a percentage?"

### GC Skew and Replication Origin
> "Compute cumulative GC skew across this circular chromosome and tell me where oriC likely is"

### Molecular Weight
> "Calculate the molecular weight of this genomic DNA as a full duplex, not a single strand"

### Melting Temperature
> "What is the Tm of this 20-mer primer using the nearest-neighbor method with Mg correction?"

### Protein Analysis
> "Give me pI, instability index, GRAVY, and aromaticity for this protein, after stripping non-standard residues"

### Sliding Window
> "Plot GC content along this sequence using 1000 bp windows"

## What the Agent Will Do

1. Import the appropriate modules (Bio.SeqUtils, Bio.SeqUtils.ProtParam)
2. Parse or accept the sequence
3. Pick the function that matches the question and its units (fraction vs percent, single vs double strand, Tm method by length)
4. Sanitize non-standard residues before protein analysis
5. Return formatted results, summarizing across records for multiple sequences

## DNA/RNA Properties Available

- **GC content**: `gc_fraction()` returns a fraction (0-1); multiply by 100 for percent
- **GC by codon position**: `GC123()` returns four percentages (0-100)
- **GC skew**: per-window (G-C)/(G+C); cumulative skew locates replication origin and terminus
- **Molecular weight**: single-strand by default; pass `double_stranded=True` for duplex DNA
- **Melting temperature**: `Tm_NN` for primers; `Tm_Wallace`/`Tm_GC` only for short or rough estimates
- **Dinucleotide and CpG composition**: observed/expected ratios

## Protein Properties Available

- **Molecular weight**: protein average mass in Daltons
- **Isoelectric point (pI)** and **charge at pH**: linear-sequence pKa only, ignores 3D and PTMs
- **Instability index**: > 40 predicts an unstable protein
- **GRAVY**: mean Kyte-Doolittle hydropathy (default scale literal is misspelled `'KyteDoolitle'`)
- **Aromaticity**: relative frequency of Phe + Trp + Tyr
- **Secondary structure fraction**: composition propensity, not a structure prediction
- **Extinction coefficient**: reduced and oxidized values at 280 nm

## Tips

- GC content returns a fraction (0-1); multiply by 100 for a percentage, and use `ambiguous='ignore'` for a faithful legacy `GC()` drop-in
- For accurate primer Tm use `Tm_NN`; `Mg` and `dNTPs` are only honored when `saltcorr` is 6 or 7
- `molecular_weight` defaults to a single strand; pass `double_stranded=True` for genomic dsDNA, and it is not exactly double
- Protein analysis requires standard amino-acid letters; strip B, Z, X, U, `*`, `-` first or it raises KeyError
- `xGC_skew` is a Tkinter graphics routine and fails headless; sum `GC_skew()` yourself for cumulative skew

## Related Skills

- seq-objects - Create and modify Seq objects before property calculation
- codon-usage - GC123 and codon-bias indices for coding-sequence analysis
- transcription-translation - Translate a CDS before protein property analysis
- sequence-io/sequence-statistics - File-level statistics (N50, totals, dataset GC)
- primer-design/primer-basics - Design primers where Tm_NN and GC content drive the choices
- restriction-analysis/restriction-sites - Locate enzyme recognition sites in the same sequence
