---
name: bio-clip-seq-clip-motif-analysis
description: Identify enriched sequence motifs at CLIP-seq binding sites for RBP binding specificity. Use when characterizing the sequence preferences of an RNA-binding protein.
tool_type: cli
primary_tool: HOMER
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# CLIP Motif Analysis

**"Find sequence motifs at my RBP binding sites"** → Discover enriched RNA sequence motifs at CLIP-seq peaks to determine the binding specificity of an RNA-binding protein.
- CLI: `findMotifs.pl peaks.fa fasta output/ -rna` (HOMER)
- CLI: `bedtools getfasta` to extract peak sequences first

## HOMER De Novo Motifs

**Goal:** Discover enriched RNA sequence motifs at CLIP-seq binding sites.

**Approach:** Extract FASTA sequences from peak regions using bedtools getfasta, then run HOMER findMotifs.pl in RNA mode to identify overrepresented motifs.

```bash
# Extract sequences from peaks
bedtools getfasta -fi genome.fa -bed peaks.bed -fo peaks.fa

# Find enriched motifs
findMotifs.pl peaks.fa fasta output_dir \
    -len 6,7,8 \
    -rna
```

## MEME-ChIP

```bash
meme-chip -oc output_dir \
    -dna \
    peaks.fa
```

## Known Motif Enrichment

```bash
# HOMER known motif scan
findMotifs.pl peaks.fa fasta output_dir \
    -rna \
    -known
```

## Related Skills

- clip-peak-calling - Get peaks
- chip-seq/motif-analysis - General motif concepts
