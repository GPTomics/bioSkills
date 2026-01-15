# sequence-manipulation

Working with sequence data programmatically using Biopython's Bio.Seq and Bio.SeqUtils modules.

**Tool type:** python
**Primary tool:** Biopython Bio.Seq, Bio.SeqUtils

## Workflow Context

```
[sequence-io] -----> Read sequences from files
    |
    v
[sequence-manipulation] <--- YOU ARE HERE
    |                        - Create/modify Seq objects
    |                        - Transcription/translation
    |                        - Motif finding
    |                        - Sequence properties
    v
[Various downstream uses]
    - Primer design (reverse-complement)
    - Expression optimization (codon-usage)
    - Feature annotation (motif-search)
    - Quality assessment (sequence-properties)
```

This category handles programmatic manipulation of sequences after they've been read from files. Used throughout bioinformatics pipelines whenever you need to analyze or transform sequence data.

## Skills

| Skill | Description | Key Functions |
|-------|-------------|---------------|
| [seq-objects](seq-objects/) | Create and modify Seq objects | `Seq()`, `MutableSeq()`, `SeqRecord()` |
| [transcription-translation](transcription-translation/) | DNA to RNA to protein | `transcribe()`, `translate()`, `back_transcribe()` |
| [reverse-complement](reverse-complement/) | Reverse and complement sequences | `reverse_complement()`, `complement()` |
| [sequence-slicing](sequence-slicing/) | Slice, extract, concatenate | Indexing, slicing, `+` operator |
| [motif-search](motif-search/) | Find patterns and motifs | `find()`, `count()`, regex, `Bio.motifs` |
| [sequence-properties](sequence-properties/) | Calculate sequence metrics | `gc_fraction()`, `GC123()`, `GC_skew()`, `molecular_weight()` |
| [codon-usage](codon-usage/) | Analyze codon bias | `CodonAdaptationIndex`, RSCU, `GC123()` |

## Core Objects

| Object | Description | Use Case |
|--------|-------------|----------|
| `Seq` | Immutable sequence | Standard sequence operations, string-like methods plus biological methods |
| `MutableSeq` | Mutable sequence | In-place modifications (editing, site-directed changes) |
| `SeqRecord` | Seq with metadata | Full records with ID, description, features, annotations |

## Example Prompts

Ask your AI agent naturally:

| Task | Example Prompt |
|------|----------------|
| Create sequences | "Create a Seq object from this DNA string and show its properties" |
| Transcription | "Transcribe this DNA sequence to RNA" |
| Translation | "Translate this coding sequence to protein" |
| Reverse complement | "Get the reverse complement of this sequence" |
| Slice sequence | "Extract positions 100-200 from this sequence" |
| Find motif | "Find all occurrences of GAATTC in my sequence" |
| GC content | "Calculate the GC content of each sequence in my FASTA" |
| GC skew | "Plot GC skew along this sequence to find the origin" |
| Codon table | "Translate using the mitochondrial codon table" |
| ORF finding | "Find all open reading frames in this sequence" |
| Concatenate | "Join these sequences together with a linker" |
| Modify sequence | "Replace all T's with U's in this sequence" |
| Molecular weight | "Calculate the molecular weight of this protein" |
| Protein analysis | "Analyze this protein: pI, stability, hydropathy" |
| Codon usage | "What is the codon usage bias in this gene?" |
| CAI calculation | "Calculate the CAI for E. coli expression" |
| Optimize codons | "Optimize this gene for yeast expression" |
| Motif files | "Parse this JASPAR motif file and search my sequence" |

## Requirements

```bash
pip install biopython
```

## Related Skills

- **sequence-io** - Read sequences from files before manipulation
- **alignment** - Align sequences for comparison
- **restriction-analysis** (planned) - Comprehensive restriction enzyme analysis using Bio.Restriction
- **database-access** - Fetch sequences from NCBI for analysis
