# sequence-manipulation

Working with sequence data programmatically using Biopython's Bio.Seq and Bio.SeqUtils modules.

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

### Seq Object
Immutable sequence object. Supports standard string operations plus biological methods.

```python
from Bio.Seq import Seq
seq = Seq('ATGCGATCGATCG')
```

### MutableSeq Object
Mutable version for in-place modifications.

```python
from Bio.Seq import MutableSeq
mut_seq = MutableSeq('ATGCGATCG')
mut_seq[0] = 'C'
```

### SeqRecord Object
Seq with metadata (ID, description, features, annotations).

```python
from Bio.SeqRecord import SeqRecord
record = SeqRecord(Seq('ATGC'), id='gene1', description='Example gene')
```

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

- **restriction-analysis** (planned) - For comprehensive restriction enzyme analysis using Bio.Restriction
- **sequence-io** - For reading/writing sequence files
