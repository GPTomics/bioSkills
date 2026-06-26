---
name: bio-reverse-complement
description: Generate reverse complements and complements of DNA/RNA sequences using Biopython, including IUPAC ambiguity codes, gapped alignments, and minus-strand features. Use when working with the opposite strand, building reverse primers, normalizing strand orientation before alignment, or extracting a coding sequence from a minus-strand feature.
tool_type: python
primary_tool: Bio.Seq
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Reverse Complement

Generate complementary and reverse complementary sequences using Biopython.

**"Get the reverse complement"** -> Produce the 5'-to-3' sequence of the opposite strand.
- Python: `seq.reverse_complement()` (BioPython `Seq`)
- CLI: `samtools faidx ref.fa region --reverse-complement` (extracts and RCs a region)

## The Governing Principle

Never hand-roll the complement table. Biopython's `reverse_complement()` already encodes the full IUPAC mapping correctly, case-insensitively, and on minus-strand features it is applied for the analyst automatically by `SeqFeature.extract()`. Every silent corruption in this domain comes from reimplementing what Biopython already does right: swapping ambiguity codes, forgetting that S/W/N are self-complementary, complementing the wrong molecule type, or reverse-complementing a second time after `extract()` already did it. Reach for the library method; reach for a guard (`molecule_type`) before it; never reach for a custom dictionary.

## Required Import

```python
from Bio.Seq import Seq
```

## Which Method for Which Question

| Question | Method | Output strand/direction |
|----------|--------|-------------------------|
| Opposite strand, conventional 5'->3' | `reverse_complement()` | 5'->3' of the complementary strand (the usual answer) |
| Base-paired sequence, same direction | `complement()` | 3'->5' of the complementary strand |
| Opposite strand of RNA, keep U | `reverse_complement_rna()` | 5'->3', emits U |
| Complement of RNA, keep U | `complement_rna()` | 3'->5', emits U |
| Coding strand from template (or vice versa) | `reverse_complement()` | the other strand, 5'->3' |
| mRNA sequence from the coding strand | `transcribe()` (NOT a complement) | same strand, T->U |

### reverse_complement()

Returns the reverse complement (5'->3' of the opposite strand). This is the most commonly used operation.

```python
seq = Seq('ATGCGATCG')
rc = seq.reverse_complement()  # Returns Seq('CGATCGCAT')
```

### complement()

Returns the complement without reversing. Less common - gives the opposite strand still written in 3'->5' order.

```python
seq = Seq('ATGCGATCG')
comp = seq.complement()  # Returns Seq('TACGCTAGC')
```

### reverse_complement_rna() and complement_rna()

For RNA, the dedicated methods emit U:

```python
rna = Seq('AUGCGAUCG')
rna.reverse_complement_rna()  # Returns Seq('CGAUCGCAU')
rna.complement_rna()          # Returns Seq('UACGCUAGC')
```

## Base Pairing and Ambiguity Codes

`reverse_complement()` complements all 15 IUPAC codes plus X correctly. The mapping is non-obvious for ambiguity codes - this is exactly why hand-rolling corrupts silently.

| Code | Bases | Complement | | Code | Bases | Complement |
|------|-------|------------|-|------|-------|------------|
| A | A | T | | M | A/C | K |
| T | T | A | | B | C/G/T | V |
| G | G | C | | V | A/C/G | B |
| C | C | G | | D | A/G/T | H |
| R | A/G | Y | | H | A/C/T | D |
| Y | C/T | R | | S | G/C | S (self) |
| K | G/T | M | | W | A/T | W (self) |
|   |       |   | | N | any | N (self) |

S, W, N, and X are SELF-complementary. The pairs that get swapped wrong by hand are B<->V and D<->H. The table is built for upper and lower case, so complementation is case-insensitive (`Seq('atRY').reverse_complement()` works).

## DNA vs RNA: the U handling rule

`reverse_complement()` runs in DNA mode: it treats any U as a T and EMITS T (docstring: "Any U in the sequence is treated as a T"). It does not raise and does not leave U.

```python
Seq('ACGU').reverse_complement()      # Returns Seq('ACGT')  -- U mapped to A, emitted as T
Seq('ACGU').reverse_complement_rna()  # Returns Seq('ACGU')  -- stays RNA
```

`transcribe()` does NOT complement. It swaps T->U on the SAME strand. Confusing "complement the template" with "transcribe the coding strand" is silent corruption. True biological transcription from the template strand is `template_dna.reverse_complement().transcribe()`.

## Gaps and Non-Table Characters

`complement` and `reverse_complement` do NOT validate the alphabet (unlike `translate()`). A gap `-` is not a table key, so it passes through unchanged and reversal preserves gap columns - the desired behavior for aligned sequences. Any other non-table character (`?`, `*`) also passes through silently.

```python
Seq('ATG-CGA--TY').reverse_complement()  # Returns Seq('RA--TCG-CAT') -- gaps preserved, Y->R
```

Because there is no alphabet check, garbage in produces garbage out without a warning (see the protein trap below).

## Code Patterns

### Visualize Double-Stranded DNA

```python
def show_dsdna(seq):
    print(f"5'-{seq}-3'")
    print(f"   {'|' * len(seq)}")
    print(f"3'-{seq.complement()}-5'")

show_dsdna(Seq('ATGCGATCG'))
```

### Check if a Sequence is Palindromic (Self-Complementary)

```python
def is_palindrome(seq):
    return seq == seq.reverse_complement()

is_palindrome(Seq('GAATTC'))  # True  -- EcoRI site
is_palindrome(Seq('ATGCGA'))  # False
```

### Reverse Complement a FASTA File

**Goal:** Produce a new FASTA file with all sequences reverse-complemented.

**Approach:** Parse records as a stream, build new SeqRecords from `.reverse_complement()`, write to output.

**Reference (BioPython 1.83+):**

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def reverse_complement_records(records):
    for record in records:
        yield SeqRecord(record.seq.reverse_complement(), id=record.id + '_rc', description=record.description + ' reverse complement')

records = SeqIO.parse('sequences.fasta', 'fasta')
SeqIO.write(reverse_complement_records(records), 'sequences_rc.fasta', 'fasta')
```

### Extract a Coding Sequence from a Minus-Strand Feature

**Goal:** Get the correct 5'->3' coding sequence for a gene annotated on the minus strand.

**Approach:** Call `feature.extract(parent.seq)`. For `strand == -1`, `extract()` ALREADY reverse-complements the slice and returns the coding sequence. Do NOT reverse-complement again.

**Reference (BioPython 1.83+):**

```python
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation

parent = Seq('AAATGGGCCCTTTAAA')
feature = SeqFeature(SimpleLocation(3, 12, strand=-1), type='CDS')
cds = feature.extract(parent)  # Already reverse-complemented; this is the coding sequence
# cds.reverse_complement()     # WRONG -- double-RC bug, valid-looking but wrong strand
```

### Search Both Strands for a Motif

**Goal:** Find a motif on both strands and report forward-strand coordinates.

**Approach:** Search the forward sequence, then search its reverse complement, mapping minus-strand hits back to forward coordinates.

**Reference (BioPython 1.83+):**

```python
def search_both_strands(seq, motif):
    motif = Seq(motif)
    results = []
    pos = seq.find(motif)
    while pos != -1:
        results.append(('+', pos))
        pos = seq.find(motif, pos + 1)
    rc = seq.reverse_complement()
    pos = rc.find(motif)
    while pos != -1:
        results.append(('-', len(seq) - pos - len(motif)))
        pos = rc.find(motif, pos + 1)
    return results

search_both_strands(Seq('ATGCGAATTCGATGAATTCGATC'), 'GAATTC')
```

## In-Place Complementation

`inplace` defaults to `False` (standardized in 1.79). On an immutable `Seq`, `inplace=True` raises `TypeError: Sequence is immutable` (a loud, useful error). In-place mutation works only on `MutableSeq`.

```python
from Bio.Seq import MutableSeq
m = MutableSeq('ATGC')
m.reverse_complement(inplace=True)  # m is now MutableSeq('GCAT')
```

## The Protein Trap

Since the 1.78 alphabet removal there is no molecule-type checking. Reverse-complementing a protein produces SILENT GARBAGE with no warning: residues that are also nucleotide codes get complemented (`Seq('MAIVMGR').reverse_complement()` -> `Seq('YCKBITK')`; M->K, V->B), while protein-only letters E, F, I, L, P, Q, Z and `*` pass through unchanged. The old `IUPAC.protein` ValueError guard is gone. Guard on the molecule type, not the Seq:

```python
if record.annotations.get('molecule_type') not in ('DNA', 'RNA'):
    raise ValueError('reverse_complement is only valid for nucleotide sequences')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| U replaced by T in result | `reverse_complement()` runs in DNA mode (U treated as T) | Use `reverse_complement_rna()` to keep RNA |
| Result is meaningless letters, no error | Reverse-complemented a protein (silent since 1.78) | Guard on `molecule_type`, not the Seq |
| Coding sequence is the wrong strand | Called `.reverse_complement()` after `extract()` on a minus-strand feature | `extract()` already RC'd it; do not RC again |
| `TypeError: Sequence is immutable` | `inplace=True` on a `Seq` | Use a `MutableSeq`, or take the returned value |
| Ambiguity codes complement wrongly | Hand-rolled complement table (B/V, D/H swapped; S/W/N not self-complementary) | Use Biopython's `reverse_complement()`; never reinvent the table |
| Same strand returned instead of complement | Used `transcribe()` thinking it complements | `transcribe()` only swaps T->U; use `reverse_complement()` for the other strand |
| `TypeError` on a plain string | Passed a `str` instead of a `Seq` | Wrap input in `Seq()` first |

## References

Cornish-Bowden A (1985) "Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984." Nucleic Acids Res 13(9):3021-3030 (PMID 2582368). Defines the IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) that Biopython's complement table implements.

## Related Skills

- seq-objects - Create and mutate Seq/MutableSeq objects to complement
- transcription-translation - transcribe() vs complement(); six-frame translation uses the reverse complement
- motif-search - Search both strands by reverse-complementing the query or sequence
- sequence-io/read-sequences - Parse FASTA/GenBank records before reverse-complementing
- primer-design/primer-basics - Reverse primers are the reverse complement of the target 3' end
- restriction-analysis/restriction-sites - Restriction sites are often palindromic (self-complementary)
- alignment-files/sam-bam-basics - BAM FLAG indicates read strand; samtools view -f 16 selects reverse reads
