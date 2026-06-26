---
name: bio-transcription-translation
description: Transcribe DNA to RNA and translate to protein using Biopython, with NCBI codon-table selection, CDS validation, and six-frame ORF finding. Use when converting a CDS or ORF to its amino-acid sequence, selecting a non-standard (mitochondrial, bacterial, ciliate) genetic code, validating a coding sequence, or scanning all reading frames.
tool_type: python
primary_tool: Bio.Seq
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Transcription and Translation

**"Translate my DNA sequence to protein"** -> Transcribe DNA to RNA and translate to protein, choosing the right genetic code and validating the reading frame.
- Python: `Seq.translate()`, `Seq.transcribe()`, `Bio.Data.CodonTable` (BioPython)

## The Governing Principle

The single most dangerous bug in translation is **silent**: a valid-but-wrong `table=` argument produces a plausible wrong protein with no error. Translating human mitochondrial DNA with the default Standard code (table 1) inserts `*` where UGA actually codes Trp and truncates at AGA/AGG (which are stops in vertebrate mito). The protein looks real and nothing complains. By contrast, an *unknown* table id or name raises `KeyError` (loud). Only valid-but-wrong tables corrupt silently.

The defense: when the input is a complete coding sequence, pass `cds=True`. It converts silent traps into loud `TranslationError` exceptions by validating start, length, and stop. Use it for ORF validation rather than trusting a clean-looking output.

Since the Biopython 1.78 alphabet removal, `transcribe()`, `back_transcribe()`, and `translate()` perform NO type checking. Transcribing a protein or translating the wrong strand returns silent garbage. Confirm the molecule and strand before converting.

## Required Import

```python
from Bio.Seq import Seq
from Bio.Data import CodonTable
```

## Transcription Is a String Operation, Not Biology

`transcribe()` is a pure T->U replacement on the coding (sense) strand; `back_transcribe()` is U->T. Neither performs splicing, intron removal, 5' capping, or poly-A addition. The Biopython tutorial states plainly that all transcribe does is replace T with U.

```python
coding_dna = Seq('ATGCGATCGATCG')
rna = coding_dna.transcribe()        # Seq('AUGCGAUCGAUCG'), T->U only
back = rna.back_transcribe()         # Seq('ATGCGATCGATCG'), U->T only
```

True biological transcription starts from the template strand, so reverse-complement first:

```python
template = Seq('CGATCGATCGCAT')
mrna = template.reverse_complement().transcribe()
```

Translation accepts DNA or RNA directly, so explicit transcription is rarely needed before `translate()`.

## Translation Basics

```python
coding_dna = Seq('ATGTTTGGT')
coding_dna.translate()               # Seq('MFG'), from DNA
Seq('AUGUUUGGU').translate()         # Seq('MFG'), from RNA
```

### Stop-Codon Behavior

`translate()` substitutes `stop_symbol` (default `'*'`) for EVERY in-frame stop, so internal stops appear as `*` mid-protein. `to_stop=True` instead halts at the first in-frame stop and does NOT append the symbol.

```python
seq = Seq('ATGTTTGGTTAAGGG')
seq.translate()                      # Seq('MFG*G'), stop shown, translation continues
seq.translate(to_stop=True)          # Seq('MFG'), halts at first stop
```

## NCBI Codon Tables: Which Code, and the Consequence of the Wrong One

Biopython exposes every NCBI genetic code by integer id or registered name. Selecting the wrong one is the #1 silent bug (see governing principle).

| ID | Name | Key reassignments vs Standard | When it matters |
|----|------|-------------------------------|-----------------|
| 1 | Standard | none (baseline) | Most nuclear genes |
| 2 | Vertebrate Mitochondrial | AGA/AGG -> STOP; AUA -> Met; UGA -> Trp | Human/vertebrate mtDNA (4 stops: UAA, UAG, AGA, AGG) |
| 3 | Yeast Mitochondrial | CUN (all four CU*) -> Thr; AUA -> Met; UGA -> Trp | **CTG -> Thr lives HERE, not table 12** |
| 4 | Mold/Protozoan Mito + Mycoplasma/Spiroplasma | UGA -> Trp (only change) | Fungal/protozoan mito; Mycoplasma |
| 5 | Invertebrate Mitochondrial | AGA/AGG -> Ser; AUA -> Met; UGA -> Trp | Insect/worm mito (AGA/AGG=Ser, not STOP as in table 2) |
| 6 | Ciliate Nuclear | UAA/UAG -> Gln; only UGA stays stop | Tetrahymena, Paramecium (single stop) |
| 11 | Bacterial/Archaeal/Plastid | same coding as Standard; expanded starts | Prokaryotes, plastids (differs from 1 mainly in initiation) |
| 12 | Alternative Yeast Nuclear | CUG -> Ser (from Leu) | *Candida* CUG-Ser clade |

**Explicit correction:** CTG -> Thr is table **3** (Yeast Mitochondrial). Table **12** is CUG -> **Ser**. Do not conflate them.

```python
seq = Seq('ATGGCCTGA')
seq.translate(table=2)                            # by NCBI integer id
seq.translate(table='Vertebrate Mitochondrial')   # by registered name

CodonTable.unambiguous_dna_by_id[2]
CodonTable.unambiguous_dna_by_name['Vertebrate Mitochondrial']
```

## Validating a Coding Sequence with cds=True

**Goal:** Translate a complete ORF and have any structural defect raise a loud error instead of producing a silent wrong protein.

**Approach:** Pass `cds=True`. It enforces four conditions, each raising `Bio.Data.CodonTable.TranslationError` on failure: (1) first codon is a start codon for the chosen table; (2) length is a multiple of 3; (3) sequence ends in a stop; (4) no internal in-frame stop. A valid alternative start (GTG/TTG/ATT) is translated as M, biologically correct for fMet initiation. The terminal stop is stripped from the output.

**Reference (BioPython 1.83+):**

```python
cds = Seq('ATGTTTGGTTAA')
cds.translate(cds=True)              # Seq('MFG'), validated, terminal stop removed

alt_start = Seq('GTGTTTGGTTAA')
alt_start.translate(table=11, cds=True)   # Seq('MFG'), GTG start -> M under bacterial code
```

Start-codon lists differ by table: table 1 = TTG/CTG/ATG; table 2 = ATT/ATC/ATA/ATG/GTG; table 11 = TTG/CTG/ATT/ATC/ATA/ATG/GTG. A start valid under one table fails under another, which is exactly the loud signal `cds=True` provides.

## translate() Parameters and Edge Cases

Signature (the `Seq.translate` METHOD): `translate(table='Standard', stop_symbol='*', to_stop=False, cds=False, gap='-')`. Note the method defaults to `gap='-'`; only the module-level `Bio.Seq.translate(sequence, ...)` function defaults to `gap=None`.

- **Partial codon** (length not a multiple of 3): emits a `BiopythonWarning` and SILENTLY drops the trailing 1-2 bases. Easy to miss in a pipeline. Under `cds=True` the same condition becomes a loud `TranslationError`.
- **Gaps:** because the method default is `gap='-'`, a full gap codon `'---'` already translates to `'-'` (e.g. `Seq('GTG---GCCATT').translate()` -> `'V-AI'`, no error). A codon mixing gap and bases (`'TT-'`) raises `TranslationError`. Pass `gap=None` (or use the module-level function) to make any `-` invalidate its codon instead.
- **Dual-coding stop tables (27 Karyorelict, 28 Condylostoma, 31 Blastocrithidia Nuclear):** these reassign a stop codon so it codes both an amino acid and stop, so `to_stop=True` raises a `ValueError` (no single truncation point).

```python
Seq('ATGTTTGG').translate()          # BiopythonWarning, trailing 'GG' dropped -> Seq('MF')
Seq('ATGTTTGG').translate(cds=True)  # TranslationError: length not a multiple of three
```

## Selenocysteine and Pyrrolysine Are Silently Lost

Selenocysteine (Sec, one-letter U) is encoded by UGA and pyrrolysine (Pyl, one-letter O) by UAG, both normally stop codons. Recoding requires a SECIS (Sec) or PYLIS (Pyl) element that Biopython does NOT detect. No NCBI table maps UGA->U or UAG->O. Naive translation therefore yields `*` mid-protein, and `to_stop=True` SILENTLY truncates the protein at that position. Real selenoproteins (GPX, TXNRD, SELENOP) come out truncated or peppered with `*`. There is no clean Biopython workaround; flag these genes and handle the recoding event manually.

## Six-Frame Translation

**Goal:** Translate a DNA sequence in all six frames (three forward, three reverse) to expose every possible protein product.

**Approach:** For each strand, offset by 0, 1, 2 bases, trim to a multiple of 3, and translate.

**Reference (BioPython 1.83+):**

```python
def six_frame_translation(seq):
    frames = []
    for strand, s in [('+', seq), ('-', seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(s) - frame) // 3)
            fragment = s[frame:frame + length]
            frames.append((strand, frame, fragment.translate()))
    return frames

seq = Seq('ATGCGATCGATCGATCGATCG')
for strand, frame, protein in six_frame_translation(seq):
    print(f'{strand}{frame}: {protein}')
```

## Find All ORFs (Start to Stop)

**Goal:** Identify all open reading frames (Met to stop) across both strands and all three frames, keeping only those above a minimum length.

**Approach:** Translate each of the six frames, then scan each translation for Met-to-stop segments meeting the threshold.

**Reference (BioPython 1.83+):**

```python
def find_orfs(seq, min_protein_length=30):
    orfs = []
    for strand, s in [('+', seq), ('-', seq.reverse_complement())]:
        for frame in range(3):
            end = frame + 3 * ((len(s) - frame) // 3)
            trans = str(s[frame:end].translate())
            aa_start = 0
            while True:
                start = trans.find('M', aa_start)
                if start == -1:
                    break
                stop = trans.find('*', start)
                if stop == -1:
                    stop = len(trans)
                orf = trans[start:stop]
                if len(orf) >= min_protein_length:
                    orfs.append((strand, frame, start * 3 + frame, orf))
                aa_start = start + 1
    return orfs

seq = Seq('ATGCGATCGATCGATCGATCGTAA')
for strand, frame, pos, orf in find_orfs(seq, min_protein_length=3):
    print(f'{strand} frame {frame} pos {pos}: {orf}')
```

## Inspect a Codon Table

```python
table = CodonTable.unambiguous_dna_by_id[2]
table.start_codons                   # ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
table.stop_codons                    # ['TAA', 'TAG', 'AGA', 'AGG']
table.forward_table['TGA']           # 'W' under vertebrate mito code
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Plausible protein, wrong residues, no error | Valid-but-wrong `table=` (e.g. mito DNA on table 1) | Select the organism's NCBI table; use `cds=True` to validate |
| `*` mid-protein or premature truncation | Selenoprotein/pyrrolysine UGA/UAG, or wrong table where UGA=Trp | Use correct mito table for UGA=Trp; Sec/Pyl recoding is not automatic |
| `TranslationError: First codon ... is not a start codon` | `cds=True` on a sequence not starting at a valid start for that table | Trim to the true start, or pick the table whose starts include it |
| `TranslationError: ... is not a multiple of three` | `cds=True` on a partial CDS | Trim to a full ORF; without `cds=True` this only warns and drops trailing bases |
| `TranslationError: Extra in frame stop codon found` | Internal stop under `cds=True` | Wrong frame, wrong table, or genuine internal stop; re-check frame/table |
| Garbage protein from a protein input | `transcribe()`/`translate()` on a non-nucleotide Seq (no type checks since 1.78) | Verify molecule type before converting |
| `KeyError` | Unknown table id or name | Use a valid NCBI id (1-6, 9-16, 21-31) or registered name |

## Decision Tree

```
Need to convert a sequence?
├── DNA <-> RNA (string-level T<->U)?
│   ├── coding strand to RNA -> seq.transcribe()
│   ├── RNA back to DNA      -> seq.back_transcribe()
│   └── template strand to mRNA -> seq.reverse_complement().transcribe()
├── DNA/RNA to protein?
│   ├── complete CDS to validate -> translate(cds=True) [loud on defects]
│   ├── stop at first stop only  -> translate(to_stop=True)
│   ├── non-standard organism    -> translate(table=N)  [pick from the table above]
│   └── show internal stops      -> translate()  [* per stop]
└── Unknown coding regions? -> six-frame translation, then scan M...* for ORFs
```

## Related Skills

- seq-objects - Create and inspect Seq objects before translation
- reverse-complement - Strand handling for six-frame translation and template-strand transcription
- codon-usage - Analyze codon bias and adaptation in coding sequences
- sequence-io/read-sequences - Parse GenBank/FASTA records and CDS features for translation

## References

The genetic-code tables and their organism assignments follow the NCBI Taxonomy "The Genetic Codes" page, compiled by Andrzej (Anjay) Elzanowski and Jim Ostell at NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). This is a maintained web resource; cite it as the NCBI page rather than as a journal article.
