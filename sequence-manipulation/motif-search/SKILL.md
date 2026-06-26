---
name: bio-motif-search
description: Find sequence motifs, degenerate IUPAC patterns, and transcription-factor binding sites in DNA/RNA using Biopython and regex, including position weight matrix (PWM/PSSM) scoring. Use when locating regulatory elements, counting overlapping motif occurrences, scanning for binding-site matches above a significance threshold, or reading motif matrices from JASPAR/MEME/TRANSFAC files. For restriction enzyme sites, use restriction-analysis/restriction-sites.
tool_type: python
primary_tool: Bio.motifs
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Motif Search

**"Search for a sequence motif or binding-site pattern"** -> Scan sequences for a fixed motif, a degenerate IUPAC consensus, or a probabilistic PWM, on one or both strands, and locate transcription-factor binding sites, regulatory elements, or custom patterns.
- Python: `Bio.SeqUtils.nt_search` (IUPAC + overlaps), `re` (regex/lookahead), `Bio.motifs` (PWM/PSSM scoring + matrix file parsing)

## The Governing Principle

Two silent failures dominate motif searching; both return a plausible-but-wrong answer with no error:

1. **Overlapping matches are dropped.** `str.count`, `str.find`, and `re.findall` consume the string left to right, so a motif that overlaps its own next occurrence is undercounted. Target `AAGCGCGCGAA`, motif `GCGC`: `str.count` returns 1, the true answer is 2 (starts at positions 2 and 4). Use a zero-width lookahead `re.finditer(r'(?=(GCGC))', target)` or `Bio.SeqUtils.nt_search`, both of which report overlaps.
2. **A PSSM score is a likelihood in bits, not a probability.** `pssm.calculate` returns log2-odds versus background. A "high-looking" threshold chosen by eye is arbitrary and non-reproducible; derive the threshold from the score distribution at a chosen false-positive rate. And a PSSM scans only the strand it is given, so scoring just the forward strand silently misses roughly half of real sites on double-stranded DNA.

## Which Approach for Which Question

| Question | Tool |
|----------|------|
| Position of first exact hit | `str.find` / `Seq.find` (returns -1 if absent) |
| All exact hits, possibly overlapping | `re.finditer(r'(?=(motif))', seq)` |
| Degenerate IUPAC consensus (e.g. `GATNNTC`), with overlaps | `Bio.SeqUtils.nt_search(seq, motif)` |
| Flexible / repeat / variable-spacer pattern | `re` with explicit character classes and quantifiers |
| Graded match to many aligned sites (binding sites) | `Bio.motifs` PWM -> PSSM, score and threshold |
| Match significance / false-positive control | `pssm.distribution(...).threshold_fpr(fpr)` |
| Restriction enzyme recognition sites | restriction-analysis/restriction-sites |

## IUPAC Degenerate Motifs

A degenerate motif expands each ambiguity code to a regex character class:

| Code | Class | Code | Class | Code | Class |
|------|-------|------|-------|------|-------|
| N | `[ACGT]` | R | `[AG]` | Y | `[CT]` |
| W | `[AT]` | S | `[GC]` | K | `[GT]` |
| M | `[AC]` | B | `[CGT]` | D | `[AGT]` |
| H | `[ACT]` | V | `[ACG]` | | |

B, D, H, V each exclude A, C, G, T respectively (the code preceding the one they drop is a mnemonic).

```python
IUPAC_DNA = {'N': '[ACGT]', 'R': '[AG]', 'Y': '[CT]', 'W': '[AT]', 'S': '[GC]',
             'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]'}

def iupac_to_regex(pattern):
    return ''.join(IUPAC_DNA.get(base, base) for base in pattern)

# 'GATNNTC' -> 'GAT[ACGT][ACGT]TC'
```

`Bio.SeqUtils.nt_search` expands IUPAC ambiguity in the query motif automatically and reports overlapping hits, so it is the shortest correct path for a degenerate consensus.

## Overlapping Matches (the count trap)

**Goal:** Report every start position of a motif, including self-overlapping occurrences.

**Approach:** Use a zero-width lookahead so the regex engine never consumes the matched text; recover the match string from the inner capture group. For IUPAC motifs, prefer `nt_search`, which both expands ambiguity and reports overlaps.

**Reference (BioPython 1.83+):**
```python
import re
from Bio.SeqUtils import nt_search

target = 'AAGCGCGCGAA'

starts = [match.start(1) for match in re.finditer(r'(?=(GCGC))', target)]  # [2, 4]
hits = [(match.start(1), match.group(1)) for match in re.finditer(r'(?=([AG]CG[CT]))', target)]

result = nt_search(target, 'GCGC')  # ['GCGC', 2, 4]
pattern, positions = result[0], result[1:]
```

`nt_search` returns a heterogeneous list: `result[0]` is the (expanded) pattern string and `result[1:]` are the 0-based start positions. When there are no hits it returns just `[pattern]` (length 1), so test `len(result) > 1` before indexing rather than truthiness.

## Bio.motifs PWM / PSSM Pipeline

**Goal:** Build a probabilistic model from a set of aligned binding sites and score a target sequence for graded matches.

**Approach:** Create a motif from instances or a matrix file, set pseudocounts and background, read the recomputed PSSM, then scan. The count matrix `m.counts['A', 0]` is indexed `[base, position]`.

**Reference (BioPython 1.83+):**
```python
from Bio import motifs
from Bio.Seq import Seq

m = motifs.create([Seq('TACAA'), Seq('TACGA'), Seq('TACTA'), Seq('TGCAA')])  # alphabet defaults to ACGT

m.pseudocounts = 0.5        # set BEFORE reading m.pssm (see trap below)
m.background = None         # None gives uniform 0.25; or pass a dict of base frequencies

pwm = m.pwm                 # normalized frequencies (property)
pssm = m.pssm               # log2-odds vs background (property; RECOMPUTED on each access)

m.consensus                 # most frequent base per column
m.degenerate_consensus      # IUPAC-degenerate consensus
```

`m.counts.normalize(pseudocounts=0.5)` returns a position weight matrix and `pwm.log_odds()` returns a PSSM; these are equivalent to reading `m.pwm` / `m.pssm` after setting `m.pseudocounts`.

### The Pseudocount / -inf Trap (silent)

A column where some base has count 0 gives that base frequency 0 and a log-odds of negative infinity; any target carrying that base at that position then scores `-inf` and is unmatchable. This is common with short motifs or few instances. Setting `m.pseudocounts` (a flat 0.5, or sqrt(N) with N the number of instances; scalar or per-base dict) makes every cell finite by shrinking toward background.

Critically, `m.pssm` is recomputed from `m.pseudocounts` and `m.background` on every access. Set both BEFORE reading `m.pssm` (or `pwm.log_odds()`), or the matrix is silently wrong.

### Score, Threshold, and P-value

`pssm.calculate(seq)` returns the log2-odds score in bits for each window (a relative likelihood, not a probability). `pssm.search(seq, threshold=...)` yields `(position, score)` pairs at or above the threshold.

To convert a bit score into a false-positive rate, build the null distribution and ask it for a threshold:

```python
dist = pssm.distribution(background=m.background, precision=10**4)
threshold = dist.threshold_fpr(0.01)        # 1% false-positive rate
threshold = dist.threshold_fnr(0.1)         # 10% false-negative rate
threshold = dist.threshold_balanced(1000)   # rate_proportion = FNR:FPR ratio (FNR = FPR x rate_proportion), NOT a sequence length; default 1.0 gives FPR=FNR
```

Choosing a threshold "because it looks high" is the classic non-reproducible error. Higher `precision` gives finer threshold resolution at the cost of memory.

### Both Strands

`pssm.calculate` scans only the strand it is handed. `pssm.search` defaults to `both=True`, scanning both strands in one call; with `both=True` a hit at negative position `p` lies on the reverse strand and its forward-coordinate start is `len(seq) + p`. To handle strands separately, set `both=False` and scan the reverse-complemented PSSM explicitly:

```python
combined = list(pssm.search(seq, threshold=3.0))  # both strands; reverse hits have NEGATIVE positions

rc_pssm = pssm.reverse_complement()
forward = list(pssm.search(seq, threshold=3.0, both=False))
reverse = list(rc_pssm.search(seq, threshold=3.0, both=False))
```

## Reading Motif Matrix Files

`motifs.read(handle, fmt)` reads exactly one motif; `motifs.parse(handle, fmt)` returns an iterator over many. The format string must match the file layout exactly.

| `fmt` string | File type / source |
|--------------|--------------------|
| `jaspar` | multi-motif JASPAR PFM collection (use `parse`) |
| `pfm` | single JASPAR-style PFM (use `read`) |
| `pfm-four-columns` | CIS-BP, HOMER, HOCOMOCO (A C G T as columns) |
| `pfm-four-rows` | ScerTF, YeTFaSCo (A C G T as rows) |
| `sites` | JASPAR sites file (use `read`) |
| `meme` | MEME program output (use `parse`) |
| `minimal` | MEME minimal text format |
| `transfac` | TRANSFAC matrices |
| `mast`, `alignace`, `clusterbuster`, `xms` | respective tool outputs |

`'cisbp'`, `'homer'`, and `'hocomoco'` are NOT valid strings; those databases use `pfm-four-columns`. The four-columns versus four-rows distinction is the most common mix-up: a 4-column matrix read as `pfm-four-rows` parses without error but produces a meaningless transposed motif.

```python
from Bio import motifs

with open('collection.jaspar') as handle:
    for m in motifs.parse(handle, 'jaspar'):
        print(m.matrix_id, m.name, m.consensus)

m.format('jaspar')      # serialize back out
m.format('transfac')
```

## Common Motif Patterns

| Motif | Pattern | Description |
|-------|---------|-------------|
| Start codon | `ATG` | Translation initiation |
| Kozak | `[AG]CCATGG` | Eukaryotic translation initiation |
| TATA box | `TATA[AT]A[AT]` | Core promoter element |
| GC box (Sp1) | `GGGCGG` | Promoter element |
| CAAT box | `CCAAT` | Promoter element |
| Poly-A signal | `AATAAA` | mRNA polyadenylation |
| E-box (bHLH) | `CA[ACGT]{2}TG` | bHLH TF binding |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Count is too low | `str.count`/`re.findall` skip overlaps | `re.finditer(r'(?=(motif))', seq)` or `nt_search` |
| `IndexError` on `nt_search` result | No hits returns `[pattern]` (length 1) | Test `len(result) > 1` before reading `result[1:]` |
| Every target scores `-inf` | Count-0 cell gives -inf log-odds | Set `m.pseudocounts` (0.5 or sqrt(N)) before reading `m.pssm` |
| PSSM scores look wrong | Pseudocounts/background set after reading `m.pssm` | Set them first; `m.pssm` is recomputed on each access |
| Roughly half of sites missed | Only forward strand scanned | Score `pssm.reverse_complement()` or pass `both=True` |
| Threshold not reproducible | Cutoff chosen by eye | `pssm.distribution(...).threshold_fpr(fpr)` |
| `ValueError` parsing matrix | Wrong `fmt` (4-column vs 4-row, `jaspar` vs `pfm`) | Match `fmt` to the actual layout |
| No matches | Case or strand mismatch | `.upper()` both; check reverse complement |

## Related Skills

- seq-objects - Create Seq objects for searching
- reverse-complement - Reverse-complement the target to search the opposite strand
- transcription-translation - ORF and codon-context motifs in coding sequences
- sequence-properties - GC content and per-sequence properties around hits
- restriction-analysis/restriction-sites - Restriction enzyme recognition sites
- chip-seq/motif-analysis - De novo motif discovery and enrichment in peak sets
- database-access/entrez-fetch - Download motif matrices from JASPAR/NCBI
