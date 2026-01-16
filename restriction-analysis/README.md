# Restriction Analysis

Restriction enzyme analysis using Biopython Bio.Restriction for finding cut sites, restriction mapping, enzyme selection, and fragment prediction.

## Overview

This category covers in silico restriction enzyme analysis: finding where enzymes cut a DNA sequence, creating restriction maps, selecting appropriate enzymes for cloning, and predicting fragment sizes. Bio.Restriction includes data for 800+ enzymes from REBASE.

**Tool type:** `python`
**Primary tools:** Biopython Bio.Restriction

## Skills

| Skill | Description |
|-------|-------------|
| [restriction-sites](restriction-sites/) | Find where enzymes cut a sequence |
| [restriction-mapping](restriction-mapping/) | Create restriction maps, visualize cut positions |
| [enzyme-selection](enzyme-selection/) | Choose enzymes by criteria (cutters, overhangs, compatibility) |
| [fragment-analysis](fragment-analysis/) | Predict fragment sizes, simulate gel electrophoresis |

## Workflow

```
DNA Sequence (FASTA/GenBank)
    |
    v
[restriction-sites] --> Find cut sites for specific enzymes
    |
    v
[restriction-mapping] --> Visualize cut positions on sequence
    |
    v
[enzyme-selection] ---> Choose enzymes for cloning strategy
    |
    v
[fragment-analysis] --> Predict and analyze fragments
```

## Enzyme Classes

| Class | Description | Example |
|-------|-------------|---------|
| Type II | Standard restriction enzymes | EcoRI, BamHI, HindIII |
| Type IIS | Cut outside recognition site | BsaI, BbsI |
| Blunt cutters | Leave blunt ends | EcoRV, SmaI |
| 5' overhang | Leave 5' sticky ends | EcoRI, BamHI |
| 3' overhang | Leave 3' sticky ends | PstI, KpnI |
| Neoschizomers | Same site, different cut | BamHI/BglII |
| Isoschizomers | Same site, same cut | HpaII/MspI |

## Common Enzymes Reference

| Enzyme | Recognition | Cut | Overhang |
|--------|-------------|-----|----------|
| EcoRI | GAATTC | G^AATTC | 5' AATT |
| BamHI | GGATCC | G^GATCC | 5' GATC |
| HindIII | AAGCTT | A^AGCTT | 5' AGCT |
| EcoRV | GATATC | GAT^ATC | Blunt |
| PstI | CTGCAG | CTGCA^G | 3' TGCA |
| NotI | GCGGCCGC | GC^GGCCGC | 5' GGCC |
| XhoI | CTCGAG | C^TCGAG | 5' TCGA |
| SalI | GTCGAC | G^TCGAC | 5' TCGA |
| BsaI | GGTCTC(N)1 | Cut outside | 5' 4nt |

## Example Prompts

### Finding Cut Sites
- "Find all EcoRI sites in this sequence"
- "Where does BamHI cut in my plasmid?"
- "Show all restriction sites for common cloning enzymes"

### Restriction Mapping
- "Create a restriction map of this sequence"
- "Map EcoRI, BamHI, and HindIII sites"
- "Show distances between cut sites"

### Enzyme Selection
- "Find enzymes that cut this sequence exactly once"
- "Which enzymes don't cut my insert?"
- "Find enzymes with compatible sticky ends"
- "List all 6-cutter enzymes that cut my sequence"
- "Is my insert compatible with Golden Gate cloning?"
- "Find enzymes not affected by Dam methylation"

### Fragment Analysis
- "What fragments will EcoRI produce?"
- "Predict the gel pattern for this digest"
- "Calculate fragment sizes for a double digest"

## Requirements

```bash
pip install biopython
```

## Key Classes and Functions

| Class/Function | Purpose |
|----------------|---------|
| RestrictionBatch | Analyze with multiple enzymes |
| Analysis | Perform restriction analysis |
| AllEnzymes | Collection of all enzymes |
| CommOnly | Commercially available enzymes |
| enzyme.search() | Find sites in sequence |
| enzyme.catalyze() | Get fragments after digestion |
| enzyme.ovhgseq | Get overhang sequence |
| enzyme.is_blunt() | Check if blunt cutter |
| enzyme.is_5overhang() | Check if 5' overhang |
| enzyme.is_3overhang() | Check if 3' overhang |
| enzyme.is_dam_methylable() | Check Dam methylation sensitivity |
| enzyme.is_dcm_methylable() | Check Dcm methylation sensitivity |

## Notes

- **Sequence must be Seq object** - not plain string
- **Linear vs circular** - Analysis handles both, specify with `linear=True/False`
- **Case insensitive** - Recognition sites match regardless of case
- **Ambiguous bases** - Some enzymes recognize ambiguous sites (N, R, Y, etc.)
- **Methylation sensitivity** - Some enzymes are blocked by methylation (check enzyme.is_methylable())

## Related Skills

- **sequence-io** - Read sequences for restriction analysis
- **sequence-manipulation** - Work with restriction fragments
- **sequence-properties** - Analyze fragment sequences

## References

- [Bio.Restriction Documentation](https://biopython.org/docs/latest/api/Bio.Restriction.html)
- [REBASE Database](http://rebase.neb.com/)
- [Biopython Tutorial - Restriction](https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:restriction)
