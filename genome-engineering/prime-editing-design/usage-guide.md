# Prime Editing Design - Usage Guide

## Overview

Design pegRNAs for prime editing to make precise insertions, deletions, or substitutions without double-strand breaks or donor templates.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Design a pegRNA to correct the G551D mutation in CFTR"
- "Create prime editing guides to insert a FLAG tag after the start codon"
- "Design pegRNA to delete a 10bp sequence at this position"

## Example Prompts

### Point Mutations

> "Design a pegRNA to change A to G at position chr7:117559593"

> "Correct the sickle cell mutation (HBB E6V) using prime editing"

### Insertions

> "Insert a 6xHis tag at the C-terminus of my gene"

> "Design pegRNA to add a loxP site at this genomic location"

### Deletions

> "Delete the 3bp causing the deltaF508 mutation in CFTR"

> "Remove a cryptic splice site at this position"

### PE3 Strategy

> "Design PE3 nicking guides to improve my pegRNA efficiency"

> "Find optimal second nick sites for my pegRNA"

## What the Agent Will Do

1. Parse target sequence and identify edit site
2. Find suitable PAM sites near the edit
3. Design spacer sequence (20nt)
4. Optimize PBS length (13-17nt) for stability
5. Generate RT template encoding the edit
6. Optionally design PE3 nicking guide
7. Assemble full pegRNA sequence

## Tips

- **PBS length** - 13nt is a good starting point; optimize between 10-17nt
- **RT length** - Minimum 10nt; extend to 15-20nt for difficult edits
- **PE3 nicking** - Second nick 40-100bp away improves efficiency
- **Edit position** - Place edit 3-10nt from the 3' end of RT template
- **GC content** - Avoid extreme GC in PBS (aim for 40-60%)
