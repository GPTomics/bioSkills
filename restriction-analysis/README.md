# restriction-analysis

## Overview

Restriction enzyme analysis using Biopython Bio.Restriction. Find cut sites, create restriction maps, select enzymes for cloning, and predict fragment sizes. Includes data for 800+ enzymes from REBASE.

**Tool type:** python | **Primary tools:** Bio.Restriction

## Skills

| Skill | Description |
|-------|-------------|
| restriction-sites | Find where enzymes cut a sequence |
| restriction-mapping | Create restriction maps, visualize cut positions |
| enzyme-selection | Choose enzymes by criteria (cutters, overhangs, compatibility) |
| fragment-analysis | Predict fragment sizes, simulate gel electrophoresis |

## Example Prompts

- "Find all EcoRI sites in this sequence"
- "Where does BamHI cut in my plasmid?"
- "Show all restriction sites for common cloning enzymes"
- "Create a restriction map of this sequence"
- "Map EcoRI, BamHI, and HindIII sites"
- "Show distances between cut sites"
- "Find enzymes that cut this sequence exactly once"
- "Which enzymes don't cut my insert?"
- "Find enzymes with compatible sticky ends"
- "List all 6-cutter enzymes that cut my sequence"
- "Is my insert compatible with Golden Gate cloning?"
- "Find enzymes not affected by Dam methylation"
- "What fragments will EcoRI produce?"
- "Predict the gel pattern for this digest"
- "Calculate fragment sizes for a double digest"

## Requirements

```bash
pip install biopython
```

## Related Skills

- **sequence-io** - Read sequences for restriction analysis
- **sequence-manipulation** - Work with restriction fragments
- **primer-design** - Design primers around restriction sites
