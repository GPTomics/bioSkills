# restriction-analysis

## Overview

Restriction enzyme analysis using Biopython Bio.Restriction. Find cut sites, build restriction maps, select enzymes for cloning and diagnostics, predict fragment sizes and gel patterns, and design or validate Type IIS scarless assembly. Includes data for the full REBASE enzyme set (1088 enzymes; 623 commercially available).

**Tool type:** python | **Primary tools:** Bio.Restriction

## Skills

| Skill | Description |
|-------|-------------|
| restriction-sites | Find where enzymes cut a sequence (linear or circular), and count cutters |
| restriction-mapping | Order cut sites, compute distances, draw text or graphical maps, overlay features |
| enzyme-selection | Choose enzymes by cut frequency, overhang, length, availability, compatibility, and methylation sensitivity |
| fragment-analysis | Predict fragment sizes and gel patterns for single and double digests |
| golden-gate-assembly | Design and validate Type IIS scarless assembly (domestication, fusion overhangs) |

## Method Selection

| Task | Skill |
|------|-------|
| "Where does this enzyme cut?" / "Does it cut at all?" | restriction-sites |
| "Order the sites / draw a map / relate sites to features" | restriction-mapping |
| "Pick enzymes to clone insert into vector / diagnostic digest / methylation-safe enzyme" | enzyme-selection |
| "What fragments / what gel bands will I get?" | fragment-analysis |
| "Design / domesticate / validate a Golden Gate or MoClo assembly" | golden-gate-assembly |

## Example Prompts

- "Find all EcoRI sites in this sequence, treating it as circular"
- "Where does BamHI cut in my plasmid?"
- "Which commercial enzymes cut my sequence exactly once?"
- "Create a restriction map of this sequence with inter-site distances"
- "Order the EcoRI, BamHI, and HindIII sites along my plasmid"
- "Find enzymes that cut the vector once but not my insert"
- "Find enzymes with sticky ends compatible with BamHI"
- "Which of these enzymes are blocked by Dam methylation in DH5-alpha DNA?"
- "What fragments will an EcoRI + BamHI double digest produce?"
- "Predict the gel pattern and tell me which bands co-migrate"
- "Scan my parts for internal BsaI sites before Golden Gate assembly"
- "Validate that my four fusion overhangs assemble uniquely"

## Requirements

```bash
pip install biopython
# matplotlib is optional, only for graphical restriction maps
pip install matplotlib
```

## Related Skills

- **sequence-io** - Read sequences for restriction analysis
- **sequence-manipulation** - Work with restriction fragments and reading frames
- **primer-design** - Add restriction sites to PCR primer tails
- **genome-engineering** - Build the constructs that Golden Gate assembly assembles
