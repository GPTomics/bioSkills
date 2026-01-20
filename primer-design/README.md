# primer-design

## Overview

Design and validate PCR primers using primer3-py, the Python binding for Primer3.

**Tool type:** python | **Primary tools:** primer3-py

## Skills

| Skill | Description |
|-------|-------------|
| primer-basics | Design PCR primers for a target sequence with primer3-py |
| qpcr-primers | Design qPCR primers and probes with Tm matching |
| primer-validation | Check primers for specificity, dimers, and secondary structures |

## Example Prompts

- "Design primers for this sequence"
- "Find primers flanking position 500-800 in my gene"
- "Design qPCR primers with a TaqMan probe"
- "Check my primers for dimers"
- "Design primers with Tm between 58-62C"
- "Find primers that avoid this SNP position"
- "Check if my primer has secondary structures"
- "Design primers for a 200bp amplicon"

## Requirements

```bash
pip install primer3-py biopython
```

## Related Skills

- **sequence-io** - Load sequences for primer design
- **sequence-manipulation** - Extract target regions
- **database-access** - BLAST primers for specificity
