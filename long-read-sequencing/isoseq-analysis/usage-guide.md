# Iso-Seq Analysis - Usage Guide

## Overview

Process PacBio Iso-Seq data for full-length transcript discovery and isoform characterization.

## Prerequisites

```bash
conda install -c bioconda pbccs lima isoseq3 sqanti3 minimap2
```

## Quick Start

- "Process my Iso-Seq data"
- "Discover novel isoforms from PacBio"
- "QC my transcript assembly with SQANTI"

## Example Prompts

### Processing

> "Run the Iso-Seq3 pipeline on my subreads"

> "Cluster my refined Iso-Seq reads"

### Quality Control

> "Classify isoforms with SQANTI3"

> "How many novel isoforms did I find?"

### Quantification

> "Quantify isoform expression from Iso-Seq"

> "Compare isoform usage between samples"

## What the Agent Will Do

1. Generate CCS reads from subreads
2. Remove primers and polyA tails
3. Cluster into high-quality isoforms
4. Classify against reference annotation
5. Report novel and known isoforms

## Tips

- **CCS quality** - min-rq 0.9 is standard for Iso-Seq
- **Primers** - Must match your library prep kit
- **SQANTI categories** - FSM=known, NIC/NNC=novel
- **Collapse** - Removes redundancy before quantification
- **HQ vs LQ** - Use high-quality (HQ) transcripts for analysis
