# Transmission Inference - Usage Guide

## Overview

Infer pathogen transmission chains and identify who-infected-whom using TransPhylo and genomic epidemiology methods.

## Prerequisites

```r
# R
install.packages('TransPhylo')
install.packages('ape')
```

```bash
pip install networkx pandas matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Identify transmission clusters in my sequences"
- "Who is the likely source of this outbreak?"
- "Find superspreaders in my outbreak data"

## Example Prompts

### Transmission Chains

> "Infer who infected whom from my dated phylogeny"

> "Build a transmission network from these outbreak sequences"

### Superspreader Analysis

> "Identify superspreading events in this outbreak"

> "Which cases were responsible for most transmissions?"

### Source Attribution

> "Who was the index case in this outbreak?"

> "Trace the transmission chain back to the source"

## What the Agent Will Do

1. Load dated phylogeny and sample metadata
2. Run TransPhylo MCMC inference
3. Extract consensus transmission tree
4. Identify transmission pairs with confidence
5. Calculate R0 and identify superspreaders
6. Visualize transmission network

## Tips

- **Input** - Requires time-scaled tree (from TreeTime/BEAST)
- **Generation time** - Pathogen-specific; affects inference
- **Unsampled cases** - TransPhylo can infer unsampled intermediates
- **Superspreaders** - Typically defined as >3 secondary cases
- **Confidence** - Lower SNP distance = higher confidence in link
