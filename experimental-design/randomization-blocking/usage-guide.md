# Randomization and Blocking Usage Guide

## Overview

This guide covers the structural design of biological experiments: deciding the experimental unit, randomizing treatments, replicating the right thing, and removing known nuisance variation by blocking. The central goal is that the analysis model mirrors how the experiment was actually run, so that the resulting inference is valid by construction. The single most common failure it prevents is pseudoreplication: counting cells, wells, or technical aliquots as independent replicates when the true sample size is the number of animals, donors, or cages.

## Prerequisites

```r
# R packages
install.packages(c('designit', 'lme4', 'lmerTest', 'pwr'))
```

Familiarity with the difference between a fixed effect (a question) and a random effect (a level of the randomization/sampling structure) is helpful but not required; the agent will explain the mapping.

## Quick Start

Tell your AI agent what you want to do:
- "What is the experimental unit in my study, and what is my real n?"
- "I have 10,000 cells from 3 mice per group; can I test on cells?"
- "Help me randomize 24 samples to 3 processing days without confounding"
- "Should this be a factorial, split-plot, or nested design?"
- "Write the mixed-model random-effects structure that matches my design"

## Example Prompts

### Experimental Unit and Pseudoreplication

> "I measured a marker in 200 cells from each of 4 control and 4 treated animals. A reviewer says my n is 4, not 1600. How should I analyze this?"

> "Treatment is delivered through the drinking water of co-housed cages. What is my experimental unit?"

### Randomization and Blocking

> "I can only process 8 of my 24 samples per day over 3 days. How do I assign and randomize them so that processing day is not confounded with condition?"

> "Set up a randomized complete block design for 2 conditions across 3 litters and give me the model formula to analyze it."

### Design Structure

> "My incubator can only hold one temperature at a time, but I want to test temperature and genotype. Is this a split-plot, and how does that change the analysis?"

> "I want to test genotype and drug together; design a factorial and tell me how to read the interaction."

## What the Agent Will Do

1. Trace the randomization to identify the experimental unit and the true sample size.
2. Flag pseudoreplication and recommend aggregation (pseudobulk) or a nested random effect.
3. Choose a design layout (CRD, RCBD, Latin square, incomplete block, factorial, split-plot, nested) for the question and constraints.
4. Generate a seeded, optionally restricted randomization including run order.
5. Write the matching mixed-model random-effects structure and note small-sample degrees-of-freedom corrections.

## Tips

- The experimental unit is the smallest entity independently assigned to a treatment; that count, not the number of measurements, is the sample size.
- Aggregate observational units to the experimental unit before testing, or model them as a nested random effect.
- A sequencing lane, chip, or incubator run is almost always a whole-plot factor; analyze it as a split-plot, not a flat factorial.
- Block only on factors with real between-block variation; blocking on noise wastes degrees of freedom.
- Always randomize run order, not just treatment assignment, and record the random seed.
- Analyze as randomized: a blocked or stratified design must include those factors in the model.

## Related Skills

- batch-design - Assigning samples to sequencing batches/lanes and batch-effect correction
- sample-size - The experimental unit defines what is replicated and counted
- power-analysis - Blocking and nesting change the effective error variance
- multiple-testing - The design fixes what counts as a family of tests
- single-cell/preprocessing - Pseudobulk aggregation to the donor for scRNA-seq
- clinical-biostatistics/power-and-sample-size - Randomization and design in the regulated-trial regime
