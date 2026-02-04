# Reaction Enumeration - Usage Guide

## Overview
Generate virtual compound libraries using reaction SMARTS transformations. Enumerate combinatorial libraries from building blocks for lead optimization.

## Prerequisites
```bash
pip install rdkit
```

## Quick Start
Tell your AI agent what you want to do:
- "Generate an amide library from my acids and amines"
- "Enumerate Suzuki coupling products from my building blocks"
- "Create a virtual library of 1000 compounds"
- "Apply R-group decoration to my core scaffold"

## Example Prompts

### Single Reactions
> "Run amide coupling between my carboxylic acids and amines."

> "Enumerate products from Buchwald amination with these building blocks."

### Combinatorial Libraries
> "Create a combinatorial library using these three sets of building blocks."

> "Generate all possible products from this multi-step synthesis."

### Validation
> "Filter the enumerated products for drug-likeness."

> "Remove any products with reactive groups or invalid valences."

## What the Agent Will Do
1. Parse reaction SMARTS and validate
2. Generate all reactant combinations
3. Run reactions to produce products
4. Sanitize and deduplicate products
5. Validate products (MW, valence, reactive groups)
6. Export enumerated library

## Tips
- Reaction SMARTS use atom mapping: [C:1] maps between reactants and products
- Always validate reaction with rxn.Validate() before running
- Products may need sanitization (SanitizeMol) after reaction
- Deduplicate products by canonical SMILES
- Common reactions: amide coupling, Suzuki, reductive amination, ester formation
- Filter products for drug-likeness after enumeration

## Related Skills
- molecular-io - Save enumerated libraries
- molecular-descriptors - Filter by properties
- admet-prediction - Screen for drug-likeness
