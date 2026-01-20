# Interactive Annotation Usage Guide

## Overview

Interactive annotation allows expert-guided cell type labeling in IMC images, essential for training classifiers, validating automated results, and rare cell type identification.

## When to Use This Skill

- Training supervised classifiers
- Validating automated phenotyping
- Identifying rare cell populations
- Ground truth generation for benchmarking
- Complex tissue microenvironments

## Annotation Workflow

### 1. Image Preparation
- Load multichannel IMC image
- Load cell segmentation mask
- Create marker overlays for guidance

### 2. Manual Annotation
- Use napari for interactive labeling
- Select cells and assign types
- Use marker visualization as guide

### 3. Propagation (Optional)
- Use annotated cells to train classifier
- Propagate to unannotated cells
- Review low-confidence predictions

### 4. Validation
- Generate gallery plots by cell type
- Review representative cells
- Correct misclassifications

## Napari Usage

### Keyboard Shortcuts
| Key | Action |
|-----|--------|
| 1-9 | Switch to annotation label 1-9 |
| P | Paint mode |
| F | Fill mode |
| E | Erase mode |

### Recommended Setup
- DNA channel always visible (nuclear reference)
- Use additive blending for overlays
- Adjust contrast per channel
- Color-code by marker function

## Annotation Guidelines

### Cell Type Definitions
Define clear criteria before starting:
- Marker positive/negative thresholds
- Morphological features
- Spatial context rules

### Quality Standards
- Annotate at least 50-100 cells per type
- Include edge cases
- Balance across tissue regions
- Document ambiguous cases

## Training Data Requirements

| Cell Type | Minimum Cells | Recommended |
|-----------|--------------|-------------|
| Common types | 100 | 500+ |
| Rare types | 20 | 100+ |
| Ambiguous types | 50 | 200+ |

## Common Issues

### Overlapping markers
- Use multiple channels simultaneously
- Check spatial context
- Define hierarchy rules

### Segmentation errors
- Annotate only well-segmented cells
- Flag segmentation issues separately
- Consider re-segmentation for problem areas

### Inconsistency
- Take breaks during long sessions
- Review previous annotations
- Use validation galleries

## Best Practices

1. Define cell types before starting
2. Create marker overlays for each type
3. Start with clear examples
4. Annotate across full tissue area
5. Validate annotations regularly
6. Save work frequently
7. Document decisions

## References

- napari: napari.org
- IMC annotation: doi:10.1038/s41467-020-20430-7
- Training data guidelines: doi:10.1038/s41592-021-01308-y
