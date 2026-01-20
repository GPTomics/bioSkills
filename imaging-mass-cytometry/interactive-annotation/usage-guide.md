# Interactive Annotation - Usage Guide

## Overview
Interactive annotation enables expert-guided cell type labeling in IMC images for training classifiers, validating automated results, and identifying rare populations.

## Prerequisites
```bash
pip install napari napari-imc scikit-learn
# For GUI support
pip install pyqt5
```

## Quick Start
Tell your AI agent what you want to do:
- "Set up napari for interactive cell annotation"
- "Create training data for a cell type classifier"
- "Validate my automated phenotyping results"

## Example Prompts

### Annotation Setup
> "Load my IMC image and segmentation mask in napari for annotation"

> "Create marker overlays for cell type annotation in napari"

### Manual Labeling
> "Set up napari to label CD8 T cells, macrophages, and tumor cells"

> "Create a labels layer for interactive cell type annotation"

### Classifier Training
> "Train a cell type classifier from my manual annotations"

> "Propagate my annotations to unlabeled cells using a random forest"

### Validation
> "Generate gallery plots showing cells from each annotated type"

> "Review low-confidence predictions from my classifier"

## What the Agent Will Do
1. Load multichannel IMC image and segmentation mask
2. Create marker overlays for visual guidance
3. Set up napari with labels layer for annotation
4. Enable interactive cell labeling with keyboard shortcuts
5. Export annotations as labeled mask or CSV
6. Optionally train classifier and propagate to unlabeled cells

## Tips
- Define clear cell type criteria before starting annotation
- Napari shortcuts: 1-9 switch labels, P for paint, F for fill, E for erase
- Keep DNA channel always visible as nuclear reference
- Annotate at least 50-100 cells per type for classifier training
- Include edge cases and ambiguous cells in training data
- Take breaks during long sessions to maintain consistency
- Save work frequently and document annotation decisions
