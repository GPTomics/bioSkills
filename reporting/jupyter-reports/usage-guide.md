# Jupyter Reports Usage Guide

## Overview

This guide covers creating reproducible Jupyter notebooks with parameterization for automated analysis pipelines.

## Prerequisites

```bash
pip install jupyter papermill nbconvert
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a parameterized Jupyter notebook for my DE analysis"
- "Run my analysis notebook on multiple samples"
- "Convert my notebook to HTML for sharing"
- "Set up an automated notebook pipeline"

## Example Prompts

### Creating Templates

> "Create a Jupyter notebook template for RNA-seq QC that accepts sample ID and input directory as parameters"

> "Make my existing analysis notebook parameterized so I can run it on different datasets"

### Batch Execution

> "Run my QC notebook on all 12 samples and generate separate reports"

> "Execute my analysis template with different FDR thresholds and compare results"

### Export

> "Convert all my analysis notebooks to HTML reports"

> "Create a PDF report from my notebook with the code hidden"

## What the Agent Will Do

1. Create or modify notebook with parameter cells
2. Set up papermill execution code
3. Configure output formats and locations
4. Handle batch processing if multiple inputs
5. Generate final reports in requested format

## Tips

- Tag parameter cells in Jupyter (View > Cell Toolbar > Tags)
- Use meaningful parameter names and default values
- Include markdown cells explaining each analysis step
- Clear outputs before version control commits
- Test templates with default parameters first
