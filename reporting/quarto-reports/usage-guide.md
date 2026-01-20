# Quarto Reports Usage Guide

## Overview

Quarto is a next-generation scientific publishing system supporting R, Python, Julia, and Observable. It replaces R Markdown with enhanced features.

## Quick Start Prompts

- "Create a Quarto document for my Python analysis"
- "Set up a multi-format report (HTML + PDF)"
- "Add interactive figures to my Quarto document"
- "Build a Quarto website for my project"

## Why Quarto?

| Feature | R Markdown | Quarto |
|---------|------------|--------|
| Languages | R primary | R, Python, Julia |
| Presentations | Limited | RevealJS, PowerPoint |
| Websites | blogdown | Built-in |
| Books | bookdown | Built-in |

## Workflow

1. **Create** .qmd file
2. **Write** code and narrative
3. **Render** with `quarto render`
4. **Preview** with `quarto preview`

## Requirements

```bash
# Install Quarto
# Download from https://quarto.org/docs/download/

# Check installation
quarto check

# For Python
pip install jupyter matplotlib

# For R
install.packages(c('knitr', 'rmarkdown'))
```

## Key Features

- **Unified syntax** across languages
- **Code cell options** with `#|` prefix
- **Cross-references** for figures/tables
- **Callouts** for notes/warnings
- **Freeze** for reproducibility

## Related Skills

- **reporting/rmarkdown-reports** - R Markdown alternative
- **data-visualization/ggplot2-fundamentals** - Visualizations
