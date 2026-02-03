# Quarto Reports - Usage Guide

## Overview
Quarto is a next-generation scientific publishing system supporting R, Python, Julia, and Observable with enhanced features over R Markdown.

## Prerequisites
```bash
# Install Quarto - download from https://quarto.org/docs/download/
quarto check

# For Python
pip install jupyter matplotlib

# For R
install.packages(c('knitr', 'rmarkdown'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a Quarto document for my Python analysis"
- "Set up a multi-format report (HTML + PDF)"
- "Build a Quarto website for my project"

## Example Prompts

### Basic Documents
> "Create a Quarto document for my scRNA-seq analysis in Python"

> "Convert my Jupyter notebook to a Quarto report with better formatting"

### Multi-Format Output
> "Set up a Quarto project that renders to both HTML and PDF"

> "Create a report that generates Word output for my collaborators"

### Interactive Features
> "Add interactive Plotly figures to my Quarto document"

> "Create a Quarto dashboard with tabsets for different samples"

### Projects and Websites
> "Set up a Quarto website to document my analysis pipeline"

> "Create a Quarto book for my lab protocol collection"

## What the Agent Will Do
1. Create a .qmd file with appropriate YAML front matter
2. Configure code cell options using `#|` syntax
3. Set up cross-references for figures, tables, and equations
4. Configure rendering for your target output formats
5. Add callouts, tabs, or other layout features as needed

## Quarto vs R Markdown

| Feature | R Markdown | Quarto |
|---------|------------|--------|
| Languages | R primary | R, Python, Julia |
| Presentations | Limited | RevealJS, PowerPoint |
| Websites | blogdown | Built-in |
| Books | bookdown | Built-in |

## Tips
- Use `#| freeze: true` to cache computed outputs across renders
- Cross-references use `@fig-label` syntax for figures and `@tbl-label` for tables
- Use callout blocks (note, warning, tip) for important information
- Set `execute: freeze: auto` in `_quarto.yml` for project-wide caching
- Preview live changes with `quarto preview` during development

## Related Skills
- reporting/rmarkdown-reports - R Markdown for R-only workflows
- data-visualization/ggplot2-fundamentals - Creating visualizations for reports
