# R Markdown Reports Usage Guide

## Overview

R Markdown combines code, results, and narrative into reproducible documents. Output to HTML, PDF, Word, or presentations.

## Quick Start Prompts

- "Create an RMarkdown template for RNA-seq analysis"
- "Add interactive tables to my report"
- "Set up parameterized reports for multiple samples"
- "Cache long-running computations"

## Document Types

| Output | Best For |
|--------|----------|
| html_document | Interactive reports, sharing |
| pdf_document | Publication, archiving |
| word_document | Collaborator editing |
| ioslides_presentation | Slides |

## Workflow

1. **Create .Rmd** file with YAML header
2. **Write** narrative and code chunks
3. **Knit** to output format
4. **Share** HTML/PDF with collaborators

## Requirements

```r
install.packages(c('rmarkdown', 'knitr', 'DT', 'kableExtra'))
# For PDF: install tinytex
tinytex::install_tinytex()
```

## Key Features

- **Code folding** - Hide/show code in HTML
- **Tabs** - Organize related content
- **Caching** - Speed up re-knitting
- **Parameters** - Reusable templates

## Related Skills

- **reporting/quarto-reports** - Next-gen alternative
- **data-visualization/ggplot2-fundamentals** - Figures
