# R Markdown Reports - Usage Guide

## Overview
R Markdown combines code, results, and narrative into reproducible documents. Output to HTML, PDF, Word, or presentations.

## Prerequisites
```r
install.packages(c('rmarkdown', 'knitr', 'DT', 'kableExtra'))
# For PDF: install tinytex
tinytex::install_tinytex()
```

## Quick Start
Tell your AI agent what you want to do:
- "Create an RMarkdown template for RNA-seq analysis"
- "Add interactive tables to my report"
- "Set up parameterized reports for multiple samples"

## Example Prompts

### Basic Reports
> "Create an RMarkdown document for my DESeq2 analysis results"

> "Generate a PDF report summarizing my variant calling pipeline"

### Interactive Elements
> "Add interactive DT tables to display my gene expression results"

> "Include collapsible code chunks in my analysis report"

### Parameterized Reports
> "Set up a parameterized RMarkdown that generates reports for each sample"

> "Create a template report that accepts a sample ID as input"

### Advanced Features
> "Cache the long-running PCA computation in my report"

> "Add tabbed sections for different comparisons in my DE analysis"

## What the Agent Will Do
1. Create an .Rmd file with appropriate YAML header for your output format
2. Structure code chunks with proper options (eval, echo, cache)
3. Add narrative sections explaining the analysis
4. Configure output-specific features (interactive tables, code folding)
5. Set up figure and table cross-references if needed

## Document Types

| Output | Best For |
|--------|----------|
| html_document | Interactive reports, sharing |
| pdf_document | Publication, archiving |
| word_document | Collaborator editing |
| ioslides_presentation | Slides |

## Tips
- Use `cache=TRUE` on expensive code chunks to speed up re-knitting
- Set `code_folding: hide` in YAML to collapse code by default
- Use DT::datatable() for interactive tables in HTML output
- Parameters in YAML header enable batch report generation
- Use `fig.width` and `fig.height` chunk options to control figure sizing

## Related Skills
- reporting/quarto-reports - Next-gen alternative with multi-language support
- data-visualization/ggplot2-fundamentals - Creating publication-quality figures
