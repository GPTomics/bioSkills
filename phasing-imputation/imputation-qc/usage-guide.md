# Imputation QC - Usage Guide

## Overview
Quality control of imputed data is essential before downstream analysis, as poor quality imputation can introduce false positives in GWAS and bias effect estimates.

## Prerequisites
```bash
# bcftools for VCF manipulation
conda install -c bioconda bcftools

# Python packages for QC analysis
pip install pandas numpy matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Run QC on my imputed VCF and generate quality plots"
- "Filter my imputed data by INFO score and MAF"
- "Check concordance between imputed and typed variants"

## Example Prompts
### Quality Assessment
> "Analyze imputation quality in my VCF and generate a report with INFO score distributions"

### Filtering
> "Filter my imputed VCF to keep only variants with R2 > 0.8 and MAF > 0.01"

### Concordance Check
> "Validate imputation accuracy by checking concordance at typed variant positions"

### MAF Analysis
> "Plot INFO score vs minor allele frequency to assess imputation quality across the frequency spectrum"

## What the Agent Will Do
1. Extract INFO scores and allele frequencies from imputed VCF
2. Calculate summary statistics (mean R2, proportion passing thresholds)
3. Generate QC plots (R2 histogram, R2 vs MAF scatter, cumulative distribution)
4. Filter VCF by specified quality thresholds
5. Write summary report

## Tips
- Always filter by INFO score before GWAS (typically R2 > 0.3)
- Rare variants impute less accurately - expect lower R2 for MAF < 0.01
- Check that reference panel ancestry matches your study
- Verify concordance at typed positions to validate imputation
- Plot INFO by chromosome to detect systematic issues
- Use stricter thresholds (R2 > 0.8) for fine-mapping and PRS

## Key Quality Metrics

### INFO Score (R2) Interpretation
| R2 Score | Quality | Typical Use |
|----------|---------|-------------|
| > 0.9 | Excellent | All analyses |
| 0.8 - 0.9 | Good | Fine-mapping, PRS |
| 0.5 - 0.8 | Moderate | GWAS |
| 0.3 - 0.5 | Low | GWAS discovery only |
| < 0.3 | Poor | Exclude |

## Complete QC Workflow

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess

class ImputationQC:
    def __init__(self, vcf_path):
        self.vcf = vcf_path
        self.info = None
        self.stats = {}

    def extract_info(self, output_file='info_scores.txt'):
        cmd = f"bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/DR2\\t%INFO/AF\\n' {self.vcf} > {output_file}"
        subprocess.run(cmd, shell=True, check=True)
        self.info = pd.read_csv(output_file, sep='\t', names=['CHR', 'POS', 'ID', 'REF', 'ALT', 'R2', 'AF'])
        self.info['MAF'] = self.info['AF'].apply(lambda x: min(x, 1-x) if pd.notna(x) else np.nan)
        return self.info

    def calculate_stats(self):
        self.stats = {
            'total_variants': len(self.info),
            'mean_r2': self.info['R2'].mean(),
            'median_r2': self.info['R2'].median(),
            'pct_r2_above_03': 100 * (self.info['R2'] >= 0.3).mean(),
            'pct_r2_above_08': 100 * (self.info['R2'] >= 0.8).mean(),
        }
        return self.stats

    def plot_qc(self, output_prefix):
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        axes[0, 0].hist(self.info['R2'].dropna(), bins=50, edgecolor='black', alpha=0.7)
        axes[0, 0].axvline(0.3, color='red', linestyle='--', label='R2=0.3')
        axes[0, 0].axvline(0.8, color='orange', linestyle='--', label='R2=0.8')
        axes[0, 0].set_xlabel('INFO Score (R2)')
        axes[0, 0].set_ylabel('Variant Count')
        axes[0, 0].legend()

        sample = self.info.dropna().sample(min(50000, len(self.info)))
        axes[0, 1].scatter(sample['MAF'], sample['R2'], alpha=0.1, s=1)
        axes[0, 1].set_xlabel('Minor Allele Frequency')
        axes[0, 1].set_ylabel('INFO Score (R2)')
        axes[0, 1].axhline(0.3, color='red', linestyle='--')

        bins = [0, 0.001, 0.01, 0.05, 0.1, 0.5]
        self.info['MAF_bin'] = pd.cut(self.info['MAF'], bins=bins)
        self.info.boxplot(column='R2', by='MAF_bin', ax=axes[1, 0])
        axes[1, 0].set_xlabel('MAF Bin')
        axes[1, 0].set_ylabel('INFO Score')
        plt.suptitle('')

        sorted_r2 = np.sort(self.info['R2'].dropna())
        axes[1, 1].plot(sorted_r2, np.arange(len(sorted_r2)) / len(sorted_r2))
        axes[1, 1].axvline(0.3, color='red', linestyle='--')
        axes[1, 1].set_xlabel('INFO Score (R2)')
        axes[1, 1].set_ylabel('Cumulative Proportion')

        plt.tight_layout()
        plt.savefig(f'{output_prefix}_qc_plots.png', dpi=150)
        plt.close()

    def filter_vcf(self, r2_threshold=0.3, maf_threshold=0.01, output=None):
        if output is None:
            output = self.vcf.replace('.vcf.gz', f'_r2{r2_threshold}_maf{maf_threshold}.vcf.gz')
        cmd = f"bcftools view -i 'INFO/DR2 >= {r2_threshold} && INFO/AF >= {maf_threshold} && INFO/AF <= {1-maf_threshold}' {self.vcf} -Oz -o {output}"
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(f'bcftools index {output}', shell=True, check=True)
        return output

    def generate_report(self, output_prefix):
        if self.info is None:
            self.extract_info()
        self.calculate_stats()
        self.plot_qc(output_prefix)
        with open(f'{output_prefix}_summary.txt', 'w') as f:
            f.write('Imputation QC Report\n')
            f.write('=' * 50 + '\n\n')
            for k, v in self.stats.items():
                f.write(f'{k}: {v:.4f}\n' if isinstance(v, float) else f'{k}: {v}\n')

# Usage
qc = ImputationQC('imputed.vcf.gz')
qc.generate_report('imputation_qc')
filtered_vcf = qc.filter_vcf(r2_threshold=0.3, maf_threshold=0.01)
```

## Concordance with Typed Variants

```bash
# Extract typed variants
bcftools view -i 'INFO/TYPED=1' imputed.vcf.gz -Oz -o typed_imputed.vcf.gz

# Compare with original
bcftools gtcheck -g original.vcf.gz typed_imputed.vcf.gz > concordance.txt

# Parse results
grep "^DC" concordance.txt  # Discordance rate
```

## Filter Commands

```bash
# Check INFO score distribution
bcftools query -f '%INFO/DR2\n' imputed.vcf.gz | \
    awk '{if($1<0.3) low++; else if($1<0.8) med++; else high++}
    END {print "Low (<0.3):", low; print "Medium (0.3-0.8):", med; print "High (>0.8):", high}'

# Filter by INFO and MAF
bcftools view -i 'INFO/DR2 > 0.3' imputed.vcf.gz -Oz -o filtered.vcf.gz
bcftools view -i 'MAF > 0.01' filtered.vcf.gz -Oz -o common.vcf.gz
```

## Troubleshooting

### Low R2 Overall
- Check reference panel ancestry match
- Verify strand alignment
- Check genotyping quality of input data

### Low R2 for Rare Variants
- Expected behavior - rare variants impute poorly
- Use larger/better matched reference panel
- Consider TOPMed for rare variants

### Systematic Chromosome Differences
- Check genetic map coverage
- Look for reference panel issues
- May indicate liftover problems
