# Imputation QC Usage Guide

## Overview

Quality control of imputed data is essential before downstream analysis. Poor quality imputation can introduce false positives in GWAS and bias effect estimates.

## Key Quality Metrics

### INFO Score (R2)
- Measures squared correlation between true and imputed genotypes
- Range: 0 (poor) to 1 (perfect)
- Based on ratio of observed to expected variance

### Interpretation
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
from pathlib import Path

class ImputationQC:
    def __init__(self, vcf_path):
        self.vcf = vcf_path
        self.info = None
        self.stats = {}

    def extract_info(self, output_file='info_scores.txt'):
        '''Extract INFO scores from VCF.'''
        cmd = f"bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/DR2\\t%INFO/AF\\n' {self.vcf} > {output_file}"
        subprocess.run(cmd, shell=True, check=True)

        self.info = pd.read_csv(output_file, sep='\t',
            names=['CHR', 'POS', 'ID', 'REF', 'ALT', 'R2', 'AF'])
        self.info['MAF'] = self.info['AF'].apply(lambda x: min(x, 1-x) if pd.notna(x) else np.nan)

        return self.info

    def calculate_stats(self):
        '''Calculate summary statistics.'''
        self.stats = {
            'total_variants': len(self.info),
            'mean_r2': self.info['R2'].mean(),
            'median_r2': self.info['R2'].median(),
            'pct_r2_above_03': 100 * (self.info['R2'] >= 0.3).mean(),
            'pct_r2_above_08': 100 * (self.info['R2'] >= 0.8).mean(),
            'mean_maf': self.info['MAF'].mean(),
        }
        return self.stats

    def plot_qc(self, output_prefix):
        '''Generate QC plots.'''
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # R2 histogram
        axes[0, 0].hist(self.info['R2'].dropna(), bins=50, edgecolor='black', alpha=0.7)
        axes[0, 0].axvline(0.3, color='red', linestyle='--', label='R2=0.3')
        axes[0, 0].axvline(0.8, color='orange', linestyle='--', label='R2=0.8')
        axes[0, 0].set_xlabel('INFO Score (R2)')
        axes[0, 0].set_ylabel('Variant Count')
        axes[0, 0].set_title('INFO Score Distribution')
        axes[0, 0].legend()

        # R2 vs MAF
        sample = self.info.dropna().sample(min(50000, len(self.info)))
        axes[0, 1].scatter(sample['MAF'], sample['R2'], alpha=0.1, s=1)
        axes[0, 1].set_xlabel('Minor Allele Frequency')
        axes[0, 1].set_ylabel('INFO Score (R2)')
        axes[0, 1].set_title('INFO vs MAF')
        axes[0, 1].axhline(0.3, color='red', linestyle='--')

        # R2 by MAF bin
        bins = [0, 0.001, 0.01, 0.05, 0.1, 0.5]
        self.info['MAF_bin'] = pd.cut(self.info['MAF'], bins=bins)
        self.info.boxplot(column='R2', by='MAF_bin', ax=axes[1, 0])
        axes[1, 0].set_xlabel('MAF Bin')
        axes[1, 0].set_ylabel('INFO Score')
        axes[1, 0].set_title('INFO by MAF Bin')
        plt.suptitle('')

        # Cumulative distribution
        sorted_r2 = np.sort(self.info['R2'].dropna())
        axes[1, 1].plot(sorted_r2, np.arange(len(sorted_r2)) / len(sorted_r2))
        axes[1, 1].axvline(0.3, color='red', linestyle='--')
        axes[1, 1].axhline((self.info['R2'] < 0.3).mean(), color='red', linestyle=':')
        axes[1, 1].set_xlabel('INFO Score (R2)')
        axes[1, 1].set_ylabel('Cumulative Proportion')
        axes[1, 1].set_title('Cumulative Distribution')

        plt.tight_layout()
        plt.savefig(f'{output_prefix}_qc_plots.png', dpi=150)
        plt.close()

    def filter_vcf(self, r2_threshold=0.3, maf_threshold=0.01, output=None):
        '''Filter VCF by quality thresholds.'''
        if output is None:
            output = self.vcf.replace('.vcf.gz', f'_r2{r2_threshold}_maf{maf_threshold}.vcf.gz')

        cmd = f"bcftools view -i 'INFO/DR2 >= {r2_threshold} && INFO/AF >= {maf_threshold} && INFO/AF <= {1-maf_threshold}' {self.vcf} -Oz -o {output}"
        subprocess.run(cmd, shell=True, check=True)

        cmd_index = f"bcftools index {output}"
        subprocess.run(cmd_index, shell=True, check=True)

        return output

    def generate_report(self, output_prefix):
        '''Generate full QC report.'''
        if self.info is None:
            self.extract_info()

        self.calculate_stats()
        self.plot_qc(output_prefix)

        # Write summary
        with open(f'{output_prefix}_summary.txt', 'w') as f:
            f.write('Imputation QC Report\n')
            f.write('=' * 50 + '\n\n')
            for k, v in self.stats.items():
                if isinstance(v, float):
                    f.write(f'{k}: {v:.4f}\n')
                else:
                    f.write(f'{k}: {v}\n')

            f.write('\n\nVariants by R2 threshold:\n')
            for thresh in [0.1, 0.3, 0.5, 0.8, 0.9]:
                n = (self.info['R2'] >= thresh).sum()
                pct = 100 * n / len(self.info)
                f.write(f'  R2 >= {thresh}: {n:,} ({pct:.1f}%)\n')

        print(f'Report written to {output_prefix}_summary.txt')
        print(f'Plots written to {output_prefix}_qc_plots.png')

# Usage
qc = ImputationQC('imputed.vcf.gz')
qc.generate_report('imputation_qc')
filtered_vcf = qc.filter_vcf(r2_threshold=0.3, maf_threshold=0.01)
```

## Concordance with Typed Variants

Validate imputation by checking accuracy at typed positions:

```bash
# Extract typed variants
bcftools view -i 'INFO/TYPED=1' imputed.vcf.gz -Oz -o typed_imputed.vcf.gz

# Compare with original
bcftools gtcheck -g original.vcf.gz typed_imputed.vcf.gz > concordance.txt

# Parse results
grep "^DC" concordance.txt  # Discordance rate
```

## Batch Effects

Check for systematic differences between batches:

```python
def check_batch_effects(vcf, sample_batch_file):
    '''Check for batch effects in imputation quality.'''
    # sample_batch_file: sample\tbatch format

    batches = pd.read_csv(sample_batch_file, sep='\t', names=['sample', 'batch'])

    # Extract sample-level metrics
    # e.g., mean heterozygosity, missing rate

    # Compare between batches
    # Flag if significant differences
    pass
```

## Troubleshooting

### Low R2 Overall
- Check reference panel ancestry match
- Verify strand alignment
- Check genotyping quality

### Low R2 for Rare Variants
- Expected: rare variants impute poorly
- Use larger/better matched reference
- Consider TOPMed for rare variants

### Systematic Chromosome Differences
- Check genetic map coverage
- Look for reference panel issues
- May indicate liftover problems
