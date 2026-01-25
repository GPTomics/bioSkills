#!/usr/bin/env python3
"""
Batch notebook execution with papermill
Demonstrates running parameterized notebooks on multiple samples
"""

import papermill as pm
from pathlib import Path
import json

# =============================================================================
# Basic Notebook Execution
# =============================================================================

def run_single_notebook(template, output, parameters):
    """Execute a notebook with given parameters."""
    pm.execute_notebook(
        template,
        output,
        parameters=parameters,
        kernel_name='python3'
    )
    print(f'Generated: {output}')


# =============================================================================
# Batch Processing Multiple Samples
# =============================================================================

def batch_process_samples(template_path, output_dir, samples):
    """Run analysis notebook on multiple samples.

    Args:
        template_path: Path to parameterized notebook template
        output_dir: Directory for output notebooks
        samples: List of dicts with sample info and parameters
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for sample in samples:
        sample_id = sample['sample_id']
        output_path = output_dir / f'{sample_id}_analysis.ipynb'

        try:
            pm.execute_notebook(
                template_path,
                str(output_path),
                parameters=sample,
                kernel_name='python3'
            )
            results.append({'sample': sample_id, 'status': 'success', 'output': str(output_path)})
        except pm.PapermillExecutionError as e:
            results.append({'sample': sample_id, 'status': 'failed', 'error': str(e)})
            print(f'Failed: {sample_id} - {e}')

    return results


# =============================================================================
# Example Usage
# =============================================================================

if __name__ == '__main__':
    # Define samples with their parameters
    samples = [
        {
            'sample_id': 'sample_A',
            'input_file': 'data/sample_A_counts.csv',
            'condition': 'treated',
            'fdr_threshold': 0.05
        },
        {
            'sample_id': 'sample_B',
            'input_file': 'data/sample_B_counts.csv',
            'condition': 'control',
            'fdr_threshold': 0.05
        },
        {
            'sample_id': 'sample_C',
            'input_file': 'data/sample_C_counts.csv',
            'condition': 'treated',
            'fdr_threshold': 0.05
        }
    ]

    # Run batch processing
    # Assumes 'analysis_template.ipynb' exists with parameter cell
    results = batch_process_samples(
        template_path='analysis_template.ipynb',
        output_dir='reports/',
        samples=samples
    )

    # Save execution summary
    with open('reports/execution_summary.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print('\n=== Execution Summary ===')
    success = sum(1 for r in results if r['status'] == 'success')
    print(f'Successful: {success}/{len(results)}')


# =============================================================================
# Template Notebook Structure
# =============================================================================

# The template notebook should have a cell tagged as "parameters":
#
# ```python
# # Parameters (tag this cell as "parameters" in Jupyter)
# sample_id = "default_sample"
# input_file = "data/default.csv"
# condition = "control"
# fdr_threshold = 0.05
# output_dir = "results/"
# ```
#
# Papermill will inject new values for these variables when executing
