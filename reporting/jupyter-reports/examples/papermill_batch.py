#!/usr/bin/env python3
"""Batch papermill execution: per-sample notebooks, scrap aggregation, code-hidden render."""
# Reference: papermill 2.6+, scrapbook 0.5+, nbconvert 7.16+ | Verify API if version differs

import json
import subprocess
from pathlib import Path

import papermill as pm
import scrapbook as sb

# Per-cell timeout: papermill's default is forever, which silently hangs a pipeline
CELL_TIMEOUT = 1800


def run_one(template, output, parameters):
    """Execute one parameterized notebook in a fresh kernel; the output .ipynb is the evidence."""
    pm.execute_notebook(template, output, parameters=parameters, kernel_name='python3',
                        log_output=True, execution_timeout=CELL_TIMEOUT)


def render_report(notebook):
    """Render a code-hidden HTML report from an executed notebook (no LaTeX needed)."""
    subprocess.run(['jupyter', 'nbconvert', '--to', 'html', '--no-input', str(notebook)], check=True)


def batch_process(template, output_dir, samples):
    """Run the template over a sample sheet; preserve failures as artifacts, render survivors."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for sample in samples:
        sample_id = sample['sample_id']
        out = output_dir / f'{sample_id}.ipynb'
        try:
            run_one(template, str(out), sample)
            render_report(out)
            results.append({'sample': sample_id, 'status': 'success', 'output': str(out)})
        except pm.PapermillExecutionError as e:
            # papermill still wrote out with the traceback captured; keep it for forensics
            results.append({'sample': sample_id, 'status': 'failed', 'error': str(e)})
            print(f'Failed: {sample_id} - {e}')
    return results


def aggregate_scraps(output_dir):
    """Collect sb.glue()'d metrics from every executed notebook into one cohort table."""
    book = sb.read_notebooks(str(output_dir))
    return book.papermill_dataframe


if __name__ == '__main__':
    samples = [
        {'sample_id': 'sample_A', 'input_file': 'data/sample_A_counts.csv', 'condition': 'treated', 'fdr_threshold': 0.05},
        {'sample_id': 'sample_B', 'input_file': 'data/sample_B_counts.csv', 'condition': 'control', 'fdr_threshold': 0.05},
        {'sample_id': 'sample_C', 'input_file': 'data/sample_C_counts.csv', 'condition': 'treated', 'fdr_threshold': 0.05},
    ]
    results = batch_process('analysis_template.ipynb', 'reports/', samples)
    Path('reports/execution_summary.json').write_text(json.dumps(results, indent=2))
    success = sum(1 for r in results if r['status'] == 'success')
    print(f'Successful: {success}/{len(results)}')

    # cohort = aggregate_scraps('reports/')   # requires the template to sb.glue() its metrics

# Template notebook: tag ONE cell 'parameters' with defaults; papermill injects overrides after it.
#   # Parameters (tag this cell 'parameters')
#   sample_id = 'default'; input_file = 'data/default.csv'; condition = 'control'; fdr_threshold = 0.05
# To harvest metrics, add in the template: import scrapbook as sb; sb.glue('n_de_genes', n_sig)
