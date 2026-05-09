'''Build a structural multiple sequence alignment with Foldmason.

Produces both an amino-acid MSA (result_aa.fa) and a 3Di structural MSA (result_3di.fa).
Use --refine-iters > 0 for iterative refinement; --report-mode 1 for an HTML LDDT report,
--report-mode 2 for a machine-readable JSON LDDT report.
'''
# Reference: foldmason 1+ | Verify CLI flags if version differs

import subprocess
from pathlib import Path

def foldmason_msa(structures, result_prefix, tmp_dir='tmp/', refine_iters=100, report_mode=1):
    structure_paths = [str(p) for p in structures]
    cmd = [
        'foldmason', 'easy-msa',
        *structure_paths,
        result_prefix,
        tmp_dir,
        '--refine-iters', str(refine_iters),
        '--report-mode', str(report_mode),
    ]
    subprocess.run(cmd, check=True)
    if report_mode == 1:
        report_path = f'{result_prefix}.html'
    elif report_mode == 2:
        report_path = f'{result_prefix}.json'
    else:
        report_path = None
    return {
        'amino_msa': f'{result_prefix}_aa.fa',
        'structural_msa': f'{result_prefix}_3di.fa',
        'guide_tree': f'{result_prefix}.nw',
        'report': report_path,
    }

def summarize(amino_msa_path):
    n_seqs = 0
    length = 0
    with open(amino_msa_path) as f:
        for line in f:
            if line.startswith('>'):
                n_seqs += 1
            elif n_seqs == 1:
                length += len(line.strip())
    return n_seqs, length

if __name__ == '__main__':
    structures = sorted(Path('structures').glob('*.pdb'))
    print(f'Aligning {len(structures)} structures with Foldmason...')

    outputs = foldmason_msa(structures, 'family_msa', refine_iters=100)
    n_seqs, length = summarize(outputs['amino_msa'])
    print(f'\nMSA: {n_seqs} sequences, {length} columns')
    print(f'Amino-acid MSA: {outputs["amino_msa"]}')
    print(f'3Di MSA: {outputs["structural_msa"]}')
    print(f'Guide tree (Newick): {outputs["guide_tree"]}')
    if outputs['report']:
        print(f'LDDT report: {outputs["report"]}')
