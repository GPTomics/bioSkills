# Reference: MAFFT 7.520+, PAL2NAL 14+ | Verify CLI flags if version differs

import subprocess
from pathlib import Path


def protein_guided_codon_alignment(protein_fasta, cds_fasta, output_fasta, output_format='fasta'):
    protein_aligned = Path(protein_fasta).stem + '_aligned.fasta'

    print('Step 1: Aligning protein sequences with MAFFT L-INS-i...')
    with open(protein_aligned, 'w') as out:
        subprocess.run(
            ['mafft', '--localpair', '--maxiterate', '1000', '--thread', '4', str(protein_fasta)],
            stdout=out, stderr=subprocess.PIPE, check=True
        )

    fmt_flag = {'fasta': 'fasta', 'paml': 'paml', 'clustal': 'clustalw', 'codon': 'codon'}
    pal2nal_fmt = fmt_flag.get(output_format, 'fasta')

    print('Step 2: Threading DNA onto protein alignment with PAL2NAL...')
    with open(output_fasta, 'w') as out:
        subprocess.run(
            ['pal2nal.pl', protein_aligned, str(cds_fasta), '-output', pal2nal_fmt],
            stdout=out, stderr=subprocess.PIPE, check=True
        )

    print(f'Codon alignment written to {output_fasta}')


if __name__ == '__main__':
    protein_guided_codon_alignment(
        protein_fasta='proteins.fasta',
        cds_fasta='coding_sequences.fasta',
        output_fasta='codon_aligned.fasta',
        output_format='fasta'
    )
