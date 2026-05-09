'''Pairwise structural alignment with TMalign / USalign and parse the structured output.

Returns TM-score (normalised by shorter chain), RMSD, alignment length, and the
superposed mobile structure as a PDB. Uses -outfmt 2 for tabular output that is
robust to TM-align release versions.
'''
# Reference: TM-align 20220412+, US-align 20231222+ | Verify CLI flags if version differs

import subprocess

def tm_align(reference_pdb, mobile_pdb, output_pdb=None, multimer=False):
    binary = 'USalign' if multimer else 'TMalign'
    cmd = [binary, mobile_pdb, reference_pdb, '-outfmt', '2']
    if multimer:
        cmd += ['-mm', '1', '-ter', '0']
    if output_pdb:
        cmd += ['-o', output_pdb]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return parse_outfmt2(result.stdout)

def parse_outfmt2(stdout):
    '''Parse TMalign / USalign -outfmt 2 tabular output.

    Header columns: PDBchain1 PDBchain2 TM1 TM2 RMSD ID1 ID2 IDali L1 L2 Lali
    '''
    for line in stdout.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split()
        if len(fields) < 11:
            continue
        return {
            'pdb1': fields[0],
            'pdb2': fields[1],
            'tm1': float(fields[2]),
            'tm2': float(fields[3]),
            'rmsd': float(fields[4]),
            'id1': float(fields[5]),
            'id2': float(fields[6]),
            'idali': float(fields[7]),
            'length1': int(fields[8]),
            'length2': int(fields[9]),
            'length_align': int(fields[10]),
        }
    return None

def interpret_tmscore(tm):
    if tm > 0.8:
        return 'equivalent topology (homologous)'
    if tm > 0.5:
        return 'same fold'
    if tm > 0.2:
        return 'weak structural similarity'
    return 'statistically random'

if __name__ == '__main__':
    result = tm_align('reference.pdb', 'mobile.pdb', output_pdb='superposed.pdb')
    tm_short = max(result['tm1'], result['tm2'])
    print(f'TM-score (normalised by shorter chain): {tm_short:.3f}')
    print(f'TM-score normalised by chain 1: {result["tm1"]:.3f}')
    print(f'TM-score normalised by chain 2: {result["tm2"]:.3f}')
    print(f'RMSD: {result["rmsd"]:.2f} A')
    print(f'Sequence identity in alignment: {result["idali"]*100:.1f}%')
    print(f'Interpretation: {interpret_tmscore(tm_short)}')
