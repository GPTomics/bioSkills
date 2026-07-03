'''Detect and rank ligand-binding pockets de novo on an apo structure.

Runs fpocket and P2Rank on a PDB file, parses their ranked pocket lists, and
prints druggability/ligandability. A ranked pocket is a HYPOTHESIS, not a
binding site: geometry finds every concavity (including crystallization-additive
sites), and druggability scores were trained on holo/druggable sets, so they
under-score apo, shallow, and cryptic pockets. Both fpocket and P2Rank are
external CLI tools; the parse functions are pure and testable without them.
'''
# Reference: fpocket 4.1+, P2Rank 2.4+ | Verify CLI flags if version differs
import csv
import subprocess
import sys
from pathlib import Path

DRUGGABILITY_TRIAGE = 0.5  # rule-of-thumb: above ~0.5 a drug-like binder is plausible
                           # (Schmidtke 2010 J Med Chem 53:5858); apo/cryptic sites fall below
                           # it while still being real, so it is triage, not a cutoff law


def run_fpocket(pdb_path):
    subprocess.run(['fpocket', '-f', pdb_path], check=True)
    stem = Path(pdb_path).stem
    return Path(pdb_path).with_name(f'{stem}_out')  # fpocket writes <stem>_out/ beside input


def parse_fpocket_info(out_dir):
    info = next(Path(out_dir).glob('*_info.txt'))
    pockets, current = [], None
    for raw in info.read_text().splitlines():
        line = raw.strip()
        if line.startswith('Pocket'):
            current = {'pocket': int(line.split()[1])}
            pockets.append(current)
        elif ':' in line and current is not None:
            key, val = (part.strip() for part in line.split(':', 1))
            current[key] = float(val)
    return pockets


def run_p2rank(pdb_path, out_dir='p2rank_out'):
    subprocess.run(['prank', 'predict', '-f', pdb_path, '-o', out_dir], check=True)
    return Path(out_dir) / f'{Path(pdb_path).name}_predictions.csv'  # keeps input extension


def parse_p2rank(csv_path):
    with open(csv_path) as fh:
        reader = csv.DictReader(fh, skipinitialspace=True)  # P2Rank pads columns with spaces
        rows = [{k.strip(): v.strip() for k, v in row.items()} for row in reader]
    return [{'rank': int(r['rank']), 'score': float(r['score']),
             'probability': float(r['probability'])} for r in rows]


def rank_fpocket_pockets(out_dir):
    pockets = sorted(parse_fpocket_info(out_dir),
                     key=lambda p: p['Druggability Score'], reverse=True)
    druggable = [p for p in pockets if p['Druggability Score'] >= DRUGGABILITY_TRIAGE]
    return pockets, druggable


def main(pdb_path):
    fpocket_out = run_fpocket(pdb_path)
    ranked, druggable = rank_fpocket_pockets(fpocket_out)
    print(f'fpocket: {len(ranked)} pockets, {len(druggable)} above triage {DRUGGABILITY_TRIAGE}')
    for p in ranked[:5]:
        print(f"  pocket {p['pocket']}: drug={p['Druggability Score']:.2f} score={p['Score']:.2f}")

    predictions = parse_p2rank(run_p2rank(pdb_path))
    print(f'P2Rank: {len(predictions)} pockets')
    for p in predictions[:5]:
        print(f"  rank {p['rank']}: probability={p['probability']:.2f} score={p['score']:.2f}")

    # Agreement between geometry (fpocket) and ML (P2Rank) on the top site raises confidence;
    # every ranked pocket remains a hypothesis to corroborate against holo/known biology.


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'protein.pdb')
