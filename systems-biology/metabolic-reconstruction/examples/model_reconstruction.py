'''Run CarveMe/gapseq and sanity-check the draft. The inspection runs on a loaded SBML draft;
the reconstruction commands are shown as the CLI calls they wrap (carve/gapseq must be installed).
'''
# Reference: CarveMe 1.6+, gapseq 1.2+, cobrapy 0.29+ | Verify API if version differs

import subprocess
import cobra

TYPICAL_RXN_RANGE = (1000, 2500)   # typical bacterial draft; far outside flags an annotation problem


def carveme_command(fasta_path, output_path, universe='bacteria', gapfill_medium='M9'):
    '''Build the CarveMe command. Gram/universe is a -u VALUE, not a --grampos flag.'''
    cmd = ['carve', str(fasta_path), '-o', str(output_path), '-u', universe]
    if gapfill_medium:
        cmd += ['--gapfill', gapfill_medium]
    return cmd


def run_carveme(fasta_path, output_path, universe='bacteria', gapfill_medium='M9'):
    subprocess.run(carveme_command(fasta_path, output_path, universe, gapfill_medium), check=True)
    return output_path


def inspect_draft(model):
    '''Inventory the draft's soft spots before curation.'''
    orphans = [r for r in model.reactions if not r.genes]   # gap-filled, spontaneous, or transport
    namespace = 'ModelSEED' if any(m.id.startswith('cpd') for m in model.metabolites[:50]) else 'BiGG-like'
    return {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes),
        'exchanges': len(model.exchanges),
        'orphan_reactions': len(orphans),
        'orphan_fraction': len(orphans) / len(model.reactions),
        'grows': model.slim_optimize() > 1e-3,   # true by construction if gap-filled
        'namespace': namespace,
    }


def main():
    # A real run: run_carveme('genome.faa', 'model.xml', universe='gramneg', gapfill_medium='M9')
    # Here the textbook model stands in for a loaded draft so the inspection is runnable.
    model = cobra.io.load_model('textbook')
    m = inspect_draft(model)
    print('=== Draft inventory (curate the soft spots, not the growth number) ===')
    for k, v in m.items():
        print(f'  {k}: {v:.3f}' if isinstance(v, float) else f'  {k}: {v}')

    lo, hi = TYPICAL_RXN_RANGE
    if not lo <= m['reactions'] <= hi:
        print(f'  NOTE: {m["reactions"]} reactions is outside the typical {lo}-{hi} draft range (small demo model).')
    print('  Reminder: "grows" is guaranteed by gap-filling and is NOT evidence the model is correct.')


if __name__ == '__main__':
    main()
