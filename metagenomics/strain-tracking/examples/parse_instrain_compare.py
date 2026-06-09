"""Apply the shared-strain definition to an inStrain genomeWide_compare table.

Read genome-level calls from inStrain's `genomeWide_compare.tsv` (breadth column `percent_compared`);
the per-scaffold `comparisonsTable.tsv` instead names the breadth column `percent_genome_compared`.

A shared strain is a THRESHOLD claim, not a discovery: a pair shares a strain only when popANI is at
or above the cutoff AND enough of the genome was co-covered to trust it. Report the threshold and the
percent_compared; report co-detection separately, because absence is not absence.
"""
# Reference: pandas 2.2+ | Verify API if version differs
import sys
import pandas as pd

POPANI_THRESHOLD = 0.99999    # Olm 2021 same-strain cutoff; IS the operational definition
MIN_BREADTH = 0.50            # percent_compared floor; below this a genome is not confidently present


def call_shared_strains(compare_df, popani=POPANI_THRESHOLD, min_breadth=MIN_BREADTH):
    """compare_df: inStrain genomeWide_compare table with popANI and percent_compared columns."""
    comparable = compare_df[compare_df['percent_compared'] >= min_breadth].copy()
    comparable['shared_strain'] = comparable['popANI'] >= popani
    dropped = len(compare_df) - len(comparable)
    return comparable, dropped


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if path:
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.DataFrame({
            'genome': ['mag_1', 'mag_2', 'mag_3', 'mag_4'],
            'name1': ['A', 'A', 'A', 'A'], 'name2': ['B', 'B', 'B', 'B'],
            'popANI': [0.999995, 0.9990, 0.999999, 0.99998],
            'percent_compared': [0.72, 0.61, 0.30, 0.55],
        })

    shared, dropped = call_shared_strains(df)
    print(f'Threshold: popANI >= {POPANI_THRESHOLD}, breadth >= {MIN_BREADTH}')
    print(f'Pairs dropped for insufficient co-coverage (absence is not absence): {dropped}')
    print(shared[['genome', 'popANI', 'percent_compared', 'shared_strain']].to_string(index=False))
    n = int(shared['shared_strain'].sum())
    print(f'\nShared-strain calls: {n} (mag_3 was excluded - only 30% co-covered, not "not shared")')
