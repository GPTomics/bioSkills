'''Parse QUAST, BUSCO/compleasm, and Merqury outputs into one three-axis verdict.

No single number is quality, so this prints contiguity, completeness, and correctness
together and flags the documented traps (N50-only, high-D haplotigs, sub-Q40 accuracy).
'''
# Reference: pandas 2.2+ | Verify API if version differs
import re
import pandas as pd

QV_FLOOR = 40            # Merqury QV40 ~ 1 error/10 kb; the EBP/VGP minimum (Rhie 2021)
COMPLETE_MIN = 95.0      # BUSCO/compleasm Complete% field convention for a well-covered lineage
FRAGMENTED_MAX = 5.0     # high Fragmented signals contiguity/base-quality problems behind a good Complete%
DUP_HAPLOTIG_FLAG = 5.0  # Duplicated% >5-8% with no known WGD suggests uncollapsed haplotigs (purge_dups)
KMER_COMPLETE_MIN = 95.0 # Merqury k-mer completeness; lower = sequence BUSCO genes cannot see is absent


def parse_quast(report_tsv):
    df = pd.read_csv(report_tsv, sep='\t', index_col=0).T
    row = df.iloc[0]
    return {k: row[k] for k in row.index if k in ('N50', 'NG50', 'auN', 'L50', '# contigs', 'Total length')}


def parse_busco_or_compleasm(summary_file):
    text = open(summary_file).read()
    m = re.search(r'C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%', text)
    if m:
        return dict(zip(('complete', 'single', 'duplicated', 'fragmented', 'missing'), map(float, m.groups())))
    sd = {cat.lower(): float(pct) for cat, pct in re.findall(r'\b([SDFIM]):([\d.]+)%', text)}
    return {'complete': sd.get('s', 0) + sd.get('d', 0), 'duplicated': sd.get('d', 0), 'fragmented': sd.get('f', 0) + sd.get('i', 0), 'missing': sd.get('m', 0)}


def parse_merqury(qv_file, completeness_stats=None):
    overall_qv = float(open(qv_file).read().split('\n')[0].split('\t')[3])
    kmer_completeness = None
    if completeness_stats:
        kmer_completeness = float(open(completeness_stats).read().split('\n')[0].split('\t')[4])
    return {'qv': overall_qv, 'kmer_completeness': kmer_completeness}


def verdict(quast, comp, merq):
    flags = []
    if merq['qv'] < QV_FLOOR:
        flags.append(f"correctness: QV {merq['qv']:.1f} below EBP floor {QV_FLOOR} (1 error/10 kb)")
    if comp['complete'] < COMPLETE_MIN:
        flags.append(f"completeness: Complete {comp['complete']:.1f}% < {COMPLETE_MIN}% -- on a good HiFi/T2T genome try compleasm before concluding incompleteness")
    if comp['fragmented'] > FRAGMENTED_MAX:
        flags.append(f"completeness: Fragmented {comp['fragmented']:.1f}% high -- contiguity/base-quality problem behind the Complete%")
    if comp['duplicated'] > DUP_HAPLOTIG_FLAG:
        flags.append(f"completeness: Duplicated {comp['duplicated']:.1f}% -- uncollapsed haplotigs until proven WGD; check size vs estimate + spectra-cn -> purge_dups")
    if merq['kmer_completeness'] is not None and merq['kmer_completeness'] < KMER_COMPLETE_MIN:
        flags.append(f"completeness: k-mer completeness {merq['kmer_completeness']:.1f}% -- whole-genome sequence absent that BUSCO cannot see")
    return flags


if __name__ == '__main__':
    import sys
    quast = parse_quast(sys.argv[1])
    comp = parse_busco_or_compleasm(sys.argv[2])
    merq = parse_merqury(sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    print('Axis A contiguity :', quast)
    print('Axis B completeness:', comp)
    print('Axis C correctness :', merq)
    issues = verdict(quast, comp, merq)
    print('\nFlags:' if issues else '\nNo threshold flags raised across the three axes.')
    for f in issues:
        print(' -', f)
