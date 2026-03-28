# Reference: pandas 2.2+ | Verify API if version differs
import gzip
import pandas as pd


def parse_gtf(gtf_path):
    '''Parse GTF into DataFrame with 0-based half-open coordinates.'''
    records = []
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attrs = {}
            for item in fields[8].strip().rstrip(';').split(';'):
                item = item.strip()
                if ' ' in item:
                    key, val = item.split(' ', 1)
                    attrs[key] = val.strip('"')
            records.append({'chrom': fields[0], 'feature': fields[2],
                            'start': int(fields[3]) - 1, 'end': int(fields[4]),
                            'strand': fields[6], **attrs})
    return pd.DataFrame(records)


def annotate_peaks(peaks_path, gtf_path, promoter_window=2000):
    '''Annotate peaks with nearest gene, signed distance to TSS, and feature.

    Promoter window defines max distance from TSS to classify as promoter.
    Feature priority: promoter > exon > intron > intergenic.
    Signed distance: negative = upstream, positive = downstream.
    '''
    gtf = parse_gtf(gtf_path)
    genes = gtf[gtf['feature'] == 'gene'].copy()
    genes['tss'] = genes.apply(lambda r: r['start'] if r['strand'] == '+' else r['end'], axis=1)
    exons = gtf[gtf['feature'] == 'exon']

    peaks = pd.read_csv(peaks_path, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'name', 'score'])
    peaks['center'] = (peaks['start'] + peaks['end']) // 2

    results = []
    for _, peak in peaks.iterrows():
        chrom_genes = genes[genes['chrom'] == peak['chr']]
        abs_dists = (chrom_genes['tss'] - peak['center']).abs()
        nearest = chrom_genes.loc[abs_dists.idxmin()]

        raw_dist = peak['center'] - nearest['tss']
        signed_dist = -raw_dist if nearest['strand'] == '-' else raw_dist

        # Feature classification: promoter > exon > intron > intergenic
        if abs(raw_dist) <= promoter_window:
            feature = 'promoter'
        else:
            chrom_exons = exons[exons['chrom'] == peak['chr']]
            in_exon = ((chrom_exons['start'] <= peak['center']) & (peak['center'] < chrom_exons['end'])).any()
            in_gene = ((chrom_genes['start'] <= peak['center']) & (peak['center'] < chrom_genes['end'])).any()
            feature = 'exon' if in_exon else ('intron' if in_gene else 'intergenic')

        results.append({'chr': peak['chr'], 'start': peak['start'], 'end': peak['end'],
                        'nearest_gene': nearest['gene_name'], 'distance_to_tss': int(signed_dist),
                        'feature': feature})

    return pd.DataFrame(results)


annotations = annotate_peaks('peaks.bed', 'genes.gtf.gz', promoter_window=2000)
annotations.to_csv('annotations.tsv', sep='\t', index=False)
