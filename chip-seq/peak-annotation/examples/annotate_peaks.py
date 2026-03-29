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
    '''Annotate peaks with gene assignment and genomic feature using
    host-gene convention: exon/intron peaks are assigned to the gene
    whose body contains the peak. Intergenic and promoter peaks use
    the nearest TSS. Feature priority: promoter > exon > intron > intergenic.
    Signed distance: negative = upstream, positive = downstream.
    '''
    gtf = parse_gtf(gtf_path)
    genes = gtf[gtf['feature'] == 'gene'].copy()
    genes['tss'] = genes.apply(lambda r: r['start'] if r['strand'] == '+' else r['end'], axis=1)
    exons = gtf[gtf['feature'] == 'exon']

    peaks = pd.read_csv(peaks_path, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'peak_id', 'score'])
    peaks['center'] = (peaks['start'] + peaks['end']) // 2

    results = []
    for _, peak in peaks.iterrows():
        chrom_genes = genes[genes['chrom'] == peak['chr']]
        chrom_exons = exons[exons['chrom'] == peak['chr']]
        abs_dists = (chrom_genes['tss'] - peak['center']).abs()
        nearest_tss_gene = chrom_genes.loc[abs_dists.idxmin()]

        if abs_dists.min() <= promoter_window:
            feature, assigned = 'promoter', nearest_tss_gene
        else:
            exon_hits = chrom_exons[(chrom_exons['start'] <= peak['center']) & (peak['center'] < chrom_exons['end'])]
            gene_hits = chrom_genes[(chrom_genes['start'] <= peak['center']) & (peak['center'] < chrom_genes['end'])]
            if len(exon_hits) > 0:
                host_gene_name = exon_hits.iloc[0].get('gene_name', '')
                host = chrom_genes[chrom_genes['gene_name'] == host_gene_name]
                feature, assigned = 'exon', host.iloc[0] if len(host) > 0 else nearest_tss_gene
            elif len(gene_hits) > 0:
                closest_host = gene_hits.loc[(gene_hits['tss'] - peak['center']).abs().idxmin()]
                feature, assigned = 'intron', closest_host
            else:
                feature, assigned = 'intergenic', nearest_tss_gene

        raw_dist = peak['center'] - assigned['tss']
        signed_dist = -raw_dist if assigned['strand'] == '-' else raw_dist

        results.append({'peak_id': peak['peak_id'], 'chr': peak['chr'],
                        'start': peak['start'], 'end': peak['end'],
                        'nearest_gene': assigned['gene_name'],
                        'distance_to_tss': int(signed_dist), 'feature': feature})

    return pd.DataFrame(results)


annotations = annotate_peaks('peaks.bed', 'genes.gtf.gz', promoter_window=2000)
annotations.to_csv('annotations.tsv', sep='\t', index=False)
