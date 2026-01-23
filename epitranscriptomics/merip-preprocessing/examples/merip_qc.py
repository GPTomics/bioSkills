'''MeRIP-seq QC: check IP enrichment and replicate correlation'''
import pysam
import numpy as np
from pathlib import Path

def calculate_enrichment(ip_bam, input_bam, regions_bed):
    '''Calculate IP/Input enrichment at known m6A regions'''
    ip = pysam.AlignmentFile(ip_bam, 'rb')
    inp = pysam.AlignmentFile(input_bam, 'rb')

    enrichments = []
    with open(regions_bed) as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')[:3]
            start, end = int(start), int(end)

            ip_count = ip.count(chrom, start, end)
            input_count = inp.count(chrom, start, end)

            # Avoid division by zero
            if input_count > 0:
                enrichments.append(ip_count / input_count)

    ip.close()
    inp.close()

    return np.median(enrichments) if enrichments else 0

def check_replicate_correlation(bam_files, bin_size=10000):
    '''Compute pairwise correlation between BAM coverage profiles'''
    from scipy import stats

    coverages = []
    for bam_file in bam_files:
        bam = pysam.AlignmentFile(bam_file, 'rb')
        coverage = []
        for chrom in bam.references[:22]:  # Autosomes only
            length = bam.get_reference_length(chrom)
            for start in range(0, length, bin_size):
                end = min(start + bin_size, length)
                coverage.append(bam.count(chrom, start, end))
        coverages.append(coverage)
        bam.close()

    # Pairwise Spearman correlation
    # Good replicates: r > 0.9 for same condition
    correlations = {}
    for i in range(len(bam_files)):
        for j in range(i + 1, len(bam_files)):
            r, p = stats.spearmanr(coverages[i], coverages[j])
            correlations[(Path(bam_files[i]).stem, Path(bam_files[j]).stem)] = r

    return correlations

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 3:
        print('Usage: merip_qc.py <ip.bam> <input.bam> [known_m6a.bed]')
        sys.exit(1)

    ip_bam = sys.argv[1]
    input_bam = sys.argv[2]

    if len(sys.argv) > 3:
        regions = sys.argv[3]
        enrich = calculate_enrichment(ip_bam, input_bam, regions)
        # Good MeRIP: enrichment > 2 at known m6A sites
        print(f'Median IP/Input enrichment at known sites: {enrich:.2f}')

    corr = check_replicate_correlation([ip_bam, input_bam])
    for pair, r in corr.items():
        print(f'{pair[0]} vs {pair[1]}: r = {r:.3f}')
