'''Calculate Tumor Mutational Burden from somatic VCF'''

from cyvcf2 import VCF
import argparse

# --- Panel sizes in megabases ---
# The capture region size determines TMB normalization
# Use the exact size from your panel's documentation
PANEL_SIZES_MB = {
    'foundation_cdx': 0.8,      # FoundationOne CDx
    'msk_impact': 1.14,         # MSK-IMPACT
    'tso500': 1.94,             # TruSight Oncology 500
    'oncomine': 1.5,            # Oncomine Comprehensive
    'wes': 30.0,                # Whole exome (approximate coding)
    'wgs': 3000.0,              # Whole genome
}

# --- Clinical thresholds ---
# FDA approval for pembrolizumab: TMB >= 10 mut/Mb
# Some studies use higher thresholds (16, 20)
TMB_THRESHOLDS = {
    'fda': 10,
    'conservative': 16,
    'strict': 20,
}


def is_coding_nonsynonymous(variant):
    '''Check if variant is coding nonsynonymous based on VEP annotation

    Adjust for your annotation format:
    - VEP: CSQ field (used here)
    - SnpEff: ANN field
    - Funcotator: FUNCOTATION field
    '''
    csq = variant.INFO.get('CSQ')
    if not csq:
        return False

    # VEP CSQ format: Allele|Consequence|IMPACT|...
    # Count nonsynonymous coding consequences
    nonsynonymous_consequences = {
        'missense_variant',
        'stop_gained',
        'stop_lost',
        'start_lost',
        'frameshift_variant',
        'inframe_insertion',
        'inframe_deletion',
        'protein_altering_variant',
    }

    for transcript in csq.split(','):
        fields = transcript.split('|')
        if len(fields) < 2:
            continue
        consequences = set(fields[1].split('&'))
        if consequences & nonsynonymous_consequences:
            return True
    return False


def get_vaf(variant):
    '''Extract variant allele frequency from Mutect2 format

    Adjust for your variant caller:
    - Mutect2: AD field (ref,alt counts)
    - Strelka: TIR/TAR fields
    - VarDict: AF field
    '''
    try:
        ad = variant.format('AD')[0]  # First sample (tumor)
        if ad is not None and sum(ad) > 0:
            return ad[1] / sum(ad)
    except (KeyError, TypeError, IndexError):
        pass
    return 0.0


def calculate_tmb(vcf_path, panel_size_mb, min_vaf=0.05, min_depth=100):
    '''Calculate TMB with quality filtering

    Args:
        vcf_path: Path to annotated somatic VCF
        panel_size_mb: Panel capture size in megabases
        min_vaf: Minimum variant allele frequency (default 5%)
        min_depth: Minimum total depth (default 100)

    Filters applied:
    - VAF >= 5%: Reduce false positives from errors
    - Depth >= 100: Ensure reliable calls
    - gnomAD AF <= 1%: Exclude germline polymorphisms
    - Coding nonsynonymous: Standard TMB definition

    Returns:
        dict with TMB value and counts
    '''
    vcf = VCF(vcf_path)

    total_variants = 0
    passing_variants = 0
    excluded_lowvaf = 0
    excluded_lowdepth = 0
    excluded_germline = 0
    excluded_noncoding = 0

    for variant in vcf:
        total_variants += 1

        # Depth filter
        depth = variant.INFO.get('DP', 0)
        if depth < min_depth:
            excluded_lowdepth += 1
            continue

        # VAF filter
        vaf = get_vaf(variant)
        if vaf < min_vaf:
            excluded_lowvaf += 1
            continue

        # Germline filter (gnomAD)
        # Field name varies by annotation version
        gnomad_af = variant.INFO.get('gnomAD_AF') or variant.INFO.get('AF_popmax') or 0
        if gnomad_af > 0.01:
            excluded_germline += 1
            continue

        # Coding nonsynonymous filter
        if not is_coding_nonsynonymous(variant):
            excluded_noncoding += 1
            continue

        passing_variants += 1

    tmb = passing_variants / panel_size_mb

    return {
        'tmb': round(tmb, 2),
        'passing_variants': passing_variants,
        'total_variants': total_variants,
        'excluded_lowvaf': excluded_lowvaf,
        'excluded_lowdepth': excluded_lowdepth,
        'excluded_germline': excluded_germline,
        'excluded_noncoding': excluded_noncoding,
        'panel_size_mb': panel_size_mb,
    }


def classify_tmb(tmb_value, threshold='fda'):
    '''Classify TMB as high or low

    FDA threshold (10 mut/Mb) used for pembrolizumab approval
    '''
    cutoff = TMB_THRESHOLDS.get(threshold, 10)
    return 'TMB-High' if tmb_value >= cutoff else 'TMB-Low'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate Tumor Mutational Burden')
    parser.add_argument('vcf', help='Input somatic VCF (annotated with VEP)')
    parser.add_argument('--panel', choices=list(PANEL_SIZES_MB.keys()),
                        default='wes', help='Panel name for size')
    parser.add_argument('--panel-size', type=float,
                        help='Custom panel size in Mb (overrides --panel)')
    parser.add_argument('--min-vaf', type=float, default=0.05,
                        help='Minimum VAF (default: 0.05)')
    parser.add_argument('--min-depth', type=int, default=100,
                        help='Minimum depth (default: 100)')
    args = parser.parse_args()

    panel_size = args.panel_size or PANEL_SIZES_MB[args.panel]

    print(f'Calculating TMB from: {args.vcf}')
    print(f'Panel size: {panel_size} Mb')
    print(f'Filters: VAF >= {args.min_vaf}, Depth >= {args.min_depth}')
    print()

    results = calculate_tmb(args.vcf, panel_size, args.min_vaf, args.min_depth)

    print(f"TMB: {results['tmb']} mutations/Mb")
    print(f"Classification: {classify_tmb(results['tmb'])}")
    print()
    print('Variant counts:')
    print(f"  Passing (nonsynonymous coding): {results['passing_variants']}")
    print(f"  Excluded - low VAF: {results['excluded_lowvaf']}")
    print(f"  Excluded - low depth: {results['excluded_lowdepth']}")
    print(f"  Excluded - germline: {results['excluded_germline']}")
    print(f"  Excluded - non-coding: {results['excluded_noncoding']}")
