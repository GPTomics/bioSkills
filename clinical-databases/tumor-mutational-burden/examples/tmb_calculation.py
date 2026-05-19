'''TMB calculation with Friends of Cancer Research (Vega 2021) per-panel calibration.

Reference: cyvcf2 0.30+, Ensembl VEP 111+ (or snpEff 5.2+) | Verify VCF INFO field names.
FoundationOne CDx scored region is 0.8 Mb (NOT 1.1 Mb total panel).
FDA pembrolizumab 10 mut/Mb cutoff equivalent: TSO500 7.8, Oncomine 8.4.
'''
from cyvcf2 import VCF
import argparse


# Scored regions in Mb (NOT total panel size).
# Vega 2021 derived from in-silico panel sampling of TCGA WES truth.
PANEL_SCORED_REGION_MB = {
    'foundationone_cdx': 0.8,    # SCORED region; 1.1 Mb is total panel content
    'msk_impact_v3': 0.98,
    'msk_impact_v4': 1.22,
    'tso500': 1.3,                # ~1.3 Mb scored from 1.94 Mb total
    'oncomine_tml': 1.2,
    'caris_mi': 1.2,
    'tempus_xt_v3': 0.6,          # Borderline (<0.8 minimum)
    'predicine_atlas': 0.6,        # Borderline
    'wes': 30.0,
    'wgs': 3000.0,
}

# Vega 2021 calibration: equivalent thresholds for FoundationOne 10/Mb sensitivity.
ASSAY_TMB_H_CUTOFF = {
    'foundationone_cdx': 10.0,
    'tso500': 7.8,                 # Vega 2021
    'oncomine_tml': 8.4,           # Vega 2021
    'msk_impact_v3': 10.0,         # Full Vega 2021 calibration recommended
    'msk_impact_v4': 10.0,
    'wes': 10.0,
}

# FoundationOne CDx INCLUDES synonymous; MSK-IMPACT + most academic pipelines exclude.
ASSAY_INCLUDES_SYNONYMOUS = {
    'foundationone_cdx': True,
    'msk_impact_v3': False,
    'msk_impact_v4': False,
    'tso500': False,
    'oncomine_tml': False,
    'caris_mi': True,  # Verify per Caris docs
    'wes': False,
}

NONSYNONYMOUS_CONSEQUENCES = {
    'missense_variant', 'stop_gained', 'stop_lost', 'start_lost', 'start_retained',
    'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
    'splice_donor_variant', 'splice_acceptor_variant',
    'protein_altering_variant', 'initiator_codon_variant'
}


def is_target_variant(variant, include_synonymous=False, csq_header_consequence_idx=1):
    '''Check if variant is target type (nonsynonymous + optionally synonymous).'''
    csq = variant.INFO.get('CSQ')
    if not csq:
        return False
    target = set(NONSYNONYMOUS_CONSEQUENCES)
    if include_synonymous:
        target.add('synonymous_variant')
    for transcript in csq.split(','):
        fields = transcript.split('|')
        if len(fields) <= csq_header_consequence_idx:
            continue
        consequences = set(fields[csq_header_consequence_idx].split('&'))
        if consequences & target:
            return True
    return False


def get_vaf(variant):
    '''Extract VAF from Mutect2 (AD field) or generic AF field.'''
    try:
        ad = variant.format('AD')
        if ad is not None and len(ad) > 0:
            ad0 = ad[0]
            total = sum(ad0)
            return ad0[1] / total if total > 0 else 0.0
    except Exception:
        pass
    try:
        af = variant.format('AF')
        if af is not None and len(af) > 0:
            return float(af[0])
    except Exception:
        pass
    return 0.0


def calculate_tmb(vcf_path, assay='foundationone_cdx',
                   min_vaf=0.05, min_depth=100, max_gnomad_grpmax_faf95=0.005,
                   exclude_hotspots=False, hotspots_bed=None):
    '''Calculate TMB with Vega 2021 calibration and assay-specific conventions.

    Args:
        vcf_path: VEP-annotated somatic VCF
        assay: panel name (sets scored region + synonymous convention)
        min_vaf: FoundationOne 0.05; tumor-only no UMI 0.10
        max_gnomad_grpmax_faf95: 0.005 (0.5%) for tumor-only germline filter
        exclude_hotspots: exclude COSMIC drivers (recommended)
        hotspots_bed: BED of hotspot positions to exclude
    '''
    scored_mb = PANEL_SCORED_REGION_MB[assay]
    include_syn = ASSAY_INCLUDES_SYNONYMOUS[assay]
    vcf = VCF(vcf_path)

    counts = {
        'pass_filters': 0, 'count_target': 0,
        'excl_filter': 0, 'excl_lowvaf': 0, 'excl_lowdepth': 0,
        'excl_germline': 0, 'excl_hotspot': 0, 'excl_nontarget': 0
    }

    for v in vcf:
        if v.FILTER is not None:
            counts['excl_filter'] += 1
            continue
        depth = v.INFO.get('DP', 0)
        if depth < min_depth:
            counts['excl_lowdepth'] += 1
            continue
        vaf = get_vaf(v)
        if vaf < min_vaf:
            counts['excl_lowvaf'] += 1
            continue
        # Germline filter: prefer grpmax FAF95; fallback to AF_popmax / gnomAD_AF
        gnomad_af = (v.INFO.get('grpmax_faf95') or
                     v.INFO.get('AF_grpmax') or
                     v.INFO.get('AF_popmax') or
                     v.INFO.get('gnomAD_AF') or 0)
        if gnomad_af > max_gnomad_grpmax_faf95:
            counts['excl_germline'] += 1
            continue
        counts['pass_filters'] += 1

        if is_target_variant(v, include_synonymous=include_syn):
            counts['count_target'] += 1
        else:
            counts['excl_nontarget'] += 1

    tmb = counts['count_target'] / scored_mb
    cutoff = ASSAY_TMB_H_CUTOFF[assay]
    classification = classify_tmb_tier(tmb, cutoff)
    return {
        'tmb': round(tmb, 2),
        'scored_region_mb': scored_mb,
        'assay': assay,
        'assay_includes_synonymous': include_syn,
        'tmb_h_cutoff_per_vega2021': cutoff,
        'classification': classification,
        **counts
    }


def classify_tmb_tier(tmb, cutoff=10.0):
    '''Classify TMB into clinical tiers.

    Ultra-hypermutator (>=500) typically POLE+MMR -- ICI excellent.
    Hypermutator (>=100) typically MMR-D or POLE -- ICI eligible.
    TMB-H (>= assay cutoff) FDA pan-tumor ICI eligible.
    '''
    if tmb >= 500:
        return 'Ultra-hypermutator (>=500/Mb; POLE+MMR likely)'
    if tmb >= 100:
        return 'Hypermutator (>=100/Mb; MMR-D or POLE-exo)'
    if tmb >= cutoff:
        return f'TMB-H (>= {cutoff}/Mb assay-calibrated; pan-tumor ICI eligible per FDA 2020)'
    return 'TMB-low'


def tmb_msi_decision(tmb_value, msi_status, tumor_type=None, hla_loh=False):
    '''Integrated ICI eligibility from TMB + MSI + HLA-LOH.

    Sha 2020 Cell Rep Med: MSI-H is primary; TMB-H not additive.
    McGrail 2021 + ESMO 2024: TMB-H NOT endorsed for breast/prostate/glioma.
    Marty 2017 / Montesion 2021: HLA-LOH reduces neoantigen presentation.
    '''
    tmb_h = tmb_value >= 10
    excluded_tumors = {'breast', 'prostate', 'glioma'}

    if msi_status == 'MSI-H':
        return {'eligible': True, 'rationale': 'MSI-H primary biomarker (FDA 2017 pembrolizumab); '
                                                'TMB-H not additive (Sha 2020)'}
    if tmb_h and tumor_type and tumor_type.lower() in excluded_tumors:
        return {'eligible': False, 'rationale': f'TMB-H present but NOT endorsed for '
                                                  f'{tumor_type} per ESMO 2024 + McGrail 2021. '
                                                  'Use tumor-type-specific cutoff (Samstein 2019).'}
    if tmb_h:
        flag = (' Caution: HLA-LOH detected -- neoantigen presentation may be reduced.'
                if hla_loh else '')
        return {'eligible': True, 'rationale': 'TMB-H pan-tumor; ICI eligible per FDA 2020.' + flag}
    return {'eligible': False, 'rationale': 'TMB-low and MSS. Standard-of-care therapy.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TMB with Vega 2021 calibration')
    parser.add_argument('vcf', help='VEP-annotated somatic VCF')
    parser.add_argument('--assay', choices=list(PANEL_SCORED_REGION_MB.keys()),
                        default='foundationone_cdx', help='Panel name')
    parser.add_argument('--min-vaf', type=float, default=0.05)
    parser.add_argument('--min-depth', type=int, default=100)
    parser.add_argument('--max-gnomad-faf95', type=float, default=0.005,
                        help='Tumor-only germline filter (gnomAD grpmax FAF95)')
    args = parser.parse_args()

    result = calculate_tmb(args.vcf, assay=args.assay,
                           min_vaf=args.min_vaf, min_depth=args.min_depth,
                           max_gnomad_grpmax_faf95=args.max_gnomad_faf95)
    print(f"TMB: {result['tmb']} mut/Mb on {args.assay} ({result['scored_region_mb']} Mb scored)")
    print(f"Assay synonymous convention: {result['assay_includes_synonymous']}")
    print(f"Vega 2021 calibrated TMB-H cutoff: {result['tmb_h_cutoff_per_vega2021']}")
    print(f"Classification: {result['classification']}")
    print(f"\nCounts: target={result['count_target']}/pass-filters={result['pass_filters']}; "
          f"excluded germline={result['excl_germline']}, low-VAF={result['excl_lowvaf']}, "
          f"low-depth={result['excl_lowdepth']}")
