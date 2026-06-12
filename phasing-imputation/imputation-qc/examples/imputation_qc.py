#!/usr/bin/env python3
'''MAF-stratified imputation quality and masked dosage-r2 accuracy.'''
# Reference: cyvcf2 0.31+, numpy 1.26+, pandas 2.2+ | Verify API if version differs

import sys
import numpy as np
import pandas as pd
from cyvcf2 import VCF

QUAL_KEYS = ('DR2', 'R2', 'INFO')   # Beagle / Minimac / GLIMPSE-IMPUTE; the field name tells the engine
MAF_BINS = [0, 0.001, 0.01, 0.05, 0.5]   # rare-to-common; the rare bins are where a flat INFO cutoff silently bites

def detect_qual_key(vcf):
    header = vcf.raw_header
    for key in QUAL_KEYS:
        if 'ID=%s,' % key in header:
            return key
    return None

def quality_by_maf(vcf_path, cutoff=0.3):
    '''Report imputation quality per MAF bin so the hidden-MAF-filter effect of a flat cutoff is visible.'''
    vcf = VCF(vcf_path)
    qual_key = detect_qual_key(vcf)
    if qual_key is None:
        sys.exit('no imputation quality field (DR2/R2/INFO) in %s' % vcf_path)
    rows = []
    for v in vcf:
        q = v.INFO.get(qual_key)
        af = v.INFO.get('AF')
        if q is None or af is None:
            continue
        af = af[0] if isinstance(af, tuple) else af
        rows.append((min(af, 1 - af), float(q)))
    df = pd.DataFrame(rows, columns=['maf', 'qual'])
    df['maf_bin'] = pd.cut(df['maf'], bins=MAF_BINS)
    summary = df.groupby('maf_bin', observed=True).agg(n=('qual', 'size'), mean_qual=('qual', 'mean'), frac_pass=('qual', lambda q: (q >= cutoff).mean()))
    return qual_key, summary

def dosage_r2_by_maf(imputed_path, truth_path):
    '''Gold-standard accuracy: squared correlation of imputed dosage to masked-then-revealed true genotype, binned by MAF.'''
    truth = {}
    for v in VCF(truth_path, gts012=True):
        truth[(v.CHROM, v.POS)] = np.array(v.gt_types, dtype=float)   # gts012: 0,1,2 = HOM_REF,HET,HOM_ALT; 3 = missing
    recs = []
    for v in VCF(imputed_path):
        key = (v.CHROM, v.POS)
        if key not in truth:
            continue
        ds = v.format('DS')
        if ds is None:
            continue
        g = truth[key]
        keep = g != 3                       # drop missing truth genotypes, else they score as hom-alt
        ds = np.asarray(ds, dtype=float).ravel()[keep]
        g = g[keep]
        if len(g) > 1 and np.std(ds) > 0 and np.std(g) > 0:
            recs.append((min(v.aaf, 1 - v.aaf), np.corrcoef(ds, g)[0, 1] ** 2))
    df = pd.DataFrame(recs, columns=['maf', 'r2'])
    df['maf_bin'] = pd.cut(df['maf'], bins=MAF_BINS)
    return df.groupby('maf_bin', observed=True).agg(n=('r2', 'size'), mean_r2=('r2', 'mean'))

if __name__ == '__main__':
    vcf_path = sys.argv[1] if len(sys.argv) > 1 else 'imputed.vcf.gz'
    key, summary = quality_by_maf(vcf_path)
    print('quality field: %s' % key)
    print(summary)
