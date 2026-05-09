#!/bin/bash
# Reference: WASP 0.3.4+, GATK 4.4+, samtools 1.19+, bcftools 1.19+, bowtie2 2.5+ | Verify API if version differs
# WASP correction + GATK ASEReadCounter + within-peak ASE aggregation.
# Mandatory: phase genotypes BEFORE WASP. Reference-allele mapping bias must be corrected.

set -euo pipefail

BAM=${1:-sample.dedup.bam}                       # WASP-input ATAC BAM
VCF=${2:-sample.phased.het.vcf.gz}               # Phased heterozygous SNPs only
GENOME_FA=${3:-hg38.fa}
BT2_INDEX=${4:-hg38_bt2}                         # bowtie2 index (must match GENOME_FA)
PEAKS=${5:-consensus_peaks.bed}
WASP_DIR=${6:-/path/to/WASP}
OUTDIR=${7:-ase_out}

mkdir -p $OUTDIR/{wasp,ase,peak_aggregation}

SAMPLE=$(basename $BAM .bam)

# 0. Pre-filter VCF to biallelic heterozygous sites (GATK ASE requirement)
bcftools view -m2 -M2 -v snps -i 'GT="het"' $VCF > $OUTDIR/wasp/${SAMPLE}.het.vcf
bgzip $OUTDIR/wasp/${SAMPLE}.het.vcf
tabix -p vcf $OUTDIR/wasp/${SAMPLE}.het.vcf.gz

# 1. WASP: find reads at SNP sites
python $WASP_DIR/mapping/find_intersecting_snps.py \
    --is_paired_end \
    --is_sorted \
    --output_dir $OUTDIR/wasp/ \
    --snp_dir $OUTDIR/wasp/snp_h5/ \
    --snp_tab $OUTDIR/wasp/snp_h5/snp_tab.h5 \
    --snp_index $OUTDIR/wasp/snp_h5/snp_index.h5 \
    --haplotype $OUTDIR/wasp/snp_h5/haps.h5 \
    --samples $SAMPLE \
    $BAM

# 2. Re-align flipped-allele reads
bowtie2 -x $BT2_INDEX \
    -1 $OUTDIR/wasp/${SAMPLE}.remap.fq1.gz \
    -2 $OUTDIR/wasp/${SAMPLE}.remap.fq2.gz \
    -p 8 -S $OUTDIR/wasp/${SAMPLE}.remap.sam

samtools view -bS $OUTDIR/wasp/${SAMPLE}.remap.sam | \
    samtools sort -o $OUTDIR/wasp/${SAMPLE}.remap.bam
samtools index $OUTDIR/wasp/${SAMPLE}.remap.bam

# 3. Keep only reads that map identically to both haplotypes
python $WASP_DIR/mapping/filter_remapped_reads.py \
    $OUTDIR/wasp/${SAMPLE}.to.remap.bam \
    $OUTDIR/wasp/${SAMPLE}.remap.bam \
    $OUTDIR/wasp/${SAMPLE}.kept.bam

# 4. Merge kept reads with non-overlapping reads
samtools merge -f $OUTDIR/wasp/${SAMPLE}.wasp.bam \
    $OUTDIR/wasp/${SAMPLE}.kept.bam \
    $OUTDIR/wasp/${SAMPLE}.keep.bam
samtools index $OUTDIR/wasp/${SAMPLE}.wasp.bam

# 5. GATK ASEReadCounter on WASP-filtered BAM
gatk ASEReadCounter \
    -I $OUTDIR/wasp/${SAMPLE}.wasp.bam \
    -V $OUTDIR/wasp/${SAMPLE}.het.vcf.gz \
    -R $GENOME_FA \
    -O $OUTDIR/ase/${SAMPLE}.ase_counts.tsv \
    --output-format TABLE \
    --min-mapping-quality 30 \
    --min-base-quality 20

# 6. Within-peak ASE aggregation (Python)
python3 - <<PYEOF
import pandas as pd, numpy as np
from scipy import stats
import pybedtools

ase = pd.read_csv("$OUTDIR/ase/${SAMPLE}.ase_counts.tsv", sep="\t")
ase = ase.rename(columns={"refCount": "REF", "altCount": "ALT", "totalCount": "total"})
ase = ase[ase["total"] >= 30].copy()
ase["ref_frac"] = ase["REF"] / ase["total"]

# Map SNPs to peaks
snps_bed = ase[["contig", "position", "position", "variantID"]].copy()
snps_bed.columns = ["chrom", "start", "end", "name"]
snps_bed["start"] -= 1
peaks_bt = pybedtools.BedTool("$PEAKS")
snps_bt = pybedtools.BedTool.from_dataframe(snps_bed)
mapped = snps_bt.intersect(peaks_bt, wa=True, wb=True).to_dataframe()
mapped["peak"] = mapped[["score", "strand", "thickStart"]].astype(str).agg("_".join, axis=1)

ase = ase.merge(mapped[["name", "peak"]], left_on="variantID", right_on="name")

# Pooled binomial per peak (>= 2 SNPs)
def peak_test(g):
    if len(g) < 2: return pd.Series({"ref_frac": np.nan, "p_value": np.nan, "snp_count": len(g)})
    r, t = g["REF"].sum(), g["total"].sum()
    return pd.Series({"ref_frac": r / t, "p_value": stats.binomtest(r, t, p=0.5).pvalue, "snp_count": len(g)})

peak_ase = ase.groupby("peak").apply(peak_test)
peak_ase["adj_p"] = stats.false_discovery_control(peak_ase["p_value"].fillna(1.0))
sig = peak_ase[(peak_ase["adj_p"] < 0.05) & (abs(peak_ase["ref_frac"] - 0.5) >= 0.2)]
sig.to_csv("$OUTDIR/peak_aggregation/${SAMPLE}.peak_ase_sig.tsv", sep="\t")
print(f"Significant ASE peaks (FDR<0.05, effect>=0.2): {len(sig)}")
PYEOF

echo "Done. Outputs in $OUTDIR/{wasp,ase,peak_aggregation}/"
