#!/bin/bash
# Reference: inStrain 1.8+, dRep 3.4+, Bowtie2 2.5+, samtools 1.19+ | Verify API if version differs
# Detect a shared strain between two metagenomes with inStrain popANI. Map to YOUR OWN dRep'd MAGs,
# not database genomes - a distant reference inflates SNVs and corrupts popANI.
set -euo pipefail

MAG_DIR=${1:-mags}          # MAGs assembled and binned from THIS dataset
R1A=${2:-sampleA_R1.fq.gz}; R2A=${3:-sampleA_R2.fq.gz}
R1B=${4:-sampleB_R1.fq.gz}; R2B=${5:-sampleB_R2.fq.gz}
OUT=${6:-strain_out}
mkdir -p "$OUT"

# 1. Dereplicate MAGs into representative genomes (97-99% ANI) and build the reference set.
dRep dereplicate "${OUT}/drep" -g "${MAG_DIR}"/*.fasta
cat "${OUT}"/drep/dereplicated_genomes/*.fasta > "${OUT}/reps.fasta"
# Build a scaffold-to-bin (.stb) so inStrain reports per-genome popANI.
parse_stb.py --reverse -f "${OUT}"/drep/dereplicated_genomes/*.fasta -o "${OUT}/reps.stb"

# 2. Map each sample to the concatenated reps and profile.
bowtie2-build "${OUT}/reps.fasta" "${OUT}/reps" >/dev/null
for s in A B; do
    r1=$([ "$s" = A ] && echo "$R1A" || echo "$R1B")
    r2=$([ "$s" = A ] && echo "$R2A" || echo "$R2B")
    bowtie2 -x "${OUT}/reps" -1 "$r1" -2 "$r2" -p 8 | samtools sort -o "${OUT}/sample${s}.bam"
    inStrain profile "${OUT}/sample${s}.bam" "${OUT}/reps.fasta" \
        -o "${OUT}/sample${s}.IS" -s "${OUT}/reps.stb" -p 8
done

# 3. Compare on popANI over the co-covered genome fraction. Same strain: popANI >= 99.999%
# AND percent_compared >= 50% (genome-level genomeWide_compare.tsv breadth column; the per-scaffold
# comparisonsTable.tsv instead names it percent_genome_compared). The threshold IS the strain definition.
inStrain compare -i "${OUT}/sampleA.IS" "${OUT}/sampleB.IS" -o "${OUT}/compare" -s "${OUT}/reps.stb" -p 8
echo "Shared-strain calls in ${OUT}/compare/output/genomeWide_compare.tsv - apply popANI >= 0.99999 and percent_compared >= 50%."
