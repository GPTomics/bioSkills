#!/usr/bin/env bash
# Reference: dwgsim 0.1.13+, samtools 1.19+, seqtk 1.4+, bcftools 1.19+ | Verify CLI if version differs
#
# Build a SHAREABLE, REPRODUCIBLE test fixture for a variant-calling workflow, plus a
# normalized-VCF comparison helper. Two reproducibility-by-design points are baked in:
#   1) the simulator seed is FIXED, so the "reference" data is itself reproducible;
#   2) the expected-output comparison strips the VCF header, so timestamps and
#      program-version lines never cause a false test failure.
#
# Usage: ./make_test_fixture.sh <reference.fa> <region e.g. chr20:1000000-1010000>
set -euo pipefail

REF="${1:?usage: make_test_fixture.sh <reference.fa> <region>}"
REGION="${2:?provide a small region, e.g. chr20:1000000-1010000}"
SEED=42                    # <- part of the fixture; changing it changes the data
OUT="tests/data"
mkdir -p "$OUT"

# 1) Carve a tiny reference slice so the fixture is CI-sized and contains no sensitive sample.
samtools faidx "$REF" "$REGION" > "$OUT/region.fa"
samtools faidx "$OUT/region.fa"

# 2) Simulate paired-end reads with a FIXED seed (-z). dwgsim injects known variants and
#    writes the truth VCF, which becomes the expected-calls oracle for the smoke test.
if command -v dwgsim >/dev/null 2>&1; then
  dwgsim -z "$SEED" -N 2000 -1 100 -2 100 "$OUT/region.fa" "$OUT/tiny"
  echo "  simulated: $OUT/tiny.bwa.read1.fastq.gz + read2 ; truth: $OUT/tiny.mutations.vcf"
else
  echo "  dwgsim not installed (conda install -c bioconda dwgsim); skipping simulation"
fi

# Alternative when shareable REAL data exists: downsample deterministically instead of simulating.
#   samtools view -s ${SEED}.01 -b big.public.bam > tests/data/tiny.bam   # 1%, seed-locked

cat > "$OUT/compare_vcf.sh" <<'EOF'
#!/usr/bin/env bash
# Normalized VCF comparison: ignore the header (timestamps, ##bcftools_command, versions)
# and compare only sorted records. Use this instead of `md5sum whole.vcf`, which false-fails.
set -euo pipefail
norm() { bcftools view --no-version "$1" | grep -v '^##' | sort; }
if diff <(norm "$1") <(norm "$2") > /dev/null; then
  echo "PASS: records match (header ignored)"
else
  echo "FAIL: record-level difference"; diff <(norm "$1") <(norm "$2"); exit 1
fi
EOF
chmod +x "$OUT/compare_vcf.sh"

echo
echo "Fixture ready under $OUT/ (seed=$SEED). Commit it with the tests."
echo "Smoke test: run the pipeline on the tiny reads, then:"
echo "  $OUT/compare_vcf.sh results/calls.vcf $OUT/tiny.mutations.vcf"
