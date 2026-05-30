#!/bin/bash
# Reference: bedtools 2.31+ | Verify API if version differs
# Validate a BED file and audit the two silent failures (chrom naming, CRLF).

BED_FILE="${1:-input.bed}"
PARTNER="${2:-}"   # optional second file to check chrom-naming compatibility against

echo "Validating: $BED_FILE"

if [[ ! -f "$BED_FILE" ]]; then
    echo "ERROR: file not found"
    exit 1
fi

# Field count must be identical on every line (positional columns).
echo "Field counts (expect one value):"
awk -F'\t' '{print NF}' "$BED_FILE" | sort | uniq -c

# Inverted (start >= end) or negative coordinates -- bedtools also errors on these.
awk -F'\t' '$2 < 0 || $2 >= $3' "$BED_FILE" | head -5 | grep -q . \
    && echo "WARNING: negative or inverted intervals present" \
    || echo "Coordinates OK (no negative/inverted)"

# CRLF audit: a trailing ^M corrupts the last field silently in tolerant parsers.
grep -q $'\r' "$BED_FILE" \
    && echo "WARNING: CRLF line endings -- fix with dos2unix or sed 's/\\r\$//'" \
    || echo "Line endings OK (no CRLF)"

# Chromosome names -- compare against a partner file before any cross-file operation.
echo "Chromosomes:"
cut -f1 "$BED_FILE" | sort -u
if [[ -n "$PARTNER" && -f "$PARTNER" ]]; then
    SHARED=$(comm -12 <(cut -f1 "$BED_FILE" | sort -u) <(cut -f1 "$PARTNER" | sort -u) | wc -l)
    [[ "$SHARED" -eq 0 ]] \
        && echo "WARNING: zero shared chromosomes with $PARTNER (chr1 vs 1 mismatch?) -- intersect will be empty" \
        || echo "Shares $SHARED chromosome name(s) with $PARTNER"
fi

# Sortedness: -sorted ops assume lexicographic order matching the partner/genome file.
if diff -q <(bedtools sort -i "$BED_FILE" 2>/dev/null) "$BED_FILE" > /dev/null 2>&1; then
    echo "File is sorted (lexicographic)"
else
    echo "File is NOT lexicographically sorted (bedtools sort -i $BED_FILE)"
fi
