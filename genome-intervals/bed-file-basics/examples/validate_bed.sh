#!/bin/bash
# Validate BED file format and content

BED_FILE="${1:-input.bed}"

echo "Validating: $BED_FILE"
echo "========================"

# Check if file exists
if [[ ! -f "$BED_FILE" ]]; then
    echo "ERROR: File not found"
    exit 1
fi

# Count lines
LINES=$(wc -l < "$BED_FILE")
echo "Total lines: $LINES"

# Check column count
echo -e "\nColumn counts:"
awk -F'\t' '{print NF}' "$BED_FILE" | sort | uniq -c

# Check for invalid coordinates (start >= end)
INVALID=$(awk -F'\t' '$2 >= $3' "$BED_FILE" | wc -l)
if [[ $INVALID -gt 0 ]]; then
    echo -e "\nWARNING: $INVALID intervals with start >= end:"
    awk -F'\t' '$2 >= $3' "$BED_FILE" | head -5
else
    echo -e "\nNo invalid intervals (start >= end): OK"
fi

# Check for negative coordinates
NEGATIVE=$(awk -F'\t' '$2 < 0 || $3 < 0' "$BED_FILE" | wc -l)
if [[ $NEGATIVE -gt 0 ]]; then
    echo -e "\nWARNING: $NEGATIVE intervals with negative coordinates"
else
    echo -e "No negative coordinates: OK"
fi

# Check if sorted
if bedtools sort -i "$BED_FILE" 2>/dev/null | diff -q - <(cat "$BED_FILE") > /dev/null 2>&1; then
    echo "File is sorted: OK"
else
    echo "File is NOT sorted (consider: bedtools sort -i $BED_FILE)"
fi

# Check chromosome naming
echo -e "\nChromosomes found:"
cut -f1 "$BED_FILE" | sort -u | head -10

# Summary statistics
echo -e "\nInterval size statistics:"
awk -F'\t' '{print $3-$2}' "$BED_FILE" | \
    awk 'BEGIN{min=999999999; max=0; sum=0; n=0}
         {if($1<min)min=$1; if($1>max)max=$1; sum+=$1; n++}
         END{print "  Min: "min"\n  Max: "max"\n  Mean: "sum/n"\n  Count: "n}'

echo -e "\nValidation complete"
