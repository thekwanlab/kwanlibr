#!/usr/bin/env bash

allDB_SITES=$1
TSS_BED=$2
SORTED_BED=$3
LOOKUP_TSV=$4

module load Bioinformatics bedops/2.4.41
bedops --version

echo "Generating TSS regions bed file $TSS_BED"

awk -vOFS='\t' '($8 == "gene") { 
        if ($6 == "+") { print $1, $2, $2 + 1, substr($15, 2, length($15)-3), $5 } 
        else { print $1, $3 - 1, $3, substr($15, 2, length($15)-3), $5 }
    }' \
    /nfs/turbo/umms-kykwan/projects/reference/gtf/gencode_annotations.bed \
    > "$TSS_BED"

echo "Sorting peak sites bed file to $SORTED_BED"
sort-bed "$allDB_SITES" > "$SORTED_BED"

echo "Finding closest features and writing to: $LOOKUP_TSV"
closest-features --closest --delim '\t' \
    "$SORTED_BED" "$TSS_BED" \
    | awk '{ print $4, $9 }' \
    > "$LOOKUP_TSV"

echo "=== Annotation complete ==="
