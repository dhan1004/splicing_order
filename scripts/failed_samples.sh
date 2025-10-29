#!/bin/bash

# Automatically find and resubmit failed samples

ORIGINAL_TSV="/users/dhan30/scratch/splicing_order/data/gsm_sra_list_paired_only.tsv"
FAILED_TSV="/users/dhan30/scratch/splicing_order/data/gsm_sra_list_failed.tsv"
RESULTS_DIR="/users/dhan30/scratch/splicing_order/results"

echo "Finding failed samples..."

# Copy header
head -n 1 "$ORIGINAL_TSV" > "$FAILED_TSV"

# Find samples without final output
tail -n +2 "$ORIGINAL_TSV" | while IFS=$'\t' read -r gsm_id sra_ids; do
    gsm_id=$(echo "$gsm_id" | xargs)
    
    if [ -z "$gsm_id" ]; then
        continue
    fi
    
    final_output="$RESULTS_DIR/$gsm_id/${gsm_id}_pairwise_splicing_order.tsv"
    
    if [ ! -f "$final_output" ] || [ ! -s "$final_output" ]; then
        echo -e "${gsm_id}\t${sra_ids}" >> "$FAILED_TSV"
    fi
done

# Count failed samples
N=$(tail -n +2 "$FAILED_TSV" | wc -l)

if [ "$N" -eq 0 ]; then
    echo "No failed samples found! All samples completed successfully."
    exit 0
fi

echo "Found $N failed samples"
echo ""
echo "Failed samples TSV: $FAILED_TSV"
echo ""
