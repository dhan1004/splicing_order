#!/bin/bash

# Script to filter TSV for only paired-end samples
# Usage: bash filter_paired_samples.sh

INPUT_TSV="/users/dhan30/scratch/splicing_order/data/gsm_sra_list_for_pipeline_hek293.tsv"
OUTPUT_TSV="/users/dhan30/scratch/splicing_order/data/gsm_sra_list_paired_only.tsv"
TEMP_DIR="/users/dhan30/scratch/splicing_order/data/temp_check"

mkdir -p "$TEMP_DIR"

# Activate conda environment for fastq-dl
module load miniconda3/23.11.0s
source /oscar/runtime/software/external/miniconda3/23.11.0/etc/profile.d/conda.sh
conda activate order_env

# Copy header
head -n 1 "$INPUT_TSV" > "$OUTPUT_TSV"

echo "Checking samples for paired-end data..."
echo "This may take a while..."

total_samples=0
paired_samples=0
single_samples=0

# Read TSV (skip header)
tail -n +2 "$INPUT_TSV" | while IFS=$'\t' read -r gsm_id sra_ids; do
    total_samples=$((total_samples + 1))
    
    # Clean up inputs
    gsm_id=$(echo "$gsm_id" | xargs)
    sra_ids=$(echo "$sra_ids" | tr -d '\r' | xargs)
    
    if [ -z "$gsm_id" ] || [ -z "$sra_ids" ]; then
        continue
    fi
    
    # Get first SRA ID to check
    first_srr=$(echo "$sra_ids" | cut -d',' -f1)
    
    echo -n "Checking $gsm_id ($first_srr)... "
    
    # Download metadata only (no actual download)
    fastq-dl --only-download-metadata --outdir "$TEMP_DIR" --accession "$first_srr" 2>/dev/null
    
    # Check if it's paired-end from the metadata
    if [ -f "$TEMP_DIR/fastq-run-info.tsv" ]; then
        # Check for paired-end indicators in metadata
        if grep -q "PAIRED" "$TEMP_DIR/fastq-run-info.tsv" || \
           grep -q "_1.fastq" "$TEMP_DIR/fastq-run-info.tsv"; then
            echo "PAIRED-END ✓"
            echo -e "${gsm_id}\t${sra_ids}" >> "$OUTPUT_TSV"
            paired_samples=$((paired_samples + 1))
        else
            echo "single-end (skipping)"
            single_samples=$((single_samples + 1))
        fi
        
        rm -f "$TEMP_DIR/fastq-run-info.tsv" "$TEMP_DIR/fastq-run-mergers.tsv"
    else
        echo "ERROR checking metadata"
    fi
done

rm -rf "$TEMP_DIR"

echo ""
echo "========================================="
echo "Summary:"
echo "  Total samples checked: $total_samples"
echo "  Paired-end samples: $paired_samples"
echo "  Single-end samples: $single_samples"
echo ""
echo "Filtered TSV saved to: $OUTPUT_TSV"
echo "========================================="