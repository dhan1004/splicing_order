#!/bin/bash


# This wrapper script reads a TSV file containing GSM IDs and their SRA IDs
# and runs the main processing script for each sample.

# --- Configuration ---
# File containing GSM IDs and their SRA IDs (Tab-Separated Values)
# Expected header: gsm_id<TAB>sra_ids
INPUT_LIST_TSV="/scratch/dhan/order_of_splicing/data/gsm_sra_list_for_pipeline_hek293.tsv" # Make sure this file exists and is correctly formatted

# Number of threads for processing each sample
THREADS=8

# Base output directory for all processed samples
BASE_OUT_DIR="/scratch/dhan/order_of_splicing/results/"

# Path to main processing script
PIPELINE_SCRIPT="/scratch/dhan/order_of_splicing/scripts/order_gse_filtering.sh"
# --- End Configuration ---

# Ensure the main pipeline script is executable
if [ ! -x "$PIPELINE_SCRIPT" ]; then
    echo "ERROR: Pipeline script '$PIPELINE_SCRIPT' not found or not executable."
    echo "Please ensure it's in the current directory and has execute permissions (chmod +x $PIPELINE_SCRIPT)."
    exit 1
fi

# Ensure the input TSV file exists
if [ ! -f "$INPUT_LIST_TSV" ]; then
    echo "ERROR: Input TSV file not found: $INPUT_LIST_TSV"
    exit 1
fi

# Create base output directory if it doesn't exist
mkdir -p "$BASE_OUT_DIR"

# Activate conda environment
# echo "Activating conda environment..."
# source /home/dhan/anaconda3/envs/order_env # Ensure this path is correct

# Read the input TSV file line by line, skipping the header row
# tail -n +2 skips the first line (header)
tail -n +2 "$INPUT_LIST_TSV" | while IFS=$'\t' read -r gsm_id_col sra_ids_col; do
    # Trim potential whitespace or carriage returns (especially if file came from Windows)
    gsm_id=$(echo "$gsm_id_col" | xargs)
    sra_ids=$(echo "$sra_ids_col" | tr -d '\r' | xargs)

    # Skip if gsm_id or sra_ids are empty after trimming
    if [ -z "$gsm_id" ] || [ -z "$sra_ids" ]; then
        echo "Skipping line with empty GSM ID ('$gsm_id_col') or SRA IDs ('$sra_ids_col')."
        continue
    fi

    echo "----------------------------------------------------------------------"
    echo "Preparing to process GSM: $gsm_id, SRAs/SRX: $sra_ids"
    echo "----------------------------------------------------------------------"

    # Define a specific output directory for this sample
    sample_out_dir="${BASE_OUT_DIR}/${gsm_id}"
    mkdir -p "$sample_out_dir" # Ensure sample-specific output dir exists

    # Run the pipeline script
    # Standard output and standard error will be redirected to a log file for each sample
    log_file="${sample_out_dir}/${gsm_id}_pipeline.log"
    echo "Logging to: $log_file"

    "$PIPELINE_SCRIPT" "$THREADS" "$sample_out_dir" "$gsm_id" "$sra_ids" > "$log_file" 2>&1

    # Check the exit status of the pipeline script
    if [ $? -eq 0 ]; then
        echo "Pipeline completed successfully for $gsm_id."
    else
        echo "Pipeline FAILED for $gsm_id. Check log: $log_file"
    fi
    echo "----------------------------------------------------------------------"
done

echo "All samples from $INPUT_LIST_TSV have been processed."
# Optionally deactivate conda environment if you activated it only for this script
# conda deactivate