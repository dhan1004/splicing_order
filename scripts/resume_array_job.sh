#!/bin/bash
#SBATCH --job-name=resume_splicing
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=out_files/resume_%A_%a.out
#SBATCH --error=out_files/resume_%A_%a.err
#SBATCH --partition=batch

# Resume pipeline for samples that have STAR alignments but no informative pairs

INPUT_LIST_TSV="/users/dhan30/scratch/splicing_order/data/gsm_sra_list_for_pipeline_hek293.tsv"
THREADS=8
BASE_OUT_DIR="/users/dhan30/scratch/splicing_order/results/"
RESUME_SCRIPT="/users/dhan30/scratch/splicing_order/scripts/resume_from_alignment.sh"

# Get sample info for this array task
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${LINE_NUM}p" "$INPUT_LIST_TSV")

gsm_id=$(echo "$LINE" | cut -f1 | xargs)

if [ -z "$gsm_id" ]; then
    echo "ERROR: Empty GSM ID for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

printf "$(date +'%d/%b/%Y %H:%M:%S') | Resuming GSM: $gsm_id\n"

sample_out_dir="${BASE_OUT_DIR}/${gsm_id}"

# Run the resume script
"$RESUME_SCRIPT" "$THREADS" "$sample_out_dir" "$gsm_id"

if [ $? -eq 0 ]; then
    printf "$(date +'%d/%b/%Y %H:%M:%S') | SUCCESS: $gsm_id\n"
else
    printf "$(date +'%d/%b/%Y %H:%M:%S') | FAILED: $gsm_id\n"
    exit 1
fi