#!/bin/bash
#SBATCH --job-name=structure_array
#SBATCH --array=1-500%5
#SBATCH --output=structure_array_%A_%a.out
#SBATCH --error=structure_array_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch

# Load environment
module load miniconda3/23.11.0s
source /oscar/runtime/software/external/miniconda3/23.11.0/etc/profile.d/conda.sh
conda activate /users/dhan30/.conda/envs/order_env
export PATH="/users/dhan30/.conda/envs/order_env/bin:$PATH"

echo "Job started: $(date)"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo ""

# Input file
INPUT_FILE="/users/dhan30/scratch/splicing_order/data/foldable_intron_pairs.tsv"
OUTPUT_DIR="/users/dhan30/scratch/splicing_order/data/structure_chunks"
mkdir -p "$OUTPUT_DIR"

# Calculate how many lines total (excluding header)
TOTAL_LINES=$(tail -n +2 "$INPUT_FILE" | wc -l)
NUM_CHUNKS=500  # Must match array size
LINES_PER_CHUNK=$(( (TOTAL_LINES + NUM_CHUNKS - 1) / NUM_CHUNKS ))

echo "Total lines: $TOTAL_LINES"
echo "Lines per chunk: $LINES_PER_CHUNK"
echo ""

# Create chunk for this task
START_LINE=$(( (SLURM_ARRAY_TASK_ID - 1) * LINES_PER_CHUNK + 2 ))  # +2 for header
END_LINE=$(( START_LINE + LINES_PER_CHUNK - 1 ))

CHUNK_FILE="$OUTPUT_DIR/chunk_${SLURM_ARRAY_TASK_ID}.tsv"
OUTPUT_FILE="$OUTPUT_DIR/structure_chunk_${SLURM_ARRAY_TASK_ID}.tsv"

# Extract chunk (with header)
head -1 "$INPUT_FILE" > "$CHUNK_FILE"
sed -n "${START_LINE},${END_LINE}p" "$INPUT_FILE" >> "$CHUNK_FILE"

echo "Processing lines $START_LINE to $END_LINE"
echo "Chunk file: $CHUNK_FILE"
echo ""

# Run structure extraction on this chunk
python3 /users/dhan30/scratch/splicing_order/scripts/extract_structure_features.py \
    --intron-pairs "$CHUNK_FILE" \
    --fasta /users/dhan30/reference/hg38.fa \
    --output "$OUTPUT_FILE" \
    --padding 50

echo ""
echo "Job finished: $(date)"

# Clean up chunk file
rm "$CHUNK_FILE"

conda deactivate