#!/bin/bash
#SBATCH --job-name=add_hbonds
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=hbond_logs/hbond_%A_%a.out
#SBATCH --error=hbond_logs/hbond_%A_%a.err

# Configuration
RESULTS_DIR="/users/dhan30/scratch/splicing_order/results"
SCRIPT="/users/dhan30/scratch/splicing_order/scripts/add_hydrogen_bonds.py"
FASTA="/users/dhan30/reference/hg38.fa"

# Create log directory
mkdir -p hbond_logs

# Load conda environment
module load miniconda3/23.11.0s
source /oscar/runtime/software/external/miniconda3/23.11.0/etc/profile.d/conda.sh
conda activate order_env
export PATH="/users/dhan30/.conda/envs/order_env/bin:$PATH"

# Get sample directory and create output filename
INPUT_FILE="/users/dhan30/scratch/splicing_order/data/structure_chunks/structure_chunk_1.tsv"
SAMPLE_DIR=$(dirname "$INPUT_FILE")
BASENAME=$(basename "$INPUT_FILE" .tsv)
OUTPUT_FILE="${SAMPLE_DIR}/${BASENAME}_with_hbonds.tsv"

echo "Processing: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"

# Run the script
python3 "$SCRIPT" \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE" \
    --fasta "$FASTA" \
    --padding 50

if [ $? -eq 0 ]; then
    echo "Successfully processed $INPUT_FILE"
else
    echo "Failed to process $INPUT_FILE"
    exit 1
fi