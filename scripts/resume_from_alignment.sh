#!/bin/bash

# Resume pipeline from existing STAR alignments
# This skips download, trimming, and alignment

# Inputs
threads=${1:-8}
out_dir=$2
gsm_id=$3

printf "$(date +'%d/%b/%Y %H:%M:%S') | Resuming processing for %s...\n" "${gsm_id}"

# STAR alignment (should already exist)
star_out_prefix="${out_dir}/${gsm_id}_STAR_"
aligned_bam_file="${star_out_prefix}Aligned.sortedByCoord.out.bam"

# Intermediate files
TEMP_PREFIX="temp_${gsm_id}"
ALL_PAIRED_SAM="${out_dir}/${TEMP_PREFIX}_all_paired_reads.sam"
PAIRED_JUNCS_SAM="${out_dir}/${TEMP_PREFIX}_junction_paired_reads.sam"
INTRON_OVERLAP_BED="${out_dir}/${TEMP_PREFIX}_intron_reads.bed"
SORTED_BAM_PREFIX="${out_dir}/${TEMP_PREFIX}_intermediate_reads_sorted"

# Final outputs
informative_pairs_bam="${out_dir}/${gsm_id}_informative_pairs.bam"
splicing_order_output="${out_dir}/${gsm_id}_pairwise_splicing_order.tsv"

# Reference files
intron_bed_file="/users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz"

# Script directory
SCRIPT_DIR="/users/dhan30/scratch/splicing_order/scripts"

################################################################################
# ENVIRONMENT SETUP
################################################################################

module load miniconda3/23.11.0s
source /oscar/runtime/software/external/miniconda3/23.11.0/etc/profile.d/conda.sh
conda activate order_env

printf "$(date +'%d/%b/%Y %H:%M:%S') | Active environment: $CONDA_DEFAULT_ENV\n"

################################################################################
# CHECK IF ALIGNED BAM EXISTS
################################################################################

if [ ! -f "$aligned_bam_file" ]; then
    printf "ERROR: STAR alignment not found: $aligned_bam_file\n"
    printf "This script requires existing STAR output. Run full pipeline instead.\n"
    exit 1
fi

printf "$(date +'%d/%b/%Y %H:%M:%S') | Found STAR alignment with $(samtools view -c $aligned_bam_file) reads\n"

################################################################################
# DELETE OLD INTERMEDIATE FILES (save space)
################################################################################

printf "$(date +'%d/%b/%Y %H:%M:%S') | Cleaning old intermediate files...\n"
rm -f "${ALL_PAIRED_SAM}" "${PAIRED_JUNCS_SAM}" "${INTRON_OVERLAP_BED}" 2>/dev/null
rm -f "${SORTED_BAM_PREFIX}.bam" "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt" 2>/dev/null
rm -f "${out_dir}"/*_all_reads*.fq.gz "${out_dir}"/*_trimmed*.fq.gz 2>/dev/null

################################################################################
# EXTRACT INFORMATIVE PAIRS
################################################################################

printf "$(date +'%d/%b/%Y %H:%M:%S') | Extracting informative pairs...\n"

# Step 1: Convert to name-sorted SAM
printf "$(date +'%d/%b/%Y %H:%M:%S') |   Converting to name-sorted SAM...\n"
samtools sort -n -O SAM -@ "$threads" -o "${ALL_PAIRED_SAM}" "${aligned_bam_file}"

# Step 2: Extract junction-containing pairs
printf "$(date +'%d/%b/%Y %H:%M:%S') |   Extracting junction pairs...\n"
awk 'BEGIN {OFS="\t"}
    /^@/ {print; next}
    $6 ~ /N/ {junction_reads[$1]=1}
    END {
        while ((getline < "'"${ALL_PAIRED_SAM}"'") > 0) {
            if ($0 ~ /^@/) continue
            if ($1 in junction_reads) print
        }
    }' "${ALL_PAIRED_SAM}" > "${PAIRED_JUNCS_SAM}"

rm "${ALL_PAIRED_SAM}"

# Step 3: Find intron-overlapping reads
printf "$(date +'%d/%b/%Y %H:%M:%S') |   Finding intron overlaps (f=0.05)...\n"
samtools view -bh "${PAIRED_JUNCS_SAM}" | \
    bedtools bamtobed -split -i stdin | \
    bedtools intersect -f 0.05 -a stdin -b "${intron_bed_file}" > "${INTRON_OVERLAP_BED}"

# Step 4: Extract pairs where one has junction AND one overlaps intron
printf "$(date +'%d/%b/%Y %H:%M:%S') |   Extracting informative pairs...\n"
awk '{print $4}' "${INTRON_OVERLAP_BED}" | \
    sed 's/\/[12]$//' | \
    sort -u > "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt"

(grep "^@" "${PAIRED_JUNCS_SAM}"; \
    grep -wFf "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt" "${PAIRED_JUNCS_SAM}") | \
    samtools view -bh - | \
    samtools fixmate -m - - | \
    samtools sort -@ "$threads" -o "${SORTED_BAM_PREFIX}.bam" - ||
    { printf "ERROR: Failed to create sorted BAM\n"; exit 1; }

# Step 5: Remove duplicates
printf "$(date +'%d/%b/%Y %H:%M:%S') |   Removing duplicates...\n"
samtools markdup -r -@ "$threads" "${SORTED_BAM_PREFIX}.bam" "${informative_pairs_bam}"
samtools index "${informative_pairs_bam}"

# Cleanup intermediate files
rm "${PAIRED_JUNCS_SAM}" "${INTRON_OVERLAP_BED}" 
rm "${SORTED_BAM_PREFIX}.bam" "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt"

printf "$(date +'%d/%b/%Y %H:%M:%S') | Informative pairs extraction complete.\n"
printf "  Total reads: $(samtools view -c ${informative_pairs_bam})\n"
printf "  Read pairs: $(($(samtools view -c ${informative_pairs_bam}) / 2))\n"

################################################################################
# ANALYZE SPLICING ORDER
################################################################################

printf "$(date +'%d/%b/%Y %H:%M:%S') | Analyzing splicing order...\n"

python3 "${SCRIPT_DIR}/analyze_splicing_order.py" \
    --bam "${informative_pairs_bam}" \
    --intron-bed "${intron_bed_file}" \
    --output "${splicing_order_output}" \
    --min-reads 5 \
    --min-mapq 10 \
    --tolerance 10 ||
    { printf "ERROR: Splicing order analysis failed\n"; exit 1; }

printf "$(date +'%d/%b/%Y %H:%M:%S') | Sample %s processing complete!\n" "${gsm_id}"

conda deactivate