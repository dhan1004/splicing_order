#!/bin/bash

# SLURM directives (if running as standalone job)
#SBATCH --job-name=order_of_splicing
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

################################################################################
# CONFIGURATION
################################################################################

# Inputs
threads=$1
out_dir=$2
gsm_id=$3
srr_id_input=$4

printf "$(date +'%d/%b/%Y %H:%M:%S') | Processing sample %s...\n" "${gsm_id}"

# Raw & trimmed files (PAIRED-END ONLY)
paired_read_one_file="${out_dir}/${gsm_id}_all_reads_1.fq.gz"
paired_read_two_file="${out_dir}/${gsm_id}_all_reads_2.fq.gz"
paired_trimmed_one_file="${out_dir}/${gsm_id}_trimmed_1.fq.gz"
paired_trimmed_two_file="${out_dir}/${gsm_id}_trimmed_2.fq.gz"

# STAR alignment
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
reference_genome_dir="/users/dhan30/reference/hg38"
intron_bed_file="/users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz"

# Script directory
SCRIPT_DIR="/users/dhan30/scratch/splicing_order/scripts"

# Ensure output directory exists
mkdir -p "${out_dir}"

################################################################################
# ENVIRONMENT SETUP
################################################################################

printf "$(date +'%d/%b/%Y %H:%M:%S') | Loading conda environment...\n"
module load miniconda3/23.11.0s
source /oscar/runtime/software/external/miniconda3/23.11.0/etc/profile.d/conda.sh
conda activate order_env

printf "$(date +'%d/%b/%Y %H:%M:%S') | Active environment: $CONDA_DEFAULT_ENV\n"

################################################################################
# CHECKPOINT 1: Check if trimmed files exist
################################################################################

skip_to_alignment=true

# if [ -f "$paired_trimmed_one_file" ] && [ -f "$paired_trimmed_two_file" ]; then
#     printf "$(date +'%d/%b/%Y %H:%M:%S') | Trimmed files found. Skipping download and trimming...\n"
#     skip_to_alignment=true
# else
#     skip_to_alignment=false
# fi

################################################################################
# DOWNLOAD & TRIM (if needed)
################################################################################

if [ "$skip_to_alignment" = false ]; then

    # --- CHECKPOINT 2: Check if raw reads exist ---
    if [ -f "$paired_read_one_file" ] && [ -f "$paired_read_two_file" ] && [ -s "$paired_read_one_file" ]; then
        printf "$(date +'%d/%b/%Y %H:%M:%S') | Raw paired reads found. Skipping download...\n"
    else
        # Switch to fastq-dl environment for downloading
        conda activate /users/dhan30/.conda/envs/fastq_dl_env
        
        printf "$(date +'%d/%b/%Y %H:%M:%S') | Downloading reads for %s...\n" "${gsm_id}"
        IFS=',' read -r -a srr_id_list <<< "$srr_id_input"

        for srr in "${srr_id_list[@]}"
        do
            printf "$(date +'%d/%b/%Y %H:%M:%S') |   Downloading %s...\n" "${srr}"
            
            fastq-dl -F --outdir "${out_dir}" --cpus "${threads}" --accession "${srr}" ||
                { printf "ERROR: fastq-dl failed for %s\n" "${srr}"; exit 1; }

            srr_output_prefix="${out_dir}/${srr}"

            # Verify paired-end files were created
            if [ -f "${srr_output_prefix}_1.fastq.gz" ] && [ -f "${srr_output_prefix}_2.fastq.gz" ]; then
                cat "${srr_output_prefix}_1.fastq.gz" >> "$paired_read_one_file"
                cat "${srr_output_prefix}_2.fastq.gz" >> "$paired_read_two_file"
                rm "${srr_output_prefix}_1.fastq.gz" "${srr_output_prefix}_2.fastq.gz"
            else
                printf "ERROR: Expected paired-end data for %s but didn't find _1 and _2 files.\n" "${srr}"
                printf "This pipeline only supports paired-end RNA-seq data.\n"
                rm -rf "${out_dir}"

                exit 1
            fi

            rm -f "${out_dir}/fastq-run-info.tsv" "${out_dir}/fastq-run-mergers.tsv"
        done

        printf "$(date +'%d/%b/%Y %H:%M:%S') | Download complete for %s.\n" "${gsm_id}"
        
        # Switch back to analysis environment
        conda deactivate
        conda activate /users/dhan30/.conda/envs/order_env
    fi

    # --- Validate we have paired-end data ---
    if [ ! -f "$paired_read_one_file" ] || [ ! -f "$paired_read_two_file" ]; then
        printf "ERROR: Paired-end read files not found for %s\n" "${gsm_id}"
        printf "This pipeline requires paired-end RNA-seq data.\n"
        exit 1
    fi

    # --- Trimming with fastp ---
    # Switch to fastq_dl_env if fastp is there
    
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Trimming reads for %s...\n" "${gsm_id}"
    fastp_report_json="${out_dir}/${gsm_id}_fastp.json"
    fastp_report_html="${out_dir}/${gsm_id}_fastp.html"

    fastp \
        --in1 "$paired_read_one_file" \
        --in2 "$paired_read_two_file" \
        --out1 "$paired_trimmed_one_file" \
        --out2 "$paired_trimmed_two_file" \
        --detect_adapter_for_pe \
        --thread "$threads" \
        --json "$fastp_report_json" \
        --html "$fastp_report_html" ||
        { printf "ERROR: fastp failed for %s\n" "${gsm_id}"; exit 1; }

    rm "$paired_read_one_file" "$paired_read_two_file"
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Trimming complete for %s.\n" "${gsm_id}"
fi

################################################################################
# CHECKPOINT 3: Check if alignment exists
################################################################################

if [ -f "$aligned_bam_file" ] && [ -s "$aligned_bam_file" ]; then
    printf "$(date +'%d/%b/%Y %H:%M:%S') | STAR alignment file found. Skipping alignment...\n"
else
    # --- STAR Alignment ---
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Aligning reads for %s using STAR...\n" "${gsm_id}"

    STAR --genomeDir "$reference_genome_dir" \
        --readFilesIn "$paired_trimmed_one_file" "$paired_trimmed_two_file" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$star_out_prefix" \
        --outSAMtype BAM SortedByCoordinate \
        --runThreadN "$threads" \
        --twopassMode Basic ||
        { printf "ERROR: STAR alignment failed for %s\n" "${gsm_id}"; exit 1; }

    if [ ! -f "$aligned_bam_file" ]; then
        printf "ERROR: STAR alignment failed for %s. BAM file not found.\n" "${gsm_id}"
        exit 1
    fi

    printf "$(date +'%d/%b/%Y %H:%M:%S') | STAR alignment complete for %s.\n" "${gsm_id}"
fi

################################################################################
# EXTRACT INFORMATIVE PAIRS (Kim et al. methodology)
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
# rm "${PAIRED_JUNCS_SAM}" "${INTRON_OVERLAP_BED}" 
# rm "${SORTED_BAM_PREFIX}.bam" "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt"

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
printf "$(date +'%d/%b/%Y %H:%M:%S') | Results: %s\n" "${splicing_order_output}"

conda deactivate