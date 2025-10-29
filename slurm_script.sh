#!/bin/bash
#SBATCH -n 10
#SBATCH --mem=128G
#SBATCH --job-name=order_splice_genome
#SBATCH --output=out_files/order_genome_match_%j.out
#SBATCH --error=out_files/order_genome_match_%j.err
# #SBATCH --nodelist=masternode

printf "$(date +'%d/%b/%Y %H:%M:%S') | Running on node: $(hostname)\n"
printf "$(date +'%d/%b/%Y %H:%M:%S') | SLURM assigned node: $SLURM_NODELIST\n"
printf "$(date +'%d/%b/%Y %H:%M:%S') | Job ID: $SLURM_JOB_ID\n"

bash /scratch/dhan/order_of_splicing/scripts/run_filtering_all.sh