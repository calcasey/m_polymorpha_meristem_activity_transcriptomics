#!/usr/bin/env bash                                                                                                                                                                                         

#SBATCH --job-name=salmon_quant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --qos=short
#SBATCH --time=0-08:00:00
#SBATCH --array=1-49
#SBATCH --output=out_salmon_%a.txt

ml load build-env/.f2021
ml load build-env/f2021
ml load salmon/1.5.2-gompi-2020b

/software/f2021/software/salmon/1.5.2-gompi-2020b/bin/salmon quant \
    -i MpTak_v6_1r1_mrna_index \
    -l A \
    -1 trimmed_${SLURM_ARRAY_TASK_ID}_R1.fastq \
    -2 trimmed_${SLURM_ARRAY_TASK_ID}_R2.fastq \
    --threads 12 \
    --validateMappings \
    -o quant_${SLURM_ARRAY_TASK_ID}
