#!/usr/bin/env bash

#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:00:00
#SBATCH --array=1-49
#SBATCH --output=out_quality_trim_%a.txt

module load build-env/2020
module load trimmomatic/0.38-java-1.8

java -jar /software/2020/software/trimmomatic/0.38-java-1.8/trimmomatic-0.38.jar \
    PE \
    -threads 20 \
    -phred33 \
    R1_${SLURM_ARRAY_TASK_ID}.fastq.gz \
    R2_${SLURM_ARRAY_TASK_ID}.fastq.gz \
    trimmed_${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
    unpaired_${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
    trimmed_${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
    unpaired_${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
    ILLUMINACLIP:all_adaptors.fasta:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:15 MINLEN:25
