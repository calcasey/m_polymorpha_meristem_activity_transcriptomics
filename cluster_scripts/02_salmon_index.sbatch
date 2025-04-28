#!/usr/bin/env bash                                                                                                                                                                                         

#SBATCH --job-name=salmon_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=index_transcriptome_out.txt

ml load build-env/.f2021
ml load build-env/f2021
ml load salmon/1.5.2-gompi-2020b

/software/f2021/software/salmon/1.5.2-gompi-2020b/bin/salmon index \
    -t MpTak_v6.1r1.mrna.fasta.gz \
    -i MpTak_v6_1r1_mrna_index \
    -k 21
