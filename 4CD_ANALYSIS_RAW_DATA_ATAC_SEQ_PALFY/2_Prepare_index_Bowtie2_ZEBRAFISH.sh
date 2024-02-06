#!/bin/bash
#SBATCH --job-name Bowtie2_index
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 48
#SBATCH --mem 300000
#SBATCH --output Bowtie2_index.log

echo "Bowtie 2 index"

GENOME="/scratch/mugolini/GRCz11_105_GENOME/Danio_rerio.GRCz11.105.dna.primary_assembly.fa"

gunzip ${GENOME}.gz

module load gcc
module load bowtie2

bowtie2-build -f --threads $SLURM_CPUS_PER_TASK ${GENOME} genome_zebrafish

echo "library created succesfully"
echo "Script finished"
