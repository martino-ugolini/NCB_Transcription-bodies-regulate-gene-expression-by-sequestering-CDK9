#!/bin/bash
#SBATCH --job-name ATAC_Seq
#SBATCH --time 00:30:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --output ATAC_Seq.%A_%a.log
#SBATCH --array 1-10

INDIR_RAW="/scratch/mugolini/ATAC/RAW_DATA/"
OUTDIR_BOWTIE2="/scratch/mugolini/ATAC/Bowtie2_Output/"

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist.txt)

module load gcc
module load trimmomatic
module load fastqc
module load bowtie2
module load samtools
module load picard
module load bedtools2

echo "Sorting by name..."
samtools sort -@ $SLURM_CPUS_PER_TASK -o ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.sortedbycoordinates.bam ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.sortedbyname.bam
samtools index -@ $SLURM_CPUS_PER_TASK ${OUTDIR_BOWTIE2}${FILE}_aligned.q30.noMT.dupremoved.sortedbycoordinates.bam

echo "Script finished"

# https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/
