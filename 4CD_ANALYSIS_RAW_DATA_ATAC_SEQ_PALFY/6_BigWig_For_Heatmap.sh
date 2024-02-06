#!/bin/bash
#SBATCH --job-name ATAC_Seq
#SBATCH --time 01:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 24
#SBATCH --output ATAC_Seq.%A_%a.log
#SBATCH --array 1-4

INDIR_RAW="/scratch/mugolini/ATAC/RAW_DATA/"
OUTDIR_BOWTIE2="/scratch/mugolini/ATAC/Bowtie2_Output/"


FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist2.txt)

module load gcc

echo "ATAC-Seq Analysis"
echo "Running analysis on $FILE..."

/users/mugolini/.local/bin/bamCoverage -b ${OUTDIR_BOWTIE2}${FILE}.sortedbycoordinates.bam -o ${OUTDIR_BOWTIE2}${FILE}.sortedbycoordinates.RPKM.bw \
--outFileFormat bigwig --normalizeUsing RPKM --binSize 20 --extendReads --numberOfProcessors 20

echo "Script finished"
