#!/bin/bash
#SBATCH --job-name eSLAM_Filtering_Third_part
#SBATCH --time 1:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 48
#SBATCH --output eSLAM_Filtering_Third_part.%A_%a.log
#SBATCH --array 1-28

echo "PART 4 (B) of SLAM-SEQ Analysis"
echo "This script will create sorted bam files with 1, 2, 3, 4, >4 point mutations"

INDIR_RAW="/scratch/mugolini/eSLAM/RAW_DATA_MERGED/"
OUTDIR="/scratch/mugolini/eSLAM/HISAT3N_Output/"

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist.txt)

echo "Processing $FILE"

module load gcc
module load samtools

for idx in '1' '2' '3' '4' 'more.4'
do
  echo "Creating filtered .bam files"
  samtools view -h -b --threads $SLURM_CPUS_PER_TASK -N ${OUTDIR}mismatch_infos/${FILE}.Yf.table.${idx}.tab ${OUTDIR}${FILE}_aligned.unique.filtered.sorted.bam > ${OUTDIR}${FILE}_aligned.unique.filtered.${idx}.bam
  echo "Sort and index filtered bam file"
  samtools sort --threads $SLURM_CPUS_PER_TASK -o ${OUTDIR}${FILE}_aligned.unique.filtered.sorted.${idx}.bam ${OUTDIR}${FILE}_aligned.unique.filtered.${idx}.bam
  samtools index -@ $SLURM_CPUS_PER_TASK ${OUTDIR}${FILE}_aligned.unique.filtered.sorted.${idx}.bam
  echo "Removing unsorted filtered .bam files"
  rm ${OUTDIR}${FILE}_aligned.unique.filtered.${idx}.bam
done

echo "script finished"
