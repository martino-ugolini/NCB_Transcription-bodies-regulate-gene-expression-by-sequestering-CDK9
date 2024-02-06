#!/bin/bash
#SBATCH --job-name eSLAM_Filter_1
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 48
#SBATCH --mem 400G
#SBATCH --output eSLAM_Filter_1.%A_%a.log
#SBATCH --array 1-28

echo "PART 4 (A) of SLAM-SEQ Analysis"
echo "This script will filter the bam file based on the number of detected mutations"
echo "The python scripts are based on the pulseR package (GitHub and see Materials and Methods)"

INDIR_RAW="/scratch/mugolini/eSLAM/RAW_DATA_MERGED/"
OUTDIR="/scratch/mugolini/eSLAM/HISAT3N_Output/"
SNP="/scratch/mugolini/eSLAM/5_Make_SNP_file/SNP_output.snp"

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist.txt)

echo "Processing $FILE"

module load gcc
module load python
module load samtools

mkdir ${OUTDIR}mismatch_infos/

#Note: The 48CPUs and 400G are also hardcoded as default values in the splbam python script
echo "Start splbam"
python splbam.py --subtract $SNP --ref-base 'T' --base-change 'C' --base-qual 20 --trim5p 0 --trim3p 0 --overwrite --num-cpus $SLURM_CPUS_PER_TASK --mem 400G ${OUTDIR}${FILE}_aligned.unique.sorted.bam ${OUTDIR}mismatch_infos/ ${FILE}

echo "Start filtering bam file"
samtools view -h -b --threads $SLURM_CPUS_PER_TASK -N ${OUTDIR}mismatch_infos/${FILE}.trueconversions_list.tab ${OUTDIR}${FILE}_aligned.unique.sorted.bam > ${OUTDIR}${FILE}_aligned.unique.filtered.bam

echo "Sort and index filtered bam file"
samtools sort --threads $SLURM_CPUS_PER_TASK -o ${OUTDIR}${FILE}_aligned.unique.filtered.sorted.bam ${OUTDIR}${FILE}_aligned.unique.filtered.bam
samtools index -@ $SLURM_CPUS_PER_TASK ${OUTDIR}${FILE}_aligned.unique.filtered.sorted.bam

echo "script finished"
