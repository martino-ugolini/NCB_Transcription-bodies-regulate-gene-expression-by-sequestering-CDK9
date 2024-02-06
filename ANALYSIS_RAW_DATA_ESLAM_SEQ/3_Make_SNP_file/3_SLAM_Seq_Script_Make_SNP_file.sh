#!/bin/bash
#SBATCH --job-name eSLAM_Make_SNP_file
#SBATCH --time 4:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --output eSLAM_Make_SNP_file.%A_%a.log
#SBATCH --array 1-4

echo "PART 3 of eSLAM-SEQ Analysis"
echo "This script will identify SNPs"
echo "Create a filelist_IAA_minus.txt file that contains the names of the 4 samples that were NOT treated with IAA"

echo "Make sure that the genome.fa is unzipped"

INDIR_RAW="/scratch/mugolini/eSLAM/RAW_DATA_MERGED/"
INDIR_REFGENOME="/scratch/mugolini/GRCz11_105_GENOME/"
GENOME="${INDIR_REFGENOME}Danio_rerio.GRCz11.105.dna.primary_assembly.fa"
OUTDIR_HISAT="/scratch/mugolini/eSLAM/HISAT3N_Output/"

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist_IAA_minus.txt)
echo "Processing $FILE"

module load gcc
module load samtools
module load varscan

echo "Varscan2 call"
samtools mpileup -f ${GENOME} ${OUTDIR_HISAT}${FILE}_aligned.unique.sorted.bam | varscan pileup2snp --min-coverage 20 --min-reads2 5 --min-avg-qual 15 --min-var-freq  0.25 --p-value 0.01 > ${FILE}.varscan.snp
echo "Script finished"

# AFTER THIS SCRIPT HAS BEEN RUN ON ALL 4 SAMPLES, RUN THE FOLLOWING COMMANDS MANUALLY. IT WILL CREATE A UNION OF ALL SNPs USING AN R SCRIPT:
# module load gcc
# module load r
# Rscript --vanilla SNP_union.R 04_MU_L8.varscan.snp 05_MU_L8.varscan.snp

# The final output is the SNP_output.snp file and a Venn diagram (.pdf). The SNP_output.snp file will be used for the next analysis steps.