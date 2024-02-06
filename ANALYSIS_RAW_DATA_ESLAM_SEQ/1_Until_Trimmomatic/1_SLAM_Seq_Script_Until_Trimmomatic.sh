#!/bin/bash
#SBATCH --job-name eSLAM_until_Trimmomatic
#SBATCH --time 12:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 48
#SBATCH --output eSLAM_until_Trimmomatic.%A_%a.log
#SBATCH --array 1-28

echo "PART 1 of eSLAM-SEQ Analysis"
echo "This script will perform pre-filtering (before mapping)"
echo "The raw .fastq.gz files need to be saved in the INDIR_RAW directory"
echo "The fw and rv files must have the same name and end with _R1.fastq.gz and sample_R2.fastq.gz"
echo "Prepare a filelist.txt file in INDIR_RAW with all the file names (without the _R1.fastq.gz, _R2.fastq.gz)"

INDIR_RAW="/scratch/mugolini/eSLAM/RAW_DATA_MERGED/"
OUTDIR_FASTQC_1="/scratch/mugolini/eSLAM/FastQC_Output/"
OUTDIR_TRIMMOMATIC="/scratch/mugolini/eSLAM/Trimmomatic_Output/"
ADAPTER_TRIMMOMATIC="/scratch/mugolini/eSLAM/Trimmomatic_adapters/"
OUTDIR_FASTQC_2="/scratch/mugolini/eSLAM/FastQC_Output_2/"

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist.txt)
echo "Running analysis on $FILE"

module load gcc
module load fastqc
module load trimmomatic
module load python
module load samtools
module load picard

chmod +x fastqc.sh
chmod +x trimmomatic.sh

echo "Version of FastQC used:"
fastqc --version
echo "Version of Trimmomatic used:"
trimmomatic -version

echo "FastQC 1 started"
./fastqc.sh ${INDIR_RAW} ${OUTDIR_FASTQC_1} ${SLURM_CPUS_PER_TASK} ${FILE}
echo "Trimmomatic started"
./trimmomatic.sh ${INDIR_RAW} ${OUTDIR_TRIMMOMATIC} ${SLURM_CPUS_PER_TASK} ${FILE} ${ADAPTER_TRIMMOMATIC}
echo "FastQC 2 started"
./fastqc.sh ${OUTDIR_TRIMMOMATIC} ${OUTDIR_FASTQC_2} ${SLURM_CPUS_PER_TASK} ${FILE}

echo "unzip fastq.gz files (Trimmomatic output)"
gunzip ${OUTDIR_TRIMMOMATIC}${FILE}_R1.fastq.gz
gunzip ${OUTDIR_TRIMMOMATIC}${FILE}_R2.fastq.gz

echo "Script finished"
