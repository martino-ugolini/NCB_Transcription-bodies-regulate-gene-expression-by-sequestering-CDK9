#!/bin/bash
#SBATCH --job-name eSLAM_HISAT3N_mapping
#SBATCH --time 12:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 48
#SBATCH --output eSLAM_HISAT3N_mapping.%A_%a.log
#SBATCH --array 1-28

echo "PART 2 of eSLAM-SEQ Analysis"
echo "This script will align paired-end RNA-Seq reads after the previous filtering steps"
echo "The raw .fastq.gz files need to be saved in the INDIR_RAW directory"
echo "The fw and rv files must have the same name and end with _R1.fastq.gz and sample_R2.fastq.gz"
echo "Prepare a filelist.txt file in INDIR_RAW with all the file names (without the _R1.fastq.gz, _R2.fastq.gz)"

INDIR_RAW="/scratch/mugolini/eSLAM/RAW_DATA_MERGED/"
OUTDIR_TRIMMOMATIC="/scratch/mugolini/eSLAM/Trimmomatic_Output/"
INDIR_HISAT3N="/users/mugolini/hisat-3n/"
INDIR_INDEX="/scratch/mugolini/eSLAM/GENOME_INDEX/Genome_3n_TC/genome_3n_TC"
OUTDIR_HISAT="/scratch/mugolini/eSLAM/HISAT3N_Output/"
mkdir -p $OUTDIR_HISAT

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${INDIR_RAW}filelist.txt)
echo "Running analysis on $FILE"

module load gcc
module load python
module load samtools
module load picard

echo "Mapping started"
${INDIR_HISAT3N}hisat-3n --base-change T,C --repeat -p $SLURM_CPUS_PER_TASK --rna-strandness RF --fr --no-discordant \
--summary-file ${OUTDIR_TRIMMOMATIC}${FILE}_summary \
-x ${INDIR_INDEX} -1 ${OUTDIR_TRIMMOMATIC}${FILE}_R1.fastq -2 ${OUTDIR_TRIMMOMATIC}${FILE}_R2.fastq -S ${OUTDIR_HISAT}${FILE}_aligned.sam
echo "Mapping completed succesfully"

echo "rezip fastq files"
gzip ${OUTDIR_TRIMMOMATIC}${FILE}_R1.fastq
gzip ${OUTDIR_TRIMMOMATIC}${FILE}_R2.fastq

echo "Keeping only uniquely mapped pairs"
samtools view -h -b -f 0x2 -F 0x4 -F 0x100 --threads $SLURM_CPUS_PER_TASK --input-fmt-option 'filter=[NH] == 1' ${OUTDIR_HISAT}${FILE}_aligned.sam > ${OUTDIR_HISAT}${FILE}_aligned.unique.bam

echo "Fixmate"
picard FixMateInformation --INPUT ${OUTDIR_HISAT}${FILE}_aligned.unique.bam --OUTPUT ${OUTDIR_HISAT}${FILE}_aligned.unique.fixmate.bam --IGNORE_MISSING_MATES FALSE

echo "Sorting file by coordinate"
samtools sort --threads $SLURM_CPUS_PER_TASK -o ${OUTDIR_HISAT}${FILE}_aligned.unique.sorted.bam ${OUTDIR_HISAT}${FILE}_aligned.unique.fixmate.bam
samtools index -@ $SLURM_CPUS_PER_TASK ${OUTDIR_HISAT}${FILE}_aligned.unique.sorted.bam
samtools flagstats -@ $SLURM_CPUS_PER_TASK ${OUTDIR_HISAT}${FILE}_aligned.unique.sorted.bam

echo "Script finished"
