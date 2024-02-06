#!/bin/bash
#SBATCH --job-name ATAC_merge_sample
#SBATCH --time 01:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --output ATAC_merge_samples.log

INDIR_RAW="/scratch/mugolini/ATAC/RAW_DATA/"

FILE1="SRR9032674"
FILE2="SRR9032675"
FILEOUT="SRR9032674_5"

c1=$(zcat ${INDIR_RAW}${FILE1}_1.fastq.gz | echo $((`wc -l`/4)))
c2=$(zcat ${INDIR_RAW}${FILE1}_2.fastq.gz | echo $((`wc -l`/4)))
c3=$(zcat ${INDIR_RAW}${FILE2}_1.fastq.gz | echo $((`wc -l`/4)))
c4=$(zcat ${INDIR_RAW}${FILE2}_2.fastq.gz | echo $((`wc -l`/4)))
echo $c1;echo $c2;echo $c3;echo $c4

cat ${INDIR_RAW}${FILE1}_1.fastq.gz ${INDIR_RAW}${FILE2}_1.fastq.gz > ${OUTDIR_RAW}${FILEOUT}_1.fastq.gz
cat ${INDIR_RAW}${FILE1}_2.fastq.gz ${INDIR_RAW}${FILE2}_2.fastq.gz > ${OUTDIR_RAW}${FILEOUT}_2.fastq.gz

c1=$(zcat ${OUTDIR_RAW}${FILEOUT}_1.fastq.gz | echo $((`wc -l`/4)))
c2=$(zcat ${OUTDIR_RAW}${FILEOUT}_2.fastq.gz | echo $((`wc -l`/4)))
echo $c1;echo $c2

echo "Script finished"
