#!/bin/bash
#SBATCH --job-name ATAC_Seq
#SBATCH --time 1:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --output ATAC_Seq.%A_%a.log

OUTDIR_BOWTIE2="/scratch/mugolini/ATAC/Bowtie2_Output/"

module load gcc
module load samtools

cd ${OUTDIR_BOWTIE2}

echo "Merge files of biological replicates"
samtools merge --threads 10 -n 256c.bam SRR9032650_aligned.q30.noMT.dupremoved.sortedbyname.bam SRR9032651_aligned.q30.noMT.dupremoved.sortedbyname.bam
samtools merge --threads 10 -n High.bam SRR9032662_aligned.q30.noMT.dupremoved.sortedbyname.bam SRR9032663_aligned.q30.noMT.dupremoved.sortedbyname.bam SRR9032664_aligned.q30.noMT.dupremoved.sortedbyname.bam
samtools merge --threads 10 -n Oblong.bam SRR9032665_aligned.q30.noMT.dupremoved.sortedbyname.bam SRR9032666_aligned.q30.noMT.dupremoved.sortedbyname.bam SRR9032667_aligned.q30.noMT.dupremoved.sortedbyname.bam
samtools merge --threads 10 -n Sphere.bam SRR9032674_5_aligned.q30.noMT.dupremoved.sortedbyname.bam SRR9032676_aligned.q30.noMT.dupremoved.sortedbyname.bam

echo "Sorting by coordinates..."
samtools sort -@ 10 -o 256c.sortedbycoordinates.bam 256c.bam
samtools sort -@ 10 -o High.sortedbycoordinates.bam High.bam
samtools sort -@ 10 -o Oblong.sortedbycoordinates.bam Oblong.bam
samtools sort -@ 10 -o Sphere.sortedbycoordinates.bam Sphere.bam

echo "Make index file..."
samtools index -@ 10 256c.sortedbycoordinates.bam
samtools index -@ 10 High.sortedbycoordinates.bam
samtools index -@ 10 Oblong.sortedbycoordinates.bam
samtools index -@ 10 Sphere.sortedbycoordinates.bam

echo "Script finished"
