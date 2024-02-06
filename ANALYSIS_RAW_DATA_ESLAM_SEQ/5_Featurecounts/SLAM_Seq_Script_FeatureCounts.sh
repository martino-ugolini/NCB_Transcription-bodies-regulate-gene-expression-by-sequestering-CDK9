#!/bin/bash
#SBATCH --job-name eSLAM_FeatureCounts
#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --output eSLAM_FeatureCounts.log

echo "PART 5 of eSLAM-SEQ Analysis"
echo "This script will calculate raw counts using FeatureCounts"

INDIR_RAW="/scratch/mugolini/eSLAM/RAW_DATA_MERGED/"
INDIR_REFGENOME="/scratch/mugolini/GRCz11_105_GENOME/"
ANNOTATION="${INDIR_REFGENOME}Danio_rerio.GRCz11.105.gtf.gz"
INDIR_HISAT3N="/users/mugolini/hisat-3n/"
OUTDIR_HISAT="/scratch/mugolini/eSLAM/HISAT3N_Output/"
PATH_FEATURECOUNTS="/users/mugolini/subread-2.0.3-Linux-x86_64/bin/"

module load gcc
module load r

echo "list of .bam files to consider (filtered and not filtered)"
bamlist=$(find ${OUTDIR_HISAT} -maxdepth 1 -type f -name "*.sorted.bam")
echo $bamlist

${PATH_FEATURECOUNTS}./featureCounts -p --countReadPairs -B --minOverlap 10 -s 2 -T $SLURM_CPUS_PER_TASK -t exon -g gene_id \
-a ${ANNOTATION} -o raw_counts_proper_orientation_filtered_and_unfiltered.txt ${bamlist}

echo "list of .bam files to consider (filtered)"
bamlist=$(find ${OUTDIR_HISAT} -maxdepth 1 -type f -name "*_aligned.unique.filtered.sorted.bam")
echo $bamlist

${PATH_FEATURECOUNTS}./featureCounts -p --countReadPairs -B --minOverlap 10 -s 2 -T $SLURM_CPUS_PER_TASK -t exon -g gene_id \
-a ${ANNOTATION} -o raw_counts_proper_orientation_filtered.txt ${bamlist}

echo "Script finished"
