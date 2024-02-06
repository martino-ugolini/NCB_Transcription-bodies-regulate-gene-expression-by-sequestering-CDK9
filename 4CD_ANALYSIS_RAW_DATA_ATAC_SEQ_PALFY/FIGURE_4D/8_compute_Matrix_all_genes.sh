#!/bin/bash
#SBATCH --job-name ATAC_Seq
#SBATCH --time 1:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 10
#SBATCH --output ATAC_Seq.%A_%a.log

module load gcc

echo "ATAC-Seq Analysis"

/users/mugolini/.local/bin/computeMatrix reference-point --referencePoint center --regionsFileName 2000_all_genes.bed \
--scoreFileName ./bigwig_replicates_merged/256c.sortedbycoordinates.RPKM.bw \
-bs 20 -a 2000 -b 2000 -p 10 --missingDataAsZero --sortRegions keep --outFileName deeptools_computeMatrix_all_genes_256c.gz --outFileNameMatrix deeptools_computeMatrix_outFile_all_genes_256c

/users/mugolini/.local/bin/computeMatrix reference-point --referencePoint center --regionsFileName 2000_all_genes.bed \
--scoreFileName ./bigwig_replicates_merged/High.sortedbycoordinates.RPKM.bw \
-bs 20 -a 2000 -b 2000 -p 10 --missingDataAsZero --sortRegions keep --outFileName deeptools_computeMatrix_all_genes_High.gz --outFileNameMatrix deeptools_computeMatrix_outFile_all_genes_High

/users/mugolini/.local/bin/computeMatrix reference-point --referencePoint center --regionsFileName 2000_all_genes.bed \
--scoreFileName ./bigwig_replicates_merged/Oblong.sortedbycoordinates.RPKM.bw \
-bs 20 -a 2000 -b 2000 -p 10 --missingDataAsZero --sortRegions keep --outFileName deeptools_computeMatrix_all_genes_Oblong.gz --outFileNameMatrix deeptools_computeMatrix_outFile_all_genes_Oblong

/users/mugolini/.local/bin/computeMatrix reference-point --referencePoint center --regionsFileName 2000_all_genes.bed \
--scoreFileName ./bigwig_replicates_merged/Sphere.sortedbycoordinates.RPKM.bw \
-bs 20 -a 2000 -b 2000 -p 10 --missingDataAsZero --sortRegions keep --outFileName deeptools_computeMatrix_all_genes_Sphere.gz --outFileNameMatrix deeptools_computeMatrix_outFile_all_genes_Sphere

echo "Script finished"
