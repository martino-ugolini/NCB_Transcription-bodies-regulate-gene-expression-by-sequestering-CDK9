#!/bin/bash
#SBATCH --job-name ATAC_Seq
#SBATCH --time 1:00:00
#SBATCH --nodes 1
#SBATCH --cpus-per-task 10
#SBATCH --output ATAC_Seq.%A_%a.log

module load gcc

echo "ATAC-Seq Analysis"

/users/mugolini/.local/bin/computeMatrix reference-point --referencePoint center --regionsFileName 2000_upregulated.bed 2000_top716.bed \
--scoreFileName ./bigwig_replicates_merged/256c.sortedbycoordinates.RPKM.bw \
-bs 20 -a 2000 -b 2000 -p 10 --missingDataAsZero --sortRegions descend --sortUsing mean --outFileName deeptools_computeMatrix.gz --outFileNameMatrix deeptools_computeMatrix_outFile
/users/mugolini/.local/bin/plotHeatmap -m deeptools_computeMatrix.gz -o Deeptools_matrix_1_only_256.pdf  --heatmapWidth 2 --heatmapHeight 8 \
--colorMap 'Greys'

echo "Script finished"
