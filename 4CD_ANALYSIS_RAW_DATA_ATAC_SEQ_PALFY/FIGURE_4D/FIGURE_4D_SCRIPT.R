
library(ggplot2)
library(wesanderson)
library("RColorBrewer")
library(ggbio)
library(rtracklayer)
library(GenomicRanges)

# Before performing the analysis, you need to use deeptools to compute the computeMatrix and computeMatrix_out files for the promoters of all annotated genes.
# For that, you first create the all_genes.bed file by using the SCRIPT_FIGURES_PAPER.R. Just define
#   ALL_GENES <- unique(GTF$gene_id)
# and then expand the code by substituting these lines
#   gene_list <- list(A = genes_up, B = not_expressed, C = top_716_genes, D = ALL_GENES)
#   output <- c("upregulated", "not_expressed", "top716", "all_genes")
# The .bed file is already provided because it takes very long to generate (2000_all_genes.bed)
# Use this .bed file to compute the computeMatrix and computeMatrix_out using the script 8_compute_Matrix_all_genes.sh.
# The script will create these files for all 4 stages of ATAC-Seq (256c, High, Oblong and Sphere), but in the script we will only use the stage 256c
# NOTE: File headers need to be corrected (add # to comment in first file and remove first col name (gene...) in the second file!

#----------------------------------------------------------
# LOAD FUNDAMENTAL DATASETS
#----------------------------------------------------------

setwd("~/switchdrive/RESEARCH/2019_PHD_PROJECT/_MANUSCRIPTS/Paper/R_SCRIPT_REPRODUCE_FIGURES")

input_data <- "./raw_counts_proper_orientation_filtered_with_genename.txt"
sample_info <- "./sample_infos.csv"
gtf_input <- "./Danio_rerio.GRCz11.105.gtf.gz"
chrom_length_file <- "./Danio_rerio.GRCz11.105.dna.primary_assembly.fa.fai" #To obtain the chromosome length, gunzip the reference genome and run samtools faidx Danio_rerio.GRCz11.105.dna.primary_assembly.fa. The obtained .fai file contains the chromosome length in the second column.


INPUT_DATA <- read.table(input_data, stringsAsFactors = TRUE, row.names = 1, header=TRUE)
INPUT_DESIGN <- read.table(sample_info, sep=',', header=TRUE, row.names = 1)
DESIGN <- INPUT_DESIGN[which(INPUT_DESIGN$IAA_processing == "yes"),]

GTF <- as.data.frame(import(gtf_input)); GTF <- GTF[which(GTF$type=="gene"),]

MT_genes <- GTF$gene_id[which(GTF$seqnames=="MT")]
miR430_genes <- GTF$gene_id[which(startsWith(as.character(GTF$gene_name), "dre-mir-430"))]




#----------------------------------------------------------
# ATAC SEQ ANALYSIS
#----------------------------------------------------------

# Choose the stage
STAGE <- c("256c", "High", "Oblong", "Sphere")[1]

# Load files
TABLE <- read.table(paste0("./ANALYSIS_RAW_DATA_ATAC_SEQ_PALFY_ET_AL_FIGURE_3CD/FIGURE_3D/deeptools_computeMatrix_all_genes_", STAGE))
TABLE_NONAME <- read.table(paste0("./ANALYSIS_RAW_DATA_ATAC_SEQ_PALFY_ET_AL_FIGURE_3CD/FIGURE_3D/deeptools_computeMatrix_outFile_all_genes_", STAGE), sep = '\t', comment.char = "#", header = FALSE)

# Determine columns names
promoter_length = 2000
bin = 20
intervals <- 2*promoter_length/bin
samples <- unlist(unique(c(TABLE_NONAME[1,])))
colnames(TABLE_NONAME) <- paste0(TABLE_NONAME[1,], "_", 1:intervals)
TABLE_NONAME <- TABLE_NONAME[-1,]

# Add columns names to TABLE
ncol(TABLE)
colnames(TABLE) <- c("seqnames", "start", "end", "promoter_id", "width", "strand", colnames(TABLE_NONAME))
#promote_id is the same as s in the main figures script

TABLE_MEAN <- TABLE[,1:6]
TABLE_RANK <- TABLE[,1:6]
for (i in 1:length(samples))
{
  print(samples[i])
  st <- 7+200*(i-1)
  en <- 6+200*(i)
  TABLE_sample <- TABLE[,st:en]
  mean <- rowMeans(TABLE_sample)
  TABLE_MEAN <- cbind(TABLE_MEAN, mean)
  rank <- rank(-mean, ties.method = "random")
  TABLE_RANK <- cbind(TABLE_RANK, rank)
}
TABLE_MEAN <- data.frame(TABLE_MEAN)
colnames(TABLE_MEAN) <- c(colnames(TABLE[,1:6]) ,samples)
TABLE_RANK <- data.frame(TABLE_RANK)
colnames(TABLE_RANK) <- c(colnames(TABLE[,1:6]) ,samples)

# UPREGULATED GENES
TABLE <- read.table("./DEG/256c_miR430_vs_WT.txt", sep = "\t")
genes_up   <- TABLE$gene_ID[which(TABLE$DE == "up-regulated" & TABLE$miR430 == FALSE & TABLE$mitochondrial == FALSE)]
genes_down   <- TABLE$gene_ID[which(TABLE$DE == "down-regulated" & TABLE$miR430 == FALSE & TABLE$mitochondrial == FALSE)]
genes_not_DE <- TABLE$gene_ID[which(TABLE$DE == "not DE" & TABLE$miR430 == FALSE & TABLE$mitochondrial == FALSE & abs(TABLE$log2FoldChange) < 1)]
genes_top716 <- read.table("top_716_genes.txt")$x

length(genes_up)
length(genes_down)
length(genes_not_DE)
length(genes_top716)

# Not expressed
DATA <- INPUT_DATA[which(!(INPUT_DATA$Chr == "MT")),]
DATA <- as.matrix(DATA[,-(1:7)])
all(rownames(DESIGN) %in% colnames(DATA)) # should be TRUE
DATA <- DATA[,rownames(DESIGN)]
DATA <- DATA[rowSums(DATA) < 10,]
genes_not_expressed <- rownames(DATA)
length(genes_not_expressed)


#Determine the promoters corresponding to the 2 lists:
#V4: promoter ID; V11: gene_id
PROMOTERS_ALL <- read.table("./BED_FASTA_FILES/2000_all_genes.bed")
PROMOTERS_UP <- PROMOTERS_ALL$V4[which(PROMOTERS_ALL$V11 %in% genes_up)]
PROMOTERS_DOWN <- PROMOTERS_ALL$V4[which(PROMOTERS_ALL$V11 %in% genes_down)]
PROMOTERS_NOT_DE <- PROMOTERS_ALL$V4[which(PROMOTERS_ALL$V11 %in% genes_not_DE)]
PROMOTERS_NOT_EXP <- PROMOTERS_ALL$V4[which(PROMOTERS_ALL$V11 %in% genes_not_expressed)]
PROMOTERS_TOP716 <- PROMOTERS_ALL$V4[which(PROMOTERS_ALL$V11 %in% genes_top716)]
#PROMOTERS_ALL <- PROMOTERS_ALL$V4

# LINE PLOT
BOXPLOT_RANK <- data.frame()
max_group <- c(length(PROMOTERS_UP),length(PROMOTERS_DOWN),length(PROMOTERS_NOT_DE),length(PROMOTERS_NOT_EXP), length(PROMOTERS_TOP716))
limits <- round(seq(0, 100, by = 10)*nrow(TABLE_RANK)/100, 0) # I will calculate the percentage rank in 10% steps

for (limit in limits)
{
  for (sample in colnames(TABLE_RANK)[-(1:6)])
  {
    TABLE_RANK_LIM <-  TABLE_RANK[which(TABLE_RANK[,sample] %in% 1:limit),]
    up <- TABLE_RANK_LIM[which(TABLE_RANK_LIM$promoter_id %in% PROMOTERS_UP),sample]
    down <- TABLE_RANK_LIM[which(TABLE_RANK_LIM$promoter_id %in% PROMOTERS_DOWN),sample]
    notde <- TABLE_RANK_LIM[which(TABLE_RANK_LIM$promoter_id %in% PROMOTERS_NOT_DE),sample]
    notexp <- TABLE_RANK_LIM[which(TABLE_RANK_LIM$promoter_id %in% PROMOTERS_NOT_EXP),sample]
    top716 <- TABLE_RANK_LIM[which(TABLE_RANK_LIM$promoter_id %in% PROMOTERS_TOP716),sample]
    
    Value <- c(length(up), length(down), length(notde), length(notexp), length(top716))
    Value_Percentage <- 100*(Value/nrow(TABLE_RANK_LIM))
    Value_Percentage2 <- 100*(Value/max_group)
    
    BOXPLOT_RANK <- rbind(BOXPLOT_RANK, cbind(Value = Value, Value_Percentage = Value_Percentage, Value_Percentage2 = Value_Percentage2, Group = c("Upregulated", "Downregulated", "Not differentially expressed", "Not expressed", "Top716"), Sample = sample, Limit = limit))
  }
}
BOXPLOT_RANK$Value <- as.double(BOXPLOT_RANK$Value)
BOXPLOT_RANK$Value_Percentage <- as.double(BOXPLOT_RANK$Value_Percentage)
BOXPLOT_RANK$Value_Percentage2 <- as.double(BOXPLOT_RANK$Value_Percentage2)
BOXPLOT_RANK$Limit <- as.double(BOXPLOT_RANK$Limit)
BOXPLOT_RANK$Limit_percentage <- 100*(BOXPLOT_RANK$Limit/nrow(TABLE_RANK))
BOXPLOT_RANK$Stage <- sapply(strsplit(BOXPLOT_RANK$Sample, "_"), `[`, 2)
BOXPLOT_RANK$Group <- factor(BOXPLOT_RANK$Group, levels = c("Upregulated", "Downregulated", "Not differentially expressed", "Not expressed", "Top716"))

p <- ggplot(data=BOXPLOT_RANK, aes(x=Limit_percentage, y=Value_Percentage2, group=Group, fill = Group, color = Group)) + geom_line() +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") + theme_bw() + 
  xlab("Top x percentage of most accessible genes") + ylab("Percentage of genes") + theme(aspect.ratio=1) + 
  ggtitle(STAGE) + theme(aspect.ratio = 1) + scale_color_manual(values=c("#377EB8", "#E41A1C", "#A65628", "#A0A0A0", "#984EA3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

pdf(file=paste0(paste0("FIGURE_4D_line_plot_percentage_group_vs_percentage_all_", STAGE, ".pdf")),width=6,height=4);print(p);dev.off()






