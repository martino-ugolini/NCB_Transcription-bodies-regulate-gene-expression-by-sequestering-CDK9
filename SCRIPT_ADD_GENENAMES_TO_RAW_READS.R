rm(list=ls())

# THIS SCRIPT TAKES THE RAW COUNTS AS GENERATED BY FEATURECOUNTS (input_raw_counts) AND ADDS GENE NAMES TO THEM FOR FUTURE USE.

#----------------------------------------------------------
# LOAD LIBRARIES
#----------------------------------------------------------
library(stringr)
library(rtracklayer)

#----------------------------------------------------------
# LOAD FUNDAMENTAL DATASETS
#----------------------------------------------------------
gtf_input <- "./Danio_rerio.GRCz11.105.gtf.gz"
GTF <- as.data.frame(import(gtf_input)); GTF <- GTF[which(GTF$type=="gene"),]

input_raw_counts <- "raw_counts_proper_orientation_filtered.txt"
output_raw_counts <- str_remove(input_raw_counts, ".txt"); output_raw_counts <- paste0(output_raw_counts, "_with_genename.txt")

#----------------------------------------------------------
# SCRIPT
#----------------------------------------------------------

INPUT_DATA <- read.table(input_raw_counts, stringsAsFactors = TRUE, row.names = 1, header=TRUE)
colnames(INPUT_DATA) <- str_remove(colnames(INPUT_DATA), paste0("X", gsub("/",".",input_raw_counts_folder)))
colnames(INPUT_DATA) <- str_remove(colnames(INPUT_DATA), "aligned.unique.")
colnames(INPUT_DATA) <- str_remove(colnames(INPUT_DATA), ".sorted.bam")
INPUT_DATA <- INPUT_DATA[,c(colnames(INPUT_DATA)[1:5], sort(colnames(INPUT_DATA)[-c(1:5)]))]


nrow(INPUT_DATA)==nrow(GTF)
gene_names <- c(); check <- c(); i <- 0
for (id in rownames(INPUT_DATA)){
  n <- GTF$gene_name[which(GTF$gene_id==id)]
  if (length(n)==1) {check <- c(check,"ok")}
  else {check <- c(check,"ERROR")}
  gene_names <- c(gene_names, n)
  i <- i+1
  if ((i %% 1000) == 0) {print(i)}}
"ERROR" %in% check

INPUT_DATA <- cbind(gene_names, gene_ID = rownames(INPUT_DATA), INPUT_DATA)

write.table(INPUT_DATA, file=output_raw_counts, sep="\t")

