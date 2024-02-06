rm(list=ls())
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library(ggplot2)
library(stringr)
library(rtracklayer)
library(gtools)
library(ggpubr)
library("gridExtra")
library(ashr)
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library("dplyr")
library(wesanderson)


#----------------------------------------------------------
# LOAD FUNDAMENTAL DATASETS
#----------------------------------------------------------

input_data <- "./REVISION_ESLAM_SEQ/raw_counts_proper_orientation_filtered_with_genename.txt"
sample_info <- "./REVISION_ESLAM_SEQ/sample_infos.csv"
gtf_input <- "./Danio_rerio.GRCz11.105.gtf.gz"
chrom_length_file <- "./Danio_rerio.GRCz11.105.dna.primary_assembly.fa.fai" #To obtain the chromosome length, gunzip the reference genome and run samtools faidx Danio_rerio.GRCz11.105.dna.primary_assembly.fa. The obtained .fai file contains the chromosome length in the second column.

INPUT_DATA <- read.table(input_data, stringsAsFactors = TRUE, row.names = 1, header=TRUE)
INPUT_DESIGN <- read.table(sample_info, sep=',', header=TRUE, row.names = 1)
GTF <- as.data.frame(import(gtf_input)); GTF <- GTF[which(GTF$type=="gene"),]

MT_genes <- GTF$gene_id[which(GTF$seqnames=="MT")]
miR430_genes <- GTF$gene_id[which(startsWith(as.character(GTF$gene_name), "dre-mir-430"))]

#----------------------------------------------------------
# SET PARAMETERS DESEQ2
#----------------------------------------------------------
alpha = 0.01
log2FC_Treshold <- 2
cutoff <- 10


#----------------------------------------------------------
# PERFORM DESEQ2 ANALYSIS
#----------------------------------------------------------
DESIGN <- INPUT_DESIGN[which(INPUT_DESIGN$IAA_processing == "yes"),]
DESIGN <- cbind(DESIGN, Group = paste0(DESIGN$Stage, "_", DESIGN$Genotype))
DATA <- as.matrix(INPUT_DATA[,-(1:7)])
all(rownames(DESIGN) %in% colnames(DATA)) # should be TRUE
DATA<-DATA[,rownames(DESIGN)]
nrow(DESIGN) == ncol(DATA) # should be TRUE

FORMULA <-  as.formula("~ Group")
dds <- DESeqDataSetFromMatrix(countData=DATA, colData=DESIGN, design=FORMULA)
dds <- dds[rowSums(counts(dds))>=cutoff,]
dds$Stage <- factor(dds$Stage, levels = c("256c", "1024c"))
dds <- DESeq(dds)


#----------------------------------------------------------
# SELECT THE STAGE FOR THE ANALYSIS
#----------------------------------------------------------

condition <- c("WT_mir430", "WT_cdk9")[2] #Choose which comparison to analyze

remove_CDK9 <- TRUE #cdk9 is highly enriched in the wt vs cdk9 comparison, and should therefore not be shown. It is most likely the injected cdk9 mRNA that is showing

for (stage in c("256c", "1024c"))
{
  filename <- paste0(stage,"_WT_vs_", condition)
  
  #----------------------------------------------------------
  # FIGURE 6E AND ED3G: VOLCANO PLOT AT THE SELECTED STAGE
  #----------------------------------------------------------
  
  # Determine DEGs (WT vs mir430-/- + Inj. MiR430 at the selected stage)
  res <- lfcShrink(dds, contrast = c("Group", paste0(stage,"_", condition), paste0(stage,"_WT")), type = "ashr")
  summary(res)
  
  res_filtering <- as.data.frame(res)
  genenames <- INPUT_DATA[which(INPUT_DATA$gene_ID %in% rownames(res_filtering)),c("gene_ID", "gene_names")]
  res_filtering <- res_filtering[order(rownames(res_filtering)),]
  genenames <- genenames[order(rownames(genenames)),]
  all(genenames$gene_ID == rownames(res_filtering))
  res_filtering <- cbind(genenames, res_filtering)
  res_filtering <- res_filtering[ order(res_filtering$log2FoldChange), ]
  
  res_filtering <- cbind(res_filtering, miR430 = (res_filtering$gene_ID %in% miR430_genes))
  res_filtering <- cbind(res_filtering, mitochondrial = (res_filtering$gene_ID %in% MT_genes))
  res_filtering <- cbind(res_filtering, DE="not DE")
  res_filtering$DE[which(res_filtering$padj < alpha & res_filtering$log2FoldChange<(-log2FC_Treshold))] <- "down-regulated"
  res_filtering$DE[which(res_filtering$padj < alpha & res_filtering$log2FoldChange>(log2FC_Treshold))] <- "up-regulated"
  res_filtering <- res_filtering[which(res_filtering$miR430 == FALSE & res_filtering$mitochondrial == FALSE),]
  
  if (condition == "WT_cdk9" & remove_CDK9 == TRUE){
    res_filtering <- res_filtering[which(!(res_filtering$gene_ID == "ENSDARG00000044811")),]
    filename <- paste0("Figure_6E_ED3G_No_CDK9_", filename)
  }
  
  color_values <- c("#E41A1C", "#A0A0A0", "#377EB8")
  if (stage == "1024c" & condition == "WT_mir430"){color_values <- c("#A0A0A0", "#377EB8")}
  
  # Volcano Plot
  res_filtering <- res_filtering[which(!(is.na(res_filtering$padj))),]
  p2 <- ggplot(res_filtering, aes(x=log2FoldChange, y=-log10(padj), color=DE)) + geom_point(size = 0.5) + theme_bw() + xlab("log2FC") + ylab("-log10(padj)") + theme(legend.position = "none")
  p2 <- p2 + geom_vline(xintercept = c(-log2FC_Treshold, log2FC_Treshold)) + geom_hline(yintercept=-log10(alpha)) + scale_color_manual(values=color_values)
  p2 <- p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio=1)

  # Numbers to manually add to the top corners of the volcanoplot
  nrow(res_filtering[which(res_filtering$miR430 == FALSE & res_filtering$mitochondrial == FALSE & res_filtering$DE == "down-regulated"),])
  nrow(res_filtering[which(res_filtering$miR430 == FALSE & res_filtering$mitochondrial == FALSE & res_filtering$DE == "up-regulated"),])
  
  pdf(file=paste0("./REVISION_ESLAM_SEQ/",condition, "/", filename, "_Volcano_plot.pdf"),width=4,height=4); print(p2); dev.off()
  

  #----------------------------------------------------------
  # WRITE DEG FILE FOR LATER USE
  #----------------------------------------------------------
  write.table(res_filtering, paste0("./REVISION_ESLAM_SEQ/",condition, "/",filename,".txt"), quote = FALSE, sep = "\t")
  
  #----------------------------------------------------------
  # KARIOGRAM AT THE SELECTED STAGE
  #----------------------------------------------------------
  
  with_miR430 <- TRUE
  de <- c("down-regulated", "up-regulated")
  color <- c("#E41A1C", "#377EB8") #Chosen colors for down (first) and up (second) genes
  miR430 <- ""
  
  for (i in 1:2) #1: down-regulated, 2: up-regulated
  {
    gene_list <- rownames(res_filtering[which(res_filtering$miR430 == FALSE & res_filtering$mitochondrial == FALSE & res_filtering$DE == de[i]),])
    GTF_sub <- GTF[which(GTF$gene_id %in% gene_list & GTF$seqnames %in% 1:25),] #Keep only the genes that are not on scaffolds
    
    GTF_sub$seqnames <- factor(GTF_sub$seqnames, levels = 1:25)
    length(gene_list) - nrow(GTF_sub) #Should be zero
    gene_list_1_25 <- unique(GTF_sub$gene_id)
    
    # BARPLOT
    DISTRIBUTION <- as.data.frame(table(GTF_sub$seqnames)); colnames(DISTRIBUTION) <- c("seqname", "x")
    DISTRIBUTION$M <- sum(DISTRIBUTION$x)
    DISTRIBUTION$n <- as.data.frame(table(GTF$seqnames))$Freq[1:25]
    DISTRIBUTION$N <- nrow(GTF[which(GTF$seqnames %in% 1:25),])
    DISTRIBUTION$Expected = (DISTRIBUTION$n / DISTRIBUTION$N) * DISTRIBUTION$M
    DISTRIBUTION$EnrichmentRatio <- DISTRIBUTION$x / DISTRIBUTION$Expected
    
    all(DISTRIBUTION$N == sum(DISTRIBUTION$n))
    all(DISTRIBUTION$M == sum(DISTRIBUTION$x))
    
    p_values <- c()
    for (j in 1:nrow(DISTRIBUTION)){p_values <- c(p_values, phyper(q = DISTRIBUTION$x[j] - 1, m = DISTRIBUTION$n[j], n = DISTRIBUTION$N[j] - DISTRIBUTION$n[j], k = DISTRIBUTION$M[j], lower.tail = FALSE))}
    DISTRIBUTION$p_value <- p_values
    DISTRIBUTION$adjusted_p_value <- p.adjust(DISTRIBUTION$p_value, method="BH")
    DISTRIBUTION$Significant <- "n.s."
    DISTRIBUTION$Significant[which(DISTRIBUTION$adjusted_p_value < 0.05)] <- "sig."
    DISTRIBUTION$seqname <- factor(DISTRIBUTION$seqname, levels = rev(1:25))
    
    p <- ggplot(data=DISTRIBUTION, aes(x=seqname, y=EnrichmentRatio, fill = Significant)) + geom_bar(stat="identity") + theme_bw() + coord_flip() + scale_fill_manual(values=c("#A0A0A0", "#377EB8"))
    pdf(file=paste0("./REVISION_ESLAM_SEQ/",condition, "/", filename, "_", de[i], "_gene_distribution_paper_barplot.pdf"),width=7,height=6); print(p); dev.off()
    
    if (with_miR430 == TRUE) #Add miR430 as a gene to show on the kariogram. It will have the same color as all other genes, so you will have to find it in Illustrator later on an change its color (green in the publication)
    {
      GTF_miR430 <- GTF[which(GTF$gene_id %in% miR430_genes & GTF$seqnames == 4),]
      s <- min(GTF_miR430$start); e <- max(GTF_miR430$end); w <- e-s+1
      GTF_miR430$start[1] <- s; GTF_miR430$end[1] <- e; GTF_miR430$width <- w
      GTF_miR430 <- GTF_miR430[1,]
      miR430 <- "with_miR430"
      GTF_sub <- rbind(GTF_sub, GTF_miR430)
      GTF_sub$seqnames <- factor(GTF_sub$seqnames, levels = 1:25)
    }
    
    chrom_length <- read.table(chrom_length_file)[,c(1:2)]; colnames(chrom_length) <- c("seqnames", "length")
    chrom_length <- chrom_length[which(chrom_length$seqnames %in% 1:25),]
    chrom_length$seqnames <- factor(chrom_length$seqnames, levels = 1:25)
    chrom_length <- chrom_length[order(chrom_length$seqnames),]
    
    chr_zebrafish <- GRanges(seqnames = chrom_length$seqnames, strand = c("*"), ranges = IRanges(start=1, width=chrom_length$length))
    
    GTF_sub <- GRanges(seqnames = GTF_sub$seqnames, strand = c("*"),ranges <- IRanges(start=GTF_sub$start, end=GTF_sub$end))
    seqlengths(GTF_sub) <- chrom_length$length
    
    # KARIOGRAM
    color_value <- color[i]
    p <- autoplot(GTF_sub, layout="karyogram", aes(color="", fill="")) + scale_color_manual(values=color_value) + scale_fill_manual(values=color_value)
    pdf(file=paste0("./REVISION_ESLAM_SEQ/",condition, "/", filename, "_", de[i], "_gene_distribution_paper_", miR430, ".pdf"),width=7,height=6); print(p@ggplot); dev.off()
  }
}


beepr::beep()






