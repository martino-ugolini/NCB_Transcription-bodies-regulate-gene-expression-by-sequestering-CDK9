rm(list=ls())

#----------------------------------------------------------
# LOAD LIBRARIES
#----------------------------------------------------------

library(ggplot2)
library(rtracklayer)
library(wesanderson)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(gtools)
library(ggpubr)
library(gridExtra)
library(ashr)
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggalluvial)


#----------------------------------------------------------
# LOAD FUNDAMENTAL DATASETS
#----------------------------------------------------------

input_data <- "./raw_counts_proper_orientation_filtered_with_genename.txt"
sample_info <- "./sample_infos.csv"
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
# FIGURE S8C (Calculate percentage of MT and Nuclear transcripts)
#----------------------------------------------------------
DESIGN <- INPUT_DESIGN[which(INPUT_DESIGN$IAA_processing == "yes"),]
total_count <- colSums(INPUT_DATA[,-c(1:7)])
MT_all <- GTF$gene_id[which(GTF$seqnames=="MT")]
MT_counts_all <- as.data.frame(colSums(INPUT_DATA[MT_all,-c(1:7)]))
MT_counts_all[,1] <- MT_counts_all[,1]/total_count * 100

Nuclear_all <- GTF$gene_id[which(!(GTF$seqnames=="MT"))]
Nuclear_counts_all <- as.data.frame(colSums(INPUT_DATA[Nuclear_all,-c(1:7)]))
Nuclear_counts_all[,1] <- Nuclear_counts_all[,1]/total_count * 100

MT_counts_all[,1] + Nuclear_counts_all[,1] #The sum should always be 100

MT_counts_all  <- cbind(MT_counts_all, Group = "Mitochondrial reads", Sample = rownames(MT_counts_all)); colnames(MT_counts_all) <- c("Value", "Group", "Sample")
Nuclear_counts_all  <- cbind(Nuclear_counts_all, Group = "Nuclear reads", Sample = rownames(Nuclear_counts_all)); colnames(Nuclear_counts_all) <- c("Value", "Group", "Sample")

TABLE <- rbind(MT_counts_all, Nuclear_counts_all)
TABLE <- TABLE[which(TABLE$Sample %in% rownames(DESIGN)),]

s <- c(); g <- c()
for (sample in TABLE$Sample){s <- c(s, DESIGN[sample,"Stage"]); g <- c(g, DESIGN[sample,"Genotype"])}
TABLE <- cbind(TABLE, Stage = s, Genotype = g)
TABLE$Stage <- factor(TABLE$Stage, levels = c("256c", "1000c", "Oblong", "Sphere"))

p <- ggplot(data=TABLE, aes(x=Stage, y=Value, fill=Group)) +
  geom_boxplot() + theme_bw() + theme(legend.position = 'right') + ylab("Percentage of total reads") +
  scale_fill_manual(values = wes_palette(n=2, name="Moonrise2")) + theme(aspect.ratio=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf(file=paste0("./FIGURES/","Figure_ED3F_Barplot_mitochondrial_nuclear.pdf"),width=6,height=4); print(p); dev.off()




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
dds$Stage <- factor(dds$Stage, levels = c("256c", "1000c", "Oblong", "Sphere"))
dds <- DESeq(dds)

#----------------------------------------------------------
# FIGURE S8B: PRINCIPAL COMPONENT ANALYSIS (PCA)
#----------------------------------------------------------
rld <- rlog(dds, blind=FALSE)

pcaData <- plotPCA(rld, intgroup=c("Genotype", "Stage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p<-ggplot(pcaData, aes(PC1, PC2, fill=Stage, shape = Genotype)) +
  geom_point(size=4) + theme_bw() + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + scale_fill_manual(values = wes_palette(n=4, name="Royal1")) + theme(aspect.ratio=1) + scale_shape_manual(values=c(21,24)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file=paste0("./FIGURES/","Figure_ED3E_PCA_all_samples.pdf"),width=10,height=5); print(p); dev.off()


#----------------------------------------------------------
# SELECT THE STAGE FOR THE ANALYSIS (FIGURE 2A (VOLCANO PLOT), 2C and ED4 (KARIOGRAM AND BARPLOT QUANTIFICATION OF ENRICHMENT))
#----------------------------------------------------------

for (stage in c("256c", "1000c", "Oblong", "Sphere"))
{
  filename <- paste0(stage,"_miR430_vs_WT") #miR430 vs WT at stage
  
  #----------------------------------------------------------
  # FIGURE 2A: VOLCANO PLOT AT THE SELECTED STAGE
  #----------------------------------------------------------
  
  # Determine DEGs (WT vs mir430-/- + Inj. MiR430 at the selected stage)
  res <- lfcShrink(dds, contrast = c("Group", paste0(stage,"_miR430"), paste0(stage,"_WT")), type = "ashr")
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
  
  #----------------------------------------------------------
  # WRITE DEG FILE FOR LATER USE AND VOLCANO PLOT
  #----------------------------------------------------------
  write.table(res_filtering, paste0("./DEG/",filename,".txt"), quote = FALSE, sep = "\t")
  
  # Volcano Plot
  p2 <- ggplot(res_filtering[which(!(is.na(res_filtering$padj))),], aes(x=log2FoldChange, y=-log10(padj), color=DE)) + geom_point(size = 0.5) + theme_bw() + xlab("log2FC") + ylab("-log10(padj)") + theme(legend.position = "none")
  p2 <- p2 + geom_vline(xintercept = c(-log2FC_Treshold, log2FC_Treshold)) + geom_hline(yintercept=-log10(alpha)) + scale_color_manual(values=c("#E41A1C", "#A0A0A0", "#377EB8"))
  p2 <- p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio=1)
  #  p <- p + scale_x_continuous(breaks = seq(-6, 10, 4)) + scale_y_continuous(breaks = seq(0, 14, 4)) (Add this line for the graph at 256c stage)
  
  # Numbers to manually add to the top corners of the volcanoplot
  nrow(res_filtering[which(res_filtering$miR430 == FALSE & res_filtering$mitochondrial == FALSE & res_filtering$DE == "down-regulated"),])
  nrow(res_filtering[which(res_filtering$miR430 == FALSE & res_filtering$mitochondrial == FALSE & res_filtering$DE == "up-regulated"),])
  
  pdf(file=paste0("./FIGURES/" ,"Figure_2A_", filename, "_Volcano_plot.pdf"),width=4,height=4); print(p2); dev.off()
  

  #----------------------------------------------------------
  # FIGURE 2C: KARIOGRAM AT THE SELECTED STAGE
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
    pdf(file=paste0("./FIGURES/","Figure_ED4B_", filename, "_", de[i], "_gene_distribution_paper_barplot.pdf"),width=7,height=6); print(p); dev.off()
    
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
    pdf(file=paste0("./FIGURES/","Figure_2C_ED4A_", filename, "_", de[i], "_gene_distribution_paper_", miR430, ".pdf"),width=7,height=6); print(p@ggplot); dev.off()
  }
}


#----------------------------------------------------------
# FIGURE 2E: DIFFERENTIAL EXPRESSION RECOVERY (BOXPLOT)
#----------------------------------------------------------

# IDENTIFY THE GENE IDs OF THE GENES THAT ARE UP- AND DOWN-REGULATED AT 256c
stage <- "256c"
res <- read.table(paste0("./DEG/" ,stage,"_miR430_vs_WT",".txt"), sep = "\t")
res <- cbind(res, Stage = stage)
gene_list_down <- res$gene_ID[res$miR430 == FALSE & res$mitochondrial == FALSE & res$DE == "down-regulated"]
gene_list_up <-   res$gene_ID[res$miR430 == FALSE & res$mitochondrial == FALSE & res$DE == "up-regulated"]

length(gene_list_down)
length(gene_list_up)

# FOR THE GENES THAT WE JUST SELECTED, DETERMINE THEIR DEGREE OF MISREGULATION AT THE 4 STAGES
res_down <- data.frame()
res_up <- data.frame()

for (stage in c("256c", "1000c", "Oblong", "Sphere"))
{
  res <- read.table(paste0("./DEG/" ,stage,"_miR430_vs_WT",".txt"), sep = "\t")
  res <- cbind(res, Stage = stage)
  
  res_down <- rbind(res_down, res[which(res$gene_ID %in% gene_list_down),])
  res_up <- rbind(res_up, res[which(res$gene_ID %in% gene_list_up),])
}

res_down$Stage <- factor(res_down$Stage, levels = c("256c", "1000c", "Oblong", "Sphere"))
res_up$Stage <- factor(res_up$Stage, levels = c("256c", "1000c", "Oblong", "Sphere"))

res_down <- cbind(res_down, Condition = "down")
res_up <- cbind(res_up, Condition = "up")

RES <- rbind(res_down, res_up)

p <- ggplot(RES, aes(x=Stage, y=log2FoldChange, fill=Condition)) +
  geom_hline(yintercept=0) + geom_boxplot(position=position_dodge(1), outlier.size = 0.1) + theme_bw() + theme(legend.position = "right") + scale_fill_manual(values=c("#E41A1C", "#377EB8")) + theme(aspect.ratio=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept=0)
pdf(file=paste0("./FIGURES/","Figure_2E_Boxplot_paper.pdf"),width=6,height=4); print(p); dev.off()

#----------------------------------------------------------
# FIGURE 4A, 4B AND ED5: TIMING OF INDUCTION OF THE UP-REGULATED GENES
#----------------------------------------------------------

# DETERMINE WHEN GENES ARE INDUCED
stage_1 <- "256c"
for (stage_2 in c("1000c", "Oblong", "Sphere"))
{
  filename <- paste0(stage_2, "_vs_", stage_1, "_WT")

  res_x <- lfcShrink(dds, contrast = c("Group", paste0(stage_2,"_WT"), paste0(stage_1,"_WT")), type = "ashr")
  res_x <- as.data.frame(res_x)
  genenames <- INPUT_DATA[which(INPUT_DATA$gene_ID %in% rownames(res_x)),c("gene_ID", "gene_names")]
  res_x <- res_x[order(rownames(res_x)),]
  genenames <- genenames[order(rownames(genenames)),]
  all(genenames$gene_ID == rownames(res_x))
  res_x <- cbind(genenames, res_x)
  res_x <- cbind(res_x, miR430 = (res_x$gene_ID %in% miR430_genes))
  res_x <- cbind(res_x, mitochondrial = (res_x$gene_ID %in% MT_genes))
  res_x <- cbind(res_x, DE="not DE")
  res_x$DE[which(res_x$padj < alpha & res_x$log2FoldChange<(-log2FC_Treshold))] <- "down-regulated"
  res_x$DE[which(res_x$padj < alpha & res_x$log2FoldChange>(log2FC_Treshold))] <- "up-regulated"
  res_x <- res_x[which(res_x$miR430 == FALSE & res_x$mitochondrial == FALSE),]

  # WRITE DEG FILE FOR LATER USE
  write.table(res_x, paste0("./DEG/",filename,".txt"), quote = FALSE, sep = "\t")
}  

# LOAD THE GENES THAT ARE MISREGULATED AT 256C (MIR430 VS WT)
res_y <- read.table(paste0("./DEG/" ,"256c","_miR430_vs_WT",".txt"), sep = "\t")
res_y <- res_y[which(res_y$miR430 == FALSE & res_y$mitochondrial == FALSE & res_y$DE == "up-regulated"),]
print(nrow(res_y))
DEG <- as.character(res_y$gene_ID)

# Determine at which stage these genes that are misregulated at 256c are induced and same them in a list for later use
for (stage_2 in c("1000c", "Oblong", "Sphere"))
{
  res_x <- read.table(paste0("./DEG/", stage_2, "_vs_", stage_1, "_WT.txt"), sep = "\t")
  res_x <- res_x[which(res_x$log2FoldChange > 2),]

  DE <- intersect(as.character(res_y$gene_ID), as.character(res_x$gene_ID))
  write.table(DE, paste0("./DEREPRESSED/", stage_2, "_vs_", stage_1,"_derepressed.tab"))
  print(length(DE))
  res_y <- res_y[which(!(res_y$gene_ID %in% DE)),]
}
write.table(res_y$gene_ID, "./DEREPRESSED/not_derepressed.tab")


# MAKE SCATTERPLOTS FIGURE 4B AND ED5 AND BOXPLOT FIGURE 4A

#Scatterplot parameters
lab_y <- c(expression(paste("256c mir430 vs. Wild type ", log[2], "(Fold Change)")))
lab_x <- list(
  A = c(expression(paste("Wild type 1000c vs. 256c ", log[2], "(Fold Change)"))),
  B = c(expression(paste("Wild type Oblong vs. 256c ", log[2], "(Fold Change)"))),
  C = c(expression(paste("Wild type Sphere vs. 256c ", log[2], "(Fold Change)"))))
lims_x <- c(-4.9, 10.3)
lims_y <- c(0, 8.8)

#Boxplot parameters
BOXPLOT <- data.frame()

# Load res_y and prepare list DEG
res_y <- read.table(paste0("./DEG/" ,"256c","_miR430_vs_WT",".txt"), sep = "\t")
res_y <- res_y[which(res_y$miR430 == FALSE & res_y$mitochondrial == FALSE & res_y$DE == "up-regulated"),]
print(nrow(res_y))
DEG <- as.character(res_y$gene_ID)

#A few lists that we will need
output_list <- c("1000c", "Oblong", "Sphere", "not_derepressed") #loop through i
group_file <- c("1000c_vs_256c_derepressed.tab", "Oblong_vs_256c_derepressed.tab", "Sphere_vs_256c_derepressed.tab", "not_derepressed.tab") #loop through i
group_boxplot <- c("Induced at 1024c", "Induced at Oblong", "Induced at Sphere", "Not induced") #loop through i


# This for-loop will create the scatterplots and prepare on-the-fly the table for the boxplot
for (i in 1:4) #I will create 4 panels, and each panel will highlight genes that are derepressed at 1000c, Oblong, Sphere or not derepressed
{
  output <- output_list[i]
  GROUP <- read.table(paste0("./DEREPRESSED/", group_file[i]))$x
  
  plot_list <- list()
  
  stage_1 <- "256c"
  stage_2 <- c("1000c", "Oblong", "Sphere") #loop through j
  for (j in 1:3) #each panel is made of 3 scatterplots, showing different WT comparisons on the x-axis (1000c vs. 256c, Oblong vs. 256c and Sphere vs 256c)
  {
    res_x <- read.table(paste0("./DEG/", stage_2[j], "_vs_", stage_1, "_WT.txt"), sep = "\t")

    SCATTERPLOT <- data.frame(x = res_x[DEG,"log2FoldChange"], y = res_y[DEG,"log2FoldChange"], Group = "Up"); rownames(SCATTERPLOT) <- DEG
    SCATTERPLOT$Group[which(rownames(SCATTERPLOT) %in% GROUP)] <- "G"
    table(SCATTERPLOT$Group)
    
    p <- ggplot(SCATTERPLOT, aes(x=x, y=y, color = Group)) + geom_point(size = 0.5) + theme_bw() + ylab(lab_y) + xlab(lab_x[[j]]) + theme(legend.position = "none") + theme(aspect.ratio=1) + 
      geom_vline(xintercept = 2) + geom_hline(yintercept=2) + 
      xlim(lims_x) + ylim(lims_y) + 
      scale_color_manual(values=c("#377EB8", "#D3D3D3", "#E41A1C", "#90cdde")) + theme(panel.grid.minor = element_blank())
    plot_list[[j]] <- p
    
    BOXPLOT_SUB <- data.frame(
      gene_ID = rownames(SCATTERPLOT[which(SCATTERPLOT$Group == "G"),]),
      log2FoldChange = SCATTERPLOT$x[which(SCATTERPLOT$Group == "G")],
      Group = group_boxplot[i],
      Comparison = stage_2[j])
    BOXPLOT <- rbind(BOXPLOT, BOXPLOT_SUB)

  }
  
  
  pdf(file=paste0("./FIGURES/","Figure_4B_ED5_induced_at_", output, ".pdf"),width=15,height=4); grid.arrange(grobs = plot_list, ncol = 3); dev.off()
}

# Make boxplot
for (g in unique(BOXPLOT$gene_ID))
{
  gr <- as.character(unique(BOXPLOT$Group[which(BOXPLOT$gene_ID == g)]))
  BOXPLOT <- rbind(BOXPLOT, cbind(gene_ID = g, log2FoldChange = 0, Group = gr, Comparison = "256c"))
}

BOXPLOT$Group <- factor(BOXPLOT$Group, levels = group_boxplot)
BOXPLOT$Comparison <- factor(BOXPLOT$Comparison, levels = c("256c", "1000c", "Oblong", "Sphere"))
BOXPLOT$log2FoldChange <- as.double(BOXPLOT$log2FoldChange)


p<-ggplot(BOXPLOT, aes(x=Comparison, y=log2FoldChange, fill = Comparison)) + theme_bw() + facet_grid(. ~ Group) +
  geom_boxplot(position=position_dodge(1), outlier.size = 0.05) + geom_hline(yintercept=2, linetype = "dashed") + ylab("log2(Expression level) relative to 256c") + theme(aspect.ratio=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
pdf(file=paste0("./FIGURES/", "Figure_4A_induced_genes_boxplot.pdf"),width=8,height=4); print(p); dev.off()


#----------------------------------------------------------
# FIGURE ED4: KARIOGRAM OF THE ZGA INDUCED GENES
#----------------------------------------------------------
stage_1 <- "256c"
stage_2 <- c("1000c", "Oblong", "Sphere")

for (i in 1:length(stage_2))
{
  
  filename <- paste0(stage_2[i], "_vs_", stage_1, "_WT")
  res_filtering <- read.table(paste0("./DEG/",filename,".txt"), sep = "\t")
  
  with_miR430 <- TRUE
  de <- "up-regulated"
  color <- "#377EB8" #Chosen colors for down (first) and up (second) genes
  miR430 <- ""
  
  gene_list <- rownames(res_filtering[which(res_filtering$miR430 == FALSE & res_filtering$mitochondrial == FALSE & res_filtering$DE == de),])
  length(gene_list)
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
  pdf(file=paste0("./FIGURES/","Figure_ED4D_", filename, "_", de, "_ZGA_gene_distribution_barplot_paper_", miR430, ".pdf"),width=7,height=6); print(p); dev.off()
  
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
  
  chr_zebrafish <- GRanges(seqnames = chrom_length$seqnames, strand = c("*"),ranges = IRanges(start=1, width=chrom_length$length))
  
  GTF_sub <- GRanges(seqnames = GTF_sub$seqnames, strand = c("*"), ranges <- IRanges(start=GTF_sub$start, end=GTF_sub$end))
  seqlengths(GTF_sub) <- chrom_length$length
  
  # KARIOGRAM
  color_value <- color
  p <- autoplot(GTF_sub, layout="karyogram", aes(color="", fill="")) + scale_color_manual(values=color_value) + scale_fill_manual(values=color_value)
  pdf(file=paste0("./FIGURES/","Figure_ED4C_", filename, "_", de, "_ZGA_gene_distribution_paper_", miR430, ".pdf"),width=7,height=6); print(p@ggplot); dev.off()
}



#----------------------------------------------------------
# CREATE GENE LISTS FOR WHICH YOU WANT TO EXTRACT PROMOTERS (FOR FIGURE 3CD AND FIGURE S13)
#----------------------------------------------------------

# Up-regulated genes at 256c (for Figure S8 and Figure 3CDE)
TABLE <- read.table("./DEG/256c_miR430_vs_WT.txt", sep = "\t")
genes_up   <- TABLE$gene_ID[which(TABLE$DE == "up-regulated" & TABLE$miR430 == FALSE & TABLE$mitochondrial == FALSE)]
length(genes_up)

# Not expressed (for Figure S8)
DESIGN <- INPUT_DESIGN[which(INPUT_DESIGN$IAA_processing == "yes"),]
DESIGN <- cbind(DESIGN, Group = paste0(DESIGN$Stage, "_", DESIGN$Genotype))
DATA <- INPUT_DATA[which(!(INPUT_DATA$Chr == "MT")),]
DATA <- as.matrix(INPUT_DATA[,-(1:7)])
all(rownames(DESIGN) %in% colnames(DATA)) # should be TRUE
DATA <- DATA[,rownames(DESIGN)]
DATA <- DATA[rowSums(DATA) < 10,]
not_expressed <- rownames(DATA)
length(not_expressed)

# Top 716 expressed genes (for figure 3CDE)
stage_TPM <- "256c"
TPM <- INPUT_DATA[,rownames(INPUT_DESIGN[which(INPUT_DESIGN$Genotype == "WT" & INPUT_DESIGN$IAA_processing == "yes" & INPUT_DESIGN$Stage == stage_TPM),])]
TPM <- TPM/INPUT_DATA$Length
TPM <- t( t(TPM) * 1e6 / colSums(TPM))
TPM_mean <- as.data.frame(rowMeans(TPM)); colnames(TPM_mean) <- "meanTPM"
TPM_mean$gene_id <- rownames(TPM_mean)
TPM_mean <- TPM_mean[which(!(TPM_mean$gene_id %in% c(MT_genes, miR430_genes))),]
TPM_mean_ranked <- cbind(TPM_mean, rank = rank(-TPM_mean$meanTPM, ties.method= "random"))
top_716_genes <- TPM_mean_ranked$gene_id[which(TPM_mean_ranked$rank %in% 1:716)]
write.table(top_716_genes, "top_716_genes.txt")
length(top_716_genes)

gene_list <- list(A = genes_up, B = not_expressed, C = top_716_genes)
output <- c("upregulated", "not_expressed", "top716")
length(gene_list[[1]])
length(gene_list[[2]])
length(gene_list[[3]])

# To create Figure 3D, you have to create the .bed file for all genes. See comments in the corresponding script
# The creation of this list takes a lot of time
# ALL_GENES <- unique(GTF$gene_id)
# gene_list <- list(A = ALL_GENES)
# output <- c("all_genes")

# Make BED Files for the promoters of these 3 lists
# Promotes are defined as plus/minus 2000bp centered around the genes TSS


promoter_length <- 2000
GTF_PROMOTERS <- as.data.frame(import(gtf_input))


for (i in 1:length(gene_list))
{
  out <- c()
  
  PROMOTERS = data.frame()
  j <- 0
  for (gene in gene_list[[i]])
  {
    if ((j %% 1000) == 0){print(j)}; j <- j+1
    
    # Do not consider things like retained introns, which we don't care about
    g <- GTF_PROMOTERS[which(GTF_PROMOTERS$gene_id == gene & GTF_PROMOTERS$type == "transcript" & !(GTF_PROMOTERS$transcript_biotype %in% c("retained_intron", "sense_intronic", "nonsense_mediated_decay"))),]

    # Choose transcripts without the mRNA_start_NF and cds_start_NF tag (incomplete transcript). Do it only if that would not remove all entries!
    g <- g[which(!(g$tag %in% c("mRNA_start_NF", "cds_start_NF"))),]

    if (nrow(g) == 0){print("No rows!"); print(gene); next}
    if (length(unique(g$strand)) > 1){print("Multiple strands"); print(gene); next}
    
    # Check how many genes have multiple TSSs
    if((g$strand[1] == "-" & !(length(unique(g$end)) == 1)) | (g$strand[1] == "+" & !(length(unique(g$start)) == 1))){
      out <- c(out, gene)
    }
    
    # Choose the most upstream promoter
    if(g$strand[1] == "-"){
      TSS <- max(g$end)
    }
    if(g$strand[1] == "+"){
      TSS <- min(g$start)
    }
    
    start_p <- TSS - promoter_length
    end_p <- TSS + promoter_length
    
    if (start_p < 0){start_p <- 0}
    
    g <- g[1,]; g$type <- "promoter";g$start <- start_p;g$end <- end_p;g$width <- g$end - g$start + 1
    PROMOTERS <- rbind(PROMOTERS, g)
  }
  
  print(length(out))
  
  #s: gene_name + promoter coordinates
  #f: just promoter coordinates
  
  PROMOTERS <- PROMOTERS[which(PROMOTERS$width == 2*promoter_length+1),]
  PROMOTERS <- cbind(PROMOTERS, s = paste0(PROMOTERS$gene_name,":", PROMOTERS$seqnames, ":", PROMOTERS$start, "-",PROMOTERS$end))
  PROMOTERS <- PROMOTERS[,c(1:3,26,4:25)]
  PROMOTERS <- cbind(PROMOTERS, f = paste0(PROMOTERS$seqnames, ":", PROMOTERS$start, "-",PROMOTERS$end))
  
  #Some genes give rise to the same promoter. We want each promoter region to be present once, so we keep only one gene (randomly) for each identical promoter region f
  PROMOTERS_CORRECTED <- data.frame()
  for (f in unique(PROMOTERS$f))
  {
    PR <- PROMOTERS[which(PROMOTERS$f == f),]
    PROMOTERS_CORRECTED <- rbind(PROMOTERS_CORRECTED, PR[1,])
  }

  options(scipen = 999) #Will avoid using scientific notation
  write.table(PROMOTERS_CORRECTED, file = paste0("./BED_FASTA_FILES/", promoter_length, "_" ,output[i], ".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  options(scipen = 0)
}

# In order to obtain fasta files from these BED files, use bedtools getfasta:
# module load gcc
# module load bedtools2

#bedtools getfasta -fi Danio_rerio.GRCz11.105.dna.primary_assembly.fa -bed 2000_upregulated.bed -fo 2000_upregulated.fa -nameOnly
#bedtools getfasta -fi Danio_rerio.GRCz11.105.dna.primary_assembly.fa -bed 2000_not_expressed.bed -fo 2000_not_expressed.fa -nameOnly
#bedtools getfasta -fi Danio_rerio.GRCz11.105.dna.primary_assembly.fa -bed 2000_top716.bed -fo 2000_top716.fa -nameOnly


#----------------------------------------------------------
# FIGURE 4E (MOTIF ENRICHMENT ANALYSIS SCATTERPLOT)
#----------------------------------------------------------

# Before generating this Figure, you must use the FASTA FILES generated before and the motif database (provided) to run SEA (MEME SUITE) (standard settings).
# You will compare the up-regulated genes with the non-expressed genes to assess enrichment, and then the opposite comparison to assess depletion.
# Save both obtained TSV files and use them to generate the following script.
# One needs to take the inverse of the depletion values
# The motif database is based on the JASPAR vertebrate database. Only the motifs of TFs that are expressed (>10 RPKM) according to Lee et al. 2013, are kept.
# Zebrafish-specific motifs for Pou5f3 and Nanog were added from the Leichsenring and Palfy papers.
# Default settings of SEA: Motifs with p-value < 0.05 and E-value < 10 are considered.

group_enriched <- "SEA_up_vs_not_expressed_enriched.tsv"
group_depleted <- "SEA_not_expressed_vs_up_depleted.tsv"

pluripotency <- c("POU5F1", "POU5F1B", "Pou5f1::Sox2", "Pou5f1_pre-MBT_Leichsenring", "POU2F1","POU2F1::SOX2", "Nanog-like_3.5hpf-zebrafish",
                  "POU1F1", "POU3F1", "POU3F3", "POU2F3", "POU3F2", "POU1F1", "POU2F2", "POU2F1::SOX2", "POU3F1", "POU3F4")

file_enriched <- paste0("./BED_FASTA_FILES/", group_enriched)
file_depleted <- paste0("./BED_FASTA_FILES/", group_depleted)

TAB_enriched <- read.table(file_enriched, header = TRUE)
TAB_depleted <- read.table(file_depleted, header = TRUE)

TAB_enriched$TP_in_up <- TAB_enriched$TP.
TAB_depleted$TP_in_up <- TAB_depleted$FP.
TAB_depleted$ENR_RATIO <- 1/TAB_depleted$ENR_RATIO

TAB <- rbind(TAB_enriched, TAB_depleted)
TAB <- TAB[which(TAB$PVALUE < 0.05),]
TAB$log2ENR_RATIO <- log2(TAB$ENR_RATIO)
TAB$Pluripotent <- "Not pluripotency factor"; TAB$Pluripotent[which(TAB$ALT_ID %in% pluripotency)] <- "Pluripotency factor"

TAB_P <- TAB[which(TAB$Pluripotent == "Pluripotency factor"),]

p <- ggplot(TAB, aes(x=TP_in_up, y=log2ENR_RATIO, color = Pluripotent)) + geom_point() + theme_bw() + ylab("log2(Enrichment Ratio)") + 
  xlab("True Positive Percentage") + scale_color_manual(values=c("#A0A0A0", "#E41A1C")) + geom_hline(yintercept=0) + 
  ggtitle(group_enriched) + theme(aspect.ratio=1) +
  geom_point(data = TAB_P, mapping = aes(x=TP_in_up, y=log2ENR_RATIO))
p
pdf(file="./FIGURES/FIGURE_4E.pdf", width=6, height=4); print(p); dev.off()



#----------------------------------------------------------
# FIGURE 3A (BOXPLOT OF MEAN TPM)
#----------------------------------------------------------

# Up-regulated, Down-regulated and Not DE genes at 256c
TABLE <- read.table("./DEG/256c_miR430_vs_WT.txt", sep = "\t")
genes_up   <- TABLE$gene_ID[which(TABLE$DE == "up-regulated" & TABLE$miR430 == FALSE & TABLE$mitochondrial == FALSE)]
genes_down   <- TABLE$gene_ID[which(TABLE$DE == "down-regulated" & TABLE$miR430 == FALSE & TABLE$mitochondrial == FALSE)]
genes_not_DE <- TABLE$gene_ID[which(TABLE$DE == "not DE" & TABLE$miR430 == FALSE & TABLE$mitochondrial == FALSE & abs(TABLE$log2FoldChange) < 1)]
length(genes_up)
length(genes_down)
length(genes_not_DE)

#All genes
all_genes <- unique(GTF$gene_id)

# Not expressed
DATA <- INPUT_DATA[which(!(INPUT_DATA$Chr == "MT")),]
DATA <- as.matrix(DATA[,-(1:7)])
all(rownames(DESIGN) %in% colnames(DATA)) # should be TRUE
DATA <- DATA[,rownames(DESIGN)]
DATA <- DATA[rowSums(DATA) < 10,]
genes_not_expressed <- rownames(DATA)
length(genes_not_expressed)


# DETERMINE TPM
stage_TPM <- "256c"
TPM <- INPUT_DATA[,rownames(INPUT_DESIGN[which(INPUT_DESIGN$Genotype == "WT" & INPUT_DESIGN$IAA_processing == "yes" & INPUT_DESIGN$Stage == stage_TPM),])]
TPM <- TPM/INPUT_DATA$Length
TPM <- t( t(TPM) * 1e6 / colSums(TPM))
TPM_mean <- as.data.frame(rowMeans(TPM)); colnames(TPM_mean) <- "meanTPM"
TPM_mean$gene_id <- rownames(TPM_mean)
TPM_mean <- TPM_mean[which(!(TPM_mean$gene_id %in% c(MT_genes, miR430_genes))),]
TPM_mean_ranked <- cbind(TPM_mean, rank = rank(-TPM_mean$meanTPM, ties.method= "random"))

all <-    cbind(meanTPM = TPM_mean_ranked[TPM_mean_ranked$gene_id %in% all_genes,"meanTPM"], Group = "All")
genes_not_expressed <- cbind(meanTPM = TPM_mean_ranked[TPM_mean_ranked$gene_id %in% genes_not_expressed,"meanTPM"], Group = "Not expressed")
genes_not_DE <- cbind(meanTPM = TPM_mean_ranked[TPM_mean_ranked$gene_id %in% genes_not_DE,"meanTPM"], Group = "Not DE")
up <-     cbind(meanTPM = TPM_mean_ranked[TPM_mean_ranked$gene_id %in% genes_up,"meanTPM"], Group = "Upregulated genes")
down <-     cbind(meanTPM = TPM_mean_ranked[TPM_mean_ranked$gene_id %in% genes_down,"meanTPM"], Group = "Downregulated genes")
TABLE <- as.data.frame(rbind(all, genes_not_expressed, genes_not_DE, up, down))
TABLE$meanTPM <- as.double(TABLE$meanTPM)
TABLE$Group <- factor(TABLE$Group, levels = c("All", "Not expressed", "Not DE", "Upregulated genes", "Downregulated genes"))

p <- ggplot(TABLE, aes(x=Group, y=meanTPM)) + geom_boxplot(outlier.shape = NA, fill = "#A0A0A0") + theme_bw() +
  scale_y_continuous(limits = c(0, 5)) + theme(aspect.ratio=1) + ylim(c(0,2.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
pdf(file=paste0(paste0("./FIGURES/","Figure_3A_mean_TPM_boxplot_no_outliers_shown_", stage_TPM, ".pdf")),width=4,height=4);print(p);dev.off()



#---------------------------------------
#--------- ALLUVIAL PLOT (2D) ---------
#---------------------------------------

#de <- "down-regulated"; color_values <- c("#E41A1C", "#A0A0A0"); levels_condition <- c("down-regulated", "not DE")
de <- "up-regulated" ; color_values <- c("#377EB8", "#A0A0A0"); levels_condition <- c("up-regulated", "not DE")

stage <- c("256c", "1000c", "Oblong", "Sphere")

DE_SLAM_256c <- read.table(paste0("./DEG/", stage[1], "_miR430_vs_WT.txt"), sep = "\t")
DE_SLAM_1024c <- read.table(paste0("./DEG/", stage[2], "_miR430_vs_WT.txt"), sep = "\t")
DE_SLAM_Oblong <- read.table(paste0("./DEG/", stage[3], "_miR430_vs_WT.txt"), sep = "\t")
DE_SLAM_Sphere <- read.table(paste0("./DEG/", stage[4], "_miR430_vs_WT.txt"), sep = "\t")

x <- list(
  S_256c =   DE_SLAM_256c$gene_ID[which(DE_SLAM_256c$DE == de)], 
  S_1024c =  DE_SLAM_1024c$gene_ID[which(DE_SLAM_1024c$DE == de)], 
  S_Oblong = DE_SLAM_Oblong$gene_ID[which(DE_SLAM_Oblong$DE == de)], 
  S_Sphere = DE_SLAM_Sphere$gene_ID[which(DE_SLAM_Sphere$DE == de)] )
all_DE <- unique(c(x[[1]], x[[2]], x[[3]], x[[4]]))


GENES <- data.frame()

for (stage_i in stage)
{
  res <- read.table(paste0("./DEG/", stage_i, "_miR430_vs_WT.txt"), sep = "\t")
  res <- res[which(res$gene_ID %in% all_DE),]
  res$stage <- stage_i
  GENES <- rbind(GENES, res)
}

GENES$stage <- paste0("Stage_", GENES$stage); stage <- paste0("Stage_", stage)
GENES$stage <- factor(GENES$stage , levels = stage)
GENES$DE <- factor(GENES$DE, levels = levels_condition)

p <- ggplot(GENES, aes(x = stage, stratum = DE, alluvium = gene_ID, fill = DE, label = DE)) +
  geom_stratum() +
  geom_text(stat = "stratum",aes(label = after_stat(stratum)), size = 3) +
  scale_fill_manual(values=color_values) + ylab("Number of genes") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") + theme(aspect.ratio=1) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "none")

pdf(file=paste0("./FIGURES/","Figure_2D_", de, "_Alluvial.pdf"),width=4,height=4); p; dev.off()
getwd()


