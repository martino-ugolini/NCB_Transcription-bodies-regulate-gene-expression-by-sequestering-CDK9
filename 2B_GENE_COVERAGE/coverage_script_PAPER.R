
library(rtracklayer)
library(IRanges)
library(GenomicRanges)
library(ggbio)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(gridExtra)
library(ggpubr)
library(scales)
library(cowplot)
library("RColorBrewer")
library(dplyr)

rm(list=ls())

#Load the GTF file
gtf_input <- "./Danio_rerio.GRCz11.105.gtf.gz"
GTF_file <- import(gtf_input)
txdb <- makeTxDbFromGRanges(GTF_file)

# SET PARAMETERS
path_coverage <- "./" # The bam files that end on _aligned.unique.filtered.sorted.1.bam, ...2.bam, ...3.bam, ...4.bam and ...more.4.bam should be in this folder
stage <- "256c"
path_output <- "./"
acc <- 0.001
color_palette_coverage <- c("#08306B", "#2171B5", "#6BAED6", "#9ECAE1", "#DEEBF7")

# CHOOSE HERE THE GENE YOU WANT TO VISUALIZE
#UP
gene_ID <- "ENSDARG00000002601"; samples <- c('miR430_A2', 'WT_A1')
#DOWN
#gene_ID <- "ENSDARG00000022978"; samples <- c('miR430_A2', 'WT_A1')



#This function plots a composite coverage
get_max_coverage <- function(path, sample_list)
{
  max_value <- c()
  for (sample in sample_list)
  {
    DATA <- data.frame()
    for (idx in c('1', '2', '3', '4', 'more.4'))
    {
      bam_file <- paste0(path,sample,'_aligned.unique.filtered.sorted.',idx,".bam") #indexed .bam file (.bai file with same name)
      xcov <- readGAlignments(bam_file, param=ScanBamParam(which=gene_range))
      xcov <- coverage(xcov); #xcov[z]; xcov$chr4[ranges(z)]
      xcov <- as.numeric(xcov[[seqnames(gene_range)]][ranges(gene_range)])
      xcov <- as.data.frame(cbind(coverage = xcov, interval = c(start(gene_range):end(gene_range)), filter = idx))
      DATA <- rbind(DATA, xcov)
    }
    DATA$coverage <- as.numeric(DATA$coverage)
    DATA$interval <- as.numeric(DATA$interval)
    DATA <- DATA %>%
      group_by(interval) %>%
      summarise(sum = sum(coverage), n = n())
    max_value <- c(max_value, max(DATA$sum))
  }
  return(max(max_value))
}


#This function plots a composite coverage
plot_coverage <- function(path, sample, acc)
{
  DATA <- data.frame()
  for (idx in c('1', '2', '3', '4', 'more.4'))
  {
    bam_file <- paste0(path,sample,'_aligned.unique.filtered.sorted.',idx,".bam") #indexed .bam file (.bai file with same name)
    xcov <- readGAlignments(bam_file, param=ScanBamParam(which=gene_range))
    xcov <- coverage(xcov); #xcov[z]; xcov$chr4[ranges(z)]
    xcov <- as.numeric(xcov[[seqnames(gene_range)]][ranges(gene_range)])
    xcov <- as.data.frame(cbind(coverage = xcov, interval = c(start(gene_range):end(gene_range)), filter = idx))
    DATA <- rbind(DATA, xcov)
  }
  DATA$coverage <- as.numeric(DATA$coverage)
  DATA$interval <- as.numeric(DATA$interval)
  DATA$filter <- factor(DATA$filter, levels=c('more.4', '4', '3', '2', '1'))
  
  plot <- ggplot(DATA, aes(y=coverage, x=interval, fill=filter)) +  geom_area() + theme_classic() + scale_fill_manual(values=color_palette_coverage)
  #plot <- plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  plot <- plot + xlim(start(gene_range),end(gene_range)) + xlab("")
  plot <- plot + ggtitle(sample) + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + scale_x_continuous(label = unit_format(unit = "Mb", sep="",accuracy = acc, scale = 1e-6)) + theme(legend.position="none")
  return(plot)
}


setwd(path_coverage)
gene_range <- GTF_file[which(GTF_file$gene_id == gene_ID),]; gene_range <- range(gene_range, ignore.strand=TRUE)
gene_name <- unique(GTF_file$gene_name[which(GTF_file$gene_id == gene_ID)])
gene_name <- paste(gene_name, collapse = ";")

plots <- list()

#Plot coverage
#The sorted and indexed .bam files for 1, 2, 3, 4 and more than 4 point mutations must be in the same path folder
#brewer.pal(n = 9, name = "YlOrRd")
max_coverage <- get_max_coverage(path_coverage, samples)
for (i in 1:length(samples))
{
  p <- plot_coverage(path_coverage, samples[i], acc) + ylim(0, 1.1*max_coverage)
  plots[[i]] <- p
}

# Plot gene model + blue rectangle
plot_model <- autoplot(txdb, which = gene_range, gap.geom = "chevron", label=FALSE, fill='black') + theme_classic2()
plot_model <- plot_model + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
plot_model <- plot_model + scale_x_continuous(label = unit_format(unit = "Mb", sep="",accuracy = acc, scale = 1e-6))
plots[[length(samples) + 1]] <- plot_model@ggplot


#Output coverage
pdf(file=paste0(path_output,gene_ID, " - (", gene_name, ") - ", stage, ".pdf"),width=6,height=8)
plot_grid(plotlist = plots, ncol = 1, align='v', axis = "l", rel_heights = c(6,6,2))
dev.off()
