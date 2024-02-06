
library(ggplot2)
library(tidyr)
library(stringr)
library("RColorBrewer")
library(DescTools)
library(ggpubr)
library(cowplot)
library(plyr)
library(gridExtra)
library("readxl")
library("RColorBrewer")

rm(list=ls())

# Create a list of all spots
SPOTS <- read_excel("filtered_spots.xlsx")
SPOTS$GROUP <- factor(SPOTS$GROUP, levels = c("WT_512", "HOM_512", "WT_1024", "HOM_1024"))
#SPOTS <- SPOTS[which((SPOTS$Freq_10 + SPOTS$Freq_11 + SPOTS$Freq_01) > 0),]

SPOTS$Replicate <- paste0(sapply(strsplit(SPOTS$NUCLEUS,"_"), `[`, 1), "_", SPOTS$GROUP)
SPOTS$sum <- SPOTS$Freq_10 + SPOTS$Freq_11 + SPOTS$Freq_01

#Correct: Cell cycle begins when first spot appears
SPOTS_CORRECTED <- data.frame()


for (nu in unique(SPOTS$NUCLEUS))
{
  SPOTS_sub <- SPOTS[which(SPOTS$NUCLEUS == nu),]
  SPOTS_sub <- SPOTS_sub[order(SPOTS_sub$FRAME),]

  while(head(SPOTS_sub$sum, n=1) == 0)
  {
    SPOTS_sub <- SPOTS_sub[-1,]
  }
  
  while(tail(SPOTS_sub$sum, n=1) == 0)
  {
    SPOTS_sub <- SPOTS_sub[-nrow(SPOTS_sub),]
  }
  
  
  if (nrow(SPOTS_sub) == 1)
  {print("Error! Only one frmae"); print(nu); next}
  
  SPOTS_sub$FRAME <- 0:(nrow(SPOTS_sub)-1)
#  SPOTS_sub$FRAME_P <- (SPOTS_sub$FRAME)/(nrow(SPOTS_sub) - 1) #Check that it is actually correct
  SPOTS_CORRECTED <- rbind(SPOTS_CORRECTED, SPOTS_sub)
}


SPOTS <- SPOTS_CORRECTED

# Calculate MEAN FREQUENCY/Nucleus
FREQ <- data.frame()
for (gr in unique(SPOTS$GROUP))
{
  for (fr in unique(SPOTS$FRAME[which(SPOTS$GROUP == gr)]))
  {
    #print(gr); print(fr)
    SPOTS_sub <- SPOTS[which(SPOTS$FRAME == fr & SPOTS$GROUP == gr),]
    tot_nuclei <- length(unique(SPOTS_sub$NUCLEUS))
    FREQ_total <- c(sum(SPOTS_sub$Freq_10)/tot_nuclei,
                    sum(SPOTS_sub$Freq_11)/tot_nuclei,
                    sum(SPOTS_sub$Freq_01)/tot_nuclei)

    FREQ_percentage <- FREQ_total/sum(FREQ_total)
    FREQ_sub <- cbind(SPOT_TYPE = c("10", "11", "01"), FREQ_total = FREQ_total, FREQ_percentage = FREQ_percentage, FRAME = fr, GROUP = gr, TOT_NUCLEI = tot_nuclei)
    FREQ <- rbind(FREQ, FREQ_sub)
  }
}
FREQ$SPOT_TYPE <- factor(FREQ$SPOT_TYPE, levels = c("01", "11", "10"))
FREQ$GROUP <- factor(FREQ$GROUP, levels = c("WT_512", "HOM_512", "WT_1024", "HOM_1024"))
FREQ$FREQ_total <- as.double(FREQ$FREQ_total)
FREQ$FREQ_percentage <- 100*as.double(FREQ$FREQ_percentage)
FREQ$FRAME <- as.double(FREQ$FRAME)
FREQ$TOT_NUCLEI <- as.integer(FREQ$TOT_NUCLEI)
FREQ$TIME <- FREQ$FRAME*120

FREQ <- FREQ[which(FREQ$TOT_NUCLEI > 2),]

p1 <- ggplot(data=FREQ, aes(x=TIME, y=FREQ_total, fill = SPOT_TYPE)) + geom_bar(stat="identity") + theme_bw() + facet_grid(cols = vars(GROUP)) + xlab("Time since appearence of first spots [sec]")  + scale_fill_manual(values = c("#4DAF4A", "#FFD92F", "#E41A1C")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2 <- ggplot(data=FREQ, aes(x=TIME, y=FREQ_percentage, fill=SPOT_TYPE)) + geom_bar(stat="identity") + theme_bw() + facet_grid(cols = vars(GROUP)) + xlab("Time since appearence of first spots [sec]") + scale_fill_manual(values = c("#4DAF4A", "#FFD92F", "#E41A1C")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file="Figure_5G_total_spots_per_nucleus_over_time_barplot.pdf",width=8,height=3);print(p1);dev.off()
pdf(file="Figure_5G_total_spots_per_nucleus_over_time_percentage_barplot.pdf",width=8,height=3);print(p2);dev.off()

