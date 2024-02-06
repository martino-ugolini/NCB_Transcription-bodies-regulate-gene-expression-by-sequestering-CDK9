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
library(reshape2)

rm(list=ls())

# Create a list of all spots
TABLE_MAIN <- read_excel("cdk9_quantification.xlsx")
TABLE_NO_PRIMARY <- read_excel("cdk9_quantification_no_primary.xlsx")

MeanIntensityNoPrimary_DAPI   <- mean(TABLE_NO_PRIMARY$Mean[which(TABLE_NO_PRIMARY$Channel_Description == "DAPI")])
MeanIntensityNoPrimary_CDK9   <- mean(TABLE_NO_PRIMARY$Mean[which(TABLE_NO_PRIMARY$Channel_Description == "CDK9")])
MeanIntensityNoPrimary_RNAPII <- mean(TABLE_NO_PRIMARY$Mean[which(TABLE_NO_PRIMARY$Channel_Description == "RNAPII_Ser2P")])


PLOT_TABLE <- data.frame()
for (samp in unique(TABLE_MAIN$Sample))
{
  SUB <- TABLE_MAIN[which(TABLE_MAIN$Sample == samp),]
  
  #if (!(nrow(SUB == 9))){print(samp)}
  #if (!(setequal(SUB$Channel, c(1,2,3,1,2,3,1,2,3)))){print(samp)}
  #if (!(setequal(SUB$Channel_Description, c("DAPI", "CDK9", "RNAPII_Ser2P", "DAPI", "CDK9", "RNAPII_Ser2P", "DAPI", "CDK9", "RNAPII_Ser2P")))){print(samp)}

  Percentage_Area <- 100 * (SUB$Area[4] + SUB$Area[7]) / SUB$Area[1]
  
  Percentage_Intensity_DAPI   <- 100 * (SUB$Area[4]*SUB$Mean[4] + SUB$Area[7]*SUB$Mean[7]) / (SUB$Area[1] * (SUB$Mean[1] - MeanIntensityNoPrimary_DAPI) )
  Percentage_Intensity_CDK9   <- 100 * (SUB$Area[5]*SUB$Mean[5] + SUB$Area[8]*SUB$Mean[8]) / (SUB$Area[2] * (SUB$Mean[2] - MeanIntensityNoPrimary_CDK9) )
  Percentage_Intensity_RNAPII <- 100 * (SUB$Area[6]*SUB$Mean[6] + SUB$Area[9]*SUB$Mean[9]) / (SUB$Area[3] * (SUB$Mean[3] - MeanIntensityNoPrimary_RNAPII) )
  
  Mean_Intensity_2_Bodies_DAPI   <- ( (SUB$Area[4]*SUB$Mean[4] + SUB$Area[7]*SUB$Mean[7]) / (SUB$Area[4] + SUB$Area[7]) )
  Mean_Intensity_2_Bodies_CDK9   <- ( (SUB$Area[5]*SUB$Mean[5] + SUB$Area[8]*SUB$Mean[8]) / (SUB$Area[5] + SUB$Area[8]) )
  Mean_Intensity_2_Bodies_RNAPII <- ( (SUB$Area[6]*SUB$Mean[6] + SUB$Area[9]*SUB$Mean[9]) / (SUB$Area[6] + SUB$Area[9]) )
  
  Mean_Intensity_Enrichment_DAPI   <- Mean_Intensity_2_Bodies_DAPI / (SUB$Mean[1] - MeanIntensityNoPrimary_DAPI)
  Mean_Intensity_Enrichment_CDK9   <- Mean_Intensity_2_Bodies_CDK9 / (SUB$Mean[2] - MeanIntensityNoPrimary_CDK9)
  Mean_Intensity_Enrichment_RNAPII <- Mean_Intensity_2_Bodies_RNAPII / (SUB$Mean[3] - MeanIntensityNoPrimary_RNAPII)
  
  Values <- c(Percentage_Area, Percentage_Intensity_DAPI, Percentage_Intensity_CDK9, Percentage_Intensity_RNAPII, Mean_Intensity_Enrichment_DAPI, Mean_Intensity_Enrichment_CDK9, Mean_Intensity_Enrichment_RNAPII)
  Description <- c("Percentage_Area", "Percentage_Intensity_DAPI", "Percentage_Intensity_CDK9", "Percentage_Intensity_RNAPII", "Mean_Intensity_Enrichment_DAPI", "Mean_Intensity_Enrichment_CDK9", "Mean_Intensity_Enrichment_RNAPII")
  Channel = c("Area", "DAPI", "CDK9", "RNAPII_Ser2P", "DAPI", "CDK9", "RNAPII_Ser2P")
  Type = c("Area", "Percentage_Intensity", "Percentage_Intensity", "Percentage_Intensity", "Mean_Intensity_Enrichment", "Mean_Intensity_Enrichment", "Mean_Intensity_Enrichment")
  
  
  PLOT_TABLE <- rbind(PLOT_TABLE,
                      cbind(Values = Values, Description = Description, Channel = Channel, Type = Type, Sample = samp, Replicate = SUB$Replicate[1]))
}

PLOT_TABLE$Values <- as.double(PLOT_TABLE$Values)
PLOT_TABLE$Channel <- factor(PLOT_TABLE$Channel, levels = c("Area", "DAPI", "CDK9", "RNAPII_Ser2P"))




# AREA
p <- ggplot(data=PLOT_TABLE[which(PLOT_TABLE$Type == "Area"),], aes(x=Channel, y=Values, fill = Channel)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("Channel")  + ylab("% Area of the 2 bodies") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 3)+ scale_fill_manual(values=c("#A0A0A0"))
p
pdf(file="Percentage_Area.pdf",width=6,height=4); print(p); dev.off()

nrow(PLOT_TABLE[which(PLOT_TABLE$Type == "Area"),])

# MEAN INTENSITY ENRICHMENT
p <- ggplot(data=PLOT_TABLE[which(PLOT_TABLE$Type == "Percentage_Intensity" & PLOT_TABLE$Channel %in% c("CDK9", "RNAPII_Ser2P")),], aes(x=Channel, y=Values, fill = Channel)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("")  + ylab("% Intensity of the 2 bodies") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1.5)+ scale_fill_manual(values=c("#008000", "#E41A1C"))
p
pdf(file="Percentage_Intensity.pdf",width=6,height=4); print(p); dev.off()

nrow(PLOT_TABLE[which(PLOT_TABLE$Type == "Percentage_Intensity" & PLOT_TABLE$Channel == "CDK9"),])
nrow(PLOT_TABLE[which(PLOT_TABLE$Type == "Percentage_Intensity" & PLOT_TABLE$Channel == "RNAPII_Ser2P"),])


# PERCENTAGE INTENSITY
p <- ggplot(data=PLOT_TABLE[which(PLOT_TABLE$Type == "Mean_Intensity_Enrichment" & PLOT_TABLE$Channel %in% c("CDK9", "RNAPII_Ser2P")),], aes(x=Channel, y=Values, fill = Channel)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("")  + ylab("Enrichment of Mean Intensity of the 2 bodies") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1.5)+ scale_fill_manual(values=c("#008000", "#E41A1C")) + geom_hline(yintercept = 1)
p
pdf(file="Mean_Intensity_Enrichment.pdf",width=6,height=4); print(p); dev.off()

nrow(PLOT_TABLE[which(PLOT_TABLE$Type == "Mean_Intensity_Enrichment" & PLOT_TABLE$Channel == "CDK9"),])
nrow(PLOT_TABLE[which(PLOT_TABLE$Type == "Mean_Intensity_Enrichment" & PLOT_TABLE$Channel == "RNAPII_Ser2P"),])

STAT <- PLOT_TABLE[which(PLOT_TABLE$Type == "Area"),]
median(STAT$Values)
table(STAT$Replicate)
wilcox.test(x = STAT$Values,  mu = 0, alternative = 'g')


STAT <- PLOT_TABLE[which(PLOT_TABLE$Type == "Percentage_Intensity" & PLOT_TABLE$Channel == "CDK9"),]
median(STAT$Values)
table(STAT$Replicate)
wilcox.test(x = STAT$Values,  mu = 0, alternative = 'g')

STAT <- PLOT_TABLE[which(PLOT_TABLE$Type == "Percentage_Intensity" & PLOT_TABLE$Channel == "RNAPII_Ser2P"),]
table(STAT$Replicate)
wilcox.test(x = STAT$Values,  mu = 0, alternative = 'g')

STAT <- PLOT_TABLE[which(PLOT_TABLE$Type == "Mean_Intensity_Enrichment" & PLOT_TABLE$Channel == "CDK9"),]
table(STAT$Replicate)
wilcox.test(x = log2(STAT$Values),  mu = 0)

STAT <- PLOT_TABLE[which(PLOT_TABLE$Type == "Mean_Intensity_Enrichment" & PLOT_TABLE$Channel == "RNAPII_Ser2P"),]
table(STAT$Replicate)
wilcox.test(x = log2(STAT$Values),  mu = 0)


p <- ggplot(data=STAT, aes(x=Channel, y=Values, fill = Channel)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("")  + ylab("Enrichment of Mean Intensity of the 2 bodies") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1.5)+ scale_fill_manual(values=c("#008000", "#E41A1C")) + geom_hline(yintercept = 1)
p





