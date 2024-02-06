library(ggplot2)
library("ggpubr")

rm(list = ls())

# LOAD DATA NO RESCUE
output <- "Summary_No_Rescue_publication"
TABLE_1 <- read.csv("No_rescue_1_GFP_Sensor_Analysis.csv", sep=";")
TABLE_2 <- read.csv("No_rescue_2_GFP_Sensor_Analysis.csv", sep=";")
TABLE_3 <- read.csv("No_rescue_3_GFP_Sensor_Analysis.csv", sep=";")
TABLE <- rbind(TABLE_1, TABLE_2, TABLE_3)

# LOAD DATA WITH RESCUE
output <- "Summary_With_Rescue_publication"
TABLE_1 <- read.csv("With_rescue_1_GFP_Sensor_Analysis.csv", sep=";")
TABLE_2 <- read.csv("With_rescue_2_GFP_Sensor_Analysis.csv", sep=";")
TABLE_3 <- read.csv("With_rescue_3_GFP_Sensor_Analysis.csv", sep=";")
TABLE_4 <- read.csv("With_rescue_4_GFP_Sensor_Analysis.csv", sep=";")
TABLE <- rbind(TABLE_1, TABLE_2, TABLE_3, TABLE_4)

# ANALYSIS GFP/RFP Signal level
TABLE <- TABLE[which(TABLE$Signal_type == "RFP"),]
TABLE$Genotype <- ordered(TABLE$Genotype, levels = c("miR430+/+", "miR430+/-", "miR430-/-"))

#p<-ggplot(TABLE, aes(x=Genotype, y=Signal_ratio)) + geom_violin(fill = "#999999") + geom_boxplot(width=0.1) + theme_bw() + theme(legend.position="none") + xlab("")
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(c(0,4))
#p <- p + ylab("Normalized GFP/RFP Signal") 

p<-ggplot(TABLE, aes(x=Genotype, y=Signal_ratio))
p <- p + theme_bw() + theme(legend.position="none") + xlab("") + geom_dotplot(binaxis='y', stackdir='center', dotsize = 1.2, binwidth = 0.05, fill = "#999999")
p <- p + stat_summary(fun=median, geom="point", shape=18, size=3, color="red")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(c(0,4)) + ylab("Normalized GFP/RFP Signal") 

pdf(paste0(output,".pdf"), width=4, height=4); p; dev.off()

kruskal.test(TABLE$Signal_ratio_log, TABLE$Genotype)
pairwise.wilcox.test(TABLE$Signal_ratio_log, TABLE$Genotype, p.adjust.method = "bonferroni")

table(TABLE$Genotype)
