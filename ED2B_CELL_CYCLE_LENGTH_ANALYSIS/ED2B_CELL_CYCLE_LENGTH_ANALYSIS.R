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
TABLE <- read_excel("Timelapse_analysis.xlsx")
TABLE$Condition <- factor(TABLE$Condition, levels = c("WT", "mir430", "mir430_rescue"))
TABLE$Replicate <- factor(TABLE$Replicate, levels = c("1", "2", "3"))

# QUALITY CHECK
all(TABLE$'256_start' == TABLE$'128_end' + 1)
all(TABLE$'512_start' == TABLE$'256_end' + 1)
all(TABLE$'1024_start' == TABLE$'512_end' + 1)
all(TABLE$Length_128 == (TABLE$'128_end' - TABLE$'128_start' + 1))
all(TABLE$Length_256 == (TABLE$'256_end' - TABLE$'256_start' + 1))
all(TABLE$Length_512 == (TABLE$'512_end' - TABLE$'512_start' + 1))
all(TABLE$Length_1024 == (TABLE$'1024_end' - TABLE$'1024_start' + 1))

TABLE <- melt(TABLE, id.vars=colnames(TABLE)[1:12])
TABLE$value <- TABLE$value * 0.5 #in minutes

TABLE$variable <- as.character(TABLE$variable)
TABLE$variable[TABLE$variable == "Length_128"] <- "128c"
TABLE$variable[TABLE$variable == "Length_256"] <- "256c"
TABLE$variable[TABLE$variable == "Length_512"] <- "512c"
TABLE$variable[TABLE$variable == "Length_1024"] <- "1024c"
TABLE$variable <- factor(TABLE$variable, levels = c("128c", "256c", "512c","1024c"))

p <- ggplot(data=TABLE, aes(x=Condition, y=value, fill = Condition)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("Stage")  + ylab("Cell cycle duration [min]") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)+ scale_fill_manual(values=c("#A0A0A0", "#377EB8", "#E41A1C")) + facet_grid(cols = vars(variable))
p
pdf(file="Timelapse_analysis.pdf",width=6,height=4); print(p); dev.off()

TABLE$Condtion2 <- paste(TABLE$Condition, TABLE$variable)
table(TABLE$Condtion2)


for (stage in c("128c", "256c", "512c", "1024c"))
{
  print(stage)
  print(kruskal.test(TABLE$value[which(TABLE$variable == stage)], TABLE$Condition[which(TABLE$variable == stage)]))
  print(pairwise.wilcox.test(TABLE$value[which(TABLE$variable == stage)], TABLE$Condition[which(TABLE$variable == stage)], p.adjust.method = "bonferroni"))
  print("------------------------------------------------")
}

#TABLE$summary <- paste0(TABLE$Condition, TABLE$variable)
#table(TABLE$summary)
#26≤n≤29


#p1 <- ggplot(data=TABLE[which(TABLE$variable == stage),], aes(x=Condition, y=value, fill = Condition)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("Stage")  +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)+ scale_fill_manual(values=c("#A0A0A0", "#377EB8", "#E41A1C")) + ggtitle(stage)
#p2 <- ggplot(data=TABLE[which(TABLE$variable == stage),], aes(x=Condition, y=value, fill = Replicate)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("Stage")  +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)+ scale_fill_manual(values=c("#A0A0A0", "#377EB8", "#E41A1C")) + ggtitle(stage)








