
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


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# Create a list of all spots
TABLE <- data.frame()

SPOTS <- read_excel("20221020_Analysis_WT_1.xlsx"); SPOTS$GROUP <- "WT_1"; SPOTS$GENOTYPE <- "Wild type"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20221020_Analysis_WT_2.xlsx"); SPOTS$GROUP <- "WT_2"; SPOTS$GENOTYPE <- "Wild type"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20221020_Analysis_WT_3.xlsx"); SPOTS$GROUP <- "WT_3"; SPOTS$GENOTYPE <- "Wild type"; TABLE <- rbind(TABLE, SPOTS)

SPOTS <- read_excel("20221020_Analysis_HOM_Rescue_1.xlsx"); SPOTS$GROUP <- "HOM_RESC_1"; SPOTS$GENOTYPE <- "Resc_HOM"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20221020_Analysis_HOM_Rescue_2.xlsx"); SPOTS$GROUP <- "HOM_RESC_2"; SPOTS$GENOTYPE <- "Resc_HOM"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20221020_Analysis_HOM_Rescue_3.xlsx"); SPOTS$GROUP <- "HOM_RESC_3"; SPOTS$GENOTYPE <- "Resc_HOM"; TABLE <- rbind(TABLE, SPOTS)

SPOTS <- read_excel("20230131_Analysis_WT_CDK9_1.xlsx"); SPOTS$GROUP <- "WT_CDK9_1"; SPOTS$GENOTYPE <- "Cdk9 overexpression"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20230202_Analysis_WT_CDK9_2.xlsx"); SPOTS$GROUP <- "WT_CDK9_2"; SPOTS$GENOTYPE <- "Cdk9 overexpression"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20230203_Analysis_WT_CDK9_3.xlsx"); SPOTS$GROUP <- "WT_CDK9_3"; SPOTS$GENOTYPE <- "Cdk9 overexpression"; TABLE <- rbind(TABLE, SPOTS)

SPOTS <- read_excel("20231030_Analysis_WT_MiR430_1.xlsx"); SPOTS$GROUP <- "WT_MiR430_1"; SPOTS$GENOTYPE <- "MiR430 injection"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20231115_Analysis_WT_MiR430_2.xlsx"); SPOTS$GROUP <- "WT_MiR430_2"; SPOTS$GENOTYPE <- "MiR430 injection"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20231213_Analysis_WT_MiR430_3.xlsx"); SPOTS$GROUP <- "WT_MiR430_3"; SPOTS$GENOTYPE <- "MiR430 injection"; TABLE <- rbind(TABLE, SPOTS)
SPOTS <- read_excel("20231216_Analysis_WT_MiR430_4.xlsx"); SPOTS$GROUP <- "WT_MiR430_4"; SPOTS$GENOTYPE <- "MiR430 injection"; TABLE <- rbind(TABLE, SPOTS)



TABLE <- melt(TABLE, id.vars=c("Condition","Date","Replicate","Series","Stage","Stage2","NucleusID", "GROUP", "GENOTYPE"))
TABLE$Stage <- factor(TABLE$Stage, levels = c("128", "256", "512", "1024"))
TABLE$GENOTYPE <- factor(TABLE$GENOTYPE, levels <- c("Wild type", "Resc_HOM", "Cdk9 overexpression", "MiR430 injection"))


#############################
##  miR430 NEGATIVE SPOTS  ##
#############################

TABLE_SUMMARY <- data_summary(TABLE, "value", c("Stage", "GENOTYPE", "variable"))
colnames(TABLE_SUMMARY) <- c("Stage", "GENOTYPE", "Type", "Mean", "SD")


View(TABLE_SUMMARY[which(TABLE_SUMMARY$Type == "mir430negative"),])


TABLE$GENOTYPE <- factor(TABLE$GENOTYPE, levels = c("Wild type", "Cdk9 overexpression", "Resc_HOM", "MiR430 injection"))


p <- ggplot(data=TABLE[which(TABLE$variable == "mir430negative"),], aes(x=GENOTYPE, y=value, fill = GENOTYPE)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("Stage")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)+ scale_fill_manual(values=c("#A0A0A0", "#377EB8", "#E41A1C", "#000000")) +
  facet_grid(cols = vars(Stage))
p

p_red <- ggplot(data=TABLE[which(TABLE$variable == "mir430negative" & TABLE$GENOTYPE %in% c("Wild type", "Resc_HOM", "Cdk9 overexpression")),], aes(x=GENOTYPE, y=value, fill = GENOTYPE)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("Stage")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)+ scale_fill_manual(values=c("#A0A0A0", "#377EB8", "#E41A1C")) +
  facet_grid(cols = vars(Stage))

pdf("Figure_ED10_mir430negative_boxplot.pdf", width=6, height=4); p; dev.off()
pdf("Figure_6G_mir430negative_boxplot.pdf", width=6, height=4); p_red; dev.off()


##################
##  STATISTICS  ##
##################

TABLE <- TABLE[which(TABLE$variable == "mir430negative"),]
TABLE$Condition2 <- paste0(TABLE$Condition, TABLE$Stage)
View(table(TABLE$Condition2))

St <- c("128", "256", "512", "1024")[1]

p2 <- ggplot(data=TABLE[which(TABLE$Stage == St),], aes(x=GENOTYPE, y=value, fill = GENOTYPE)) + geom_boxplot(position=position_dodge(1), outlier.size=0.5) + theme_bw() + xlab("Stage")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(aspect.ratio = 1)+ scale_fill_manual(values=c("#A0A0A0", "#377EB8", "#E41A1C", "#000000")) +
  facet_grid(cols = vars(Stage))
p2

kruskal.test(TABLE$value[which(TABLE$Stage == St)], TABLE$GENOTYPE[which(TABLE$Stage == St)])
pairwise.wilcox.test(TABLE$value[which(TABLE$Stage == St)], TABLE$GENOTYPE[which(TABLE$Stage == St)], p.adjust.method = "bonferroni")

