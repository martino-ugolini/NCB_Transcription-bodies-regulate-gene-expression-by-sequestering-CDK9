library(ggplot2)
library("ggpubr")

#http://www.sthda.com/english/wiki/two-way-anova-test-in-r

rm(list = ls())

min_per_frame <- 10

# LOAD DATA NO RESCUE
output <- "Summary_No_Rescue"
TABLE_1 <- read.csv("20200903_no_rescue_1.csv", sep=";")
TABLE_2 <- read.csv("20200907_no_rescue_2.csv", sep=";")
TABLE_6 <- read.csv("20201006_no_rescue_6.csv", sep=";")
TABLE <- rbind(TABLE_1, TABLE_2, TABLE_6)

# LOAD DATA WITH RESCUE
output <- "Summary_With_Rescue"
TABLE_1 <- read.csv("20200909_with_rescue_1.csv", sep=";")
TABLE_4 <- read.csv("20201013_with_rescue_4.csv", sep=";")
TABLE_5 <- read.csv("20201021_with_rescue_5.csv", sep=";")
TABLE_7 <- read.csv("20201107_with_rescue_7.csv", sep=";")
TABLE <- rbind(TABLE_1, TABLE_4, TABLE_5, TABLE_7)

TABLE$Frame <- min_per_frame * TABLE$Frame
TABLE$Genotype <- ordered(TABLE$Genotype, levels = c("miR430 +/+", "miR430 +/-", "miR430 -/-"))

#TABLE <- TABLE[which(TABLE$Feature == "100% Epiboly"),]; output2 <- "100_epi_publication"; lim <- c(450, 800)
TABLE <- TABLE[which(TABLE$Feature == "Kupffer's vesicle"),]; output2 <- "kupf_ves_publicationi"; lim <- c(550, 900)

summary(TABLE$Genotype)

# MAKE GRAPH

#p<-ggplot(TABLE, aes(x=Genotype, y=Frame)) + geom_violin(fill = "#999999") + geom_boxplot(width=0.1) + theme_bw() + theme(legend.position="none") + xlab("")
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + ylim(lim) + ylab("Time [min]")

p<-ggplot(TABLE, aes(x=Genotype, y=Frame))
p <- p + theme_bw() + theme(legend.position="none") + xlab("") + geom_dotplot(binaxis='y', stackdir='center', dotsize = 1.2, binwidth = 5, fill = "#999999")
p <- p + stat_summary(fun=median, geom="point", shape=18, size=3, color="red")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(c(lim)) + ylab("Time [min]") 

pdf(paste0(output,"_", output2,".pdf"), width=4, height=4); p; dev.off()




kruskal.test(TABLE$Frame, TABLE$Genotype)
pairwise.wilcox.test(TABLE$Frame, TABLE$Genotype, p.adjust.method = "bonferroni")
p






