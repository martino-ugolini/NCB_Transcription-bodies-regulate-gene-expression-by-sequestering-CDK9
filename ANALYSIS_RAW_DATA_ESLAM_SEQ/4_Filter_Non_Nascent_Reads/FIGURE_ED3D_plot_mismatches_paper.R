# plot mismatches 
# we follow same nomenclature/structure as GRAND-SLAM for comparison

# Script adapted from PulseR package (see Materials and Methods)
# This Script can be run after part 4A has been executed (and mismatch_infos folder has been created)


library(ggplot2)
library(plyr)
library(purrr)
library(data.table)
library(wesanderson)

rm(list=ls())
loc <-  "./PART_A/mismatch_infos/"

print("Setting the path of the folder in which the *.mismatches.tab.gz files are located.")
which <- TRUE
basename <- "mismatches.pdf"

# first construct full dataframe
mismatches <- map(list.files(loc, pattern = 'mismatches.tab.gz$', full.names = T), read.table, head = TRUE)
names <- vapply(strsplit(list.files(loc, pattern = 'mismatches.tab.gz$'), ".", fixed = TRUE), "[", "", 1)
names(mismatches) <- names
for(i in seq_along(mismatches)) {
    mismatches[[i]]$Condition <- names[i]
}
all.mismatches <- do.call("rbind", mismatches)

if (which) {
    used <- map(list.files(loc, pattern = 'mismatches-used.tab.gz$', full.names = T), read.table, head = TRUE)
    names <- vapply(strsplit(list.files(loc, pattern = 'mismatches-used.tab.gz$'), ".", fixed = TRUE), "[", "", 1)
    names(used) <- names
    for(i in seq_along(used)) {
        used[[i]]$Condition <- names[i]
        used[[i]] <- used[[i]] %>% dplyr::rename(Used=Mismatches)
    }
    all.used <- do.call("rbind", used)


    fields = colnames(all.used)
    fields = fields[fields!='Used']

    all.mismatches <- merge(all.mismatches, all.used, by=fields, all = T)
    # reorder
    o <- order(all.mismatches$Condition, all.mismatches$Orientation, all.mismatches$Genomic, all.mismatches$Read)
    all.mismatches <- all.mismatches[o,]
    
    all.mismatches <- all.mismatches %>%
                        dplyr::mutate(Used = dplyr::coalesce(Used ,Mismatches)) %>%
                        dplyr::select(-Mismatches) %>% dplyr::rename(Mismatches=Used)
    basename <- "mismatches-used.pdf"
}

# prep dataframe
all.mismatches$Mismatch <- paste0(all.mismatches$Genomic,"->",all.mismatches$Read)
all.mismatches$Mismatch <- factor(all.mismatches$Mismatch,levels=unique(all.mismatches$Mismatch))
all.mismatches$Rate <- all.mismatches$Mismatches/all.mismatches$Coverage
all.mismatches$se <- sqrt(all.mismatches$Rate*(1-all.mismatches$Rate)/all.mismatches$Coverage)
all.mismatches$Condition <- factor(all.mismatches$Condition,as.character(unique(all.mismatches$Condition)))

TABLE <- all.mismatches
TABLE$Genotype <- map(strsplit(as.character(TABLE$Condition), split = "_"), 1)
TABLE$Stage2 <- map(strsplit(as.character(TABLE$Condition), split = "_"), 2)
TABLE$Stage[TABLE$Stage2 %in% c("A1", "A2", "A3", "A4")] <- "256c"
TABLE$Stage[TABLE$Stage2 %in% c("B1", "B2", "B3", "B4")] <- "1000c"
TABLE$Stage[TABLE$Stage2 %in% c("C1", "C2", "C3", "C4")] <- "Oblong"
TABLE$Stage[TABLE$Stage2 %in% c("D1", "D2", "D3", "D4")] <- "Sphere"
TABLE$IAA <- "Treated"
TABLE$IAA[TABLE$Stage2 %in% c("A4", "B4", "C4", "D4")] <- "Not treated"
TABLE$Stage <- as.factor(unlist(TABLE$Stage))
TABLE$Genotype <- as.factor(unlist(TABLE$Genotype))
TABLE$IAA <- as.factor(TABLE$IAA)

for (or in c("First", "Second"))
{
  p <- ggplot(TABLE[TABLE$Orientation==or,], aes(x=Mismatch, y=Rate, fill=IAA)) + geom_boxplot(outlier.size = 0.05) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=wes_palette(n=2, name="Moonrise2")) + theme(aspect.ratio=1)
  pdf(file=paste0("ED3D_", or, "_mismatch_plot.pdf"),width=6,height=4); print(p); dev.off()
}



