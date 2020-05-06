#set directory to project folder/ shared drive for N2vsWYA-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R files/N2vsW-DE")
#read in N2vsWYA YA data "R_vs_N_DESeq2" csv file and assign data the label "N2vsWYA" 
#increase max print if over limit
N2vsWYA<-read.csv("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R Files/Data/W_vs_N_DESeq2.csv", header = TRUE, stringsAsFactors = FALSE)

#check data for possible loss of rows and columns
dim(N2vsWYA)
head(N2vsWYA)
tail(N2vsWYA)
nrow(N2vsWYA)
names(N2vsWYA)
sapply(N2vsWYA, class)

#load tidyr package to Use drop_na function from tidyverse package to format [] filtered data
library(tidyr)

#check data for possible loss of rows and columns then export
dim(N2vsWYAnew1)
head(N2vsWYAnew1)
tail(N2vsWYAnew1)
nrow(N2vsWYAnew1)
names(N2vsWYAnew1)
sapply(N2vsWYAnew1, class)

write.table(N2vsWYAnew1, file="N2vsWYAnew1.csv", row.names = F, sep=",")

#subset DE data from N2vsWYA according to log2foldchange>=+2/-2, and p adjusted value<=0.01 and name results "N2vsWYAsigDE"
N2vsWYAsigDEUpDown<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1 & N2vsWYA$pvalue<=0.01|N2vsWYA$log2FoldChange<=-1 & N2vsWYA$pvalue<=0.01,)
#check data
head(N2vsWYAsigDEUpDown)
tail(N2vsWYAsigDEUpDown)
nrow(N2vsWYAsigDEUpDown)
names(N2vsWYAsigDEUpDown)
sapply(N2vsWYAsigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsWYAsigDEup<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1 & N2vsWYA$pvalue<=0.01,)
#check data
dim(N2vsWYAsigDEup)
head(N2vsWYAsigDEup)
tail(N2vsWYAsigDEup)
nrow(N2vsWYAsigDEup)
names(N2vsWYAsigDEup)
sapply(N2vsWYAsigDEup, class)

N2vsWYAsigDEdown<-subset(N2vsWYA, N2vsWYA$log2FoldChange<=-1 & N2vsWYA$pvalue<=0.01,)
#check data 
dim(N2vsWYAsigDEdown)
head(N2vsWYAsigDEdown)
tail(N2vsWYAsigDEdown)
nrow(N2vsWYAsigDEdown)
names(N2vsWYAsigDEdown)
sapply(N2vsWYAsigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsWYAsigDEupdownAlt<- N2vsWYA[abs(N2vsWYA$log2FoldChange)>=1 & (N2vsWYA$pvalue<=0.01),]
N2vsWYAsigDEupdownAltNAfree<- N2vsWYAsigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYAsigDEupdownAltNAfree)
head(N2vsWYAsigDEupdownAltNAfree)
tail(N2vsWYAsigDEupdownAltNAfree)
nrow(N2vsWYAsigDEupdownAltNAfree)
names(N2vsWYAsigDEupdownAltNAfree)
sapply(N2vsWYAsigDEupdownAltNAfree, class)

#Export if numbers same for both methods
write.table(N2vsWYAsigDEupdownAlt, file="N2vsWYAsigDEupdownAlt.csv", row.names = F, sep=",")

N2vsWYAsigDEupAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)>1 & (N2vsWYA$pvalue<=0.01),]
N2vsWYAsigDEupAltNAfree<-N2vsWYAsigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYAsigDEupAltNAfree)
head(N2vsWYAsigDEupAltNAfree)
tail(N2vsWYAsigDEupAltNAfree)
nrow(N2vsWYAsigDEupAltNAfree)
names(N2vsWYAsigDEupAltNAfree)
sapply(N2vsWYAsigDEupAltNAfree, class)

N2vsWYAsigDEdownAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)<= -1 & (N2vsWYA$pvalue<=0.01),]
N2vsWYAsigDEdownAltNAfree<-N2vsWYAsigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYAsigDEdownAltNAfree)
head(N2vsWYAsigDEdownAltNAfree)
tail(N2vsWYAsigDEdownAltNAfree)
nrow(N2vsWYAsigDEdownAltNAfree)
names(N2vsWYAsigDEdownAltNAfree)
sapply(N2vsWYAsigDEdownAltNAfree, class)

#check if sum of upreg dtaset and downreg dataset are equivalent to the updownreg dataset



#export 3 subsets as csv files without row names
write.table(N2vsWYAsigDEup, file="N2vsWYAsigDEup.csv", row.names = F, sep=",")
write.table(N2vsWYAsigDEUpDown, file="N2vsWYAsigDEUpDown.csv", row.names = F, sep=",")
write.table(N2vsWYAsigDEdown, file="N2vsWYAsigDEdown.csv", row.names = F, sep=",")

write.table(N2vsWYAsigDEup, file="N2vsWYAsigDEup.txt", row.names = F, sep=" ")
write.table(N2vsWYAsigDEUpDown, file="N2vsWYAsigDEUpDown.txt", row.names = F, sep=",")
write.table(N2vsWYAsigDEdown, file="N2vsWYAsigDEdown.txt", row.names = F, sep=" ")

##possible directions: H. sapiens-C elegans homologue under certain conditions e.g. ("N2vsW-DE")disease , life stage etc 
#enhancedvolcano plot
#enhancedvolcano 
#Turn on ggplot2 and enhancedvolcano package
library(ggplot2)
library(EnhancedVolcano)

#create volc plot
EnhancedVolcano(N2vsWYA,
                lab = rownames(N2vsWYA), #required
                x = "log2FoldChange", 
                y = "pvalue", 
                xlim = c(-7, 7),
                ylim = c(0,12),
                FCcutoff = 1,
                pCutoff = 10e-2,
                pointSize = 1,
                labSize = 2.5, 
                selectLab = 0,
                title = "Differentialy expressed genes - WDR-5 mutant vs N2 (Young Adult stage)",
                titleLabSize = 14,
                axisLabSize = 12,
                subtitle= "Fold change (FC) against p-value (P)",
                legendLabels=c('NS', '+2 or -2 FC', 'p-value<0.01', '+2/-2 FC and p-value<0.01'),
                legendPosition = "bottom",
                legendLabSize = 12,
                caption = 'Total genes=18163, mis regulated genes=1175 (UP =873, DOWN=302)',
                captionLabSize = 10,
                col = c("grey30", "forestgreen", "royalblue", "red2"), 
                colAlpha = 1)
