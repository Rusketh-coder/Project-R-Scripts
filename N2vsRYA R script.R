#set directory to that Dropbox file for N2vsRYA-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R files/N2vsR-DE")
#read in N2vsRYA YA data "R_vs_N_DESeq2" csv file and assign data the label "N2vsRYA" 
#increase max print if over limit
N2vsRYA<-read.csv("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R Files/Data/R_vs_N_DESeq2.csv", header = TRUE, stringsAsFactors = FALSE)

#check data for possible loss of rows and columns
dim(N2vsRYA)
head(N2vsRYA)
tail(N2vsRYA)
nrow(N2vsRYA)
names(N2vsRYA)
sapply(N2vsRYA, class)

#load tidyr package to Use drop_na function from tidyverse package to format [] filtered data
library(tidyr)

#check data for possible loss of rows and columns then export
dim(N2vsRYAnew1)
head(N2vsRYAnew1)
tail(N2vsRYAnew1)
nrow(N2vsRYAnew1)
names(N2vsRYAnew1)
sapply(N2vsRYAnew1, class)

write.table(N2vsRYAnew1, file="N2vsRYAnew1.csv", row.names = F, sep=",")

#subset DE data from N2vsRYA according to log2foldchange>=+2/-2, and p adjusted value<=0.01 and name results "N2vsRYAsigDE"
N2vsRYAsigDEUpDown<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1 & N2vsRYA$pvalue<=0.01|N2vsRYA$log2FoldChange<=-1 & N2vsRYA$pvalue<=0.01,)
#check data
head(N2vsRYAsigDEUpDown)
tail(N2vsRYAsigDEUpDown)
nrow(N2vsRYAsigDEUpDown)
names(N2vsRYAsigDEUpDown)
sapply(N2vsRYAsigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsRYAsigDEup<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1 & N2vsRYA$pvalue<=0.01,)
#check data
dim(N2vsRYAsigDEup)
head(N2vsRYAsigDEup)
tail(N2vsRYAsigDEup)
nrow(N2vsRYAsigDEup)
names(N2vsRYAsigDEup)
sapply(N2vsRYAsigDEup, class)

N2vsRYAsigDEdown<-subset(N2vsRYA, N2vsRYA$log2FoldChange<=-1 & N2vsRYA$pvalue<=0.01,)
#check data 
dim(N2vsRYAsigDEdown)
head(N2vsRYAsigDEdown)
tail(N2vsRYAsigDEdown)
nrow(N2vsRYAsigDEdown)
names(N2vsRYAsigDEdown)
sapply(N2vsRYAsigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsRYAsigDEupdownAlt<- N2vsRYA[abs(N2vsRYA$log2FoldChange)>=1 & (N2vsRYA$pvalue<=0.01),]
N2vsRYAsigDEupdownAltNAfree<- N2vsRYAsigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYAsigDEupdownAltNAfree)
head(N2vsRYAsigDEupdownAltNAfree)
tail(N2vsRYAsigDEupdownAltNAfree)
nrow(N2vsRYAsigDEupdownAltNAfree)
names(N2vsRYAsigDEupdownAltNAfree)
sapply(N2vsRYAsigDEupdownAltNAfree, class)

#Export if numbers same for both methods
write.table(N2vsRYAsigDEupdownAlt, file="N2vsRYAsigDEupdownAlt.csv", row.names = F, sep=",")

N2vsRYAsigDEupAlt<-N2vsRYA[(N2vsRYA$log2FoldChange)>1 & (N2vsRYA$pvalue<=0.01),]
N2vsRYAsigDEupAltNAfree<-N2vsRYAsigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYAsigDEupAltNAfree)
head(N2vsRYAsigDEupAltNAfree)
tail(N2vsRYAsigDEupAltNAfree)
nrow(N2vsRYAsigDEupAltNAfree)
names(N2vsRYAsigDEupAltNAfree)
sapply(N2vsRYAsigDEupAltNAfree, class)

N2vsRYAsigDEdownAlt<-N2vsRYA[(N2vsRYA$log2FoldChange)<= -1 & (N2vsRYA$pvalue<=0.01),]
N2vsRYAsigDEdownAltNAfree<-N2vsRYAsigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYAsigDEdownAltNAfree)
head(N2vsRYAsigDEdownAltNAfree)
tail(N2vsRYAsigDEdownAltNAfree)
nrow(N2vsRYAsigDEdownAltNAfree)
names(N2vsRYAsigDEdownAltNAfree)
sapply(N2vsRYAsigDEdownAltNAfree, class)

#check if sum of upreg dtaset and downreg dataset are equivalent to the updownreg dataset
2256+692


#export 3 subsets as csv files without row names
write.table(N2vsRYAsigDEup, file="N2vsRYAsigDEup.csv", row.names = F, sep=",")
write.table(N2vsRYAsigDEUpDown, file="N2vsRYAsigDEUpDown.csv", row.names = F, sep=",")
write.table(N2vsRYAsigDEdown, file="N2vsRYAsigDEdown.csv", row.names = F, sep=",")

write.table(N2vsRYAsigDEup, file="N2vsRYAsigDEup.txt", row.names = F, sep=" ")
write.table(N2vsRYAsigDEUpDown, file="N2vsRYAsigDEUpDown.txt", row.names = F, sep=",")
write.table(N2vsRYAsigDEdown, file="N2vsRYAsigDEdown.txt", row.names = F, sep=" ")

##possible directions: H. sapiens-C elegans homologue under certain conditions e.g. ("N2vsR-DE")disease , life stage etc 
#enhancedvolcano plot
#enhancedvolcano 
#Turn on ggplot2 and enhancedvolcano package
library(ggplot2)
library(EnhancedVolcano)

#create volc plot
EnhancedVolcano(N2vsRYA,
                lab = rownames(N2vsRYA), #required
                x = "log2FoldChange", 
                y = "pvalue", 
                xlim = c(-7, 7),
                ylim = c(0,12),
                FCcutoff = 1,
                pCutoff = 10e-2,
                pointSize = 1,
                labSize = 2.5, 
                selectLab = 0,
                title = "Differentialy expressed genes - RBBP-5 mutant vs N2 (Young Adult stage)",
                titleLabSize = 14,
                axisLabSize = 12,
                subtitle= "Fold change (FC) against p-value (P)",
                legendLabels=c('NS', '+2 or -2 FC', 'p-value<0.01', '+2/-2 FC and p-value<0.01'),
                legendPosition = "bottom",
                legendLabSize = 12,
                caption = 'Total genes=18163, mis regulated genes=2971 (UP =2276, DOWN=695)',
                captionLabSize = 10,
                col = c("grey30", "forestgreen", "royalblue", "red2"), 
                colAlpha = 1)