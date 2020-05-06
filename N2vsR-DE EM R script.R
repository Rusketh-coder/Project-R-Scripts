#set directory to that Dropbox file for N2vsR-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R files/N2vsR-DE")
#read in N2vsR "rbbp5_vs_N2" csv file and assign data the label "N2vsR" 
#increase max print if over limit
N2vsR<-read.csv("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R Files/Data/rbbp5_vs_N2.csv", header = TRUE, stringsAsFactors = FALSE)
options(max.print=1e6)
#check data for possible loss of rows and columns
dim(N2vsR)
head(N2vsR)
tail(N2vsR)
nrow(N2vsR)
names(N2vsR)
sapply(N2vsR, class)
#subset DE data from N2vsR according to log2foldchange>=+1/-1, and p val<0.01 and name results "N2vsRsigDE"
N2vsRsigDEUpDown<-subset(N2vsR, N2vsR$log2FoldChange>=1 & N2vsR$pval<=0.01|N2vsR$log2FoldChange<=-1 & N2vsR$pval<=0.01, )
#check data
head(N2vsRsigDEUpDown)
tail(N2vsRsigDEUpDown)
nrow(N2vsRsigDEUpDown)
names(N2vsRsigDEUpDown)
sapply(N2vsRsigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsRsigDEup<-N2vsR[(N2vsR$log2FoldChange)>=1 & (N2vsR$pval<=0.01),]
#check data
dim(N2vsRsigDEup)
head(N2vsRsigDEup)
tail(N2vsRsigDEup)
nrow(N2vsRsigDEup)
names(N2vsRsigDEup)
sapply(N2vsRsigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsRsigDEupdownGP<-subset(N2vsR, abs(N2vsR$log2FoldChange)>=1 & N2vsR$pval<=0.01,)
#check data
dim(N2vsRsigDEupdownGP)
head(N2vsRsigDEupdownGP)
tail(N2vsRsigDEupdownGP)
nrow(N2vsRsigDEupdownGP)
names(N2vsRsigDEupdownGP)
sapply(N2vsRsigDEupdownGP, class)

N2vsRsigDEdown<-N2vsR[(N2vsR$log2FoldChange)<=-1 & (N2vsR$pval<0.01),]
#check data 
dim(N2vsRsigDEdown)
head(N2vsRsigDEdown)
tail(N2vsRsigDEdown)
nrow(N2vsRsigDEdown)
names(N2vsRsigDEdown)
sapply(N2vsRsigDEdown, class)
#check if sum of upreg dtaset and downreg dataset are equivalent to the updownreg dataset
1453+108

#export 3 subsets as csv files without row names
write.table(N2vsRsigDEup, file="N2vsRsigDEup.csv", row.names = F, sep=",")
write.table(N2vsRsigDEupdownGP, file="N2vsRsigDEupdownGP.csv", row.names = F, sep=",")
write.table(N2vsRsigDEdown, file="N2vsRsigDEdown.csv", row.names = F, sep=",")

write.table(N2vsRsigDEup, file="N2vsRsigDEup.txt", row.names = F, sep=" ")
write.table(N2vsRsigDEupdownGP, file="N2vsRsigDEupdownGP.txt", row.names = F, sep=" ")
write.table(N2vsRsigDEdown, file="N2vsRsigDEdown.txt", row.names = F, sep=" ")

#activate ggplot2 and import dataset for volcano plot creation
library(ggplot2)
library(EnhancedVolcano)

##possible directions: H. sapiens-C elegans homologue under certain conditions e.g. dis2vsR-DE")ease , life stage etc 
#enhancedvolcano plot
#enhancedvolcano 
#Turn on ggplot2 and enhancedvolcano package
library(ggplot2)
library(EnhancedVolcano)

#create volc plot
EnhancedVolcano(N2vsR,
                lab = rownames(N2vsR), #required
                x = "log2FoldChange", 
                y = "pval", 
                xlim = c(-7, 7),
                ylim = c(0,12),
                FCcutoff = 1,
                pCutoff = 10e-2,
                pointSize = 1,
                labSize = 2.5, 
                selectLab = 0,
                title = "Differentialy expressed genes - RBBP-5 mutant vs N2 (Embryonic Stage)",
                titleLabSize = 14,
                axisLabSize = 12,
                subtitle= "Fold change (FC) against p-value (P)",
                legendLabels=c('NS', '+2 or -2 FC', 'p-val<0.01', '+2/-2 FC and p-val<0.01'),
                legendPosition = "bottom",
                legendLabSize = 12,
                caption = 'Total genes=17009, mis regulated genes=1561 (UP =1453, DOWN=108)',
                captionLabSize = 10,
                col = c("grey30", "forestgreen", "royalblue", "red2"), 
                colAlpha = 1)



