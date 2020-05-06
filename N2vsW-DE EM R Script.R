#set directory to that Dropbox file for N2vsW-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R files/N2vsW-DE")
#read in N2vsW "wdr5_vs_N2" csv file and assign data the label "N2vsW" 
#increase max print if over limit
N2vsW<-read.csv("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R Files/Data/wdr5_vs_N2.csv", header = TRUE, stringsAsFactors = FALSE)
options(max.print=1e6)
#check data for possible loss of rows and columns
dim(N2vsW)
head(N2vsW)
tail(N2vsW)
nrow(N2vsW)
names(N2vsW)
sapply(N2vsW, class)
#subset DE data from N2vsW according to log2foldchange>=+2/-2, and p adjusted value<0.01 and name results "N2vsWsigDE"
N2vsWsigDEUpDown<-N2vsW[abs(N2vsW$log2FoldChange)>=1 & (N2vsW$pval<0.01),]
#check data
head(N2vsWsigDEUpDown)
tail(N2vsWsigDEUpDown)
nrow(N2vsWsigDEUpDown)
names(N2vsWsigDEUpDown)
sapply(N2vsWsigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsWsigDEup<-N2vsW[(N2vsW$log2FoldChange)>=1 & (N2vsW$pval<0.01),]
#check data
dim(N2vsWsigDEup)
head(N2vsWsigDEup)
tail(N2vsWsigDEup)
nrow(N2vsWsigDEup)
names(N2vsWsigDEup)
sapply(N2vsWsigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsWsigDEupdown<- N2vsW[abs(N2vsW$log2FoldChange)>=1 & (N2vsW$pval<0.01),]
#check data
dim(N2vsWsigDEupdown)
head(N2vsWsigDEupdown)
tail(N2vsWsigDEupdown)
nrow(N2vsWsigDEupdown)
names(N2vsWsigDEupdown)
sapply(N2vsWsigDEupdown, class)

N2vsWsigDEdown<-N2vsW[(N2vsW$log2FoldChange) <=-1 & (N2vsW$pval<0.01),]
#check data 
dim(N2vsWsigDEdown)
head(N2vsWsigDEdown)
tail(N2vsWsigDEdown)
nrow(N2vsWsigDEdown)
names(N2vsWsigDEdown)
sapply(N2vsWsigDEdown, class)
#check if sum of upreg dtaset and downreg dataset are equivalent to the updownreg dataset
2973+1147

#export 3 subsets as csv files without row names
write.table(N2vsWsigDEup, file="N2vsWsigDEup.csv", row.names = F, sep=",")
write.table(N2vsWsigDEupdown, file="N2vsWsigDEupdown.csv", row.names = F, sep=",")
write.table(N2vsWsigDEdown, file="N2vsWsigDEdown.csv", row.names = F, sep=",")

write.table(N2vsWsigDEup, file="N2vsWsigDEup.txt", row.names = F, sep=" ")
write.table(N2vsWsigDEupdown, file="N2vsWsigDEupdown.txt", row.names = F, sep=" ")
write.table(N2vsWsigDEdown, file="N2vsWsigDEdown.txt", row.names = F, sep=" ")

#activate ggplot2 and enhanced volcano import dataset for first volcano plot creation
library(ggplot2)
library(EnhancedVolcano)
#create volc plot
EnhancedVolcano(N2vsW,
                lab = rownames(N2vsW),
                x = "log2FoldChange", 
                y = "pval", 
                xlim = c(-7, 7),
                ylim = c(0,12),
                FCcutoff = 1,
                pCutoff = 10e-2,
                pointSize = 1,
                labSize = 2.5, 
                selectLab = 0,
                title = "Differentialy expressed genes - WDR-5 mutant vs N2 (Embryonic Stage)",
                titleLabSize = 14,
                axisLabSize = 12,
                subtitle= "Fold change (FC) against p-value (P)",
                legendLabels=c('NS', '+2 or -2 FC', 'p-value<0.01', '+2/-2 FC and p-value<0.01'),
                legendPosition = "bottom",
                legendLabSize = 12,
                caption = 'Total genes=16927, mis regulated genes=4120 (UP =2973, DOWN=1147)',
                captionLabSize = 10,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 1)


