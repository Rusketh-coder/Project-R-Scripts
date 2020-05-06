#set directory to project folder/ shared drive for N2 vs WDR-5 YA-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/N2vsW-DE")
#read in N2vsWYA data "R_vs_N_DEsSeq2.csv" csv file and assign data the label "N2vsWYA" 
#increase max print if over limit
N2vsWYA<-read.csv("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R Files/Data/W_vs_N_DESeq2.csv", header = TRUE, stringsAsFactors = FALSE, )
      

#check data for possible loss of rows and columns
dim(N2vsWYA)
head(N2vsWYA)
tail(N2vsWYA)
nrow(N2vsWYA)
names(N2vsWYA)
sapply(N2vsWYA, class)

#load tidyr package to Use drop_na function from tidyverse package to format [] filtered data
library(tidyr)

#export this data for futher analysis on DAVID (GO), WebGestalt (GO/GSEA)
write.table(N2vsWYA, file="N2vsWYA.csv", row.names = F, sep=",")

#subset DE data from N2vsWYA according to log2foldchange>=+1.5/-1.5, and p adjusted value<=0.05 and name results "N2vsWYA1sigDE"
N2vsWYA1sigDEUpDown<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=0.5849625007 & N2vsWYA$pvalue<=0.05|N2vsWYA$log2FoldChange<=-0.5849625007 & N2vsWYA$pvalue<=0.05,)
#check data
head(N2vsWYA1sigDEUpDown)
tail(N2vsWYA1sigDEUpDown)
nrow(N2vsWYA1sigDEUpDown)
names(N2vsWYA1sigDEUpDown)
sapply(N2vsWYA1sigDEUpDown, class)

write.table(N2vsWYA1sigDEUpDown, file="N2vsWYA1sigDEUpDown.csv", row.names = F, sep=",")

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsWYA1sigDEup<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=0.5849625007 & N2vsWYA$pvalue<=0.05,)
#check data
dim(N2vsWYA1sigDEup)
head(N2vsWYA1sigDEup)
tail(N2vsWYA1sigDEup)
nrow(N2vsWYA1sigDEup)
names(N2vsWYA1sigDEup)
sapply(N2vsWYA1sigDEup, class)

N2vsWYA1sigDEdown<-subset(N2vsWYA, N2vsWYA$log2FoldChange<=-0.5849625007 & N2vsWYA$pvalue<=0.05,)
#check data 
dim(N2vsWYA1sigDEdown)
head(N2vsWYA1sigDEdown)
tail(N2vsWYA1sigDEdown)
nrow(N2vsWYA1sigDEdown)
names(N2vsWYA1sigDEdown)
sapply(N2vsWYA1sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsWYA1sigDEupdownAlt<- N2vsWYA[abs(N2vsWYA$log2FoldChange)>=0.5849625007 & (N2vsWYA$pvalue<=0.05),]
N2vsWYA1sigDEupdownAltNAfree<- N2vsWYA1sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA1sigDEupdownAltNAfree)
head(N2vsWYA1sigDEupdownAltNAfree)
tail(N2vsWYA1sigDEupdownAltNAfree)
nrow(N2vsWYA1sigDEupdownAltNAfree)
names(N2vsWYA1sigDEupdownAltNAfree)
sapply(N2vsWYA1sigDEupdownAltNAfree, class)

write.table(N2vsWYA1sigDEupdownAlt, file="N2vsWYA1sigDEupdownAlt.csv", row.names = F, sep=",")

N2vsWYA1sigDEupAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)>=0.5849625007 & (N2vsWYA$pvalue<=0.05),]
N2vsWYA1sigDEupAltNAfree<-N2vsWYA1sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA1sigDEupAltNAfree)
head(N2vsWYA1sigDEupAltNAfree)
tail(N2vsWYA1sigDEupAltNAfree)
nrow(N2vsWYA1sigDEupAltNAfree)
names(N2vsWYA1sigDEupAltNAfree)
sapply(N2vsWYA1sigDEupAltNAfree, class)

N2vsWYA1sigDEdownAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)<= -0.5849625007 & (N2vsWYA$pvalue<=0.05),]
N2vsWYA1sigDEdownAltNAfree<-N2vsWYA1sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA1sigDEdownAltNAfree)
head(N2vsWYA1sigDEdownAltNAfree)
tail(N2vsWYA1sigDEdownAltNAfree)
nrow(N2vsWYA1sigDEdownAltNAfree)
names(N2vsWYA1sigDEdownAltNAfree)
sapply(N2vsWYA1sigDEdownAltNAfree, class)

#subset DE data from N2vsWYA according to log2foldchange>=+1.5/-1.5, and p adjusted value<=0.01 and name results "N2vsWYA2vsigDE"
N2vsWYA2vsigDEUpDown<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=0.5849625007 & N2vsWYA$pvalue<=0.01|N2vsWYA$log2FoldChange<=-0.5849625007 & N2vsWYA$pvalue<=0.01,)
#check data
head(N2vsWYA2vsigDEUpDown)
tail(N2vsWYA2vsigDEUpDown)
nrow(N2vsWYA2vsigDEUpDown)
names(N2vsWYA2vsigDEUpDown)
sapply(N2vsWYA2vsigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsWYA2vsigDEup<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=0.5849625007 & N2vsWYA$pvalue<=0.01,)
#check data
dim(N2vsWYA2vsigDEup)
head(N2vsWYA2vsigDEup)
tail(N2vsWYA2vsigDEup)
nrow(N2vsWYA2vsigDEup)
names(N2vsWYA2vsigDEup)
sapply(N2vsWYA2vsigDEup, class)

N2vsWYA2vsigDEdown<-subset(N2vsWYA, N2vsWYA$log2FoldChange<=-0.5849625007 & N2vsWYA$pvalue<=0.01,)
#check data 
dim(N2vsWYA2vsigDEdown)
head(N2vsWYA2vsigDEdown)
tail(N2vsWYA2vsigDEdown)
nrow(N2vsWYA2vsigDEdown)
names(N2vsWYA2vsigDEdown)
sapply(N2vsWYA2vsigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsWYA2vsigDEupdownAlt<- N2vsWYA[abs(N2vsWYA$log2FoldChange)>=0.5849625007 & (N2vsWYA$pvalue<=0.01),]
N2vsWYA2vsigDEupdownAltNAfree<- N2vsWYA2vsigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA2vsigDEupdownAltNAfree)
head(N2vsWYA2vsigDEupdownAltNAfree)
tail(N2vsWYA2vsigDEupdownAltNAfree)
nrow(N2vsWYA2vsigDEupdownAltNAfree)
names(N2vsWYA2vsigDEupdownAltNAfree)
sapply(N2vsWYA2vsigDEupdownAltNAfree, class)

N2vsWYA2vsigDEupAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)>=0.5849625007 & (N2vsWYA$pvalue<=0.01),]
N2vsWYA2vsigDEupAltNAfree<-N2vsWYA2vsigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA2vsigDEupAltNAfree)
head(N2vsWYA2vsigDEupAltNAfree)
tail(N2vsWYA2vsigDEupAltNAfree)
nrow(N2vsWYA2vsigDEupAltNAfree)
names(N2vsWYA2vsigDEupAltNAfree)
sapply(N2vsWYA2vsigDEupAltNAfree, class)

N2vsWYA2vsigDEdownAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)<= -0.5849625007 & (N2vsWYA$pvalue<=0.01),]
N2vsWYA2vsigDEdownAltNAfree<-N2vsWYA2vsigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA2vsigDEdownAltNAfree)
head(N2vsWYA2vsigDEdownAltNAfree)
tail(N2vsWYA2vsigDEdownAltNAfree)
nrow(N2vsWYA2vsigDEdownAltNAfree)
names(N2vsWYA2vsigDEdownAltNAfree)
sapply(N2vsWYA2vsigDEdownAltNAfree, class)


#subset DE data from N2vsWYA according to log2foldchange>=+2/-2, and p adjusted value<=0.05 and name results "N2vsWYA3sigDE"
N2vsWYA3sigDEUpDown<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1 & N2vsWYA$pvalue<=0.05|N2vsWYA$log2FoldChange<=-1 & N2vsWYA$pvalue<=0.05,)
#check data
head(N2vsWYA3sigDEUpDown)
tail(N2vsWYA3sigDEUpDown)
nrow(N2vsWYA3sigDEUpDown)
names(N2vsWYA3sigDEUpDown)
sapply(N2vsWYA3sigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsWYA3sigDEup<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1 & N2vsWYA$pvalue<=0.05,)
#check data
dim(N2vsWYA3sigDEup)
head(N2vsWYA3sigDEup)
tail(N2vsWYA3sigDEup)
nrow(N2vsWYA3sigDEup)
names(N2vsWYA3sigDEup)
sapply(N2vsWYA3sigDEup, class)

N2vsWYA3sigDEdown<-subset(N2vsWYA, N2vsWYA$log2FoldChange<=-1 & N2vsWYA$pvalue<=0.05,)
#check data 
dim(N2vsWYA3sigDEdown)
head(N2vsWYA3sigDEdown)
tail(N2vsWYA3sigDEdown)
nrow(N2vsWYA3sigDEdown)
names(N2vsWYA3sigDEdown)
sapply(N2vsWYA3sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsWYA3sigDEupdownAlt<- N2vsWYA[abs(N2vsWYA$log2FoldChange)>=1 & (N2vsWYA$pvalue<=0.05),]
N2vsWYA3sigDEupdownAltNAfree<- N2vsWYA3sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA3sigDEupdownAltNAfree)
head(N2vsWYA3sigDEupdownAltNAfree)
tail(N2vsWYA3sigDEupdownAltNAfree)
nrow(N2vsWYA3sigDEupdownAltNAfree)
names(N2vsWYA3sigDEupdownAltNAfree)
sapply(N2vsWYA3sigDEupdownAltNAfree, class)

N2vsWYA3sigDEupAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)>1 & (N2vsWYA$pvalue<=0.05),]
N2vsWYA3sigDEupAltNAfree<-N2vsWYA3sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA3sigDEupAltNAfree)
head(N2vsWYA3sigDEupAltNAfree)
tail(N2vsWYA3sigDEupAltNAfree)
nrow(N2vsWYA3sigDEupAltNAfree)
names(N2vsWYA3sigDEupAltNAfree)
sapply(N2vsWYA3sigDEupAltNAfree, class)

N2vsWYA3sigDEdownAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)<= -1 & (N2vsWYA$pvalue<=0.05),]
N2vsWYA3sigDEdownAltNAfree<-N2vsWYA3sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA3sigDEdownAltNAfree)
head(N2vsWYA3sigDEdownAltNAfree)
tail(N2vsWYA3sigDEdownAltNAfree)
nrow(N2vsWYA3sigDEdownAltNAfree)
names(N2vsWYA3sigDEdownAltNAfree)
sapply(N2vsWYA3sigDEdownAltNAfree, class)


#subset DE data from N2vsWYA according to log2foldchange>=+2/-2, and p adjusted value<=0.01 and name results "N2vsWYA4sigDE"
N2vsWYA4sigDEUpDown<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1 & N2vsWYA$pvalue<=0.01|N2vsWYA$log2FoldChange<=-1 & N2vsWYA$pvalue<=0.01,)
#check data
head(N2vsWYA4sigDEUpDown)
tail(N2vsWYA4sigDEUpDown)
nrow(N2vsWYA4sigDEUpDown)
names(N2vsWYA4sigDEUpDown)
sapply(N2vsWYA4sigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsWYA4sigDEup<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1 & N2vsWYA$pvalue<=0.01,)
#check data
dim(N2vsWYA4sigDEup)
head(N2vsWYA4sigDEup)
tail(N2vsWYA4sigDEup)
nrow(N2vsWYA4sigDEup)
names(N2vsWYA4sigDEup)
sapply(N2vsWYA4sigDEup, class)

N2vsWYA4sigDEdown<-subset(N2vsWYA, N2vsWYA$log2FoldChange<=-1 & N2vsWYA$pvalue<=0.01,)
#check data 
dim(N2vsWYA4sigDEdown)
head(N2vsWYA4sigDEdown)
tail(N2vsWYA4sigDEdown)
nrow(N2vsWYA4sigDEdown)
names(N2vsWYA4sigDEdown)
sapply(N2vsWYA4sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsWYA4sigDEupdownAlt<- N2vsWYA[abs(N2vsWYA$log2FoldChange)>=1 & (N2vsWYA$pvalue<=0.01),]
N2vsWYA4sigDEupdownAltNAfree<- N2vsWYA4sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA4sigDEupdownAltNAfree)
head(N2vsWYA4sigDEupdownAltNAfree)
tail(N2vsWYA4sigDEupdownAltNAfree)
nrow(N2vsWYA4sigDEupdownAltNAfree)
names(N2vsWYA4sigDEupdownAltNAfree)
sapply(N2vsWYA4sigDEupdownAltNAfree, class)

#Export if numbers same for both methods
write.table(N2vsWYA4sigDEupdownAlt, file="N2vsWYA4sigDEupdownAlt.csv", row.names = F, sep=",")

N2vsWYA4sigDEupAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)>1 & (N2vsWYA$pvalue<=0.01),]
N2vsWYA4sigDEupAltNAfree<-N2vsWYA4sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA4sigDEupAltNAfree)
head(N2vsWYA4sigDEupAltNAfree)
tail(N2vsWYA4sigDEupAltNAfree)
nrow(N2vsWYA4sigDEupAltNAfree)
names(N2vsWYA4sigDEupAltNAfree)
sapply(N2vsWYA4sigDEupAltNAfree, class)

N2vsWYA4sigDEdownAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)<= -1 & (N2vsWYA$pvalue<=0.01),]
N2vsWYA4sigDEdownAltNAfree<-N2vsWYA4sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA4sigDEdownAltNAfree)
head(N2vsWYA4sigDEdownAltNAfree)
tail(N2vsWYA4sigDEdownAltNAfree)
nrow(N2vsWYA4sigDEdownAltNAfree)
names(N2vsWYA4sigDEdownAltNAfree)
sapply(N2vsWYA4sigDEdownAltNAfree, class)

#subset DE data from N2vsWYA according to log2foldchange>=+3/-3, and p adjusted value<=0.001 and name results "N2vsWYA5sigDE"
N2vsWYA5sigDEUpDown<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1.584962501 & N2vsWYA$pvalue<=0.001|N2vsWYA$log2FoldChange<=-1.584962501 & N2vsWYA$pvalue<=0.001,)
#check data
head(N2vsWYA5sigDEUpDown)
tail(N2vsWYA5sigDEUpDown)
nrow(N2vsWYA5sigDEUpDown)
names(N2vsWYA5sigDEUpDown)
sapply(N2vsWYA5sigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsWYA5sigDEup<-subset(N2vsWYA, N2vsWYA$log2FoldChange>=1.584962501 & N2vsWYA$pvalue<=0.001,)
#check data
dim(N2vsWYA5sigDEup)
head(N2vsWYA5sigDEup)
tail(N2vsWYA5sigDEup)
nrow(N2vsWYA5sigDEup)
names(N2vsWYA5sigDEup)
sapply(N2vsWYA5sigDEup, class)

N2vsWYA5sigDEdown<-subset(N2vsWYA, N2vsWYA$log2FoldChange<=-1.584962501 & N2vsWYA$pvalue<=0.001,)
#check data 
dim(N2vsWYA5sigDEdown)
head(N2vsWYA5sigDEdown)
tail(N2vsWYA5sigDEdown)
nrow(N2vsWYA5sigDEdown)
names(N2vsWYA5sigDEdown)
sapply(N2vsWYA5sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsWYA5sigDEupdownAlt<- N2vsWYA[abs(N2vsWYA$log2FoldChange)>=1.584962501 & (N2vsWYA$pvalue<=0.001),]
N2vsWYA5sigDEupdownAltNAfree<- N2vsWYA5sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA5sigDEupdownAltNAfree)
head(N2vsWYA5sigDEupdownAltNAfree)
tail(N2vsWYA5sigDEupdownAltNAfree)
nrow(N2vsWYA5sigDEupdownAltNAfree)
names(N2vsWYA5sigDEupdownAltNAfree)
sapply(N2vsWYA5sigDEupdownAltNAfree, class)

N2vsWYA5sigDEupAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)>1.584962501 & (N2vsWYA$pvalue<=0.001),]
N2vsWYA5sigDEupAltNAfree<-N2vsWYA5sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA5sigDEupAltNAfree)
head(N2vsWYA5sigDEupAltNAfree)
tail(N2vsWYA5sigDEupAltNAfree)
nrow(N2vsWYA5sigDEupAltNAfree)
names(N2vsWYA5sigDEupAltNAfree)
sapply(N2vsWYA5sigDEupAltNAfree, class)

N2vsWYA5sigDEdownAlt<-N2vsWYA[(N2vsWYA$log2FoldChange)<= -1.584962501 & (N2vsWYA$pvalue<=0.001),]
N2vsWYA5sigDEdownAltNAfree<-N2vsWYA5sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsWYA5sigDEdownAltNAfree)
head(N2vsWYA5sigDEdownAltNAfree)
tail(N2vsWYA5sigDEdownAltNAfree)
nrow(N2vsWYA5sigDEdownAltNAfree)
names(N2vsWYA5sigDEdownAltNAfree)
sapply(N2vsWYA5sigDEdownAltNAfree, class)





                
#################################################################################################################################################

#set directory to that Dropbox file for N2vsR-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R files/N2vsR-DE")
#read in N2vsR YA data "R_vs_N_DEsSeq2.csv" csv file and assign data the label "N2vsR" 
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
      
#export this data for futher analysis on DAVID (GO), WebGestalt (GO/GSEA)
write.table(N2vsRYA, file="N2vsRYA.csv", row.names = F, sep=",")
      
#subset DE data from N2vsRYA according to log2foldchange>=+1.5/-1.5, and p adjusted value<=0.05 and name results "N2vsRYA1sigDE"
N2vsRYA1sigDEUpDown<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=0.5849625007 & N2vsRYA$pvalue<=0.05|N2vsRYA$log2FoldChange<=-0.5849625007 & N2vsRYA$pvalue<=0.05,)
#check data
head(N2vsRYA1sigDEUpDown)
tail(N2vsRYA1sigDEUpDown)
nrow(N2vsRYA1sigDEUpDown)
names(N2vsRYA1sigDEUpDown)
sapply(N2vsRYA1sigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsRYA1sigDEup<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=0.5849625007 & N2vsRYA$pvalue<=0.05,)
#check data
dim(N2vsRYA1sigDEup)
head(N2vsRYA1sigDEup)
tail(N2vsRYA1sigDEup)
nrow(N2vsRYA1sigDEup)
names(N2vsRYA1sigDEup)
sapply(N2vsRYA1sigDEup, class)

N2vsRYA1sigDEdown<-subset(N2vsRYA, N2vsRYA$log2FoldChange<=-0.5849625007 & N2vsRYA$pvalue<=0.05,)
#check data 
dim(N2vsRYA1sigDEdown)
head(N2vsRYA1sigDEdown)
tail(N2vsRYA1sigDEdown)
nrow(N2vsRYA1sigDEdown)
names(N2vsRYA1sigDEdown)
sapply(N2vsRYA1sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsRYA1sigDEupdownAlt<- N2vsRYA[abs(N2vsRYA$log2FoldChange)>=0.5849625007 & (N2vsRYA$pvalue<=0.05),]
N2vsRYA1sigDEupdownAltNAfree<- N2vsRYA1sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA1sigDEupdownAltNAfree)
head(N2vsRYA1sigDEupdownAltNAfree)
tail(N2vsRYA1sigDEupdownAltNAfree)
nrow(N2vsRYA1sigDEupdownAltNAfree)
names(N2vsRYA1sigDEupdownAltNAfree)
sapply(N2vsRYA1sigDEupdownAltNAfree, class)

N2vsRYA1sigDEupAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)>=0.5849625007 & (N2vsRYA$pvalue<=0.05),]
N2vsRYA1sigDEupAltNAfree<-N2vsRYA1sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA1sigDEupAltNAfree)
head(N2vsRYA1sigDEupAltNAfree)
tail(N2vsRYA1sigDEupAltNAfree)
nrow(N2vsRYA1sigDEupAltNAfree)
names(N2vsRYA1sigDEupAltNAfree)
sapply(N2vsRYA1sigDEupAltNAfree, class)

N2vsRYA1sigDEdownAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)<= -0.5849625007 & (N2vsRYA$pvalue<=0.05),]
N2vsRYA1sigDEdownAltNAfree<-N2vsRYA1sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA1sigDEdownAltNAfree)
head(N2vsRYA1sigDEdownAltNAfree)
tail(N2vsRYA1sigDEdownAltNAfree)
nrow(N2vsRYA1sigDEdownAltNAfree)
names(N2vsRYA1sigDEdownAltNAfree)
sapply(N2vsRYA1sigDEdownAltNAfree, class)

#subset DE data from N2vsRYA according to log2foldchange>=+1.5/-1.5, and p adjusted value<=0.01 and name results "N2vsRYA2vsigDE"
N2vsRYA2vsigDEUpDown<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=0.5849625007 & N2vsRYA$pvalue<=0.01|N2vsRYA$log2FoldChange<=-0.5849625007 & N2vsRYA$pvalue<=0.01,)
#check data
head(N2vsRYA2vsigDEUpDown)
tail(N2vsRYA2vsigDEUpDown)
nrow(N2vsRYA2vsigDEUpDown)
names(N2vsRYA2vsigDEUpDown)
sapply(N2vsRYA2vsigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsRYA2vsigDEup<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=0.5849625007 & N2vsRYA$pvalue<=0.01,)
#check data
dim(N2vsRYA2vsigDEup)
head(N2vsRYA2vsigDEup)
tail(N2vsRYA2vsigDEup)
nrow(N2vsRYA2vsigDEup)
names(N2vsRYA2vsigDEup)
sapply(N2vsRYA2vsigDEup, class)

N2vsRYA2vsigDEdown<-subset(N2vsRYA, N2vsRYA$log2FoldChange<=-0.5849625007 & N2vsRYA$pvalue<=0.01,)
#check data 
dim(N2vsRYA2vsigDEdown)
head(N2vsRYA2vsigDEdown)
tail(N2vsRYA2vsigDEdown)
nrow(N2vsRYA2vsigDEdown)
names(N2vsRYA2vsigDEdown)
sapply(N2vsRYA2vsigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsRYA2vsigDEupdownAlt<- N2vsRYA[abs(N2vsRYA$log2FoldChange)>=0.5849625007 & (N2vsRYA$pvalue<=0.01),]
N2vsRYA2vsigDEupdownAltNAfree<- N2vsRYA2vsigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA2vsigDEupdownAltNAfree)
head(N2vsRYA2vsigDEupdownAltNAfree)
tail(N2vsRYA2vsigDEupdownAltNAfree)
nrow(N2vsRYA2vsigDEupdownAltNAfree)
names(N2vsRYA2vsigDEupdownAltNAfree)
sapply(N2vsRYA2vsigDEupdownAltNAfree, class)

N2vsRYA2vsigDEupAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)>=0.5849625007 & (N2vsRYA$pvalue<=0.01),]
N2vsRYA2vsigDEupAltNAfree<-N2vsRYA2vsigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA2vsigDEupAltNAfree)
head(N2vsRYA2vsigDEupAltNAfree)
tail(N2vsRYA2vsigDEupAltNAfree)
nrow(N2vsRYA2vsigDEupAltNAfree)
names(N2vsRYA2vsigDEupAltNAfree)
sapply(N2vsRYA2vsigDEupAltNAfree, class)

N2vsRYA2vsigDEdownAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)<= -0.5849625007 & (N2vsRYA$pvalue<=0.01),]
N2vsRYA2vsigDEdownAltNAfree<-N2vsRYA2vsigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA2vsigDEdownAltNAfree)
head(N2vsRYA2vsigDEdownAltNAfree)
tail(N2vsRYA2vsigDEdownAltNAfree)
nrow(N2vsRYA2vsigDEdownAltNAfree)
names(N2vsRYA2vsigDEdownAltNAfree)
sapply(N2vsRYA2vsigDEdownAltNAfree, class)


#subset DE data from N2vsRYA according to log2foldchange>=+2/-2, and p adjusted value<=0.05 and name results "N2vsRYA3sigDE"
N2vsRYA3sigDEUpDown<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1 & N2vsRYA$pvalue<=0.05|N2vsRYA$log2FoldChange<=-1 & N2vsRYA$pvalue<=0.05,)
#check data
head(N2vsRYA3sigDEUpDown)
tail(N2vsRYA3sigDEUpDown)
nrow(N2vsRYA3sigDEUpDown)
names(N2vsRYA3sigDEUpDown)
sapply(N2vsRYA3sigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsRYA3sigDEup<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1 & N2vsRYA$pvalue<=0.05,)
#check data
dim(N2vsRYA3sigDEup)
head(N2vsRYA3sigDEup)
tail(N2vsRYA3sigDEup)
nrow(N2vsRYA3sigDEup)
names(N2vsRYA3sigDEup)
sapply(N2vsRYA3sigDEup, class)

N2vsRYA3sigDEdown<-subset(N2vsRYA, N2vsRYA$log2FoldChange<=-1 & N2vsRYA$pvalue<=0.05,)
#check data 
dim(N2vsRYA3sigDEdown)
head(N2vsRYA3sigDEdown)
tail(N2vsRYA3sigDEdown)
nrow(N2vsRYA3sigDEdown)
names(N2vsRYA3sigDEdown)
sapply(N2vsRYA3sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsRYA3sigDEupdownAlt<- N2vsRYA[abs(N2vsRYA$log2FoldChange)>=1 & (N2vsRYA$pvalue<=0.05),]
N2vsRYA3sigDEupdownAltNAfree<- N2vsRYA3sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA3sigDEupdownAltNAfree)
head(N2vsRYA3sigDEupdownAltNAfree)
tail(N2vsRYA3sigDEupdownAltNAfree)
nrow(N2vsRYA3sigDEupdownAltNAfree)
names(N2vsRYA3sigDEupdownAltNAfree)
sapply(N2vsRYA3sigDEupdownAltNAfree, class)

N2vsRYA3sigDEupAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)>1 & (N2vsRYA$pvalue<=0.05),]
N2vsRYA3sigDEupAltNAfree<-N2vsRYA3sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA3sigDEupAltNAfree)
head(N2vsRYA3sigDEupAltNAfree)
tail(N2vsRYA3sigDEupAltNAfree)
nrow(N2vsRYA3sigDEupAltNAfree)
names(N2vsRYA3sigDEupAltNAfree)
sapply(N2vsRYA3sigDEupAltNAfree, class)

N2vsRYA3sigDEdownAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)<= -1 & (N2vsRYA$pvalue<=0.05),]
N2vsRYA3sigDEdownAltNAfree<-N2vsRYA3sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA3sigDEdownAltNAfree)
head(N2vsRYA3sigDEdownAltNAfree)
tail(N2vsRYA3sigDEdownAltNAfree)
nrow(N2vsRYA3sigDEdownAltNAfree)
names(N2vsRYA3sigDEdownAltNAfree)
sapply(N2vsRYA3sigDEdownAltNAfree, class)


#subset DE data from N2vsRYA according to log2foldchange>=+2/-2, and p adjusted value<=0.01 and name results "N2vsRYA4sigDE"
N2vsRYA4sigDEUpDown<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1 & N2vsRYA$pvalue<=0.01|N2vsRYA$log2FoldChange<=-1 & N2vsRYA$pvalue<=0.01,)
#check data
head(N2vsRYA4sigDEUpDown)
tail(N2vsRYA4sigDEUpDown)
nrow(N2vsRYA4sigDEUpDown)
names(N2vsRYA4sigDEUpDown)
sapply(N2vsRYA4sigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsRYA4sigDEup<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1 & N2vsRYA$pvalue<=0.01,)
#check data
dim(N2vsRYA4sigDEup)
head(N2vsRYA4sigDEup)
tail(N2vsRYA4sigDEup)
nrow(N2vsRYA4sigDEup)
names(N2vsRYA4sigDEup)
sapply(N2vsRYA4sigDEup, class)

N2vsRYA4sigDEdown<-subset(N2vsRYA, N2vsRYA$log2FoldChange<=-1 & N2vsRYA$pvalue<=0.01,)
#check data 
dim(N2vsRYA4sigDEdown)
head(N2vsRYA4sigDEdown)
tail(N2vsRYA4sigDEdown)
nrow(N2vsRYA4sigDEdown)
names(N2vsRYA4sigDEdown)
sapply(N2vsRYA4sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsRYA4sigDEupdownAlt<- N2vsRYA[abs(N2vsRYA$log2FoldChange)>=1 & (N2vsRYA$pvalue<=0.01),]
N2vsRYA4sigDEupdownAltNAfree<- N2vsRYA4sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA4sigDEupdownAltNAfree)
head(N2vsRYA4sigDEupdownAltNAfree)
tail(N2vsRYA4sigDEupdownAltNAfree)
nrow(N2vsRYA4sigDEupdownAltNAfree)
names(N2vsRYA4sigDEupdownAltNAfree)
sapply(N2vsRYA4sigDEupdownAltNAfree, class)

#Export if numbers same for both methods
write.table(N2vsRYA4sigDEupdownAltNAfree, file="N2vsRYA4sigDEupdownAltNAfree.csv", row.names = F, sep=",")


N2vsRYA4sigDEupAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)>1 & (N2vsRYA$pvalue<=0.01),]
N2vsRYA4sigDEupAltNAfree<-N2vsRYA4sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA4sigDEupAltNAfree)
head(N2vsRYA4sigDEupAltNAfree)
tail(N2vsRYA4sigDEupAltNAfree)
nrow(N2vsRYA4sigDEupAltNAfree)
names(N2vsRYA4sigDEupAltNAfree)
sapply(N2vsRYA4sigDEupAltNAfree, class)

N2vsRYA4sigDEdownAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)<= -1 & (N2vsRYA$pvalue<=0.01),]
N2vsRYA4sigDEdownAltNAfree<-N2vsRYA4sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA4sigDEdownAltNAfree)
head(N2vsRYA4sigDEdownAltNAfree)
tail(N2vsRYA4sigDEdownAltNAfree)
nrow(N2vsRYA4sigDEdownAltNAfree)
names(N2vsRYA4sigDEdownAltNAfree)
sapply(N2vsRYA4sigDEdownAltNAfree, class)

#subset DE data from N2vsRYA according to log2foldchange>=+3/-3, and p adjusted value<=0.001 and name results "N2vsRYA5sigDE"
N2vsRYA5sigDEUpDown<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1.584962501 & N2vsRYA$pvalue<=0.001|N2vsRYA$log2FoldChange<=-1.584962501 & N2vsRYA$pvalue<=0.001,)
#check data
head(N2vsRYA5sigDEUpDown)
tail(N2vsRYA5sigDEUpDown)
nrow(N2vsRYA5sigDEUpDown)
names(N2vsRYA5sigDEUpDown)
sapply(N2vsRYA5sigDEUpDown, class)

#The approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsRYA5sigDEup<-subset(N2vsRYA, N2vsRYA$log2FoldChange>=1.584962501 & N2vsRYA$pvalue<=0.001,)
#check data
dim(N2vsRYA5sigDEup)
head(N2vsRYA5sigDEup)
tail(N2vsRYA5sigDEup)
nrow(N2vsRYA5sigDEup)
names(N2vsRYA5sigDEup)
sapply(N2vsRYA5sigDEup, class)

N2vsRYA5sigDEdown<-subset(N2vsRYA, N2vsRYA$log2FoldChange<=-1.584962501 & N2vsRYA$pvalue<=0.001,)
#check data 
dim(N2vsRYA5sigDEdown)
head(N2vsRYA5sigDEdown)
tail(N2vsRYA5sigDEdown)
nrow(N2vsRYA5sigDEdown)
names(N2vsRYA5sigDEdown)
sapply(N2vsRYA5sigDEdown, class)

#My approach to achieving the same results as Dr Poulin in filtering according to the selcted confidence levels with [] function and dplyr drop_na function
N2vsRYA5sigDEupdownAlt<- N2vsRYA[abs(N2vsRYA$log2FoldChange)>=1.584962501 & (N2vsRYA$pvalue<=0.001),]
N2vsRYA5sigDEupdownAltNAfree<- N2vsRYA5sigDEupdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA5sigDEupdownAltNAfree)
head(N2vsRYA5sigDEupdownAltNAfree)
tail(N2vsRYA5sigDEupdownAltNAfree)
nrow(N2vsRYA5sigDEupdownAltNAfree)
names(N2vsRYA5sigDEupdownAltNAfree)
sapply(N2vsRYA5sigDEupdownAltNAfree, class)

N2vsRYA5sigDEupAlt<- N2vsRYA[(N2vsRYA$log2FoldChange)>1.584962501 & (N2vsRYA$pvalue<=0.001),]
N2vsRYA5sigDEupAltNAfree<-N2vsRYA5sigDEupAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA5sigDEupAltNAfree)
head(N2vsRYA5sigDEupAltNAfree)
tail(N2vsRYA5sigDEupAltNAfree)
nrow(N2vsRYA5sigDEupAltNAfree)
names(N2vsRYA5sigDEupAltNAfree)
sapply(N2vsRYA5sigDEupAltNAfree, class)

N2vsRYA5sigDEdownAlt<-N2vsRYA[(N2vsRYA$log2FoldChange)<= -1.584962501 & (N2vsRYA$pvalue<=0.001),]
N2vsRYA5sigDEdownAltNAfree<-N2vsRYA5sigDEdownAlt %>% drop_na(log2FoldChange)
#check data
dim(N2vsRYA5sigDEdownAltNAfree)
head(N2vsRYA5sigDEdownAltNAfree)
tail(N2vsRYA5sigDEdownAltNAfree)
nrow(N2vsRYA5sigDEdownAltNAfree)
names(N2vsRYA5sigDEdownAltNAfree)
sapply(N2vsRYA5sigDEdownAltNAfree, class)

      
      
    