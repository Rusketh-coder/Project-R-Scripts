#set directory to project folder/ shared drive for N2 vs WDR-5 YA-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/N2vsW-DE")
#read in N2vsW YA data "wdr5_vs_N2" csv file and assign data the label "N2vsW" 
#increase max print if over limit
N2vsW<-read.csv("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R Files/Data/wdr5_vs_N2.csv", header = TRUE, stringsAsFactors = FALSE, )


#check data for possible loss of rows and columns
dim(N2vsW)
head(N2vsW)
tail(N2vsW)
nrow(N2vsW)
names(N2vsW)
sapply(N2vsW, class)

#use complete.cases function to remove na value rows and double check with na.omit function
N2vsWnew1 <- N2vsW[complete.cases(N2vsW),]
N2vsWnew2 <- na.omit(N2vsWnew1)

#check data for possible loss of rows and columns then export
dim(N2vsWnew1)
head(N2vsWnew1)
tail(N2vsWnew1)
nrow(N2vsWnew1)
names(N2vsWnew1)
sapply(N2vsWnew1, class)

#export this data for futher analysis on DAVID (GO), WebGestalt (GO/GSEA)
write.table(N2vsWnew1, file="N2vsWnew2.csv", row.names = F, sep=",")

#subset DE data from N2vsW according to log2foldchange>=+2/-2, and p adjusted value<0.05 and name results "N2vsR1sigDE"
N2vsW1sigDEUpDown<-subset(N2vsW, N2vsW$log2FoldChange>=0.5849625007 & N2vsW$pval<=0.05|N2vsW$log2FoldChange<=-0.5849625007 & N2vsW$pval<=0.05, )
#check data
head(N2vsW1sigDEUpDown)
tail(N2vsW1sigDEUpDown)
nrow(N2vsW1sigDEUpDown)
names(N2vsW1sigDEUpDown)
sapply(N2vsW1sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsW1sigDEup<-subset(N2vsW, N2vsW$log2FoldChange>=0.5849625007 & N2vsW$pval<=0.05,)
#check data
dim(N2vsW1sigDEup)
head(N2vsW1sigDEup)
tail(N2vsW1sigDEup)
nrow(N2vsW1sigDEup)
names(N2vsW1sigDEup)
sapply(N2vsW1sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsW1sigDEupdownAlt<- N2vsW[abs(N2vsW$log2FoldChange)>=0.5849625007 & (N2vsW$pval<=0.05),]
#check data
dim(N2vs1sigDEupdown)
head(N2vsW1sigDEupdown)
tail(N2vsW1sigDEupdown)
nrow(N2vsW1sigDEupdown)
names(N2vsW1sigDEupdown)
sapply(N2vsW1sigDEupdown, class)

N2vsW1sigDEdown<-subset(N2vsW, N2vsW$log2FoldChange<=-0.5849625007 & N2vsW$pval<=0.05,)
#check data 
dim(N2vsW1sigDEdown)
head(N2vsW1sigDEdown)
tail(N2vsW1sigDEdown)
nrow(N2vsW1sigDEdown)
names(N2vsW1sigDEdown)
sapply(N2vsW1sigDEdown, class)



#subset DE data from N2vsW according to log2foldchange>=+2/-2, and p adjusted value<0.01 and name results "N2vsW2sigDE"
N2vsW2sigDEUpDown<-subset(N2vsW, N2vsW$log2FoldChange>=0.5849625007 & N2vsW$pval<=0.01|N2vsW$log2FoldChange<=-0.5849625007 & N2vsW$pval<=0.01, )
#check data
head(N2vsW2sigDEUpDown)
tail(N2vsW2sigDEUpDown)
nrow(N2vsW2sigDEUpDown)
names(N2vsW2sigDEUpDown)
sapply(N2vsW2sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsW2sigDEup<-subset(N2vsW, N2vsW$log2FoldChange>=0.5849625007 & N2vsW$pval<=0.01,)
#check data
dim(N2vsW2sigDEup)
head(N2vsW2sigDEup)
tail(N2vsW2sigDEup)
nrow(N2vsW2sigDEup)
names(N2vsW2sigDEup)
sapply(N2vsW2sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsW2sigDEupdownAlt<- N2vsW[abs(N2vsW$log2FoldChange)>=0.5849625007 & (N2vsW$pval<=0.01),]
#check data
dim(N2vsW2sigDEupdown)
head(N2vsW2sigDEupdown)
tail(N2vsW2sigDEupdown)
nrow(N2vsW2sigDEupdown)
names(N2vsW2sigDEupdown)
sapply(N2vsW2sigDEupdown, class)

N2vsW2sigDEdown<-subset(N2vsW, N2vsW$log2FoldChange<=-0.5849625007 & N2vsW$pval<=0.01,)
#check data 
dim(N2vsW2sigDEdown)
head(N2vsW2sigDEdown)
tail(N2vsW2sigDEdown)
nrow(N2vsW2sigDEdown)
names(N2vsW2sigDEdown)
sapply(N2vsW2sigDEdown, class)


#subset DE data from N2vsW according to log2foldchange>=+2/-2, and p adjusted value<0.05 and name results "N2vsW3sigDE"
N2vsW3sigDEUpDown<-subset(N2vsW, N2vsW$log2FoldChange>=1 & N2vsW$pval<=0.05|N2vsW$log2FoldChange<=-1 & N2vsW$pval<=0.05, )
#check data
head(N2vsW3sigDEUpDown)
tail(N2vsW3sigDEUpDown)
nrow(N2vsW3sigDEUpDown)
names(N2vsW3sigDEUpDown)
sapply(N2vsW3sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsW3sigDEup<-subset(N2vsW, N2vsW$log2FoldChange>=1 & N2vsW$pval<=0.05,)
#check data
dim(N2vsW3sigDEup)
head(N2vsW3sigDEup)
tail(N2vsW3sigDEup)
nrow(N2vsW3sigDEup)
names(N2vsW3sigDEup)
sapply(N2vsW3sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsW3sigDEupdownAlt<- N2vsW[abs(N2vsW$log2FoldChange)>=1 & (N2vsW$pval<=0.05),]
#check data
dim(N2vsW3sigDEupdown)
head(N2vsW3sigDEupdown)
tail(N2vsW3sigDEupdown)
nrow(N2vsW3sigDEupdown)
names(N2vsW3sigDEupdown)
sapply(N2vsW3sigDEupdown, class)

N2vsW3sigDEdown<-subset(N2vsW, N2vsW$log2FoldChange<=-1 & N2vsW$pval<=0.05,)
#check data 
dim(N2vsW3sigDEdown)
head(N2vsW3sigDEdown)
tail(N2vsW3sigDEdown)
nrow(N2vsW3sigDEdown)
names(N2vsW3sigDEdown)
sapply(N2vsW3sigDEdown, class)



#subset DE data from N2vsW according to log2foldchange>=+2/-2, and p adjusted value<0.01 and name results "N2vsW4sigDE"
N2vsW4sigDEUpDown<-subset(N2vsW, N2vsW$log2FoldChange>=1 & N2vsW$pval<=0.01|N2vsW$log2FoldChange<=-1 & N2vsW$pval<=0.01, )
#check data
head(N2vsW4sigDEUpDown)
tail(N2vsW4sigDEUpDown)
nrow(N2vsW4sigDEUpDown)
names(N2vsW4sigDEUpDown)
sapply(N2vsW4sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsW4sigDEup<-subset(N2vsW, N2vsW$log2FoldChange>=1 & N2vsW$pval<=0.01,)
#check data
dim(N2vsW4sigDEup)
head(N2vsW4sigDEup)
tail(N2vsW4sigDEup)
nrow(N2vsW4sigDEup)
names(N2vsW4sigDEup)
sapply(N2vsW4sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsW4sigDEupdownAlt<- N2vsW[abs(N2vsW$log2FoldChange)>=1 & (N2vsW$pval<=0.01),]
#check data
dim(N2vsW4sigDEupdown)
head(N2vsW4sigDEupdown)
tail(N2vsW4sigDEupdown)
nrow(N2vsW4sigDEupdown)
names(N2vsW4sigDEupdown)
sapply(N2vsW4sigDEupdown, class)

N2vsW4sigDEdown<-subset(N2vsW, N2vsW$log2FoldChange<=-1 & N2vsW$pval<=0.01,)
#check data 
dim(N2vsW4sigDEdown)
head(N2vsW4sigDEdown)
tail(N2vsW4sigDEdown)
nrow(N2vsW4sigDEdown)
names(N2vsW4sigDEdown)
sapply(N2vsW4sigDEdown, class)


#subset DE data from N2vsW according to log2foldchange>=+2/-2, and p adjusted value<0.001 and name results "N2vsW5sigDE"
N2vsW5sigDEUpDown<-subset(N2vsW, N2vsW$log2FoldChange>=1.584962501 & N2vsW$pval<=0.001|N2vsW$log2FoldChange<=-1.584962501 & N2vsW$pval<=0.001, )
#check data
head(N2vsW5sigDEUpDown)
tail(N2vsW5sigDEUpDown)
nrow(N2vsW5sigDEUpDown)
names(N2vsW5sigDEUpDown)
sapply(N2vsW5sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsW5sigDEup<-subset(N2vsW, N2vsW$log2FoldChange>=1.584962501 & N2vsW$pval<=0.001,)
#check data
dim(N2vsW5sigDEup)
head(N2vsW5sigDEup)
tail(N2vsW5sigDEup)
nrow(N2vsW5sigDEup)
names(N2vsW5sigDEup)
sapply(N2vsW5sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsW5sigDEupdownAlt<- N2vsW[abs(N2vsW$log2FoldChange)>=1.584962501 & (N2vsW$pval<=0.001),]
#check data
dim(N2vsW5sigDEupdown)
head(N2vsW5sigDEupdown)
tail(N2vsW5sigDEupdown)
nrow(N2vsW5sigDEupdown)
names(N2vsW5sigDEupdown)
sapply(N2vsW5sigDEupdown, class)

N2vsW5sigDEdown<-subset(N2vsW, N2vsW$log2FoldChange<=-1.584962501 & N2vsW$pval<=0.001,)
#check data 
dim(N2vsW5sigDEdown)
head(N2vsW5sigDEdown)
tail(N2vsW5sigDEdown)
nrow(N2vsW5sigDEdown)
names(N2vsW5sigDEdown)
sapply(N2vsW5sigDEdown, class)


#################################################################################################################################################

#set directory to that Dropbox file for N2vsR-DE Project
setwd("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R files/N2vsR-DE")
#read in N2vsR YA data "rbbp5_vs_N2" csv file and assign data the label "N2vsR" 
#increase max print if over limit
N2vsR<-read.csv("C:/Users/96nor/Dropbox/UG BioInfo 2019-20/Charles R Files/Data/rbbp5_vs_N2.csv", header = TRUE, stringsAsFactors = FALSE)

#check data for possible loss of rows and columns
dim(N2vsR)
head(N2vsR)
tail(N2vsR)
nrow(N2vsR)
names(N2vsR)
sapply(N2vsR, class)

#use complete.cases function to remove na value rows and double check with na.omit function
N2vsRnew1 <- N2vsR[complete.cases(N2vsR),]
N2vsRnew2 <- na.omit(N2vsRnew1)

#check data for possible loss of rows and columns then export
dim(N2vsRnew1)
head(N2vsRnew1)
tail(N2vsRnew1)
nrow(N2vsRnew1)
names(N2vsRnew1)
sapply(N2vsRnew1, class)

#export this data for futher analysis on DAVID (GO), WebGestalt (GO/GSEA)
write.table(N2vsRnew1, file="N2vsRnew2.csv", row.names = F, sep=",")

#subset DE data from N2vsR according to log2foldchange>=+2/-2, and p adjusted value<0.05 and name results "N2vsR1sigDE"
N2vsR1sigDEUpDown<-subset(N2vsR, N2vsR$log2FoldChange>=0.5849625007 & N2vsR$pval<=0.05|N2vsR$log2FoldChange<=-0.5849625007 & N2vsR$pval<=0.05, )
#check data
head(N2vsR1sigDEUpDown)
tail(N2vsR1sigDEUpDown)
nrow(N2vsR1sigDEUpDown)
names(N2vsR1sigDEUpDown)
sapply(N2vsR1sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsR1sigDEup<-subset(N2vsR, N2vsR$log2FoldChange>=0.5849625007 & N2vsR$pval<=0.05,)
#check data
dim(N2vsR1sigDEup)
head(N2vsR1sigDEup)
tail(N2vsR1sigDEup)
nrow(N2vsR1sigDEup)
names(N2vsR1sigDEup)
sapply(N2vsR1sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsR1sigDEupdownAlt<- N2vsR[abs(N2vsR$log2FoldChange)>=0.5849625007 & (N2vsR$pval<=0.05),]
#check data
dim(N2vsR1sigDEupdown)
head(N2vsR1sigDEupdown)
tail(N2vsR1sigDEupdown)
nrow(N2vsR1sigDEupdown)
names(N2vsR1sigDEupdown)
sapply(N2vsR1sigDEupdown, class)

N2vsR1sigDEdown<-subset(N2vsR, N2vsR$log2FoldChange<=-0.5849625007 & N2vsR$pval<=0.05,)
#check data 
dim(N2vsR1sigDEdown)
head(N2vsR1sigDEdown)
tail(N2vsR1sigDEdown)
nrow(N2vsR1sigDEdown)
names(N2vsR1sigDEdown)
sapply(N2vsR1sigDEdown, class)



#subset DE data from N2vsR according to log2foldchange>=+2/-2, and p adjusted value<0.01 and name results "N2vsR2sigDE"
N2vsR2sigDEUpDown<-subset(N2vsR, N2vsR$log2FoldChange>=0.5849625007 & N2vsR$pval<=0.01|N2vsR$log2FoldChange<=-0.5849625007 & N2vsR$pval<=0.01, )
#check data
head(N2vsR2sigDEUpDown)
tail(N2vsR2sigDEUpDown)
nrow(N2vsR2sigDEUpDown)
names(N2vsR2sigDEUpDown)
sapply(N2vsR2sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsR2sigDEup<-subset(N2vsR, N2vsR$log2FoldChange>=0.5849625007 & N2vsR$pval<=0.01,)
#check data
dim(N2vsR2sigDEup)
head(N2vsR2sigDEup)
tail(N2vsR2sigDEup)
nrow(N2vsR2sigDEup)
names(N2vsR2sigDEup)
sapply(N2vsR2sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsR2sigDEupdownAlt<- N2vsR[abs(N2vsR$log2FoldChange)>=0.5849625007 & (N2vsR$pval<=0.01),]
#check data
dim(N2vsR2sigDEupdown)
head(N2vsR2sigDEupdown)
tail(N2vsR2sigDEupdown)
nrow(N2vsR2sigDEupdown)
names(N2vsR2sigDEupdown)
sapply(N2vsR2sigDEupdown, class)

N2vsR2sigDEdown<-subset(N2vsR, N2vsR$log2FoldChange<=-0.5849625007 & N2vsR$pval<=0.01,)
#check data 
dim(N2vsR2sigDEdown)
head(N2vsR2sigDEdown)
tail(N2vsR2sigDEdown)
nrow(N2vsR2sigDEdown)
names(N2vsR2sigDEdown)
sapply(N2vsR2sigDEdown, class)


#subset DE data from N2vsR according to log2foldchange>=+2/-2, and p adjusted value<0.05 and name results "N2vsR3sigDE"
N2vsR3sigDEUpDown<-subset(N2vsR, N2vsR$log2FoldChange>=1 & N2vsR$pval<=0.05|N2vsR$log2FoldChange<=-1 & N2vsR$pval<=0.05, )
#check data
head(N2vsR3sigDEUpDown)
tail(N2vsR3sigDEUpDown)
nrow(N2vsR3sigDEUpDown)
names(N2vsR3sigDEUpDown)
sapply(N2vsR3sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsR3sigDEup<-subset(N2vsR, N2vsR$log2FoldChange>=1 & N2vsR$pval<=0.05,)
#check data
dim(N2vsR3sigDEup)
head(N2vsR3sigDEup)
tail(N2vsR3sigDEup)
nrow(N2vsR3sigDEup)
names(N2vsR3sigDEup)
sapply(N2vsR3sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsR3sigDEupdownAlt<- N2vsR[abs(N2vsR$log2FoldChange)>=1 & (N2vsR$pval<=0.05),]
#check data
dim(N2vsR3sigDEupdown)
head(N2vsR3sigDEupdown)
tail(N2vsR3sigDEupdown)
nrow(N2vsR3sigDEupdown)
names(N2vsR3sigDEupdown)
sapply(N2vsR3sigDEupdown, class)

N2vsR3sigDEdown<-subset(N2vsR, N2vsR$log2FoldChange<=-1 & N2vsR$pval<=0.05,)
#check data 
dim(N2vsR3sigDEdown)
head(N2vsR3sigDEdown)
tail(N2vsR3sigDEdown)
nrow(N2vsR3sigDEdown)
names(N2vsR3sigDEdown)
sapply(N2vsR3sigDEdown, class)



#subset DE data from N2vsR according to log2foldchange>=+2/-2, and p adjusted value<0.01 and name results "N2vsR4sigDE"
N2vsR4sigDEUpDown<-subset(N2vsR, N2vsR$log2FoldChange>=1 & N2vsR$pval<=0.01|N2vsR$log2FoldChange<=-1 & N2vsR$pval<=0.01, )
#check data
head(N2vsR4sigDEUpDown)
tail(N2vsR4sigDEUpDown)
nrow(N2vsR4sigDEUpDown)
names(N2vsR4sigDEUpDown)
sapply(N2vsR4sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsR4sigDEup<-subset(N2vsR, N2vsR$log2FoldChange>=1 & N2vsR$pval<=0.01,)
#check data
dim(N2vsR4sigDEup)
head(N2vsR4sigDEup)
tail(N2vsR4sigDEup)
nrow(N2vsR4sigDEup)
names(N2vsR4sigDEup)
sapply(N2vsR4sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsR4sigDEupdownAlt<- N2vsR[abs(N2vsR$log2FoldChange)>=1 & (N2vsR$pval<=0.01),]
#check data
dim(N2vsR4sigDEupdown)
head(N2vsR4sigDEupdown)
tail(N2vsR4sigDEupdown)
nrow(N2vsR4sigDEupdown)
names(N2vsR4sigDEupdown)
sapply(N2vsR4sigDEupdown, class)

N2vsR4sigDEdown<-subset(N2vsR, N2vsR$log2FoldChange<=-1 & N2vsR$pval<=0.01,)
#check data 
dim(N2vsR4sigDEdown)
head(N2vsR4sigDEdown)
tail(N2vsR4sigDEdown)
nrow(N2vsR4sigDEdown)
names(N2vsR4sigDEdown)
sapply(N2vsR4sigDEdown, class)




#subset DE data from N2vsR according to log2foldchange>=+2/-2, and p adjusted value<0.001 and name results "N2vsR5sigDE"
N2vsR5sigDEUpDown<-subset(N2vsR, N2vsR$log2FoldChange>=1.584962501 & N2vsR$pval<=0.001|N2vsR$log2FoldChange<=-1.584962501 & N2vsR$pval<=0.001, )
#check data
head(N2vsR5sigDEUpDown)
tail(N2vsR5sigDEUpDown)
nrow(N2vsR5sigDEUpDown)
names(N2vsR5sigDEUpDown)
sapply(N2vsR5sigDEUpDown, class)

#Alternate approach to subsetting according to Dr Poulin to differentiate sigDE up and down regulated genes
N2vsR5sigDEup<-subset(N2vsR, N2vsR$log2FoldChange>=1.584962501 & N2vsR$pval<=0.001,)
#check data
dim(N2vsR5sigDEup)
head(N2vsR5sigDEup)
tail(N2vsR5sigDEup)
nrow(N2vsR5sigDEup)
names(N2vsR5sigDEup)
sapply(N2vsR5sigDEup, class)

#Sorting up and down regulated genes seperately and together
N2vsR5sigDEupdownAlt<- N2vsR[abs(N2vsR$log2FoldChange)>=1.584962501 & (N2vsR$pval<=0.001),]
#check data
dim(N2vsR5sigDEupdown)
head(N2vsR5sigDEupdown)
tail(N2vsR5sigDEupdown)
nrow(N2vsR5sigDEupdown)
names(N2vsR5sigDEupdown)
sapply(N2vsR5sigDEupdown, class)

N2vsR5sigDEdown<-subset(N2vsR, N2vsR$log2FoldChange<=-1.584962501 & N2vsR$pval<=0.001,)
#check data 
dim(N2vsR5sigDEdown)
head(N2vsR5sigDEdown)
tail(N2vsR5sigDEdown)
nrow(N2vsR5sigDEdown)
names(N2vsR5sigDEdown)
sapply(N2vsR5sigDEdown, class)

#################################################################################################################################################


