####################
# PloidyVisualizer
###################
# Authors: Arthur Melo & Iago Hale
# Department of Agriculture, Nutrition, and Food Systems
# University of New Hampshire, Durham, NH, 03824, USA.
# September 2018
#
# This script parses the output of the downstream tool of GBS-SNP-CROP
# (Step 8, HetFreq format option selected) and generates a table of
# allele depth ratios at all heterozygous loci for all genotypes in the
# population. Based on this table, the script creates a histogram of
# allele depth ratios for each genotype, thus allowing a visual 
# inspection of the ploidy status of each individual.

setwd("/YOURPATHHERE") # SET THE PATH TO YOUR WORKING DIRECTORY 
data<-read.table(file="GSC.HetFreq.txt",sep="\t",h=T) # OPEN THE FILE GSC.HetFreq.txt from GBS-SNP-CROP Step 8

pdf("PloidyVisualizerOut.pdf",width=8,height=6) # Output can be found in file PloidyVisualizerOut.pdf
par(mar=c(5,5,3,0),cex.lab=1,cex.main=1,cex.axis=1)
for (i in 3:ncol(data)) {
  k<-data[,i]
  name<-colnames(data[i])
  k<-as.data.frame(na.omit(k))
  for (i in 1:nrow(k)) {
    if ( (k[i,1] > 0.95) || (k[i,1] < 0.05)) {
      k[i,1] <- NA
    } else {
      next
    }
  }
  colnames(k)<-"Freq"
  k<-na.omit(k)
  hist(k$Freq, xlab="Het allele depth ratio",main=name,probability=T,col="white",
       border="white",xlim=c(0,1))
  lines(density (k$Freq), col="tomato",lwd = 2)
}
dev.off()
