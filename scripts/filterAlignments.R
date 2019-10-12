## -------------------------------------------------------
## filter alignments for too much missing data

## Example: Rscript filterAlignments.R path/to/MAFFT
## -------------------------------------------------------

args<-commandArgs(trailingOnly = TRUE)

mywd<-args[1]

al<-read.table(paste(mywd,"alignInfo.txt",sep="/"),header=TRUE,as.is=1)

## note that proportions computed in bash are not rounded properly
## -> recalculate proportion of missing data in alignment ('-' or 'X')
propMiss<-(al$nGap+al$nX)/(al$len*al$nSp)
##hist(propMiss)

toKeep<-which(propMiss<=0.2)
##length(toKeep)
##hist(al$len[toKeep])

OGs<-paste0(al$OG[toKeep],".phy")

write.table(OGs,file=paste(mywd,"filteredAlignments.txt",sep="/"),
            col.names=FALSE,row.names=FALSE,quote=FALSE)
