## -------------------------------------------------------
## get orthologous groups that have genes for all species

## Example: Rscript getOGs.R path/to/OMA/Output
## -------------------------------------------------------

args<-commandArgs(trailingOnly = TRUE)

mywd<-args[1]

## read file with presence/absence indicators
dat<-read.table(paste(mywd,"PhyleticProfileOMAGroups.txt",sep="/"),
                header=TRUE,as.is=1)

fullgroups<-which(rowSums(dat[,-1])==ncol(dat)-1)

##length(fullgroups)

OGs<-dat[fullgroups,1]

OGs<-gsub("OMA0*","",OGs,perl=TRUE)

write.table(OGs,file=paste(mywd,"fullOGs.txt",sep="/"),
            col.names=FALSE,row.names=FALSE,quote=FALSE)

## -------------------------------------------------------
