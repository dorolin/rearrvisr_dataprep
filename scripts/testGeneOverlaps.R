#!/usr/bin/env Rscript

## Example: Rscript testGeneOverlaps.R SPE_info.txt SPE_overlapping_any.txt

## -----------------------------------

args = commandArgs(trailingOnly=TRUE)

if(length(args)!=2){
    stop("Input and output files must be specified. Example: Rscript testGeneOverlaps.R SPE_info.txt SPE_overlapping_any.txt")
}

infile<-args[1]
outfile<-args[2]

## -----------------------------------

genes<-read.table(infile,header=TRUE,
                  colClasses=c("character","character","integer","integer",
                      "integer","character"))


genes<-genes[order(genes$chr,genes$strand,genes$start),]

## assuming overlap is allowed if strand is different
overlap<-character()

i<-1
pchr<-genes$chr[i]
pstrand<-genes$strand[i]
pend<-genes$end[i]
pname<-genes$name[i]
tmpoverlap<-character()

for(i in 2:nrow(genes)){
    if(genes$chr[i]==pchr & genes$strand[i]==pstrand){
        if(genes$start[i]<=pend){ ## overlap
            if(length(tmpoverlap)==0){
                tmpoverlap<-paste(pname,genes$name[i],sep="; ")
            }else{
                tmpoverlap<-paste(tmpoverlap,genes$name[i],sep="; ")
            }
            pend<-max(pend,genes$end[i])
        }else{ ## no overlap
            pend<-genes$end[i]
            pname<-genes$name[i]
            if(length(tmpoverlap)>0){
                overlap<-c(overlap,tmpoverlap)
                tmpoverlap<-character()
            }
        }
    }else{ ## other strand or chromosome
        pchr<-genes$chr[i]
        pstrand<-genes$strand[i]
        pend<-genes$end[i]
        pname<-genes$name[i]
        if(length(tmpoverlap)>0){
            overlap<-c(overlap,tmpoverlap)
            tmpoverlap<-character()
        }
    }
}
if(length(tmpoverlap)>0){
    overlap<-c(overlap,tmpoverlap)
}


write.table(overlap,file=outfile,quote=FALSE,col.names=FALSE,row.names=FALSE)
