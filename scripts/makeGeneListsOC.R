## ------------------------------------------------------------------

## generate gene lists after running 'gff3Info3.pl'

## Example: Rscript makeGeneListsOC.R MEL_info.txt MEL

## ------------------------------------------------------------------

args<-commandArgs(trailingOnly = TRUE)

genefile<-args[1] ## gene info
id<-args[2] ## species ID
prfx<-args[3] ## scaffold/chromosome prefix (to be removed)
pad<-args[4] ## pad scaffold ids (use only if only numbers follow prefix)


## avoid having gene positions in output as 3e+05 etc.
options("scipen"=100)

ge<-read.table(genefile,as.is=c(1,2,6),header=T)


## change info for strand from '+' to '1' and from '-' to '-1'
ge$strand<-gsub("+","1",ge$strand,fixed=TRUE)
ge$strand<-gsub("-","-1",ge$strand,fixed=TRUE)

## get scaffold IDs and order them (should work for most cases)
strReverse<-function(x){
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

scafsraw<-unique(ge$chr)
scafs<-sort(unique(gsub(prfx,"",scafsraw)))

if(length(scafsraw) != length(scafs)){
    ## better keep original names
    scafstab<-gsub("[0-9]+$","",scafsraw,perl=TRUE)
    scafsnum<-strReverse(gsub("^([0-9]+)[^0-9]+.*$","\\1",
                              strReverse(scafsraw),perl=TRUE))
    scafsnum[scafsnum==scafstab]<-NA
    scafstab<-cbind(scafstab,scafsnum)
    scafs<-scafsraw[order(scafstab[,1],as.numeric(scafstab[,2]),na.last=TRUE)]
} else {
    ## also change names in genefile
    ge$chr<-gsub(prfx,"",ge$chr)
    ## with padding
    if(pad == 1){
        ndigits<-floor(log10(max(as.numeric(scafs))))+1
        for(i in 1:length(ge$chr)){ ## change names in genefile
            curdigits<-floor(log10(as.numeric(ge$chr[i])))+1
            padding<-""
            if(curdigits < ndigits){
                for(j in 1:(ndigits-curdigits)){
                    padding<-paste(padding,"0",sep="")
                }
            }
            ge$chr[i]<-paste(padding,ge$chr[i],sep="")
        }
        ## change scaffold names
        scafs<-sort(unique(ge$chr))
    }
}


geneinfo<-data.frame(GENE_NAME=character(),CONTIG=character(),
                     START=numeric(),END=numeric(),STRAND=character())

for(i in 1:length(scafs)){
    tmp<-which(ge$chr==scafs[i])
    myord<-order(ge$mean[tmp]) ## sort by mean pos
    genes<-paste(id,"_",ge$name[tmp[myord]],sep="")
    tmpinfo<-data.frame(GENE_NAME=genes,
                        CONTIG=rep(scafs[i],length(tmp)),
                        START=ge$mean[tmp[myord]]-1,
                        END=ge$mean[tmp[myord]]+1,
                        STRAND=ge$strand[tmp[myord]])
    geneinfo<-rbind(geneinfo,tmpinfo)
}

tmp<-which(diff(geneinfo$START)<3 &
               diff(as.numeric(as.factor(geneinfo$CONTIG)))==0)
if(length(tmp)>0){
    ##head(geneinfo[sort(unique(c(tmp,tmp+1))),])

    overl<-geneinfo[sort(unique(c(tmp,tmp+1))),]
    overlgenes<-character()
    genes<-overl$GENE_NAME[1]
    pos1<-overl$START[1]
    cont1<-overl$CONTIG[1]
    for(i in 2:nrow(overl)){
        pos2<-overl$START[i]
        cont2<-overl$CONTIG[i]
        if(cont1==cont2 & (pos2-pos1)<3){
            genes<-paste(genes,overl$GENE_NAME[i])
        }else{
            overlgenes<-c(overlgenes,genes)
            genes<-overl$GENE_NAME[i]
        }
        pos1<-pos2
        cont1<-cont2
    }
    overlgenes<-c(overlgenes,genes)

    write.table(overlgenes,file=paste(id,"_overlapping.txt",sep=""),quote=F,
            col.names=F,row.names=F)
    print(paste("Genes with overlapping midpoints written to '",id,"_overlapping.txt'; these need to be checked and potentially be excluded from further ananlyses",sep=""),quote=FALSE)
}


write.table(geneinfo,file=paste(id,"_genome.txt",sep=""),quote=F,
            col.names=F,row.names=F)

