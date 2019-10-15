## --- beta version: tweak overlapping genome locations ---

## generate gene lists after running 'gff3Info3.pl'

## Example: Rscript makeGeneListsOC_beta.R MEL_info.txt MEL

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

midpt<-ge$start+(ge$end-ge$start)/2

for(i in 1:length(scafs)){
    tmp<-which(ge$chr==scafs[i])
    ## sort by mean pos, then by midpoint, then by start
    myord<-order(ge$mean[tmp],midpt[tmp],ge$start[tmp])
    genes<-paste(id,"_",ge$name[tmp[myord]],sep="")
    tmpinfo<-data.frame(GENE_NAME=genes,
                        CONTIG=rep(scafs[i],length(tmp)),
                        START=ge$mean[tmp[myord]]-1,
                        END=ge$mean[tmp[myord]]+1,
                        STRAND=ge$strand[tmp[myord]])
    ## check that no starts/ends are overlapping
    ovr<-which(diff(tmpinfo$START)<3)
    if(length(ovr)>0){
        ## adjust pos for ovr[k]+1
        ## entries in ovr might be adjacent, so use potentially
        ##  changed value for tmpinfo$START[ovr[k]] as reference
        ## if next entry is not adjacent, make sure that no new
        ##  conflict was introduced by adjustment
        for(k in 1:length(ovr)){
            g1<-tmpinfo[ovr[k],c(3,4)] ## START END gene 1
            tmpinfo[ovr[k]+1,c(3,4)]<-g1+3 ## adjust START END gene 2
            if(k<length(ovr)){
                if(ovr[k+1]-ovr[k]==1){
                    next
                }
            }
            nxt<-ovr[k]+1
            while(length(nxt)>0){
                if((nxt+1)>nrow(tmpinfo)){
                    break
                }
                if(diff(tmpinfo[nxt:(nxt+1),3])<3){
                    ## adjust START END next gene
                    tmpinfo[nxt+1,c(3,4)]<-tmpinfo[nxt,c(3,4)]+3
                    nxt<-nxt+1
                }else{
                    nxt<-numeric()
                }
            }
        }
    }
    ## tmpinfo[sort(unique(c(ovr,ovr+1))),]
    if(sum(diff(tmpinfo$START)<3)>0){
        stop(paste("Still overlapping positions for scaffold id",i))
    }
    geneinfo<-rbind(geneinfo,tmpinfo)
}

tmp<-which(diff(geneinfo$START)<3 &
               diff(as.numeric(as.factor(geneinfo$CONTIG)))==0)
if(length(tmp)>0){
    head(geneinfo[sort(unique(c(tmp,tmp+1))),])
}


write.table(geneinfo,file=paste(id,"_genome.txt",sep=""),quote=F,
            col.names=F,row.names=F)

