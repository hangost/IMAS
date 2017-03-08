ReadBamToPosition <- function(bamfile=NULL,SNPchr=NULL,positions=NULL,interval=NULL){
    names(positions) <- c("start","end")
    which.ragnes <- GRanges(seqnames = SNPchr,ranges=IRanges(start=as.integer(positions["start"])-as.integer(interval),end=as.integer(positions["end"])+as.integer(interval)))
    what.param <- c("qname","rname","pos","mapq","cigar","seq","qual")
    param <- ScanBamParam(which = which.ragnes, what = what.param, tag = "MD")
    readed.bam <- do.call("DataFrame",scanBam(bamfile,param=param))
    if (nrow(readed.bam) == 0){
        return (NULL)
    }
    cn <- do.call(rbind,strsplit(colnames(readed.bam),"[.]"))
    colnames(readed.bam) <- cn[,ncol(cn)]
    readed.bam.mat <- as.matrix(readed.bam)
    readed.bam.mat[,"rname"] <- as.character(readed.bam[,"rname"])
    readed.bam.mat[,"seq"] <- as.character(readed.bam[,"seq"])
    readed.bam.mat[,"qual"] <- as.character(readed.bam[,"qual"])
    readed.bam.mat <- matrix(unlist(readed.bam.mat),ncol=ncol(readed.bam))
    colnames(readed.bam.mat) <- colnames(readed.bam)
    return (readed.bam.mat)
}

