\name{SplicingReads}
\alias{SplicingReads}
\title{
    Count a junction and paired-end reads
}
\description{
    This function counts the reads that are mapped to two separate exons, mapped to either splice site of two exons (called junction reads) or within each of two exons (paired end reads).
}
\usage{
    SplicingReads(samfile = NULL, test.exon = NULL, 
        spli.jun = NULL, readsinfo = "paired", inse = 40)
}
\arguments{
    \item{samfile}{
        A data frame containing read information from the \code{\link{ReadBamToPosition}} function.
    }
    \item{test.exon}{
        A data frame containing an alternative target exon and their neighboring exons.
    }
    \item{spli.jun}{
        A data frame containing spliced junction information.
    }
    \item{readsinfo}{
        Information of RNA-seq types (single- or paired- end reads).
    }
    \item{inse}{
        An insert size
    }
}

\value{
    This function returns the list object providing counts the reads that are mapped to two separate exons, mapped to either splice site of two exons (called junction reads) or within each of two exons (paired end reads).
}
\author{
    Seonggyun Han, Younghee Lee
}
\examples{
    data(bamfilestest)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    #ASdb <- Splicingfinder(transdb,Ncor=1)
    #ASdb <- ExonsCluster(ASdb,transdb)
    #bamfiles <- rbind(samplebamfiles[,"path"])
    #splicingInfo <- rbind(samplebamfiles[,"names"])
    #Total.splicingInfo <- ASdb@SplicingModel$"ES"
    #each.ES.re <- rbind(ES.fi.result[ES.fi.result[,"Index"] == "ES3",])
    #each.ranges <- rbind(unique(cbind(do.call(rbind,strsplit(each.ES.re[,"DownEX"],"-"))[,1],
    #    do.call(rbind,strsplit(each.ES.re[,"UpEX"],"-"))[,2])))
    #tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.ES.re[,"Nchr"]),each.ranges,0)
    #group.1.spl <- c(split.splice(each.ES.re[,"Do_des"],each.ES.re[,"1st_des"]),
    #    split.splice(each.ES.re[,"1st_des"],each.ES.re[,"Up_des"]))
    #group.2.spl <- split.splice(each.ES.re[,"Do_des"],each.ES.re[,"Up_des"])
    #total.reads <- SplicingReads(tested.reads,each.ES.re[,c("DownEX","1stEX","UpEX")],
    #    c(group.1.spl,group.2.spl),"paired")
}