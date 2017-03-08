\name{ReadBamToPosition}
\alias{ReadBamToPosition}
\title{
    Extract reads information from a bamfile
}
\description{
    This function extracts reads information from a bamfile using the ScanBamParam function in the \pkg{Rsamtools} package.
}
\usage{
    ReadBamToPosition(bamfile = NULL, SNPchr = NULL, positions = NULL, interval=NULL)
}
\arguments{
    \item{bamfile}{
        A path of a bamfile.
    }
    \item{SNPchr}{
        A chromosome information.
    }
    \item{positions}{
        A data frame containing genomic coordinates of exon start and end that are used to extract reads aligned to a given exon.
    }
    \item{interval}{
        An length from start and end that are used extract reads aligned to a given exon.
    }
}

\value{
    This function returns the data frame including name, chromosome number, genomic coordinate, mapping quality score, cigar string, genomic sequence, and quality score for each read.
}
\author{
    Seonggyun Han, Younghee Lee
}
\seealso{
    \code{\link[Rsamtools]{ScanBamParam}}
}
\examples{
    data(bamfilestest)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    ASdb <- Splicingfinder(transdb,Ncor=1)
    ASdb <- ExonsCluster(ASdb,transdb)
    bamfiles <- rbind(samplebamfiles[,"path"])
    splicingInfo <- rbind(samplebamfiles[,"names"])
    #Total.splicingInfo <- ASdb@SplicingModel$"ES"
    #each.ES.re <- rbind(ES.fi.result[ES.fi.result[,"Index"] == "ES3",])
    #each.ranges <- rbind(unique(cbind(do.call(rbind,strsplit(each.ES.re[,"DownEX"],"-"))[,1],
    #    do.call(rbind,strsplit(each.ES.re[,"UpEX"],"-"))[,2])))
    #tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.ES.re[,"Nchr"]),each.ranges,0)
}