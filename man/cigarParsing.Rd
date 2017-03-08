\name{cigarParsing}
\alias{cigarParsing}
\title{Parse a cigar string}
\description{
    Parse a cigar string in a data frame containing read information from the \code{\link{ReadBamToPosition}} function.
}
\usage{
    cigarParsing(samfile.reads = NULL)
}
\arguments{
    \item{samfile.reads}{
        A Data frame containing read information from the \code{\link{ReadBamToPosition}} function.
    }
}

\value{
    This function returns a data frame with an alignment length of each cigar string.
}
\references{
    Heng Li, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16): 2078?2079.
}
\author{
    Seonggyun Han, Younghee Lee
}
\seealso{
    \code{\link{ReadBamToPosition}}
}
\examples{
    data(bamfilestest)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    ASdb <- Splicingfinder(transdb)
    bamfiles <- rbind(samplebamfiles[,"path"])
    splicingInfo <- rbind(samplebamfiles[,"names"])
    #Total.splicingInfo <- ASdb@SplicingModel$"ES"
    #each.ES.re <- rbind(ES.fi.result[ES.fi.result[,"Index"] == "ES3",])
    #each.ranges <- rbind(unique(cbind(do.call(rbind,strsplit(each.ES.re[,"DownEX"],"-"))[,1],
    #        do.call(rbind,strsplit(each.ES.re[,"UpEX"],"-"))[,2])))
    #tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.ES.re[,"Nchr"]),each.ranges,0)
    #cigarParsing(rbind(tested.reads[tested.reads[,"qname"]=="test.reads2-14",]))
}