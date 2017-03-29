\name{RatioFromReads}
\alias{RatioFromReads}
\title{
    Calculate expression ratio (PSI) from bamfiles
}
\description{
    This function extracts reads information from bamfile using \pkg{Rsamtools} and calculates expression ratio (denoted as Percent Splice-In, PSI) of each alternatively spliced exon (i.e., exon skipping, intro retention, and 5- and 3- prime splice sites).
}
\usage{
    RatioFromReads(ASdb=NULL,Total.bamfiles=NULL,readsInfo=c("paired","single"),
        readLen=NULL,inserSize=NULL,minr=3,CalIndex=NULL,Ncor=1,out.dir=NULL)
}
\arguments{
    \item{ASdb}{
        An ASdb object including "SplicingModel" slot from the \code{\link{Splicingfinder}} function.
    }
    \item{Total.bamfiles}{
        A data frame containing the path and name of a bamfile from RNA-seq
    }
    \item{readsInfo}{
        Information of RNA-seq types (single- or paired-end reads)
    }
    \item{readLen}{
        The read length
    }
    \item{inserSize}{
        The insert size between paired-end reads.
    }
    \item{minr}{
        A minimum number of testable reads mapping to a given exon.
    }
    \item{CalIndex}{
        An index number in the ASdb object which will be tested in this function. 
    }
    \item{Ncor}{
        The number of cores for multi-threads.
    }
    \item{out.dir}{
        An output directory.
    }
}

\value{
    ASdb with the slot (labeled by "Ratio") containing results from the the \code{\link{RatioFromReads}} function. The "Ratio" slot contains a list object and each element of the list object returns the results assigned to three elements, which is of each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). Three elements are as follows;
    \item{ES}{A data frame for the result of Exon skipping, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), 1stEX (alternatively spliced target exon), 2ndEX (second alternatively spliced target exon which is the other one of the mutually exclusive spliced exons), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), and names of individuals.}
    \item{ASS}{A data frame for the result of Alternative splice sites, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), ShortEX (shorter spliced target exon), LongEX (longer spliced target exon), NeighborEX (neighboring down or upstream exons), Types (splicing type), and names of individuals.}
    \item{IR}{A data frame for the result of Intron retention, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), RetainEX (retained intron exon), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), and names of individuals.}
}
\author{
    Seonggyun Han, Younghee Lee
}
\seealso{
    \code{\link[IMAS]{SplicingReads}}
}
\examples{
    data(bamfilestest)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    \dontrun{
    ASdb <- Splicingfinder(transdb,Ncor=1)
    ASdb <- ExonsCluster(ASdb,transdb)
    ASdb <- RatioFromReads(ASdb,samplebamfiles,"paired",50,40,3,CalIndex="ES3")
    }
}