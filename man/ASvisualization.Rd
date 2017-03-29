\name{ASvisualization}
\alias{ASvisualization}
\title{
    Visualize the results of the ASdb object.
}
\description{
    This function makes a pdf file consisting of plots for results in the ASdb object.
}
\usage{
    ASvisualization(ASdb,CalIndex=NULL,txTable=NULL,exon.range=NULL,snpdata=NULL,
    snplocus=NULL,methyldata=NULL,methyllocus=NULL,GroupSam=NULL,
    ClinicalInfo=NULL,out.dir=NULL)
}
\arguments{
    \item{ASdb}{
        A ASdb object.
    }
    \item{CalIndex}{
        An index number in the ASdb object which will be tested in this function.
    }
    \item{txTable}{
        A data frame of transcripts including transcript IDs, Ensembl gene names, Ensembl transcript names, transcript start sites, and transcript end sites.
    }
    \item{exon.range}{
        A list of GRanges objects including total exon ranges in each transcript resulted from the \code{\link{exonsBy}} function in \pkg{GenomicFeatures}.
    }
    \item{snpdata}{
        A data frame of genotype data.
    }
    \item{snplocus}{
        A data frame consisting of locus information of SNP markers in the snpdata.
    }
    \item{methyldata}{
        A data frame consisting of methylation levels.
    }
    \item{methyllocus}{
        A data frame consisting of methylation locus.
    }
    \item{GroupSam}{
        A list object of a group of each sample.
    }
    \item{ClinicalInfo}{
        A data frame consisting of a path of bam file and identifier of each sample. 
    }
    \item{out.dir}{
        An output directory
    }
}

\value{
    This function makes pdf for plots.
}

\author{
    Seonggyun Han, Younghee Lee
}

\examples{
    data(sampleGroups)
    data(samplemethyl)
    data(samplemethyllocus)
    data(samplesnp)
    data(samplesnplocus)
    data(sampleclinical)
    data(bamfilestest)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    ASdb <- Splicingfinder(transdb,Ncor=1)
    ASdb <- ExonsCluster(ASdb,transdb)
    ASdb <- RatioFromReads(ASdb,samplebamfiles,"paired",50,40,3,CalIndex="ES3")
    ASdb <- sQTLsFinder(ASdb,samplesnp,samplesnplocus,method="lm")
    ASdb <- CompGroupAlt(ASdb,GroupSam,CalIndex="ES3")
    ASdb <- MEsQTLFinder(ASdb,sampleMedata,sampleMelocus,CalIndex="ES3",GroupSam=GroupSam,out.dir=NULL)
    Sdb <- ClinicAnalysis(ASdb,Clinical.data,CalIndex="ES3",out.dir=NULL)
    exon.range <- exonsBy(transdb,by="tx")
    sel.cn <- c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND")
    txTable <- select(transdb, keys=names(exon.range),columns=sel.cn,keytype="TXID")
    ASvisualization(ASdb,CalIndex="ES3",txTable,exon.range,samplesnp,samplesnplocus,
        sampleMedata,sampleMelocus,GroupSam,Clinical.data,out.dir="./")
}