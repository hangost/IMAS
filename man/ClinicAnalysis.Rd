\name{ClinicAnalysis}
\alias{ClinicAnalysis}
\title{Analysis for differential clinical outcomes across PSI values}
\description{
    This function separate a set of samples into two groups (low and high PSI values) using K-means clustering and perform a statistical test to identify differential survival outcomes between the groups. Internally, this function calls the \code{\link{kmeans}} and \code{\link{survdiff}} functions in the \pkg{stats} and \pkg{survival} packages, respectively.
}
\usage{
    ClinicAnalysis(ASdb, ClinicalInfo = NULL, CalIndex = NULL, 
        display = FALSE, Ncor = 1, out.dir = NULL)
}
\arguments{
    \item{ASdb}{
        An ASdb object containing "SplicingModel" and "Ratio" slots from the \code{\link{Splicingfinder}} and \code{\link{RatioFromFPKM}} functions, respectively.
    }
    \item{ClinicalInfo}{
        A data frame consisting of a path of bam file and identifier of each sample. 
    }
    \item{CalIndex}{
        An index number in the ASdb object which will be tested in this function. 
    }
    \item{display}{
        The option returns the survival Kaplan-Meier plot. (TRUE = it will return the list object with a ggplot object and table showing the result of this function, FALSE = it will return p-value.)
    }
    \item{Ncor}{
        The number of cores for multi-threads function.
    }
    \item{out.dir}{
        An output directory.
    }
}

\value{
    ASdb with the slot (labeled by "Clinical") containing results from the \code{\link{ClinicAnalysis}} function.    The "Clinical" slot contains a list object and each element of the list object returns the results assigned to three elements, which is of each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). Three elements are as follows;
    \item{ES}{A data frame for the result of Exon skipping, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), 1stEX (alternatively spliced target exon), 2ndEX (second alternatively spliced target exon which is the other one of the mutually exclusive spliced exons), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), Pvalue (p-value of Kaplan-Meier test for differential survival outcomes between low and high PSI groups), and Fdr.p (FDR values).}
    \item{ASS}{A data frame for the result of Alternative splice sites, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), ShortEX (shorter spliced target exon), LongEX (longer spliced target exon), NeighborEX (neighboring down or upstream exons), Types (splicing type), Pvalue (p-value of Kaplan-Meier test for differential survival outcomes between low and high PSI groups), and Fdr.p (FDR values).}
    \item{IR}{A data frame for the result of Intron retention, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), RetainEX (retained intron exon), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), Pvalue (p-values of Kaplan- Meier test for differential survival outcomes between low and high PSI groups), and Fdr.p (FDR values).}
}

\author{
    Seonggyun Han, Younghee Lee
}
\seealso{
    \code{\link[stats]{kmeans}},
    \code{\link[survival]{survdiff}},
    \code{\link[survival]{survfit}}
}
\examples{
    data(bamfilestest)
    data(sampleclinical)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    #ASdb <- Splicingfinder(transdb,Ncor=1)
    #ASdb <- ExonsCluster(ASdb,transdb)
    #ASdb <- RatioFromReads(ASdb,samplebamfiles,"paired",50,40,3,CalIndex="ES3")
    #ASdb <- ClinicAnalysis(ASdb,Clinical.data,CalIndex="ES3",out.dir=NULL)
}