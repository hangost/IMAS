\name{CompGroupAlt}
\alias{CompGroupAlt}
\title{
    Identify alternatively spliced exons with a differential PSIs between the groups
}
\description{
    This function performs a regression test to identify alternatively spliced exons that are differentially expressed between two groups. It will call the \code{\link{lm}} function to test a linear regression model.
}
\usage{
    CompGroupAlt(ASdb, GroupSam = NULL, Ncor = 1, CalIndex = NULL, out.dir = NULL)
}
\arguments{
    \item{ASdb}{
        An ASdb object containing "SplicingModel" and "Ratio" slots from the \code{\link{Splicingfinder}} and \code{\link{RatioFromFPKM}} functions, respectively.
    }
    \item{GroupSam}{
        A list object of a group of each sample.
    }
    \item{Ncor}{
        The number of cores for multi-threads function.
    }
    \item{CalIndex}{
        An index number in the ASdb object which will be tested in this function. 
    }
    \item{out.dir}{
        An output directory.
    }
}
\value{
    ASdb with the slot (labeled by "GroupDiff") containing results from the \code{\link{CompGroupAlt}} function. The "GroupDiff" slot consists of a list object and each element of the list object returns the results assigned to three elements, which is of each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). Three elements are as follows;
    \item{ES}{A data frame for the result of Exon skipping, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), 1stEX (alternatively spliced target exon), 2ndEX (second alternatively spliced target exon which is the other one of the mutually exclusive spliced exons), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), Diff.P (P-value of linear regression test for differential expression between groups), and Fdr.p (FDR values).}
    \item{ASS}{A data frame for the result of Alternative splice sites, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome nam), ShortEX (shorter spliced target exon), LongEX (longer spliced target exon), NeighborEX (neighboring down or upstream exons), Types (splicing type), Diff.P (P-value of linear regression test for differential expression between groups), and Fdr.p (FDR values).}
    \item{IR}{A data frame for the result of Intron retention, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), RetainEX (retained intron exon), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), Diff.P (P-value of linear regression test for differential expression between groups), and Fdr.p (FDR values).}
}
\references{
    Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{
    Seonggyun Han, Younghee Lee
}
\seealso{
    \code{\link[stats]{lm}}
}
\examples{
    data(bamfilestest)
    data(sampleGroups)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    \dontrun{
    ASdb <- Splicingfinder(transdb,Ncor=1)
    ASdb <- ExonsCluster(ASdb,transdb)
    ASdb <- RatioFromReads(ASdb,samplebamfiles,"paired",50,40,3,CalIndex="ES3")
    ASdb <- CompGroupAlt(ASdb,GroupSam,CalIndex="ES3")
    }
}