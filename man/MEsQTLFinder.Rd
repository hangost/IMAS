\name{MEsQTLFinder}
\alias{MEsQTLFinder}
\title{
    Identify methylation loci that are significantly associated with alternatively spliced exons
}
\description{
    This function performs a regression test to identify significant association between methylation levels and PSI values using a linear regression model of \code{\link{lm}} function.
}
\usage{
    MEsQTLFinder(ASdb, Total.Medata = NULL, Total.Melocus = NULL, GroupSam = NULL, 
        Ncor = 1, CalIndex = NULL, out.dir = NULL)
}
\arguments{
    \item{ASdb}{
        An ASdb object including "SplicingModel" and "Ratio" slots from the \code{\link{Splicingfinder}} and \code{\link{RatioFromFPKM}} functions, respectively.
    }
    \item{Total.Medata}{
        A data frame consisting of methylation levels.
    }
    \item{Total.Melocus}{
        A data frame consisting of methylation locus.
    }
    \item{GroupSam}{
        A list object of a group of each sample.
    }
    \item{Ncor}{
        The number of cores for multi-threads.
    }
    \item{CalIndex}{
        An index number in the ASdb object which will be tested in this function. 
    }
    \item{out.dir}{
        An output directory.
    }
}

\value{
    ASdb with the slot (labeled by "Me.sQTLs") containing the results from the \code{\link{MEsQTLFinder}} function. The "Me.sQTLs" slot is consists of a list object and each element of the list object returns the results assigned to three elements, which is of each alternative splicing type (i.e. Exon skipping, Alternative splice site, Intron retention). Three elements are as follows;
    \item{ES}{A data frame for the result of Exon skipping, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), 1stEX (alternatively spliced target exon), 2ndEX (second alternatively spliced target exon which is the other one of the mutually exclusive spliced exons), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), pByMet (P-values of linear regression test for association between methylation levels and PSI values), fdrByMet (FDR values for the pByMet column), pByGroups (P-values of t-test for differential methylation levels between two groups, and fdrByGroups ( FDR values for the pByGroups column).}
    \item{ASS}{A data frame for the result of Alternative splice sites, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromo- some nam), ShortEX (shorter spliced target exon), LongEX (longer spliced target exon), NeighborEX (neighboring down or upstream exons), Types (splicing type), pByMet (P-values of linear regression test for association between methylation levels and PSI values), fdrByMet (FDR values for the pByMet column), pByGroups (P-values of t-test for differential methylation levels between groups, and fdrByGroups (adjust FDR values for the pByGroups column.}
    \item{IR}{A data frame for the result of Intron retention, consisting of the columns named as follows; Index (index number), EnsID (gene name), Nchr (chromosome name), RetainEX (retained intron exon), DownEX (downstream exon range), UpEX (upstream exon range), Types (splicing type), pByMet (P-values of linear regression test for association between methylation levels and PSI values), fdrByMet (adjust FDR values for the pByMet column), pByGroups (P-values of t-test for differential methylation levels between the groups, and fdrByGroups (adjust FDR values for the pByGroups column.}
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
    data(samplemethyl)
    data(samplemethyllocus)
    data(sampleGroups)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    sampleDB <- system.file("extdata", "sampleDB", package="IMAS")
    transdb <- loadDb(sampleDB)
    #ASdb <- Splicingfinder(transdb,Ncor=1)
    #ASdb <- ExonsCluster(ASdb,transdb)
    #ASdb <- RatioFromReads(ASdb,samplebamfiles,"paired",50,40,3,CalIndex="ES3")
    #ASdb <- MEsQTLFinder(ASdb,sampleMedata,sampleMelocus,CalIndex="ES3",GroupSam=GroupSam,out.dir=NULL)
}