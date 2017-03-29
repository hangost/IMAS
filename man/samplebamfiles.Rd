\name{samplebamfiles}
\alias{samplebamfiles}
\docType{data}
\title{
    A data frame for example expression bam files.
    }
\description{
    A path and identifier of bam files for 50 samples. For each bam file, mapped reads were randomly generated that came from the genomic region of chr11: 100,933,178 - 100,996,889. With each simulated bam file of 50 samples, PSI level is calculated for the exon that is located in chr11: 100,962,491-100,962,607. The simulated PSI values are in the range of 0.6 to 1.0. The range of 0.9 to 1.0 of PSI values are assigned to PR-positive group and 0.5 to 0.6 to PR-negative group. The detailed overview of the data is described in the vignette.
    }
\usage{data(bamfilestest)}
\format{
    A data frame with paths and identifiers on the 50 samples
    }
\value{
    A data frame with paths and identifiers on the 50 samples
    }
\source{
    The data was provided from \pkg{IMAS}
    }
\examples{
    data(bamfilestest)
    }
\keyword{datasets}