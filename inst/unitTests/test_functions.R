test_All_functions <- function(){
    #load data sets
    data(bamfilestest)
    ext.dir <- system.file("extdata", package="IMAS")
    samplebamfiles[,"path"] <- paste(ext.dir,"/samplebam/",samplebamfiles[,"path"],".bam",sep="")
    data(sampleGroups)
    data(samplesnp)
    data(samplesnplocus)
    data(samplemethyl)
    data(samplemethyllocus)
    data(sampleclinical)
    data(bamfilestest)
    sampledb <- list.files(ext.dir,pattern="DB",full.names=TRUE)
    transdb <- loadDb(sampledb)

    # sampleGroups is a list, and each of samplesnp, samplesnplocus, sampleMedata, sampleMelocus, Clinical.data, and samplebamfiles is a matrix, and transdb is a s4 class.
    checkTrue(is.list(GroupSam))
    checkTrue(is.matrix(samplesnp))
    checkTrue(is.matrix(samplesnplocus))
    checkTrue(is.matrix(sampleMedata))
    checkTrue(is.matrix(sampleMelocus))
    checkTrue(is.matrix(Clinical.data))
    checkTrue(is.matrix(samplebamfiles))
    checkTrue(identical(typeof(transdb),"S4"))

    # The example GroupSam is consisting of two elements
    checkTrue(identical(names(GroupSam),c("GroupA","GroupB")))

    # The example samplesnp is consisting of  genotypes of SNP1, 2, 3, 4, and 5 on the 50 samples
    checkTrue(identical(rownames(samplesnp),c("SNP1","SNP2","SNP3","SNP4","SNP5")))
    checkEquals(ncol(samplesnp),50)

    # The example samplesnplocus is consisting of  genomic coordinates of SNPs matched with those of samplesnp object
    checkTrue(identical(colnames(samplesnplocus),c("SNP","CHR","locus")))
    checkTrue(identical(samplesnplocus[,"SNP"],rownames(samplesnp)))

    # The example sampleMedata is consisting of methylation status of the five methylation loci on the 50 samples
    checkEquals(nrow(sampleMedata),5)
    checkEquals(ncol(sampleMedata),50)

    # The example sampleMedata is consisting of loci information of five methylations matched with those of sampleMedata object
    checkTrue(identical(colnames(sampleMelocus),c("Methyl","CHR","locus")))
    checkTrue(identical(sampleMelocus[,"Methyl"],rownames(sampleMedata)))

    # The example Clinical.data is consisting of information of status and survivel time on the 50 samples matched with 
    checkTrue(identical(colnames(Clinical.data),c("status","sur_time")))
    checkEquals(nrow(Clinical.data),50)

    # The example samplebamfiles is consisting of information of paths and identifier for each of 50 sample
    checkTrue(identical(colnames(samplebamfiles),c("path","names")))
    checkEquals(nrow(samplebamfiles),50)

    ## Tabulate and cartegorize alternative splicing exons into the four alternative splicing patterns (As an example, search for sQTLs in the ENSG00000082175 gene)
    ASdb <- Splicingfinder(GTFdb=transdb,calGene="ENSG00000082175",Ncor=1)
    ASdb <- ExonsCluster(ASdb,transdb)
    ## the result is saved in "SplicingModel" added slot in the ASdb object. ASdb is a s4 class. The data in the slot is a list object consisting of three elements ("ES","ASS", and "IR")
    checkTrue(identical(typeof(ASdb),"S4"))
    checkTrue(is.list(slot(ASdb,"SplicingModel")))
    checkTrue(identical(names(slot(ASdb,"SplicingModel")),c("ES","ASS","IR")))

    ## Estimate expression ratio (PSI) in the alternative splicing exons included in "SplicingModel" slot in the ASdb object. 
    ASdb <- RatioFromReads(ASdb=ASdb,samplebamfiles,readsInfo="paired",readLen=50,inserSize=40,minr=3,CalIndex="ES3",Ncor=1)
    ## The result is saved in "Ratio" added slot in the ASdb object. The data in the slot is a list object consisting of three elements ("ES","ASS", and "IR")
    checkTrue(is.list(slot(ASdb,"Ratio")))
    checkTrue(identical(names(slot(ASdb,"Ratio")),c("ES","ASS","IR")))

    ## Calculates significant differences in PSIs between groups using linear regression with PSI values included in "Ratio" slot in the ASdb object. 
    ASdb <- CompGroupAlt(ASdb=ASdb,GroupSam=GroupSam,Ncor=1,CalIndex="ES3")
    ## The result is saved in "GroupDiff" added slot in the ASdb object. The data in the slot is a list object consisting of three elements ("ES","ASS", and "IR")
    checkTrue(is.list(slot(ASdb,"GroupDiff")))
    checkTrue(identical(names(slot(ASdb,"GroupDiff")),c("ES","ASS","IR")))

    ## Identify significant SNPs that are associated with AS ratio (PSI)
    ASdb <- sQTLsFinder(ASdb=ASdb,Total.snpdata=samplesnp,Total.snplocus=samplesnplocus,GroupSam=GroupSam,Ncor=1,CalIndex="ES3",method="lm")
    ## The result is saved in "sQTLs" added slot in the ASdb object. The data in the slot is a list object consisting of three elements ("ES","ASS", and "IR")
    checkTrue(is.list(slot(ASdb,"sQTLs")))
    checkTrue(identical(names(slot(ASdb,"sQTLs")),c("ES","ASS","IR")))

    ## Identifies methylation locus that may affect AS events
    ASdb <- MEsQTLFinder(ASdb=ASdb,Total.Medata=sampleMedata,Total.Melocus=sampleMelocus,GroupSam=GroupSam,Ncor=1,CalIndex="ES3")
    ## The result is saved in "Me.sQTLs" added slot in the ASdb object. The data in the slot is a list object consisting of three elements ("ES","ASS", and "IR")
    checkTrue(is.list(slot(ASdb,"Me.sQTLs")))
    checkTrue(identical(names(slot(ASdb,"Me.sQTLs")),c("ES","ASS","IR")))

    ## Identifies which exons ,that are differentially expressed in groups, are associated with clinical outcomes
    ASdb <- ClinicAnalysis(ASdb=ASdb,ClinicalInfo=Clinical.data,Ncor=1,CalIndex="ES3")
    ## The result is saved in "Clinical" added slot in the ASdb object. The data in the slot is a list object consisting of three elements ("ES","ASS", and "IR")
    checkTrue(is.list(slot(ASdb,"Clinical")))
    checkTrue(identical(names(slot(ASdb,"Clinical")),c("ES","ASS","IR")))
    
    ## All functions of this package provide to use multi-core
    ASdb_1 <- sQTLsFinder(ASdb=ASdb,Total.snpdata=samplesnp,Total.snplocus=samplesnplocus,GroupSam=GroupSam,Ncor=1,CalIndex="ES3",method="lm")
    ASdb_4 <- sQTLsFinder(ASdb=ASdb,Total.snpdata=samplesnp,Total.snplocus=samplesnplocus,GroupSam=GroupSam,Ncor=4,CalIndex="ES3",method="lm")
    identical(slot(ASdb_1,"sQTLs"),slot(ASdb_4,"sQTLs"))
}








