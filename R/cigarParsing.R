cigarParsing <- function(samfile.reads = NULL){
    total.cigar <- c("M","N","D","I","S","H","P","X","=")
    if (is.na(samfile.reads[,"cigar"]) | length(samfile.reads[,"cigar"]) ==0 ){
        return (NULL)
    }
    cigar.string <- as.matrix(samfile.reads[,"cigar"])
    cigar.sepa <- strsplit(as.matrix(cigar.string),"")[[1]]
    cigar.num <- which(is.element(cigar.sepa,total.cigar)=="TRUE")
    start.num <- 1
    final.parce.cigar <- NULL
    for (j in 1:length(cigar.num)){
        mat.count <- cigar.sepa[start.num:as.integer(cigar.num[j]-1)]
        final.parce.cigar <- rbind(final.parce.cigar,as.integer(paste(mat.count,collapse="")))
        start.num <- as.integer(cigar.num[j])+1
    }
    rownames(final.parce.cigar) <- cigar.sepa[cigar.num]
    colnames(final.parce.cigar) <- "matching"
    return (final.parce.cigar)
}