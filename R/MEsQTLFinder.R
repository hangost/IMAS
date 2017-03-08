MEsQTLFinder <- function(ASdb=NULL,Total.Medata=NULL,Total.Melocus=NULL,GroupSam=NULL,Ncor=1,CalIndex=NULL,out.dir=NULL){
    CalSigMe <- function(ratio.mat,Me.mat,overlapMe,each.Melocus,chr,each.gene){
        options(warn = -1)
        colnames(overlapMe) <- c("snp","locus")
        realNA <- colnames(ratio.mat)[ratio.mat != "NA" & ratio.mat != "NaN"]
        test.exp <- rbind(as.double(ratio.mat[,realNA]))
        colnames(test.exp) <- realNA
        Mern <- rownames(Me.mat)
        test.Me <- rbind(Me.mat[,realNA])
        rownames(test.Me) <- Mern
        each.row.ratio <- test.exp
        numsamp <- 10
        if (realNA > numsamp){
            stac.result <- lapply(1:nrow(test.Me),function(test.each.num){
                pre.result.lm <- NULL
                Meid <- rownames(test.Me)[test.each.num]
                test.each.Me <- rbind(Me.mat[test.each.num,])
                realNA <- colnames(test.each.Me)[test.each.Me != "NA" & !is.na(test.each.Me)]
                if (realNA > as.integer(ncol(test.Me)/2)){
                    test.each.Me <- rbind(test.each.Me[,realNA])
                    test.exp <- rbind(test.exp[,realNA])
                    total.mat <- rbind(test.exp,test.each.Me)
                    total.mat <- data.frame(t(total.mat))
                    colnames(total.mat) <- c("exp","me")
                    lm.auo.pvalue <- "NaN"
                    pByttest <- NULL
                    if (length(GroupSam) != 0){
                        A.methyl <- as.double(test.each.Me[,is.element(colnames(test.each.Me),GroupSam$"GroupA")])
                        B.methyl <- as.double(test.each.Me[,is.element(colnames(test.each.Me),GroupSam$"GroupB")])
                        pByttest <- t.test(A.methyl,B.methyl)$"p.value"
                    }
                    if (TRUE){
                        auo <- summary(lm(formula = exp ~ me, data = total.mat))
                        if (length(auo$fstatistic[1])>0 &length(auo$fstatistic[2])>0 & length(auo$fstatistic[3])>0){
                            lm.auo.pvalue <- pf(auo$fstatistic[1],auo$fstatistic[2],auo$fstatistic[3],lower.tail=FALSE)
                            if (lm.auo.pvalue != "NaN" & lm.auo.pvalue != "NA"){
                                pre.result.lm <- cbind(Meid,lm.auo.pvalue,pByttest,"lm")
                            }
                        }
                    }
                    pre.result.lm
                }
            })
        }
        stac.result <- do.call(rbind,stac.result)
        if (is.element("pByttest",colnames(stac.result))){
            colnames(stac.result) <- c("Meid","pByMet","pByGroups","met")
        }
        else    colnames(stac.result) <- c("Meid","pByMet","met")
        return    (stac.result)
    }
    Exon.ratio.mat <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
    names(Exon.ratio.mat) <- c("ES","ASS","IR")
    total.list <- Exon.ratio.mat
    if (ncol(ASdb@"Ratio"[["ES"]]) != 1){
        Exon.ratio.mat$"ES" <- ASdb@"Ratio"[["ES"]]
        if (length(CalIndex) != 0){
            each.result <- Exon.ratio.mat$"ES"[is.element(Exon.ratio.mat$"ES"[,"Index"],CalIndex),]
            Exon.ratio.mat$"ES" <- rbind(each.result)
        }
    }
    if (ncol(ASdb@"Ratio"[["ASS"]]) != 1){
        Exon.ratio.mat$"ASS" <- ASdb@"Ratio"[["ASS"]]
        if (length(CalIndex) != 0){
            each.result <- Exon.ratio.mat$"ASS"[is.element(Exon.ratio.mat$"ASS"[,"Index"],CalIndex),]
            Exon.ratio.mat$"ASS" <- rbind(each.result)
        }
    }
    if (ncol(ASdb@"Ratio"[["IR"]]) != 1){
        Exon.ratio.mat$"IR" <- ASdb@"Ratio"[["IR"]]
        if (length(CalIndex) != 0){
            each.result <- Exon.ratio.mat$"IR"[is.element(Exon.ratio.mat$"IR"[,"Index"],CalIndex),]
            Exon.ratio.mat$"IR" <- rbind(each.result)
        }
    }
    subtypes <- names(Exon.ratio.mat[Exon.ratio.mat != "NA"])
    registerDoParallel(cores=Ncor)
    Total.Medata <- as.matrix(Total.Medata)
    len.chr <- grep("chr",(Exon.ratio.mat[[subtypes[1]]][,"Nchr"]))
    if (length(len.chr) != 0){
        Total.Melocus[,"CHR"] <- paste("chr",gsub("chr","",Total.Melocus[,"CHR"]),sep="")
    }
    else if (length(len.chr) == 0){
        Total.Melocus[,"CHR"] <- gsub("chr","",Total.Melocus[,"CHR"])
    }
    Total.Melocus <- gsub(" ","",as.matrix(Total.Melocus))
    total.result <- NULL
    final.result <- lapply(subtypes,function(each.type){
        total.result <- NULL
        sub.exon.ratio.mat <- gsub(" ","",as.matrix(Exon.ratio.mat[[each.type]]))
        sub.Melocus <- rbind(Total.Melocus[is.element(Total.Melocus[,"CHR"],unique(sub.exon.ratio.mat[,"Nchr"])),])
        inter.Me <- intersect(rownames(Total.Medata),sub.Melocus[,"Methyl"])
        sub.Medata <- rbind(Total.Medata[inter.Me,])
        rownames(sub.Medata) <- inter.Me
        if (length(sub.exon.ratio.mat) != 0 & length(sub.Melocus) != 0 & length(sub.Medata) != 0){
            over.samples <- intersect(colnames(sub.exon.ratio.mat),colnames(sub.Medata))
            called.packages <- c("lme4","GenomicRanges","GenomicFeatures")
            i=1
            Total.chr <- unique(as.matrix(sub.Melocus[,"CHR"]))
            Total.chr <- Total.chr[order(as.integer(Total.chr))]
            for (j in 1:length(Total.chr)){
                print (paste("-------------------Processing : chr",Total.chr[j]," (",each.type,") -------------------",sep=""))
                ch.sub.exon.ratio <- rbind(sub.exon.ratio.mat[sub.exon.ratio.mat[,"Nchr"] == Total.chr[j],])
                ch.Me.locus <- rbind(sub.Melocus[sub.Melocus[,"CHR"] == Total.chr[j],])
                inter.Me <- intersect(rownames(sub.Medata),ch.Me.locus[,"Methyl"])
                ch.Me.data <- rbind(sub.Medata[inter.Me,])
                rownames(ch.Me.data) <- inter.Me
                pa.result <- foreach(i=1:nrow(ch.sub.exon.ratio),.packages=called.packages,.combine=rbind) %dopar% {
                    each.sub.exon.ratio <- rbind(ch.sub.exon.ratio[i,])
                    test.expdata <- rbind(each.sub.exon.ratio[,over.samples])
                    inter.cns <- c("DownEX","UpEX","ShortEX","LongEX","NeighborEX","ShortNeighborEX","LongNeighborEX")
                    ex.regions <- each.sub.exon.ratio[,is.element(colnames(each.sub.exon.ratio),inter.cns)]
                    ex.regions <- do.call(rbind,strsplit(ex.regions,"-"))
                    ex.regions <- unlist(strsplit(ex.regions,","))
                    ex.regions <- ex.regions[ex.regions != "NA" & ex.regions != "NaN" & !is.na(ex.regions)]
                    ex.regions <- cbind(min(as.integer(ex.regions)),max(as.integer(ex.regions)))
                    colnames(ex.regions) <- c("start","end")
                    each.ranges <- IRanges(start=as.integer(ex.regions[,"start"]),end=as.integer(ex.regions[,"end"]))
                    EX.region.range <- GRanges(seqnames=Rle(each.sub.exon.ratio[,"Nchr"]),ranges=each.ranges)
                    EX.region.range <- list(EX.region.range)
                    names(EX.region.range) <- "alterIntron"
                    each.ranges <- IRanges(start=as.integer(ch.Me.locus[,"locus"]),end=as.integer(ch.Me.locus[,"locus"]))
                    MeRagne <- GRanges(seqnames=Rle(ch.Me.locus[,"CHR"]),ranges=each.ranges,metadata=ch.Me.locus[,"Methyl"])
                    overlapMe <- findOversnp(EX.region.range,MeRagne)
                    if (length(overlapMe) != 0 & length(test.expdata[test.expdata!="NA"]) != 0){
                        inter.Me <- intersect(rownames(ch.Me.data),overlapMe[,"snp"])
                        if (length(inter.Me) != 0){
                            test.Medata <- rbind(ch.Me.data[inter.Me,over.samples])
                            rownames(test.Medata) <- inter.Me
                            test.Melocus <- rbind(ch.Me.locus[is.element(ch.Me.locus[,"Methyl"],overlapMe[,"snp"]),])
                            sig.result <- CalSigMe(test.expdata,test.Medata,overlapMe,test.Melocus,Total.chr[j],each.sub.exon.ratio[,"EnsID"])
                            inter.cn <- c("Index","EnsID","Strand","Nchr","1stEX","2ndEX","DownEX","UpEX","Types",
                                "Diff.P","ShortEX","LongEX","NeighborEX","ShortNeighborEX","LongNeighborEX","RetainEX")
                            inter.cn <- inter.cn[is.element(inter.cn,colnames(each.sub.exon.ratio))]
                            pre.inf <- rep(rbind(each.sub.exon.ratio[,inter.cn]),nrow(sig.result))
                            pre.inf <- matrix(pre.inf,nrow=nrow(sig.result),byrow=TRUE)
                            colnames(pre.inf) <- inter.cn
                            cbind(sig.result[,"Meid"],pre.inf,sig.result[,is.element(colnames(sig.result),c("pByMet","pByGroups","met"))])
                        }
                        else {NULL}
                    }
                    else {NULL}
                }
                total.result <- rbind(total.result,pa.result)
                if (length(total.result) != 0)    colnames(total.result)[1] <- "MeID"
            }
            total.result
        }
        else NULL
    })
    names(final.result) <- subtypes
    if (length(final.result$"ES") != 0){
        each.result <- unique(final.result$"ES")
        rownames(each.result) <- 1:nrow(each.result)
        p.num <- which(colnames(each.result)=="pByMet")
        fdr.met <- p.adjust(as.double(each.result[,"pByMet"]),"fdr")
        fdr.group <-    p.adjust(as.double(each.result[,"pByGroups"]),"fdr")
        each.result <- cbind(each.result[,1:as.integer(p.num-1)],pByMet=each.result[,"pByMet"],
            fdrByMet=fdr.met,pByGroups=each.result[,"pByGroups"],fdrByGroups=fdr.group)
        total.list$"ES" <- each.result
    }
    if (length(final.result$"ASS") != 0){
        each.result <- unique(final.result$"ASS")
        rownames(each.result) <- 1:nrow(each.result)
        p.num <- which(colnames(each.result)=="pByMet")
        fdr.met <- p.adjust(as.double(each.result[,"pByMet"]),"fdr")
        fdr.group <-    p.adjust(as.double(each.result[,"pByGroups"]),"fdr")
        each.result <- cbind(each.result[,1:as.integer(p.num-1)],pByMet=each.result[,"pByMet"],
            fdrByMet=fdr.met,pByGroups=each.result[,"pByGroups"],fdrByGroups=fdr.group)
        total.list$"ASS" <- each.result
    }
    if (length(final.result$"IR") != 0){
        each.result <- unique(final.result$"IR")
        rownames(each.result) <- 1:nrow(each.result)
        p.num <- which(colnames(each.result)=="pByMet")
        fdr.met <- p.adjust(as.double(each.result[,"pByMet"]),"fdr")
        fdr.group <-    p.adjust(as.double(each.result[,"pByGroups"]),"fdr")
        each.result <- cbind(each.result[,1:as.integer(p.num-1)],pByMet=each.result[,"pByMet"],
            fdrByMet=fdr.met,pByGroups=each.result[,"pByGroups"],fdrByGroups=fdr.group)
        total.list$"IR" <- each.result
    }
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",GroupDiff=ASdb@"GroupDiff",
        sQTLs=ASdb@"sQTLs",Me.sQTLs=total.list,Clinical=ASdb@"Clinical") 
    if (length(out.dir) != 0){
        system(paste("mkdir -p ",out.dir,"/AS_Me-sQTLs",sep=""))
        write.table(final.result[["ES"]],paste(out.dir,"/AS_Me-sQTLs/ES_Me-sQTLs.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.result[["ASS"]],paste(out.dir,"/AS_Me-sQTLs/ASS_Me-sQTLs.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.result[["IR"]],paste(out.dir,"/AS_Me-sQTLs/IR_Me-sQTLs.txt",sep=""),sep='\t',quote=FALSE)
    }
    return (ASdb)
}

