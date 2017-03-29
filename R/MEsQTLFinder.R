MEsQTLFinder <- function(ASdb=NULL,Total.Medata=NULL,Total.Melocus=NULL,
    GroupSam=NULL,Ncor=1,CalIndex=NULL,out.dir=NULL){
    CalSigMe <- function(ratio.mat,Me.mat,overlapMe,each.Melocus,chr){
        options(warn = -1)
        colnames(overlapMe) <- c("snp","locus")
        realnums <- ratio.mat != "NA" & ratio.mat != "NaN"
        realNA <- colnames(ratio.mat)[realnums]
        test.exp <- rbind(as.double(ratio.mat[,realNA]))
        colnames(test.exp) <- realNA
        Mern <- rownames(Me.mat)
        test.Me <- rbind(Me.mat[,realNA])
        rownames(test.Me) <- Mern
        each.row.ratio <- test.exp
        numsamp <- 10
        if (length(realNA) < numsamp)    return (NULL)
        pre.result.lm <- NULL
        stac.result <- lapply(seq_len(nrow(test.Me)),function(test.each.num){
            Meid <- rownames(test.Me)[test.each.num]
            te.Me <- rbind(Me.mat[test.each.num,])
            realNA <- colnames(te.Me)[te.Me != "NA" & !is.na(te.Me)]
            if (realNA > as.integer(ncol(test.Me)/2)){
                te.Me <- rbind(te.Me[,realNA])
                test.exp <- rbind(test.exp[,realNA])
                total.mat <- rbind(test.exp,te.Me)
                total.mat <- data.frame(t(total.mat))
                colnames(total.mat) <- c("exp","me")
                l.p <- "NaN"
                pByttest <- NULL
                if (any(seq_along(GroupSam))){
                    A.te <- is.element(colnames(te.Me),GroupSam$"GroupA")
                    B.te <- is.element(colnames(te.Me),GroupSam$"GroupB")
                    A.methyl <- as.double(te.Me[,A.te])
                    B.methyl <- as.double(te.Me[,B.te])
                    pByttest <- t.test(A.methyl,B.methyl)$"p.value"
                }
                if (TRUE){
                    auo <- summary(lm(formula = exp ~ me, data = total.mat))
                    au.f <- auo$fstatistic
                    te1 <- any(seq_along(au.f[1]))
                    te2 <- any(seq_along(au.f[2]))
                    te3 <- any(seq_along(au.f[3]))
                    if (te1 & te2 & te3){
                        l.p <- pf(au.f[1],au.f[2],au.f[3],lower.tail=FALSE)
                        if (l.p != "NaN" & l.p != "NA"){
                            pre.result.lm <- cbind(Meid,l.p,pByttest,"lm")
                        }
                    }
                }
                pre.result.lm
            }
        })
        stac.result <- do.call(rbind,stac.result)
        if (is.element("pByttest",colnames(stac.result))){
            colnames(stac.result) <- c("Meid","pByMet","pByGroups","met")
        }
        else    colnames(stac.result) <- c("Meid","pByMet","met")
        return    (stac.result)
    }
    TestMe <- function(Each.mat){
        if (ncol(Each.mat) == 1)    return (NULL)
        subn <- is.element(Total.Melocus[,"CHR"],unique(Each.mat[,"Nchr"]))
        sub.Melo <- rbind(Total.Melocus[subn,])
        inter.Me <- intersect(rownames(Total.Medata),sub.Melo[,"Methyl"])
        sub.Meda <- rbind(Total.Medata[inter.Me,])
        rownames(sub.Meda) <- inter.Me
        te1 <- any(seq_along(Each.mat))
        te2 <- any(seq_along(sub.Melo))
        te3 <- any(seq_along(sub.Meda))
        te <- te1 & te2 & te3
        int.me.lo <- as.integer(sub.Melo[,"locus"])
        Me.ran <- IRanges(start=int.me.lo,end=int.me.lo)
        Me.ran <- GRanges(seqnames=Rle(sub.Melo[,"CHR"]),
            ranges=Me.ran,metadata=sub.Melo[,"Methyl"])
        if (any(te)){
            over.sam <- intersect(colnames(Each.mat),colnames(sub.Meda))
            i=1
            pa.result <- foreach(i=seq_len(nrow(Each.mat)),
                .packages=called.packages,.combine=rbind) %dopar% {
                test.mat <- rbind(Each.mat[i,])
                test.exp <- rbind(test.mat[,over.sam])
                ex.re <- test.mat[,is.element(colnames(test.mat),inter.cns)]
                ex.re <- do.call(rbind,strsplit(ex.re,"-"))
                ex.re <- unlist(strsplit(ex.re,","))
                test.ex.re <- ex.re != "NA" & ex.re != "NaN" & !is.na(ex.re)
                ex.re <- as.integer(ex.re[test.ex.re])
                ex.re <- cbind(min(ex.re),max(ex.re))
                colnames(ex.re) <- c("start","end")
                each.ran <- IRanges(start=ex.re[,"start"],end=ex.re[,"end"])
                chr.rle <- Rle(test.mat[,"Nchr"])
                EX.ran <- GRanges(seqnames=chr.rle,ranges=each.ran)
                EX.ran <- list(EX.ran)
                names(EX.ran) <- "alterIntron"
                overMe <- findOversnp(EX.ran,Me.ran)
                if (any(seq_along(overMe))&any(which(test.exp!="NA"))){
                    int.Me <- intersect(rownames(sub.Meda),overMe[,"snp"])
                    if (any(seq_along(int.Me))){
                        te.Meda <- rbind(sub.Meda[int.Me,over.sam])
                        rownames(te.Meda) <- int.Me
                        on <- is.element(sub.Melo[,"Methyl"],overMe[,"snp"])
                        te.Melo <- rbind(sub.Melo[on,])
                        sig.re <- CalSigMe(test.exp,te.Meda,overMe,
                            te.Melo,test.mat[,"Nchr"])
                        if (any(seq_len(length(sig.re)))){
                            nsi <- nrow(sig.re)
                            o.cn.n <- is.element(inter.cn,colnames(test.mat))
                            o.in.cn <- inter.cn[o.cn.n]
                            pre.inf <- rep(rbind(test.mat[,o.in.cn]),nsi)
                            pre.inf <- matrix(pre.inf,nsi,byrow=TRUE)
                            colnames(pre.inf) <- o.in.cn
                            p.ma <- sig.re[,is.element(colnames(sig.re),
                                c("pByMet","pByGroups"))]
                            cbind(sig.re[,"Meid"],pre.inf,p.ma)
                        }
                    }
                    else {NULL}
                }
                else {NULL}
            }
            if (any(length(pa.result))){
                colnames(pa.result)[1] <- "MeID"
                }
        }
        return (pa.result)
    }
    fdr.cal <- function(each.result){
        if (!any(nrow(each.result)))    return (na.mat)
        rownames(each.result) <- 1:nrow(each.result)
        p.num <- which(colnames(each.result)=="pByMet")
        Gp.n <- which(colnames(each.result)=="pByGroups")
        fdr.met <- p.adjust(as.double(each.result[,p.num]),"fdr")
        p.num.ra <- 1:as.integer(p.num-1)
        fd.mat <- cbind(pByMet=each.result[,"pByMet"],
            fdrByMet=fdr.met)
        if (any(seq_along(Gp.n))){
            fdr.group <- p.adjust(as.double(each.result[,Gp.n]),"fdr")
            fdg.mat <- cbind(pByGroups=each.result[,"pByGroups"],
                fdrByGroups=fdr.group)
            fd.mat <- cbind(fd.mat,fdg.mat)
        }
        each.result <- cbind(each.result[,p.num.ra],fd.mat)
        return (each.result)
    }
    na.mat <- as.matrix("NA")
    registerDoParallel(cores=Ncor)
    called.packages <- c("lme4","GenomicRanges","GenomicFeatures")
    inter.cns <- c("DownEX","UpEX","ShortEX","LongEX","NeighborEX",
        "ShortNeighborEX","LongNeighborEX")
    inter.cn <- c("Index","EnsID","Strand","Nchr","1stEX","2ndEX",
        "DownEX","UpEX","Types","Diff.P","ShortEX","LongEX","NeighborEX",
        "ShortNeighborEX","LongNeighborEX","RetainEX")
    ra.mat <- ASdb@Ratio
    T.ra <- list(na.mat,na.mat,na.mat)
    names(T.ra) <- c("ES","ASS","IR")
    total.list <- T.ra
    if (any(seq_along(CalIndex))){
        ES.n <- grep("ES",CalIndex)
        ASS.n <- grep("ASS",CalIndex)
        IR.n <- grep("IR",CalIndex)
        if (any(seq_along(ES.n))){
            ea.mat <- ra.mat$ES
            T.ra$ES <- rbind(ea.mat[is.element(ea.mat[,"Index"],CalIndex),])
        }
        if (any(seq_along(ASS.n))){
            ea.mat <- ra.mat$ASS
            T.ra$ASS <- rbind(ea.mat[is.element(ea.mat[,"Index"],CalIndex),])
        }
        if (any(seq_along(IR.n))){
            ea.mat <- ra.mat$IR
            T.ra$IR <- rbind(ea.mat[is.element(ea.mat[,"Index"],CalIndex),])
        }
    }
    else    T.ra <- ra.mat
    Total.Medata <- as.matrix(Total.Medata)
    Total.Melocus <- gsub(" ","",as.matrix(Total.Melocus))
    total.result <- NULL
    ES.re <- TestMe(T.ra$ES)
    ASS.re <- TestMe(T.ra$ASS)
    IR.re <- TestMe(T.ra$IR)
    total.list$"ES" <- fdr.cal(unique(ES.re))
    total.list$"ASS" <- fdr.cal(unique(ASS.re))
    total.list$"IR" <- fdr.cal(unique(IR.re))
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",
        GroupDiff=ASdb@"GroupDiff",sQTLs=ASdb@"sQTLs",Me.sQTLs=total.list,
        Clinical=ASdb@"Clinical") 
    if (length(out.dir)){
        p.out <- paste(out.dir,"/AS_Me-sQTLs/",sep="")
        system(paste("mkdir -p ",p.out,sep=""))
        write.table(total.list[["ES"]],
            paste(p.out,"/ES_Me-sQTLs.txt",sep=""),sep='\t',quote=FALSE)
        write.table(total.list[["ASS"]],
            paste(p.out,"/ASS_Me-sQTLs.txt",sep=""),sep='\t',quote=FALSE)
        write.table(total.list[["IR"]],
            paste(p.out,"/IR_Me-sQTLs.txt",sep=""),sep='\t',quote=FALSE)
    }
    return (ASdb)
}

