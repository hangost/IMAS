CompGroupAlt <- function(ASdb=NULL,GroupSam=NULL,Ncor=1,
    CalIndex=NULL,out.dir=NULL){
    CalsigGroup <- function(ratio.mat=NULL,A.groups,B.groups){
        p.glm <- NULL
        p.result <- foreach(each.nums=seq_len(nrow(ratio.mat)),
            .packages=called.packages,.combine=rbind) %dopar% {
            each.mat <- rbind(ratio.mat[each.nums,])
            A.nums <- is.element(colnames(each.mat),A.groups)
            B.nums <- is.element(colnames(each.mat),B.groups)
            A.ra <- try(rbind(each.mat[,A.nums]),silent=TRUE)
            B.ra <- try(rbind(each.mat[,B.nums]),silent=TRUE)
            A.ra <- A.ra[,A.ra != "NA" & A.ra != "NaN" & !is.na(A.ra)]
            B.ra <- B.ra[,B.ra != "NA" & B.ra != "NaN" & !is.na(B.ra)]
            A.er <- !any(seq_along(grep("Err",A.ra))) 
            B.er <- !any(seq_along(grep("Err",B.ra)))
            A.len <- length(A.ra) > 4
            B.len <- length(B.ra) > 4
            if (A.er & B.er & A.len & B.len){
                t.A.ra <- t(rbind(A.ra,"A"))
                t.B.ra <- t(rbind(B.ra,"B"))
                total.mat <- rbind(t.A.ra,t.B.ra)
                ra.100 <- round(as.double(total.mat[,1])*100)
                total.mat <- data.frame(ra.100,total.mat[,2])
                colnames(total.mat) <- c("y","x")
                lm.result <- lm(y~x,data=total.mat)
                s.lm <- summary(lm.result)
                p.lm <- pf(s.lm$fstatistic[1],s.lm$fstatistic[2],
                    s.lm$fstatistic[3],lower.tail=FALSE)
                p.lm
            }
            else as.matrix("NA")
        }
        ratio.mat <- cbind(ratio.mat,p.result)
        colnames(ratio.mat)[ncol(ratio.mat)] <- "Diff.P"
        return (ratio.mat)
    }
    sigEnv <- environment(CalsigGroup)
    A.groups <- GroupSam$"GroupA"
    B.groups <- GroupSam$"GroupB"
    each.nums <- NULL
    each.p.result <- NULL
    called.packages <- c("lme4","GenomicRanges","GenomicFeatures")
    registerDoParallel(cores=Ncor)
    T.ra <- ASdb@"Ratio"
    not.values <- any(length(T.ra)) & T.ra != "NA" & !is.na(T.ra)
    tested.types <- names(T.ra)[not.values]
    totalTypes <- names(T.ra)
    total.p <- lapply(totalTypes,function(each.type){
        if (is.element(each.type,tested.types)){
            each.ratio.mat <- T.ra[[each.type]]
            if (any(length(CalIndex))){
                cal.num <- is.element(each.ratio.mat[,"Index"],CalIndex)
                each.ratio.mat <- rbind(each.ratio.mat[cal.num,])
            }
            if (length(each.ratio.mat)){
                each.p.result <- sigEnv$CalsigGroup(each.ratio.mat,
                    A.groups,B.groups)
                real.num <- each.p.result[,"Diff.P"] != "NA"
                each.p.result <- rbind(each.p.result[real.num,])
                inter.cn <- c("Index","EnsID","Strand","Nchr",
                    "1stEX","2ndEX","DownEX","UpEX","Types",
                    "Diff.P","ShortEX","LongEX","NeighborEX",
                    "ShortNeighborEX","LongNeighborEX","RetainEX")
                ov.nm <- is.element(colnames(each.p.result),inter.cn)
                each.p.result <- rbind(each.p.result[,ov.nm])
                fdr.p <- p.adjust(each.p.result[,"Diff.P"],"fdr")
                each.p.result <- cbind(each.p.result,Fdr.p=fdr.p)
                if (!any(seq_along(each.p.result))){
                    as.matrix("NA")
                    }
                else    each.p.result
            }
            else    as.matrix("NA")
        }
        else    as.matrix("NA")
    })
    names(total.p) <- totalTypes
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",
        GroupDiff=total.p,sQTLs=ASdb@"sQTLs",Me.sQTLs=ASdb@"Me.sQTLs",
        Clinical=ASdb@"Clinical")
    if (length(out.dir)){
        p.out <- paste(out.dir,"/AS_CompGroup/",sep="")
        system(paste("mkdir -p ",p.out,sep=""))
        write.table(total.p$"ES",
            paste(p.out,"/ES_CompGroup.txt",sep=""),sep='\t',quote=FALSE)
        write.table(total.p$"ASS",
            paste(p.out,"/ASS_CompGroup.txt",sep=""),sep='\t',quote=FALSE)
        write.table(total.p$"IR",
            paste(p.out,"/IR_CompGroup.txt",sep=""),sep='\t',quote=FALSE)
    }
    return (ASdb)
    
}

