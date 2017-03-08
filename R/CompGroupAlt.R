CompGroupAlt <- function(ASdb=NULL,GroupSam=NULL,Ncor=1,CalIndex=NULL,out.dir=NULL){
    CalsigGroup <- function(ratio.mat=NULL){
        A.groups <- GroupSam$"GroupA"
        B.groups <- GroupSam$"GroupB"
        p.result <- foreach(each.mat.nums=1:nrow(ratio.mat),.packages=called.packages,.combine=rbind) %dopar% {
            each.mat <- rbind(ratio.mat[each.mat.nums,])
            A.exp.ratio <- try(rbind(each.mat[,is.element(colnames(each.mat),A.groups)]),silent=TRUE)
            B.exp.ratio <- try(rbind(each.mat[,is.element(colnames(each.mat),B.groups)]),silent=TRUE)
            A.exp.ratio <- A.exp.ratio[,A.exp.ratio != "NA" & A.exp.ratio != "NaN" & !is.na(A.exp.ratio)]
            B.exp.ratio <- B.exp.ratio[,B.exp.ratio != "NA" & B.exp.ratio != "NaN" & !is.na(B.exp.ratio)]
            not.errors <- length(grep("Error",A.exp.ratio)) == 0 & length(grep("Error",B.exp.ratio)) == 0 
            if (not.errors & length(A.exp.ratio) != 0 & length(B.exp.ratio) != 0){
                total.mat <- rbind(t(rbind(A.exp.ratio,"A")),t(rbind(B.exp.ratio,"B")))
                total.mat <- data.frame(round(as.double(total.mat[,1])*100),total.mat[,2])
                colnames(total.mat) <- c("y","x")
                glm.result <- lm(y~x,data=total.mat)
                s.glm <- summary(glm.result)
                #p.glm <- s.glm$coefficients[nrow(s.glm$coefficients),"Pr(>|z|)"]
                p.glm <- pf(s.glm$fstatistic[1],s.glm$fstatistic[2],s.glm$fstatistic[3],lower.tail=FALSE)
                p.glm
            }
            else as.matrix("NA")
        }
        ratio.mat <- cbind(ratio.mat,p.result)
        colnames(ratio.mat)[ncol(ratio.mat)] <- "Diff.P"
        return (ratio.mat)
    }
    each.mat.nums <- NULL
    called.packages <- c("lme4","GenomicRanges","GenomicFeatures")
    registerDoParallel(cores=Ncor)
    Total.ratio.mat <- ASdb@"Ratio"
    not.values <- length(Total.ratio.mat) != 0 & Total.ratio.mat != "NA" & !is.na(Total.ratio.mat)
    tested.types <- names(Total.ratio.mat)[not.values]
    totalTypes <- names(Total.ratio.mat)
    total.p.result <- lapply(totalTypes,function(each.type){
        if (is.element(each.type,tested.types)){
            each.ratio.mat <- Total.ratio.mat[[each.type]]
            if (length(CalIndex) != 0)    each.ratio.mat <- rbind(each.ratio.mat[is.element(each.ratio.mat[,"Index"],CalIndex),])
            if (length(each.ratio.mat) != 0){
                each.p.result <- CalsigGroup(each.ratio.mat)
                inter.cn <- c("Index","EnsID","Strand","Nchr","1stEX","2ndEX","DownEX","UpEX","Types",
                    "Diff.P","ShortEX","LongEX","NeighborEX","ShortNeighborEX","LongNeighborEX","RetainEX")
                each.p.result <- rbind(each.p.result[,is.element(colnames(each.p.result),inter.cn)])
                each.p.result <- cbind(each.p.result,Fdr.p=p.adjust(each.p.result[,"Diff.P"],"fdr"))
                each.p.result
            }
            else    as.matrix("NA")
        }
        else    as.matrix("NA")
    })
    names(total.p.result) <- totalTypes
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",GroupDiff=total.p.result,
        sQTLs=ASdb@"sQTLs",Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
    if (length(out.dir) != 0){
        system(paste("mkdir -p ",out.dir,"/AS_CompGroup",sep=""))
        write.table(total.p.result[["ES"]],paste(out.dir,"/AS_CompGroup/ES_CompGroup.txt",sep=""),sep='\t',quote=FALSE)
        write.table(total.p.result[["ASS"]],paste(out.dir,"/AS_CompGroup/ASS_CompGroup.txt",sep=""),sep='\t',quote=FALSE)
        write.table(total.p.result[["IR"]],paste(out.dir,"/AS_CompGroup/IR_CompGroup.txt",sep=""),sep='\t',quote=FALSE)
    }
    return (ASdb)
    
}

