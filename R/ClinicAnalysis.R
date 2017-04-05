ClinicAnalysis <- function(ASdb,ClinicalInfo=NULL,CalIndex=NULL,
    display=FALSE,Ncor=1,out.dir=NULL){
    kmsur <- function(ClinicalInfo,txs.ratio){
        pv <- NULL
        sta.na <- ClinicalInfo[,"status"] != "NA" & 
            !is.na(ClinicalInfo[,"status"])
        sur.na <- ClinicalInfo[,"sur_time"] != "NA" & 
            !is.na(ClinicalInfo[,"sur_time"])
        over.samples <- rownames(ClinicalInfo)[sta.na & sur.na]
        over.samples <- intersect(over.samples,names(txs.ratio))
        txs.ratio <- as.double(txs.ratio[over.samples])
        names(txs.ratio) <- over.samples
        kmeans.re <- kmeans(txs.ratio,2)
        exp.groups <- kmeans.re$"cluster"
        if (kmeans.re$"centers"[1] > kmeans.re$"centers"[2]){
            low.high <- cbind(c(1,2),c("High PSI","Low PSI"))
        }
        else    low.high <- cbind(c(1,2),c("Low PSI","High PSI"))
        ClinicalInfo <- ClinicalInfo[over.samples,]
        Total.info <- cbind(ClinicalInfo,exp.groups)
        colnames(Total.info) <- c("status","sur_time","groups")
        Total.info <- data.frame(Total.info)
        Surv <- Surv(c(Total.info[,"sur_time"]),c(Total.info[,"status"]))
        Conditions <- c(Total.info[,"groups"])
        numSize = length(levels(factor(unique(Total.info[,"groups"]))))
        if(numSize < 2)    return (NULL)
        rank.p <- survdiff(Surv ~ Conditions)
        pv <- pchisq(rank.p$chisq,1, lower.tail=FALSE)
        summary_coxph <- summary(coxph(Surv ~ Conditions))
        ci <- summary_coxph$conf.int
        if (display) {
            Conditions <- gsub(1,low.high[1,2],Conditions)
            Conditions <- gsub(2,low.high[2,2],Conditions)
            surfit = survfit(Surv ~ Conditions)
            t.Conditions <- table(Conditions)
            ncondi <- paste(paste("Condition",1:length(t.Conditions),
                " : ",t.Conditions,sep=""),collapse=", ")
            plot.re <- autoplot(surfit,surv.linetype = 'dashed',
                conf.int = FALSE,censor.shape = '*', censor.size = 5)
            stat = paste("n = ", length(Conditions),"( ",ncondi,
                " )", ", "," p = ",round(pv,3), sep="")
            vdim = dim(ci)
            x=min(Surv[,1])+0.4*(max(Surv[,1])-min(Surv[,1]))
            plot.re -> pv
        }
        return(pv)
    }
    kmEnv <- environment(kmsur)
    p.cal <- function(test.mat){
        ea.re <- test.mat
        total.p <- foreach(each.num=seq_len(nrow(ea.re)),
            .packages=called.packages,.combine=rbind) %dopar% {
            Pv <- NULL
            ea.ra <- rbind(ea.re[each.num,])
            p.r <- ea.ra[!is.element(colnames(ea.ra),rownames(ClinicalInfo))]
            ea.ra <- ea.ra[,!is.na(ea.ra) & ea.ra != "NA"]
            ov.sam <- intersect(rownames(ClinicalInfo),names(ea.ra))
            ea.ra <- ea.ra[ov.sam]
            sub.Cl <- rbind(ClinicalInfo[ov.sam,])
            if (length(ov.sam) > as.integer(t.sam/3)){
                Pv <- kmEnv$kmsur(sub.Cl,ea.ra)
                if (display)    Pv
                else    rbind(c(p.r,Pv))
            }
            else    Pv
        }
        if (!display & length(total.p)){
            cn.test <- !is.element(colnames(ea.re),rownames(ClinicalInfo))
            colnames(total.p) <- c(colnames(ea.re)[cn.test],"Pvalue")
        }
        return (total.p)
    }
    each.num <- NULL
    registerDoParallel(cores=Ncor) 
    Exon.ratio.mat <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
    names(Exon.ratio.mat) <- c("ES","ASS","IR")
    clinical.mat <- Exon.ratio.mat
    if (any(length(CalIndex))){
        ES.n <- grep("ES",CalIndex)
        ASS.n <- grep("ASS",CalIndex)
        IR.n <- grep("IR",CalIndex)
        if (ncol(ASdb@Ratio$ES) != 1 & any(length(ES.n))){
            ea.mat <- ASdb@Ratio$ES
            ea.re <- ea.mat[is.element(ea.mat[,"Index"],CalIndex),]
            Exon.ratio.mat$ES <- rbind(ea.re)
        }
        if (ncol(ASdb@Ratio$ASS) != 1 & any(length(ASS.n))){
            ea.mat <- ASdb@Ratio$ASS
            ea.re <- ea.mat[is.element(ea.mat[,"Index"],CalIndex),]
            Exon.ratio.mat$ASS <- rbind(ea.re)
        }
        if (ncol(ASdb@Ratio$IR) != 1 & any(length(IR.n))){
            ea.mat <- ASdb@Ratio$IR
            ea.re <- ea.mat[is.element(ea.mat[,"Index"],CalIndex),]
            Exon.ratio.mat$IR <- rbind(ea.re)
        }
    }
    else    Exon.ratio.mat <- ASdb@Ratio
    t.sam <- nrow(ClinicalInfo)
    called.packages <- c("lme4","GenomicRanges","GenomicFeatures")
    ES.re <- p.cal(Exon.ratio.mat$"ES")
    ASS.re <- p.cal(Exon.ratio.mat$"ASS")
    IR.re <- p.cal(Exon.ratio.mat$"IR")
    pre.result <- list(ES.re,ASS.re,IR.re)
    names(pre.result) <- c("ES","ASS","IR")
    total.types <- names(pre.result)
    total.types <- total.types[lengths(pre.result) > 1]
    if (display){
        if (total.types == "ES")    return (pre.result$"ES")
        if (total.types == "ASS")    return (pre.result$"ASS")
        if (total.types == "IR")    return (pre.result$"IR")
    }
    if (is.element("ES",total.types)){
        each.re <- unique(pre.result$"ES")
        rownames(each.re) <- 1:nrow(each.re)
        fp <- p.adjust(as.double(each.re[,"Pvalue"]),"fdr")
        each.re <- cbind(each.re,Fdr.p=fp)
        clinical.mat$"ES" <- each.re
    }
    if (is.element("ASS",total.types)){
        each.re <- unique(pre.result$"ASS")
        rownames(each.re) <- 1:nrow(each.re)
        fp <- p.adjust(as.double(each.re[,"Pvalue"]),"fdr")
        each.re <- cbind(each.re,Fdr.p=fp)
        clinical.mat$"ASS" <- each.re
    }
    if (is.element("IR",total.types)){
        each.re <- unique(pre.result$"IR")
        rownames(each.re) <- 1:nrow(each.re)
        fp <- p.adjust(as.double(each.re[,"Pvalue"]),"fdr")
        each.re <- cbind(each.re,Fdr.p=fp)
        clinical.mat$"IR" <- each.re
    }
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",
        GroupDiff=ASdb@"GroupDiff",sQTLs=ASdb@"sQTLs",
        Me.sQTLs=ASdb@"Me.sQTLs",Clinical=clinical.mat)
    if (length(out.dir)){
        p.out <- paste(out.dir,"/AS_Survival/",sep="")
        system(paste("mkdir -p ",p.out,sep=""))
        write.table(clinical.mat[["ES"]],
            paste(p.out,"ES_Survival.txt",sep=""),sep='\t',quote=FALSE)
        write.table(clinical.mat[["ASS"]],
            paste(p.out,"ASS_Survival.txt",sep=""),sep='\t',quote=FALSE)
        write.table(clinical.mat[["IR"]],
            paste(p.out,"IR_Survival.txt",sep=""),sep='\t',quote=FALSE)
    }
    return (ASdb)
}

