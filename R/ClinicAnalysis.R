ClinicAnalysis <- function(ASdb,ClinicalInfo=NULL,CalIndex=NULL,display=FALSE,Ncor=1,out.dir=NULL){
    km_sur <- function(ClinicalInfo,txs.ratio){
        pv <- NULL
        sta.na <- ClinicalInfo[,"status"] != "NA" & !is.na(ClinicalInfo[,"status"])
        sur.na <- ClinicalInfo[,"sur_time"] != "NA" & !is.na(ClinicalInfo[,"sur_time"])
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
            Conditions <- gsub(2,low.high[2,2],gsub(1,low.high[1,2],Conditions))
            surfit = survfit(Surv ~ Conditions)
            t.Conditions <- table(Conditions)
            ncondi <- paste(paste("Condition",1:length(t.Conditions)," : ",t.Conditions,sep=""),collapse=", ")
            plot.re <- autoplot(surfit,surv.linetype = 'dashed', conf.int = FALSE,
                censor.shape = '*', censor.size = 5)
            stat = paste("n = ", length(Conditions),"( ",ncondi," )", ", "," p = ",round(pv,3), sep="")
            vdim = dim(ci)
            x=min(Surv[,1])+0.4*(max(Surv[,1])-min(Surv[,1]))
            plot.re -> pv
        }
        return(pv)
    }
    each.num <- NULL
    registerDoParallel(cores=Ncor) 
    Exon.ratio.mat <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
    names(Exon.ratio.mat) <- c("ES","ASS","IR")
    clinical.mat <- Exon.ratio.mat
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
    total.samples <- nrow(ClinicalInfo)
    total.types <- names(Exon.ratio.mat)
    total.types <- total.types[Exon.ratio.mat != "NA" & lengths(Exon.ratio.mat) != 0]
    called.packages <- c("lme4","GenomicRanges","GenomicFeatures")
    pre.result <- lapply(total.types,function(each.type){
        each.type.ratio.mat <- Exon.ratio.mat[[each.type]]
        total.p <- foreach(each.num=1:nrow(each.type.ratio.mat),.packages=called.packages,.combine=rbind) %dopar% {
            Pvalues <- NULL
            each.ratio.mat <- rbind(each.type.ratio.mat[each.num,])
            pre.ratio.mat <- each.ratio.mat[!is.element(colnames(each.ratio.mat),rownames(ClinicalInfo))]
            each.ratio.mat <- each.ratio.mat[,!is.na(each.ratio.mat) & each.ratio.mat != "NA"]
            inter.over.samples <- intersect(rownames(ClinicalInfo),names(each.ratio.mat))
            each.ratio.mat <- each.ratio.mat[inter.over.samples]
            sub.ClinicalInfo <- rbind(ClinicalInfo[inter.over.samples,])
            if (length(inter.over.samples) > as.integer(total.samples/3)){
                Pvalues <- km_sur(sub.ClinicalInfo,each.ratio.mat)
                if (display)    Pvalues
                else    rbind(c(pre.ratio.mat,Pvalues))
            }
            else    Pvalues
        }
        if (!display & length(total.p) != 0){
            cn.test <- !is.element(colnames(Exon.ratio.mat[[each.type]]),rownames(ClinicalInfo))
            colnames(total.p) <- c(colnames(Exon.ratio.mat[[each.type]])[cn.test],"Pvalue")
        }
        total.p
    })
    names(pre.result) <- total.types
    if (display)    return (pre.result[[total.types]])
    if (length(pre.result$"ES") != 0){
        each.result <- unique(pre.result$"ES")
        rownames(each.result) <- 1:nrow(each.result)
        each.result <- cbind(each.result,Fdr.p=p.adjust(as.double(each.result[,"Pvalue"]),"fdr"))
        clinical.mat$"ES" <- each.result
    }
    if (length(pre.result$"ASS") != 0){
        each.result <- unique(pre.result$"ASS")
        rownames(each.result) <- 1:nrow(each.result)
        each.result <- cbind(each.result,Fdr.p=p.adjust(as.double(each.result[,"Pvalue"]),"fdr"))
        clinical.mat$"ASS" <- each.result
    }
    if (length(pre.result$"IR") != 0){
        each.result <- unique(pre.result$"IR")
        rownames(each.result) <- 1:nrow(each.result)
        each.result <- cbind(each.result,Fdr.p=p.adjust(as.double(each.result[,"Pvalue"]),"fdr"))
        clinical.mat$"IR" <- each.result
    }
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",GroupDiff=ASdb@"GroupDiff",
        sQTLs=ASdb@"sQTLs",Me.sQTLs=ASdb@"Me.sQTLs",Clinical=clinical.mat)
    if (length(out.dir) != 0){
        system(paste("mkdir -p ",out.dir,"/AS_Survival",sep=""))
        write.table(clinical.mat[["ES"]],paste(out.dir,"/AS_Survival/ES_Survival.txt",sep=""),sep='\t',quote=FALSE)
        write.table(clinical.mat[["ASS"]],paste(out.dir,"/AS_Survival/ASS_Survival.txt",sep=""),sep='\t',quote=FALSE)
        write.table(clinical.mat[["IR"]],paste(out.dir,"/AS_Survival/IR_Survival.txt",sep=""),sep='\t',quote=FALSE)
    }
    return (ASdb)
}

