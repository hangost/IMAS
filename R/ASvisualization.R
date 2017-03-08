ASvisualization <- function(ASdb,CalIndex=NULL,txTable=NULL,exon.range=NULL,snpdata=NULL,snplocus=NULL,methyldata=NULL,methyllocus=NULL,GroupSam=NULL,ClinicalInfo=NULL,out.dir=NULL){
    BoxforGroup <- function(test.Group.mat,test.ratio){
        GroupA.exp <- test.ratio[,is.element(colnames(test.ratio),GroupSam$"GroupA")]
        GroupB.exp <- test.ratio[,is.element(colnames(test.ratio),GroupSam$"GroupB")]
        A.B.exp <- c(as.double(GroupA.exp),as.double(GroupB.exp))
        A.B.groups <- c(rep("GroupA",length(GroupA.exp)),
            rep("GroupB",length(GroupB.exp)))
        total.exp <- data.frame(A.B.exp,A.B.groups)
        p.nums <- is.element(colnames(test.Group.mat),c("Diff.P","Fdr.p"))
        p.val.mat <- test.Group.mat[,p.nums]
        test.Group.mat[,p.nums] <- p.round(p.val.mat)
        text.box <- rbind(test.Group.mat[,c("Index","EnsID","Diff.P","Fdr.p")])
        colnames(total.exp) <- c("Ratio","Groups")
        gplot.result <- ggplot(data=total.exp,aes(x=Groups,y=Ratio,fill=Groups))+
            geom_boxplot(width = 0.3)+geom_jitter(width = 0.03)+ylim(c(0,1))
        gplot.result <- list(plot=gplot.result,text=text.box)
        return (gplot.result)
    }
    BoxforsQTLs <- function(snp.result=NULL,test.snp.mat=NULL){
        snpplots.re <- lapply(rownames(test.snp.mat),function(snpid){
            each.result <- rbind(snp.result[snp.result[,"SNP"] == snpid])
            if (length(each.result) != 0){
                colnames(each.result) <- colnames(snp.result)
                each.snp.mat <- rbind(test.snp.mat[snpid,])
                rownames(each.snp.mat) <- snpid
                box.list <- sQTLsFinder(ASdb,each.snp.mat,snplocus,method="boxplot",GroupSam=GroupSam,CalIndex=CalIndex)
                box.exp <- box.list$"exp"
                box.OR <- box.list$"Ratios"
                tb.cn <- c("Index","SNP","EnsID","pByGeno","FdrByGeno","pByGroups","FdrByGroups","OR","lowCI","highCI")
                text.box <- rbind(each.result[,is.element(colnames(each.result),tb.cn)])
                p.val.nums <- is.element(colnames(text.box),c("pByGeno","FdrByGeno","pByGroups","FdrByGroups"))
                p.val.mat <- text.box[,p.val.nums]
                p.val.mat <- p.round(p.val.mat)
                text.box[,p.val.nums] <- p.val.mat
                if (length(GroupSam) == 0){
                    TotalRatio <- lapply(names(box.exp),function(each.geno){
                        each.genoBox <- cbind(exp=box.exp[[each.geno]],each.geno)
                        each.genoBox
                    })
                    TotalRatio <- do.call(rbind,TotalRatio)
                    TotalRatio <- data.frame(as.double(TotalRatio[,"exp"]),TotalRatio[,"each.geno"])
                    colnames(TotalRatio) <- c("Ratio","Genotype")
                    gplot.result <- ggplot(data=TotalRatio,aes(x=Genotype,y=Ratio,fill=Genotype))+
                        geom_boxplot(width=0.3)+geom_jitter(width = 0.03)
                    gplot.result <- list(plot=gplot.result,text=text.box)
                    gplot.result
                }
                else if (length(GroupSam) != 0){
                    TotalRatio <- lapply(names(box.exp),function(each.geno){
                        each.genoBox <- cbind(exp=box.exp[[each.geno]],each.geno,group="Group A")
                        each.genoBox[is.element(rownames(each.genoBox),GroupSam$"GroupB"),"group"] <- "Group B"
                        each.genoBox
                    })
                    TotalRatio <- do.call(rbind,TotalRatio)
                    TotalRatio <- data.frame(Ratio=as.double(TotalRatio[,"exp"]),Genotype=TotalRatio[,"each.geno"],Groups=TotalRatio[,"group"])
                    TotalRatio.A <- TotalRatio[TotalRatio[,"Groups"] == "Group A",]
                    TotalRatio.B <- TotalRatio[TotalRatio[,"Groups"] == "Group B",]
                    A.ggplot <- ggplot(data=TotalRatio.A,aes(x=Genotype,y=Ratio,fill=Genotype))+
                        geom_boxplot(width=0.5)+geom_jitter(width = 0.03) + ggtitle("Group A")+
                        theme(plot.title = element_text(vjust = 0)) + ylim(c(0,1))
                    B.ggplot <- ggplot(data=TotalRatio.B,aes(x=Genotype,y=Ratio,fill=Genotype))+
                        geom_boxplot(width=0.5)+geom_jitter(width = 0.03) + ggtitle("Group B")+
                        theme(plot.title = element_text(vjust = 0)) + ylim(c(0,1))
                    gplot.result <- arrangeGrob(A.ggplot,B.ggplot,ncol=2)
                    text.box[,c("OR","lowCI","highCI")] <- round(as.double(text.box[,c("OR","lowCI","highCI")]),3)
                    gplot.result <- list(plot=gplot.result,text=text.box)
                    gplot.result
                }
            }
            else    NULL
        })
        names(snpplots.re) <- rownames(test.snp.mat)
        return (snpplots.re)
    }
    
    
    BoxforMe <- function(test.Me.mat,test.ratio){
        meplots.re <- lapply(test.Me.mat[,"MeID"],function(each.id){
            over.samples <- intersect(colnames(test.ratio),colnames(methyldata))
            each.me.mat <- rbind(test.Me.mat[test.Me.mat[,"MeID"] == each.id,])
            each.me.ratio <- rbind(methyldata[rownames(methyldata) == each.id,over.samples])
            if (length(each.me.ratio) != 0){
                me.exp.ratio <- rbind(test.ratio[,over.samples],each.me.ratio[,over.samples])
                tb.cn <- c("MeID","Index","EnsID","pByMet","fdrByMet","pByGroups","fdrByGroups")
                text.box <- each.me.mat[,is.element(colnames(each.me.mat),tb.cn)]
                text.box[grep("pBy|fdrBy",names(text.box))] <- p.round(text.box[grep("pBy|fdrBy",names(text.box))])
                rownames(text.box) <- NULL
                if (length(GroupSam) == 0){
                    me.exp.ratio <- cbind(t(me.exp.ratio))
                    me.exp.ratio <- data.frame(as.double(me.exp.ratio[,1]),as.double(me.exp.ratio[,2]))
                    colnames(me.exp.ratio) <- c("Ratio","MethylationExp")
                    gplot.result <- ggplot(data=me.exp.ratio,aes(x=MethylationExp,y=Ratio))+
                        geom_point()+ylim(c(0,1))
                    gplot.result <- list(plot=gplot.result,text=text.box)
                }
                else {
                    me.exp.ratio <- cbind(t(me.exp.ratio),group="group")
                    me.exp.ratio[is.element(rownames(me.exp.ratio),GroupSam$"GroupA"),"group"] <- "GroupA"
                    me.exp.ratio[is.element(rownames(me.exp.ratio),GroupSam$"GroupB"),"group"] <- "GroupB"
                    me.exp.ratio <- data.frame(as.double(me.exp.ratio[,1]),as.double(me.exp.ratio[,2]),me.exp.ratio[,"group"])
                    colnames(me.exp.ratio) <- c("Ratio","MethylationExp","Groups")
                    lm.objec1 <- coef(lm(Ratio ~ MethylationExp, data=me.exp.ratio[me.exp.ratio[,"Groups"] == "GroupA",]))
                    lm.objec2 <- coef(lm(Ratio ~ MethylationExp, data=me.exp.ratio[me.exp.ratio[,"Groups"] == "GroupB",]))
                    lm.total.obj <- coef(lm(Ratio ~ MethylationExp, data=me.exp.ratio))
                    gplot.result.exp.me <- ggplot(data=me.exp.ratio,aes(x=MethylationExp,y=Ratio,color=Groups))+
                        geom_point()+ylim(c(0,1))+
                        geom_abline(intercept = as.double(lm.objec1["(Intercept)"]), 
                            slope = as.double(lm.objec1["MethylationExp"]),color="#F8766D")+
                        geom_abline(intercept = as.double(lm.objec2["(Intercept)"]), 
                            slope = as.double(lm.objec2["MethylationExp"]),color="#00BFC4")+
                        geom_abline(intercept = as.double(lm.total.obj["(Intercept)"]), 
                            slope = as.double(lm.total.obj["MethylationExp"]),color="grey")
                    gplot.result.me.group <- ggplot(data=me.exp.ratio,aes(x=Groups,y=MethylationExp,fill=Groups))+
                        geom_boxplot(width = 0.3)+geom_jitter(width = 0.03)
                    gplot.result <- arrangeGrob(gplot.result.exp.me,gplot.result.me.group,ncol=2)
                    gplot.result <- list(plot=gplot.result,text=rbind(text.box))
                    gplot.result
                }
                gplot.result
            }
            else NULL
        })
        names(meplots.re) <- test.Me.mat[,"MeID"]
        return (meplots.re)
    }
    
    BoxforClinical <- function(each.index){
        text.box <- ClinicalAnal[,is.element(colnames(ClinicalAnal),c("Index","EnsID","Pvalue","Fdr.p"))]
        text.box[c("Pvalue","Fdr.p")] <- p.round(text.box[c("Pvalue","Fdr.p")])
        clinic.plots <- ClinicAnalysis(ASdb,ClinicalInfo,display=TRUE,CalIndex=each.index,out.dir=NULL)
        clinic.plots <- list(plot=clinic.plots,text=rbind(text.box))
        return (clinic.plots)
    }
    
    p.round <- function(p.val.mat){
        p.val.mat <- unlist(lapply(p.val.mat,function(each.p){
            if (length(grep("e",each.p)) != 0){
                each.p <- unlist(strsplit(as.character(each.p),"e"))
                paste(round(as.double(each.p[1]),4),"e",each.p[2],sep="")
            }
            else round(as.double(each.p),4)
        }))
        return (p.val.mat)
    }
    
    Ratio <- NULL
    Groups <- NULL
    Genotype <- NULL
    MethylationExp <- NULL
    Conditions <- NULL
    v1 <- NULL
    v2 <- NULL
    pos <- NULL
    locus <- NULL
    Txs <- NULL
    Types <- NULL
    adj.pos <- NULL
    ori.pos <- NULL
    adj <- NULL
    types <- NULL
    id <- NULL
    
    if (length(grep("ES",CalIndex)) != 0)    testType <- "ES"
    if (length(grep("ASS",CalIndex)) != 0)    testType <- "ASS"
    if (length(grep("IR",CalIndex)) != 0)    testType <- "IR"
    
    testASmodel <- ASdb@"SplicingModel"[[testType]]
    if (ncol(testASmodel) == 1) return (NULL)
    short.long.ex <- is.element(colnames(testASmodel),c("ShortNeighborEX","LongNeighborEX"))
    short.long.des <- is.element(colnames(testASmodel),c("ShortNeighbor_des","LongNeighbor_des"))
    if (length(which(short.long.ex == TRUE)) != 0){
        exs <- as.integer(unlist(strsplit(testASmodel[,short.long.ex],"-")))
        exs.des <- paste(testASmodel[,short.long.des],collapse=",")
        testASmodel <- rbind(testASmodel[,!is.element(colnames(testASmodel),c("LongNeighborEX","LongNeighbor_des"))])
        testASmodel[,"ShortNeighborEX"] <- paste(min(exs),max(exs),sep="-")
        testASmodel[,"ShortNeighbor_des"] <- exs.des
        colnames(testASmodel)[which(colnames(testASmodel) == "ShortNeighborEX")] <- "NeighborEX"
        colnames(testASmodel)[which(colnames(testASmodel) == "ShortNeighbor_des")] <- "Neighbor_des"
    }
    
    
    exonRatio <- ASdb@"Ratio"[[testType]]
    sQTLs <- ASdb@"sQTLs"[[testType]]
    MesQTLs <- ASdb@"Me.sQTLs"[[testType]]
    Diffgroups <- ASdb@"GroupDiff"[[testType]]
    ClinicalAnal <- ASdb@"Clinical"[[testType]]
    Num.samples <- NULL
    if (length(GroupSam) != 0){
        Num.samples <- lengths(GroupSam)
        A.nums <- paste("A:",Num.samples["GroupA"],sep="")
        B.nums <- paste("B:",Num.samples["GroupB"],sep="")
        Num.samples <- paste("Total:",sum(Num.samples),"\n",A.nums,",",B.nums,sep="")
    }
    testRatio <- NULL
    testsQTLs <- NULL
    testMesQTLs <- NULL
    testDiffgroups <- NULL
    testClinical <- NULL
    
    testASmodel <- rbind(testASmodel[testASmodel[,"Index"] == CalIndex,])
    if (length(testASmodel) == 0) return (NULL)
    if (length(which(exonRatio != "NA")) != 0 & length(exonRatio) != 0){
        testRatio <- rbind(exonRatio[exonRatio[,"Index"] == CalIndex,])
    }
    if (length(which(sQTLs != "NA")) != 0    & length(sQTLs) != 0){
        testsQTLs <- rbind(sQTLs[sQTLs[,"Index"] == CalIndex,])
    }
    if (length(which(MesQTLs != "NA")) != 0    & length(MesQTLs) != 0){
        testMesQTLs <- rbind(MesQTLs[MesQTLs[,"Index"] == CalIndex,])
    }
    if (length(which(Diffgroups != "NA")) != 0    & length(Diffgroups) != 0){
        testDiffgroups <- rbind(Diffgroups[Diffgroups[,"Index"] == CalIndex,])
    }
    if (length(which(ClinicalAnal != "NA")) != 0    & length(ClinicalAnal) != 0){
        testClinical <- rbind(ClinicalAnal[ClinicalAnal[,"Index"] == CalIndex,])
    }
    
    total.inter.cn <- c("Do_des","1st_des","2nd_des","Up_des","Short_des","Long_des","Neighbor_des","Retain_des")
    total.exs <- unlist(strsplit(testASmodel[,is.element(colnames(testASmodel),total.inter.cn)],","))
    total.exs <- do.call(rbind,strsplit(total.exs,"-"))
    colnames(total.exs) <- c("start","end")
    total.exs <- total.exs[total.exs[,"start"] != "NA",]
    total.exs <- total.exs[order(as.integer(total.exs[,"start"])),]
    total.exs <- cbind(start=as.integer(total.exs[,"start"]),end=as.integer(total.exs[,"end"]))
    total.exs.range <- GRanges(seqnames=Rle("*"),ranges=IRanges(start=total.exs[,"start"],end=total.exs[,"end"]))
    re.total.exs.range <- reduce(total.exs.range)
    re.total.exs <- cbind(start=start(re.total.exs.range),end=end(re.total.exs.range))
    re.total.exs <- rbind(re.total.exs[order(re.total.exs[,"start"]),])
    num.rows <- nrow(re.total.exs)
    adj.mat <- cbind(re.total.exs,adj=0)
    adj.mat <- unique(rbind(adj.mat,cbind(total.exs,adj=0)))
    if (num.rows > 1){
        int.range <- cbind(start=as.integer(re.total.exs[1:as.integer(nrow(re.total.exs)-1),"end"]),
            end=as.integer(re.total.exs[2:nrow(re.total.exs),"start"]))
        int.ran <- int.range[,"end"] - int.range[,"start"]
        adj.value <- 10^as.integer(nchar(as.character(as.integer(mean(int.ran))))-3)
        adj.diff.value <- int.ran - as.integer(int.ran/adj.value)
        adj.mat <- lapply(1:num.rows,function(re.nums){
            each.mat <- rbind(re.total.exs[re.nums,])
            dw.nums <- c(1:num.rows)[1:num.rows < re.nums]
            over.ex <- rbind(total.exs[each.mat[,"start"] <= total.exs[,"start"] & each.mat[,"end"] >= total.exs[,"end"],])
            over.ex <- over.ex[!is.element(paste(over.ex,collapse=""),paste(each.mat,collapse="")),]
            if (length(dw.nums) != 0){
                sum.adj <- sum(adj.diff.value[dw.nums])
                if (length(over.ex) != 0){
                    cbind(over.ex,sum.adj)
                }
                else    c(each.mat,sum.adj)
            }
            else{
                if (length(over.ex) != 0){
                    cbind(over.ex,0)
                }
                else    c(each.mat,0)
            }
        })
        adj.mat <- do.call(rbind,adj.mat)
        colnames(adj.mat) <- c("start","end","adj")
    }
    adj.loci <- c(adj.mat[,"start"],adj.mat[,"end"])
    adj.loci <- cbind(ori.pos=adj.loci,adj.pos=adj.loci-adj.mat[,"adj"])
    re.total.exs <- cbind(re.total.exs,adj=adj.mat[is.element(paste(adj.mat[,"start"],adj.mat[,"end"]),
        paste(re.total.exs[,"start"],re.total.exs[,"end"])),"adj"])
    total.exs.adj <- merge(total.exs,adj.loci,by.x="start",by.y="ori.pos")
    total.exs.adj <- merge(total.exs.adj,adj.loci,by.x="end",by.y="ori.pos")
    total.exs.adj <- unique(total.exs.adj[,c("adj.pos.x","adj.pos.y")])
    total.exs.adj.range <- GRanges(seqnames=Rle("*"),ranges=IRanges(start=total.exs.adj[,"adj.pos.x"],
        end=total.exs.adj[,"adj.pos.y"]))
    dis.repre.mat <- disjoin(total.exs.adj.range)
    dis.repre.mat <- cbind(start=start(dis.repre.mat),end=end(dis.repre.mat))
    snp.granges <- NULL
    me.granges <- NULL
    omic.mat.snp <- NULL
    omic.mat.me <- NULL
    testSNPlocus <- NULL
    testMelocus <- NULL
    if (length(testsQTLs) != 0){
        snplocus <- as.matrix(snplocus)
        testSNPlocus <- rbind(snplocus[is.element(snplocus[,"SNP"],testsQTLs[,"SNP"]),])
    }
    if (length(testMesQTLs) != 0){
        methyllocus <- as.matrix(methyllocus)
        testMelocus <- rbind(methyllocus[is.element(methyllocus[,"Methyl"],testMesQTLs[,"MeID"]),])
    }
    if (length(testSNPlocus) != 0){
        adj.snp <- lapply(1:nrow(testSNPlocus),function(snp.nums){
            SNPinEx <- which(re.total.exs[,"start"] < testSNPlocus[snp.nums,"locus"] & 
                re.total.exs[,"end"] > testSNPlocus[snp.nums,"locus"])
            SNPinInt <- which(int.range[,"start"] < testSNPlocus[snp.nums,"locus"] & 
                int.range[,"end"] > testSNPlocus[snp.nums,"locus"])
            if (length(SNPinEx) != 0){
                c(testSNPlocus[snp.nums,],as.integer(testSNPlocus[snp.nums,"locus"]) - adj.mat[SNPinEx,"adj"])
            }
            else{
                adj.ex.end <- re.total.exs[SNPinInt,"end"]-re.total.exs[SNPinInt,"adj"]
                lo.per <- (as.integer(testSNPlocus[snp.nums,"locus"]) - int.range[SNPinInt,"start"]) / int.ran[SNPinInt]
                adj.lo.per <- as.integer(lo.per * as.integer(int.ran/adj.value)[SNPinInt])
                c(testSNPlocus[snp.nums,],adj.ex.end + adj.lo.per)
            }
        })
        adj.snp <- do.call(rbind,adj.snp)
        colnames(adj.snp) <- c("SNP","CHR","locus","adj")
        snp.granges <- cbind(rbind(adj.snp[,c("SNP","adj")]),types="SNP",color="black")
        omic.mat.snp <- data.frame(adj.pos=as.integer(adj.snp[,"adj"]),id=adj.snp[,"SNP"])
    }
    if (length(testMelocus) != 0){
        adj.me <- lapply(1:nrow(testMelocus),function(me.nums){
            MeinEx <- which(re.total.exs[,"start"] < testMelocus[me.nums,"locus"] & 
                re.total.exs[,"end"] > testMelocus[me.nums,"locus"])
            MeinInt <- which(int.range[,"start"] < testMelocus[me.nums,"locus"] & 
                int.range[,"end"] > testMelocus[me.nums,"locus"])
            if (length(MeinEx) != 0){
                c(as.matrix(testMelocus[me.nums,]),as.integer(testMelocus[me.nums,"locus"]) - adj.mat[MeinEx,"adj"])
            }
            else{
                adj.ex.end <- re.total.exs[MeinInt,"end"] - re.total.exs[MeinInt,"adj"]
                lo.per <- (as.integer(testMelocus[me.nums,"locus"]) - int.range[MeinInt,"start"]) / int.ran[MeinInt]
                adj.lo.per <- as.integer(lo.per * as.integer(int.ran/adj.value)[MeinInt])
                c(as.matrix(testMelocus[me.nums,]),adj.ex.end + adj.lo.per)
            }
        })
        adj.me <- do.call(rbind,adj.me)
        colnames(adj.me) <- c("Methyl","CHR","locus","adj")
        me.granges <- cbind(rbind(adj.me[,c("Methyl","adj")]),types="Methylation",color="black")
        omic.mat.me <- data.frame(adj.pos=as.integer(adj.me[,"adj"]),id=adj.me[,"Methyl"])
    }
    each.mat <- testASmodel
    inter.des <- c("1st_des","2nd_des","Do_des","Up_des","Short_des","Long_des","Neighbor_des","Retain_des")
    text.gplot<- rbind(each.mat[,!is.element(colnames(each.mat),inter.des)])
    text.gplot <- rbind(text.gplot[,text.gplot != "NA"])
    txids <- txTable[txTable[,"GENEID"] == each.mat[,"EnsID"],"TXID"]
    u.exon.range <- unlist(exon.range)
    tx.exon.range <- u.exon.range[is.element(names(u.exon.range),txids),]
    tx.exon.locus <- paste(start(tx.exon.range),end(tx.exon.range),sep="-")
    names(tx.exon.locus) <- names(tx.exon.range)
    total.txs <- length(unique(names(tx.exon.locus)))
    each.mat.exs <- strsplit(each.mat[,is.element(colnames(each.mat),inter.des)],",")
    if (testType == "ASS"){
        each.mat.exs <- list(c(each.mat.exs$Short_des,each.mat.exs$Long_des),each.mat.exs$Neighbor_des)
        names(each.mat.exs) <- c("Do_des","1st_des")
    }
    each.mat.exs <- each.mat.exs[each.mat.exs != "NA"]
    each.mat.dw <- NULL
    each.mat.1st <- NULL
    each.mat.2nd <- NULL
    each.mat.up <- NULL
    if (length(intersect(names(each.mat.exs),"Do_des")) != 0){
        each.mat.dw <- each.mat.exs[[intersect(names(each.mat.exs),"Do_des")]]
    }
    if (length(intersect(names(each.mat.exs),"1st_des")) != 0){
        each.mat.1st <- each.mat.exs[[intersect(names(each.mat.exs),"1st_des")]]
    }
    if (length(intersect(names(each.mat.exs),"2nd_des")) != 0){
        each.mat.2nd <- each.mat.exs[[intersect(names(each.mat.exs),"2nd_des")]]
    }
    if (length(intersect(names(each.mat.exs),"Up_des")) != 0) {
        each.mat.up <- each.mat.exs[[intersect(names(each.mat.exs),"Up_des")]]
    }
    IR.exs <- each.mat.exs$Retain_des
    
    total.names <- c("Do_des","1st_des","Up_des","2nd_des")
    allcomb <- do.call(rbind,lapply(1:length(each.mat.dw),function(fi.num){
        fi.re <- do.call(rbind,lapply(1:length(each.mat.1st),function(se.num){
            if (length(each.mat.up) != 0){
                se.re <- do.call(rbind,lapply(1:length(each.mat.up),function(th.num){
                    if (length(each.mat.2nd) != 0){
                        th.re <- do.call(rbind,lapply(1:length(each.mat.2nd),function(fo.num){
                            c(th.num,fo.num)
                        }))
                        c(se.num,th.re)
                    }
                    else    c(se.num,th.num)
                }))
                cbind(fi.num,se.re)
            }
            else    cbind(fi.num,se.num)
        }))
        fi.re
    }))
    allcomb <- rbind(allcomb)
    allcomb <- rbind(allcomb[,apply(allcomb,2,min) != 0])
    colnames(allcomb) <- total.names[is.element(total.names,names(each.mat.exs))]
    count.num <- 1
    final.re <- NULL
    for (each.num in 1:nrow(allcomb)){
        pre.re <- NULL
        adj.dw <- NULL
        adj.1st <- NULL
        adj.up <- NULL
        adj.2nd <- NULL
        hit.tx.dw <- NULL
        hit.tx.1st <- NULL
        hit.tx.up <- NULL
        hit.tx.2nd <- NULL
        if (is.element("Do_des",colnames(allcomb))){
            s.dw <- each.mat.dw[allcomb[each.num,"Do_des"]]
            adj.dw <- adj.mat[paste(adj.mat[,"start"],adj.mat[,"end"],sep="-") == s.dw,]
            adj.dw <- paste(adj.dw[c("start","end")] - adj.dw["adj"],collapse="-")
            hit.tx.dw <- names(tx.exon.locus[tx.exon.locus == s.dw])
        }
        if (is.element("1st_des",colnames(allcomb))){
            s.1st <- each.mat.1st[allcomb[each.num,"1st_des"]]
            adj.1st <- adj.mat[paste(adj.mat[,"start"],adj.mat[,"end"],sep="-") == s.1st,]
            adj.1st <- paste(adj.1st[c("start","end")] - adj.1st["adj"],collapse="-")
            hit.tx.1st <- names(tx.exon.locus[tx.exon.locus == s.1st])
        }
        if (is.element("Up_des",colnames(allcomb))){
            s.up <- each.mat.up[allcomb[each.num,"Up_des"]]
            adj.up <- adj.mat[paste(adj.mat[,"start"],adj.mat[,"end"],sep="-") == s.up,]
            adj.up <- paste(adj.up[c("start","end")] - adj.up["adj"],collapse="-")
            hit.tx.up <- names(tx.exon.locus[tx.exon.locus == s.up])
        }
        if (is.element("2nd_des",colnames(allcomb))){
            s.2nd <- each.mat.up[allcomb[each.num,"2nd_des"]]
            adj.2nd <- adj.mat[paste(adj.mat[,"start"],adj.mat[,"end"],sep="-") == s.2nd,]
            adj.2nd <- paste(adj.2nd[c("start","end")] - adj.2nd["adj"],collapse="-")
            hit.tx.2nd <- names(tx.exon.locus[tx.exon.locus == s.2nd])
        }
        adj.ex.loci <- c(adj.dw,adj.1st,adj.2nd,adj.up)
        over.txs.in <- NULL
        over.txs.sk <- NULL
        over.txs.d.1 <- NULL
        over.txs.d.2 <- NULL
        inter.up.dw <- intersect(hit.tx.up,hit.tx.dw)
        if (length(hit.tx.up) == 0 & length(hit.tx.2nd) == 0 & length(IR.exs) == 0){
            over.txs.in <- intersect(hit.tx.dw,hit.tx.1st)
        }
        else if (length(hit.tx.up) != 0 & length(hit.tx.2nd) == 0 & length(IR.exs) == 0){
            over.txs.d.u <- intersect(hit.tx.dw,hit.tx.up)
            over.txs.in <- intersect(hit.tx.1st,over.txs.d.u)
            over.txs.sk <- over.txs.d.u[!is.element(over.txs.d.u,names(tx.exon.locus)[is.element(tx.exon.locus,each.mat.1st)])]
        }
        else if (length(hit.tx.up) != 0 & length(hit.tx.2nd) != 0 & length(IR.exs) == 0){
            over.txs.d.1 <- intersect(hit.tx.dw,hit.tx.1st)
            over.txs.d.2 <- intersect(hit.tx.dw,hit.tx.2nd)
            over.txs.in <- intersect(over.txs.d.1,hit.tx.up)
            over.txs.sk <- intersect(over.txs.d.2,hit.tx.up)
        }
        if (length(over.txs.in) != 0){
            over.txs <- txTable[is.element(txTable[,"TXID"],over.txs.in),"TXNAME"]
            over.exs <- paste(sort(adj.ex.loci),collapse="~")
            if (length(over.txs.d.1) != 0)    over.exs <- paste(adj.dw,adj.1st,adj.up,sep="~")
            pre.re <- cbind(rep(over.exs,length(over.exs)),over.txs)
            colnames(pre.re) <- c("Ex","Txs")
            count.num <- count.num + 1
        }
        if (length(over.txs.sk) != 0){
            over.txs <- txTable[is.element(txTable[,"TXID"],over.txs.sk),"TXNAME"]
            over.exs <- paste(sort(c(adj.ex.loci[1],adj.ex.loci[3])),sep="~")
            if (length(over.txs.d.2) != 0)    over.exs <- paste(adj.dw,adj.2nd,adj.up,sep="~")
            pre.re <- rbind(pre.re,cbind(rep(over.exs,length(over.exs)),over.txs))
        }
        if (length(inter.up.dw) != 0 & length(IR.exs) != 0){
            over.exs.in <- paste(adj.ex.loci,collapse="~")
            over.exs.sk <- tx.exon.locus[is.element(tx.exon.locus,IR.exs)]
            over.txs.in <- txTable[is.element(txTable[,"TXID"],inter.up.dw),"TXNAME"]
            over.txs.sk <- txTable[is.element(txTable[,"TXID"],names(over.exs.sk)),"TXNAME"]
            pre.re <- rbind(Ex=cbind(rep(over.exs.in,length(over.txs.in)),over.txs.in),
                                            Txs=cbind(rep(over.exs.sk,length(over.txs.sk)),over.txs.sk))
            colnames(pre.re) <- c("Ex","Txs")
        }
        final.re <- rbind(final.re,pre.re)
    }
    final.re <- unique(final.re)
    pre.re <- tapply(paste(final.re[,"Ex"],final.re[,"Txs"]),final.re[,"Txs"],function(x){
        pre.each.mat <- do.call(rbind,strsplit(x," "))
        s.each.ex <- paste(unique(unlist(strsplit(pre.each.mat[,1],"~"))),collapse="~")
        c(s.each.ex,unique(pre.each.mat[,2]))
    })
    pre.re <- do.call(rbind,pre.re)
    colnames(pre.re) <- c("Ex","Txs")
    final.re <- tapply(paste(pre.re[,"Ex"],pre.re[,"Txs"]),pre.re[,"Ex"],function(x){
        pre.each.mat <- do.call(rbind,strsplit(x," "))
        s.each.ex <- unique(unlist(strsplit(pre.each.mat[,1],"~")))
        cbind(s.each.ex,paste(pre.each.mat[,2],collapse=","))
    })
    final.re <- do.call(rbind,final.re)
    colnames(final.re) <- c("Ex","Txs")
    repre.ex <- cbind(dis.repre.mat,"RepresentativeExon")
    final.re <- cbind(do.call(rbind,strsplit(final.re[,"Ex"],"-")),final.re[,-which(colnames(final.re) == "Ex")])
    final.re <- rbind(final.re,repre.ex)
    colnames(final.re) <- c("Exstart","Exend","Txs")
    final.re <- tapply(paste(final.re[,"Exstart"],final.re[,"Exend"],final.re[,"Txs"]),final.re[,"Txs"],function(x){
        pre.each.mat <- do.call(rbind,strsplit(x," "))
        unique(pre.each.mat)
    })
    final.re <- do.call(rbind,final.re)
    colnames(final.re) <- c("Exstart","Exend","Txs")
    u.final.txs <- unique(final.re[,"Txs"])
    adj.txs.name <- lapply(u.final.txs,function(each.u){
        txs <- NULL
        s.u.tx <- unlist(strsplit(each.u,","))
        num.s.u.tx <- length(s.u.tx)
        odd.even <- num.s.u.tx%%2
        if (num.s.u.tx > 2 & odd.even == 1){
            num.arr <- seq(2,num.s.u.tx,2)
            txs <- lapply(num.arr,function(x)    paste(paste(s.u.tx[c(as.integer(x-1),x)],collapse=","),",","\n",sep=""))
            txs <- paste(unlist(txs),s.u.tx[length(s.u.tx)],sep="")
        }
        else if (num.s.u.tx > 2 & odd.even == 0){
            num.arr <- seq(2,num.s.u.tx,2)
            txs <- lapply(num.arr,function(x)    paste(s.u.tx[c(as.integer(x-1),x)],collapse=","))
            txs <- paste(unlist(txs),collapse=",\n")
        }
        if (length(txs) != 0){
            unlist(txs)
            c(each.u,txs)
        }
        else NULL
    })
    adj.txs.name <- do.call(rbind,adj.txs.name)
    final.re[is.element(final.re[,"Txs"],adj.txs.name[,1]),"Txs"] <- adj.txs.name[,2]
    final.re <- data.frame(as.integer(final.re[,"Exstart"]),as.integer(final.re[,"Exend"]),final.re[,"Txs"])
    colnames(final.re) <- c("Exstart","Exend","Txs")
    min.ex <- as.integer(min(final.re[,"Exstart"])) - 200
    max.ex <- as.integer(max(final.re[,"Exend"])) + 200
    vline.mat <- cbind(unique(c(repre.ex[,c("start","end")])),"exons")
    vline.mat <- data.frame(as.integer(vline.mat[,1]),vline.mat[,2])
    colnames(vline.mat) <- c("pos","types")
    adj.loci <- rbind(adj.loci[is.element(adj.loci[,"adj.pos"],vline.mat[,"pos"]),])
    adj.loci <- data.frame(ori.pos=as.integer(adj.loci[,"ori.pos"]),adj.pos=adj.loci[,"adj.pos"])
    omics.names <- NULL
    if (length(snp.granges) != 0)    omics.names <- c("SNPid","SNP")
    if (length(me.granges) != 0)    omics.names <- c("Meid","Methylation",omics.names)
    ori.exs.mat <- each.mat[,each.mat != "NA" & is.element(colnames(each.mat),inter.des)]
    ori.exs.mat <- unlist(strsplit(ori.exs.mat,","))
    names(ori.exs.mat) <- paste(do.call(rbind,strsplit(names(ori.exs.mat),"_"))[,1],"Ex",sep="")
    ori.exs.mat <- cbind(ori.exs.mat,names(ori.exs.mat))
    colnames(ori.exs.mat) <- c("ori.Exs","ExTypes")
    adj.exs.mat <- cbind(adj.mat[,c("start","end")],adj.mat[,c("start","end")] - adj.mat[,"adj"])
    colnames(adj.exs.mat) <- c("ori.start","ori.end","adj.start","adj.end")
    adj.exs.mat <- cbind(paste(adj.exs.mat[,"ori.start"],adj.exs.mat[,"ori.end"],sep="-"),
        paste(adj.exs.mat[,"adj.start"],adj.exs.mat[,"adj.end"],sep="-"))
    colnames(adj.exs.mat) <- c("ori.Exs","adj.Exs")
    adj.exs.mat <- merge(ori.exs.mat,adj.exs.mat,by.x="ori.Exs",by.y="ori.Exs")
    mean.ranges <- do.call(rbind,strsplit(as.matrix(adj.exs.mat[,"adj.Exs"]),"-"))
    mean.ranges <- cbind(as.matrix(adj.exs.mat[,"ExTypes"]),apply(mean.ranges,1,function(x) as.integer(median(as.integer(x)))))
    adj.exs.mat <- cbind(mean.ranges,as.matrix(adj.exs.mat[,"adj.Exs"]))
    colnames(adj.exs.mat) <- c("Types","median.EX","Exs")
    text.exs <- cbind(paste(final.re[,"Exstart"],final.re[,"Exend"],sep="-"),as.matrix(final.re[,"Txs"]))
    colnames(text.exs) <- c("Exs","Txs")
    text.exs <- rbind(text.exs[text.exs[,"Txs"] != "RepresentativeExon",])
    final.text <- lapply(1:nrow(adj.exs.mat),function(each.nums){
        each.exs.mat <- adj.exs.mat[each.nums,]
        each.exs.mat <- cbind(rbind(text.exs[is.element(text.exs[,"Exs"],each.exs.mat["Exs"]),]),
        each.exs.mat["Types"],each.exs.mat["median.EX"])
        each.exs.mat
    })
    final.text <- do.call(rbind,final.text)
    colnames(final.text) <- c("Exs","Txs","Types","locus")
    final.text <- data.frame(final.text[,"Exs"],final.text[,"Txs"],final.text[,"Types"],as.integer(final.text[,"locus"]))
    colnames(final.text) <- c("Exs","Txs","Types","locus")
    total.txs <- length(unique(unlist(strsplit(as.matrix(final.re[,"Txs"]),","))))-1
    y.orders <- c("",omics.names,sort(unique(as.matrix(final.re[,"Txs"])),decreasing=TRUE))
    int.lines <- rbind(cbind(min.ex,y.orders),cbind(max.ex,y.orders))
    int.lines <- int.lines[order(int.lines[,"y.orders"]),]
    int.lines <- data.frame(as.integer(int.lines[,"min.ex"]),int.lines[,"y.orders"])
    colnames(int.lines) <- c("v1","v2")
    int.lines <- int.lines[!is.element(int.lines[,"v2"],c("SNP","SNPid","Meid","Methylation","")),]
    ggplot.mat <- ggplot() + geom_line(data=int.lines,aes(x=v1,y=v2),color="#363636")
    for (i in 1:nrow(final.re)){
        adj.tar <- rbind(adj.mat[is.element(paste(adj.mat[,"start"],adj.mat[,"end"],sep="-"),
        each.mat[,intersect(colnames(each.mat),c("1st_des","2nd_des","Short_des","Long_des"))]),])
        adj.tar.r <- rbind(adj.tar[,c("start","end")] - adj.tar[,"adj"])
        adj.tar.r <- paste(adj.tar.r[,"start"],adj.tar.r[,"end"],sep="-")
        test.tar <- is.element(paste(final.re[i,"Exstart"],final.re[i,"Exend"],sep="-"),adj.tar.r)
        ex.col <- "#999999"
        if (test.tar)    ex.col <- "#1d8b95"
        if (final.re[i,"Txs"] == "RepresentativeExon"){
            over.repre <- final.re[final.re[,"Exstart"] <= final.re[i,"Exstart"] & 
                final.re[,"Exend"] >= final.re[i,"Exend"],"Txs"]
            over.repre <- unique(unlist(strsplit(as.matrix(over.repre),",")))
            over.repre <- length(which(over.repre != "RepresentativeExon")) / as.integer(total.txs)
            ex.col <- colorRampPalette(c("red","grey"))(101)[as.integer(100-over.repre*100+1)]
        }
        exs <- data.frame(cbind(rbind(final.re[i,"Exstart"],final.re[i,"Exend"]),
            rbind(as.matrix(final.re[i,"Txs"]),as.matrix(final.re[i,"Txs"]))))
        colnames(exs) <- c("v1","v2")
        ggplot.mat <- ggplot.mat + geom_line(data=exs,aes(x=as.integer(as.matrix(v1)),y=v2,group=1),size=15,color=ex.col)
    }
    ggplot.mat <- ggplot.mat + scale_y_discrete(limits=y.orders) + 
        geom_vline(data=vline.mat,aes(xintercept=pos),color="lightgrey",linetype=2) +
        labs(x = "Genomic Coordinates", y ="", title=paste("Alternative Splicing Model for ","Index ",CalIndex,sep="")) +
        geom_text(data=final.text,aes(x=locus,y=Txs,label=Types),check_overlap=TRUE,color="white",size=3) +
        geom_text(data=adj.loci,aes(x=adj.pos,y="",label=paste(" ",ori.pos,sep="")),size=2.5,angle=90,check_overlap=TRUE) + 
        theme(axis.text.x = element_blank(),axis.text.y = element_text(size=10),axis.ticks = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "white"))
    if (length(snp.granges) != 0){
        sig.SNPs <- sQTLs[as.double(sQTLs[,"FdrByGeno"]) < 0.05,"SNP"]
        snp.granges[is.element(snp.granges[,"SNP"],sig.SNPs),"color"] <- "red"
        snp.granges <- data.frame(SNP=snp.granges[,"SNP"],adj=as.integer(snp.granges[,"adj"]),
            types=snp.granges[,"types"],color=snp.granges[,"color"])
        ggplot.mat <- ggplot.mat + geom_point(data=snp.granges,aes(x=adj,y=types,shape=types),
            color=snp.granges[,"color"],pch=17,cex=2)+
            geom_text(data=omic.mat.snp,aes(x=adj.pos,y="SNPid",label=id),angle=90,size=2)
    }
    if (length(me.granges) != 0){
        sig.MEs <- MesQTLs[as.double(MesQTLs[,"fdrByMet"]) < 0.05,"MeID"]
        me.granges[is.element(me.granges[,"Methyl"],sig.MEs),"color"] <- "red"
        me.granges <- data.frame(Methyl=me.granges[,"Methyl"],adj=as.integer(me.granges[,"adj"]),
            types=me.granges[,"types"],color=me.granges[,"color"])
        ggplot.mat <- ggplot.mat + geom_point(data=me.granges,aes(x=adj,y=types),
            color=me.granges[,"color"],shape="|",cex=2)+
            geom_text(data=omic.mat.me,aes(x=adj.pos,y="Meid",label=id),angle=90,size=2)
    }
    ggplot.mat <- list(plot=ggplot.mat,text=gsub("-","-\n",text.gplot))
    mytheme <- ttheme_default(base_size = 8, base_colour = "black",colhead=list(fg_params = list(parse=TRUE)))
    print ("Writing AS model figure into PDF format")
    pdf(file=paste(out.dir,"/",CalIndex,".pdf",sep=""),width = 2000, height = 2000,paper="a4r")#,8.27, height=11.69
    text.plot.box <- tableGrob(ggplot.mat$"text",rows=NULL,theme=ttheme_default(7))
    grid.arrange(ggplot.mat$"plot",text.plot.box,heights=c(1,0.3))
    if (length(testDiffgroups) != 0){
        print ("Writing Differential PSI level figure into PDF format")
        diff.plot <- BoxforGroup(testDiffgroups,testRatio)
        text.plot.box <- tableGrob(cbind(diff.plot$"text",Nofsam=Num.samples),theme = mytheme)
        grid.arrange(diff.plot$"plot",text.plot.box,heights=c(1,0.3),
            top=textGrob(paste("Differential PSI Levels between Conditions for Index ",CalIndex,sep=""),
            gp=gpar(fontsize=15,fontface = "bold")))
    }
    if (length(testsQTLs) != 0){
        print ("Writing sQTLs figure into PDF format")
        sig.testsQTLs <- rbind(testsQTLs[as.double(testsQTLs[,"FdrByGeno"]) < 0.05,])
        sQTLs.plot <- BoxforsQTLs(sig.testsQTLs,snpdata)
        processing <- lapply(sQTLs.plot,function(each.plot){
            if (length(each.plot) != 0){
                text.plot.box <- tableGrob(cbind(each.plot$"text",Nofsam=Num.samples),theme = mytheme)
                text.plot.box$widths <- unit(rep(1/ncol(text.plot.box), ncol(text.plot.box)), "npc")
                grid.arrange(each.plot$"plot",text.plot.box,heights=c(1,0.3),
                    top=textGrob(paste("sQTLs Results of ",each.plot$"text"[,"SNP"]," for Index ",CalIndex,sep=""),
                    gp=gpar(fontsize=15,fontface = "bold")))
            }
        })
    }
    if (length(testMesQTLs) != 0){
        print ("Writing Me-sQTLs figure into PDF format")
        sig.testMesQTLs <- rbind(testMesQTLs[as.double(testMesQTLs[,"pByMet"]) < 0.05,])
        ME.plots <- BoxforMe(sig.testMesQTLs,testRatio)
        processing <- lapply(ME.plots,function(each.plot){
            if (length(each.plot) != 0){
                text.plot.box <- tableGrob(cbind(each.plot$"text",Nofsam=Num.samples),rows = NULL,theme = mytheme)
                text.plot.box$widths <- unit(rep(1/ncol(text.plot.box), ncol(text.plot.box)), "npc")
                grid.arrange(each.plot$"plot",text.plot.box,heights=c(1,0.3),
                    top=textGrob(paste("Me-sQTLs Results of ",each.plot$"text"[,"MeID"]," for Index ",CalIndex,sep=""),
                    gp=gpar(fontsize=15,fontface = "bold")))
            }
        })
    }
    if (length(testClinical) != 0){
        print ("Writing Clinical analysis figure into PDF format")
        Cli.plots <- BoxforClinical(CalIndex)
        text.plot.box <- tableGrob(cbind(Cli.plots$"text",Nofsam=Num.samples),rows=NULL,theme = mytheme)
        text.plot.box$widths <- unit(rep(1/ncol(text.plot.box), ncol(text.plot.box)), "npc")
        grid.arrange(Cli.plots$"plot",text.plot.box,heights=c(1,0.3),
            top=textGrob(paste("Clinical Outcomes for Index ",CalIndex,sep=""),
            gp=gpar(fontsize=15,fontface = "bold")))
    }
    dev.off()
}
