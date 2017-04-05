ASvisualization <- function(ASdb,CalIndex=NULL,txTable=NULL,exon.range=NULL,
    snpdata=NULL,snplocus=NULL,methyldata=NULL,methyllocus=NULL,GroupSam=NULL,
    ClinicalInfo=NULL,out.dir=NULL){
    BoxforGroup <- function(te.gro,te.ra){
        A.nums <- is.element(colnames(te.ra),GroupSam$"GroupA")
        B.nums <- is.element(colnames(te.ra),GroupSam$"GroupB")
        GroupA.exp <- rbind(te.ra[,A.nums])
        GroupB.exp <- rbind(te.ra[,B.nums])
        A.B.exp <- c(as.double(GroupA.exp),as.double(GroupB.exp))
        A.B.groups <- c(rep("GroupA",length(GroupA.exp)),
            rep("GroupB",length(GroupB.exp)))
        T.exp <- data.frame(A.B.exp,A.B.groups)
        p.nums <- is.element(colnames(te.gro),c("Diff.P","Fdr.p"))
        p.val.mat <- te.gro[,p.nums]
        te.gro[,p.nums] <- p.round(p.val.mat)
        text.box <- rbind(te.gro[,c("Index","EnsID","Diff.P","Fdr.p")])
        colnames(T.exp) <- c("Ratio","Groups")
        gplot.result <- ggplot(data=T.exp,aes(x=Groups,y=Ratio,fill=Groups))+
            geom_boxplot(width = 0.3)+geom_jitter(width = 0.03)+ylim(c(0,1))
        gplot.result <- list(plot=gplot.result,text=text.box)
        return (gplot.result)
    }
    BoxforsQTLs <- function(snp.result=NULL,test.snp.mat=NULL){
        tb.cn <- c("Index","SNP","EnsID","pByGeno","FdrByGeno",
            "pByGroups","FdrByGroups","OR","lowCI","highCI")
        ICOR <- c("OR","lowCI","highCI")
        p.cn <- c("pByGeno","FdrByGeno","pByGroups","FdrByGroups")
        sQTLs.gplot <- function(each.result,snpid){
            colnames(each.result) <- colnames(snp.result)
            each.snp.mat <- rbind(test.snp.mat[snpid,])
            rownames(each.snp.mat) <- snpid
            box.list <- sQTLsFinder(ASdb,each.snp.mat,snplocus,
                method="boxplot",GroupSam=GroupSam,CalIndex=CalIndex)
            box.exp <- box.list$"exp"
            box.OR <- box.list$"Ratios"
            ov.cn <- is.element(colnames(each.result),tb.cn)
            text.box <- rbind(each.result[,ov.cn])
            p.val.nums <- is.element(colnames(text.box),p.cn)
            p.val.mat <- text.box[,p.val.nums]
            p.val.mat <- p.round(p.val.mat)
            text.box[,p.val.nums] <- p.val.mat
            G.test <- which(names(text.box) == "OR")
            if (!any(G.test)){
                T.Ra <- lapply(names(box.exp),function(each.geno){
                    ea.genoBox <- cbind(exp=box.exp[[each.geno]],each.geno)
                    ea.genoBox
                })
                T.Ra <- do.call(rbind,T.Ra)
                T.Ra <- data.frame(as.double(T.Ra[,"exp"]),
                    T.Ra[,"each.geno"])
                colnames(T.Ra) <- c("Ratio","Genotype")
                gplot.result <- ggplot(T.Ra,aes(x=Genotype,y=Ratio,
                    fill=Genotype))+geom_boxplot(width=0.3)+
                    geom_jitter(width = 0.03)
                gplot.result <- list(plot=gplot.result,text=text.box)
                gplot.result
            }
            else if (any(G.test)){
                T.Ra <- lapply(names(box.exp),function(each.geno){
                    ea.genoBox <- cbind(exp=box.exp[[each.geno]],each.geno)
                    ea.genoBox <- cbind(ea.genoBox,group="Group A")
                    Bn <- is.element(rownames(ea.genoBox),GroupSam$"GroupB")
                    ea.genoBox[Bn,"group"] <- "Group B"
                    ea.genoBox
                })
                T.Ra <- do.call(rbind,T.Ra)
                T.Ra <- data.frame(Ratio=as.double(T.Ra[,"exp"]),
                    Genotype=T.Ra[,"each.geno"],Groups=T.Ra[,"group"])
                T.RaA <- T.Ra[T.Ra[,"Groups"] == "Group A",]
                T.RaB <- T.Ra[T.Ra[,"Groups"] == "Group B",]
                AG <- ggplot(data=T.RaA,aes(x=Genotype,y=Ratio,fill=Genotype
                    ))+geom_boxplot(width=0.5)+geom_jitter(width=0.03)+
                    ylim(c(0,1))+ggtitle("Group A")+
                    theme(plot.title=element_text(vjust = 0))
                BG <- ggplot(data=T.RaB,aes(x=Genotype,y=Ratio,fill=Genotype
                    ))+geom_boxplot(width=0.5)+geom_jitter(width = 0.03)+
                    ylim(c(0,1))+ggtitle("Group B")+
                    theme(plot.title = element_text(vjust = 0))
                gplot.result <- arrangeGrob(AG,BG,ncol=2)
                text.box[,ICOR] <- round(as.double(text.box[,ICOR]),3)
                gplot.result <- list(plot=gplot.result,text=text.box)
                gplot.result
            }
        return (gplot.result)
        }
        ea.genoBox <- NULL
        gplot.result <- NULL
        if (!length(test.snp.mat))    return(NULL)
        snpplots.re <- lapply(rownames(test.snp.mat),function(snpid){
            each.result <- rbind(snp.result[snp.result[,"SNP"] == snpid])
            if (any(seq_along(each.result))){
                sQTLs.gplot(each.result,snpid)
            }
            else    NULL
        })
        names(snpplots.re) <- rownames(test.snp.mat)
        return (snpplots.re)
    }
    
    
    BoxforMe <- function(test.Me.mat,test.ratio){
        MeExp <- NULL
        cn.ratio <- colnames(test.ratio)
        cn.me <- colnames(methyldata)
        ov.sam <- intersect(cn.ratio,cn.me)
        tb.cn <- c("MeID","Index","EnsID","pByMet",
            "fdrByMet","pByGroups","fdrByGroups")
        Ac <- "#F8766D"
        Bc <- "#00BFC4"
        MesQTLs.gplot <- function(each.mat,each.id){
            ea.me.id <- rownames(methyldata) == each.id
            ea.ra <- rbind(methyldata[ea.me.id,ov.sam])
            me.ratio <- cbind(test.ratio[,ov.sam],ea.ra[,ov.sam])
            if (!any(seq_along(ea.ra)))    return (NULL)
            colnames(me.ratio) <- c("Ra","Me")
            text.box <- each.mat[,is.element(colnames(each.mat),tb.cn)]
            p.nums <- grep("pBy|fdrBy",names(text.box))
            text.box[p.nums] <- p.round(text.box[p.nums])
            rownames(text.box) <- NULL
            G.test <- which(names(text.box) == "pByGroups")
            if (!any(G.test)){
                me.ratio <- data.frame(as.double(me.ratio[,"Ra"]),
                    as.double(me.ratio[,"Me"]))
                colnames(me.ratio) <- c("Ratio","MeExp")
                gplot.result <- ggplot(data=me.ratio,aes(x=MeExp,
                    y=Ratio))+geom_point()+ylim(c(0,1))
                gplot.result <- list(plot=gplot.result,text=text.box)
            }
            else if (any(G.test)) {
                me.ratio <- cbind(me.ratio,group="group")
                A.nums <- is.element(rownames(me.ratio),GroupSam$"GroupA")
                B.nums <- is.element(rownames(me.ratio),GroupSam$"GroupB")
                me.ratio[A.nums,"group"] <- "GroupA"
                me.ratio[B.nums,"group"] <- "GroupB"
                colnames(me.ratio) <- c("Ra","Me","Gr")
                me.ratio <- data.frame(as.double(me.ratio[,"Ra"]),
                    as.double(me.ratio[,"Me"]),me.ratio[,"Gr"])
                colnames(me.ratio) <- c("Ratio","MeExp","Groups")
                A.me <- me.ratio[me.ratio[,"Groups"] == "GroupA",]
                B.me <- me.ratio[me.ratio[,"Groups"] == "GroupB",]
                ob1 <- coef(lm(Ratio ~ MeExp, data=A.me))
                ob2 <- coef(lm(Ratio ~ MeExp, data=B.me))
                lm.T <- coef(lm(Ratio ~ MeExp, data=me.ratio))
                gplot.resultE <- ggplot(data=me.ratio,aes(x=MeExp,
                    y=Ratio,color=Groups))+geom_point()+ylim(c(0,1))+
                    geom_abline(intercept = as.double(ob1["(Intercept)"]),
                        slope = as.double(ob1["MeExp"]),color=Ac)+
                    geom_abline(intercept = as.double(ob2["(Intercept)"]),
                        slope = as.double(ob2["MeExp"]),color=Bc)+
                    geom_abline(intercept = as.double(lm.T["(Intercept)"]),
                        slope = as.double(lm.T["MeExp"]),color="grey")
                gplot.resultG <- ggplot(data=me.ratio,aes(x=Groups,y=MeExp,
                    fill=Groups))+geom_boxplot(width = 0.3)+
                    geom_jitter(width = 0.03)
                gp.re <- arrangeGrob(gplot.resultE,gplot.resultG,ncol=2)
                gplot.result <- list(plot=gp.re,text=rbind(text.box))
                gplot.result
            }
            else    gplot.result <- NULL
            return (gplot.result)
        }
        if (!length(test.Me.mat))    return (NULL)
        meplots.re <- lapply(test.Me.mat[,"MeID"],function(each.id){
            te.id <- test.Me.mat[,"MeID"] == each.id
            ea.me <- rbind(test.Me.mat[te.id,])
            MesQTLs.gplot(ea.me,each.id)
            })
        names(meplots.re) <- test.Me.mat[,"MeID"]
        return (meplots.re)
        }
    BoxforClinical <- function(each.index){
        iep.m <- c("Index","EnsID","Pvalue","Fdr.p")
        p.m <- c("Pvalue","Fdr.p")
        ov.nm <- is.element(colnames(ClinicalAnal),iep.m)
        text.box <- ClinicalAnal[,ov.nm]
        text.box[p.m] <- p.round(text.box[p.m])
        clinic.plots <- ClinicAnalysis(ASdb,ClinicalInfo,display=TRUE,
            CalIndex=each.index,out.dir=NULL)
        clinic.plots <- list(plot=clinic.plots,text=rbind(text.box))
        return (clinic.plots)
    }
    
    p.round <- function(p.val.mat){
        p.val.mat <- unlist(lapply(p.val.mat,function(each.p){
            if (length(grep("e",each.p))){
                each.p <- unlist(strsplit(as.character(each.p),"e"))
                paste(round(as.double(each.p[1]),4),"e",each.p[2],sep="")
            }
            else round(as.double(each.p),4)
        }))
        return (p.val.mat)
    }
    omics.test <- function(t.sQTLs,t.lo,types){
        adj.gran <- NULL
        omic.mat <- NULL
        t.mat <- NULL
        f.mat <- list(adj.gran,omic.mat)
        if (!any(seq_along(t.sQTLs)))    return (f.mat)
        t.id <- gsub("Methylation","MeID",types)
        t.tp <- gsub("Methylation","Methyl",types)
        t.lo <- as.matrix(t.lo)
        o.nm <- is.element(t.lo[,t.tp],t.sQTLs[,t.id])
        t.mat <- rbind(t.lo[o.nm,])
        if (!any(seq_along(t.mat)))    return (f.mat)
        adj.val <- lapply(seq_len(nrow(t.mat)),function(e.nm){
            TinEx <- which(re.T.ex[,"start"] < t.mat[e.nm,"locus"]
                &re.T.ex[,"end"] > t.mat[e.nm,"locus"])
            TinInt <- which(int.range[,"start"] < t.mat[e.nm,"locus"]
                &int.range[,"end"] > t.mat[e.nm,"locus"])
            if (any(seq_along(TinEx))){
                int.lo <- as.integer(t.mat[e.nm,"locus"])
                dif.le <- int.lo - adj.mat[TinEx,"adj"]
                c(t.mat[e.nm,],dif.le)
            }
            else {
                adj.en <- re.T.ex[TinInt,"end"] - re.T.ex[TinInt,"adj"]
                lo.per <- (as.integer(t.mat[e.nm,"locus"]) - 
                    int.range[TinInt,"start"]) / int.ran[TinInt]
                adj.lo <- lo.per*as.integer(int.ran/adj.value)[TinInt]
                c(t.mat[e.nm,],adj.en+as.integer(adj.lo))
            }
        })
        adj.val <- do.call(rbind,adj.val)
        colnames(adj.val) <- c(t.tp,"CHR","locus","adj")
        adj.gran <- cbind(rbind(adj.val[,c(t.tp,"adj")]),
            types=types,color="black")
        omic.mat <- data.frame(as.integer(adj.val[,"adj"]),adj.val[,t.tp])
        colnames(omic.mat) <- c("adj.pos","id")
        f.mat <- list(adj.gran,omic.mat)
        return (f.mat)
    }
    test.na <- function(t.mat){
        te1 <- any(which(t.mat != "NA"))
        te2 <- any(length(t.mat))
        TRa <- NULL
        if (te1 & te2){
            TRa <- rbind(t.mat[t.mat[,"Index"] == CalIndex,])
        }
        return (TRa)
    }
    adjFun <- function(ea.nm,allcomb){
        pre.re <- NULL
        adj.dw <- NULL
        adj.1st <- NULL
        adj.up <- NULL
        adj.2nd <- NULL
        hit.tx.dw <- NULL
        hit.tx.1st <- NULL
        hit.tx.up <- NULL
        hit.tx.2nd <- NULL
        over.txs.in <- NULL
        over.txs.sk <- NULL
        over.txs.d.1 <- NULL
        over.txs.d.2 <- NULL
        p.mat <- paste(adj.mat[,"start"],adj.mat[,"end"],sep="-")
        if (is.element("Do_des",colnames(allcomb))){
            s.dw <- each.mat.dw[allcomb[ea.nm,"Do_des"]]
            adj.dw <- adj.mat[p.mat == s.dw,]
            adj.dw <- paste(adj.dw[s.e] - adj.dw["adj"],collapse="-")
            hit.tx.dw <- names(tx.exon.locus[tx.exon.locus == s.dw])
        }
        if (is.element("1st_des",colnames(allcomb))){
            s.1st <- each.mat.1st[allcomb[ea.nm,"1st_des"]]
            adj.1st <- adj.mat[p.mat == s.1st,]
            adj.1st <- paste(adj.1st[s.e] - adj.1st["adj"],collapse="-")
            hit.tx.1st <- names(tx.exon.locus[tx.exon.locus == s.1st])
        }
        if (is.element("Up_des",colnames(allcomb))){
            s.up <- each.mat.up[allcomb[ea.nm,"Up_des"]]
            adj.up <- adj.mat[p.mat == s.up,]
            adj.up <- paste(adj.up[s.e] - adj.up["adj"],collapse="-")
            hit.tx.up <- names(tx.exon.locus[tx.exon.locus == s.up])
        }
        if (is.element("2nd_des",colnames(allcomb))){
            s.2nd <- each.mat.up[allcomb[ea.nm,"2nd_des"]]
            adj.2nd <- adj.mat[p.mat == s.2nd,]
            adj.2nd <- paste(adj.2nd[s.e] - adj.2nd["adj"],collapse="-")
            hit.tx.2nd <- names(tx.exon.locus[tx.exon.locus == s.2nd])
        }
        adj.ex.loci <- c(adj.dw,adj.1st,adj.2nd,adj.up)
        inter.up.dw <- intersect(hit.tx.up,hit.tx.dw)
        if (!length(hit.tx.up) & !length(hit.tx.2nd) & !length(IR.exs)){
            over.txs.in <- intersect(hit.tx.dw,hit.tx.1st)
        }
        else if (length(hit.tx.up) & !length(hit.tx.2nd) & !length(IR.exs)){
            over.txs.d.u <- intersect(hit.tx.dw,hit.tx.up)
            over.txs.in <- intersect(hit.tx.1st,over.txs.d.u)
            over.lo.1 <- is.element(tx.exon.locus,each.mat.1st)
            o.f <- names(tx.exon.locus)[over.lo.1]
            over.txs.sk <- over.txs.d.u[!is.element(over.txs.d.u,o.f)]
        }
        else if (length(hit.tx.up) & length(hit.tx.2nd) & !length(IR.exs)){
            over.txs.d.1 <- intersect(hit.tx.dw,hit.tx.1st)
            over.txs.d.2 <- intersect(hit.tx.dw,hit.tx.2nd)
            over.txs.in <- intersect(over.txs.d.1,hit.tx.up)
            over.txs.sk <- intersect(over.txs.d.2,hit.tx.up)
        }
        if (length(over.txs.in)){
            ov.txid <- is.element(txTable[,"TXID"],over.txs.in)
            over.txs <- txTable[ov.txid,"TXNAME"]
            over.exs <- paste(sort(adj.ex.loci),collapse="~")
            if (length(over.txs.d.1)){
                over.exs <- paste(adj.dw,adj.1st,adj.up,sep="~")
            }
            pre.re <- cbind(rep(over.exs,length(over.exs)),over.txs)
            colnames(pre.re) <- c("Ex","Txs")
            count.num <- count.num + 1
        }
        if (length(over.txs.sk)){
            ov.txid <- is.element(txTable[,"TXID"],over.txs.sk)
            over.txs <- txTable[ov.txid,"TXNAME"]
            over.exs <- paste(sort(c(adj.ex.loci[1],adj.ex.loci[3])),sep="~")
            if (length(over.txs.d.2)){
                over.exs <- paste(adj.dw,adj.2nd,adj.up,sep="~")
            }
            ov.ex.l <- length(over.exs)
            pre.re <- rbind(pre.re,cbind(rep(over.exs,ov.ex.l),over.txs))
        }
        if (length(inter.up.dw) & length(IR.exs)){
            over.e.in <- paste(adj.ex.loci,collapse="~")
            over.e.sk <- tx.exon.locus[is.element(tx.exon.locus,IR.exs)]
            ov.tx.in <- is.element(txTable[,"TXID"],inter.up.dw)
            ov.tx.sk <- is.element(txTable[,"TXID"],names(over.e.sk))
            over.txs.in <- txTable[ov.tx.in,"TXNAME"]
            over.txs.sk <- txTable[ov.tx.sk,"TXNAME"]
            ov.t.in <- cbind(rep(over.e.in,length(over.txs.in)),over.txs.in)
            ov.t.sk <- cbind(rep(over.e.sk,length(over.txs.sk)),over.txs.sk)
            pre.re <- rbind(ov.t.in,ov.t.sk)
            colnames(pre.re) <- c("Ex","Txs")
        }
        return (pre.re)
    }

    PEnv <- environment(p.round)
    p.round <- PEnv$p.round
    adjEnv <- environment(adjFun)
    adjFun <- adjEnv$adjFun
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
    nei.ex <- c("ShortNeighborEX","LongNeighborEX")
    nei.de <- c("ShortNeighbor_des","LongNeighbor_des")
    s.e <- c("start","end")

    if (any(grep("ES",CalIndex)))    testType <- "ES"
    if (any(grep("ASS",CalIndex)))    testType <- "ASS"
    if (any(grep("IR",CalIndex)))    testType <- "IR"
    
    testASmodel <- ASdb@"SplicingModel"[[testType]]
    if (ncol(testASmodel) == 1) return (NULL)
    short.long.ex <- is.element(colnames(testASmodel),nei.ex)
    short.long.des <- is.element(colnames(testASmodel),nei.de)
    if (any(which(short.long.ex == TRUE))){
        exs <- as.integer(unlist(strsplit(testASmodel[,short.long.ex],"-")))
        exs.des <- paste(testASmodel[,short.long.des],collapse=",")
        longex <- c("LongNeighborEX","LongNeighbor_des")
        lo.de <- is.element(colnames(testASmodel),longex)
        testASmodel <- rbind(testASmodel[,!lo.de])
        testASmodel[,"ShortNeighborEX"] <- paste(min(exs),max(exs),sep="-")
        testASmodel[,"ShortNeighbor_des"] <- exs.des
        Ac <- colnames(testASmodel)
        colnames(testASmodel)[Ac == "ShortNeighborEX"] <- "NeighborEX"
        colnames(testASmodel)[Ac == "ShortNeighbor_des"] <- "Neighbor_des"
    }
    
    
    exonRatio <- ASdb@"Ratio"[[testType]]
    sQTLs <- ASdb@"sQTLs"[[testType]]
    MesQTLs <- ASdb@"Me.sQTLs"[[testType]]
    Diffgroups <- ASdb@"GroupDiff"[[testType]]
    ClinicalAnal <- ASdb@"Clinical"[[testType]]
    Num.samples <- NULL
    if (any(seq_along(GroupSam))){
        Num.samples <- lengths(GroupSam)
        A.nums <- paste("A:",Num.samples["GroupA"],sep="")
        B.nums <- paste("B:",Num.samples["GroupB"],sep="")
        Num.samples <- paste("Total:",sum(Num.samples),
            "\n",A.nums,",",B.nums,sep="")
    }
    testRatio <- NULL
    testsQTLs <- NULL
    testMesQTLs <- NULL
    testDiffgroups <- NULL
    testClinical <- NULL
    
    testASmodel <- rbind(testASmodel[testASmodel[,"Index"] == CalIndex,])
    if (!any(seq_along(testASmodel))) return (NULL)
    testRatio <- test.na(exonRatio)
    testsQTLs <- test.na(sQTLs)
    testMesQTLs <- test.na(MesQTLs)
    testDiffgroups <- test.na(Diffgroups)
    testClinical <- test.na(ClinicalAnal)
    
    total.inter.cn <- c("Do_des","1st_des","2nd_des","Up_des",
        "Short_des","Long_des","Neighbor_des","Retain_des")
    ov.cn <- is.element(colnames(testASmodel),total.inter.cn)
    total.exs <- unlist(strsplit(testASmodel[,ov.cn],","))
    total.exs <- do.call(rbind,strsplit(total.exs,"-"))
    colnames(total.exs) <- s.e
    total.exs <- total.exs[total.exs[,"start"] != "NA",]
    total.exs <- total.exs[order(as.integer(total.exs[,"start"])),]
    exs.st <- as.integer(total.exs[,"start"])
    exs.en <- as.integer(total.exs[,"end"])
    total.exs <- cbind(start=exs.st,end=exs.en)
    to.ran <- IRanges(start=exs.st,end=exs.en)
    total.exs.range <- GRanges(seqnames=Rle("*"),ranges=to.ran)
    T.exs.ran <- reduce(total.exs.range)
    re.T.ex <- cbind(start=start(T.exs.ran),end=end(T.exs.ran))
    re.T.ex <- rbind(re.T.ex[order(re.T.ex[,"start"]),])
    num.rows <- nrow(re.T.ex)
    adj.mat <- cbind(re.T.ex,adj=0)
    adj.mat <- unique(rbind(adj.mat,cbind(total.exs,adj=0)))
    if (num.rows > 1){
        st.int <- as.integer(re.T.ex[1:as.integer(nrow(re.T.ex)-1),"end"])
        en.int <- as.integer(re.T.ex[2:nrow(re.T.ex),"start"])
        int.range <- cbind(start=st.int,end=en.int)
        int.ran <- int.range[,"end"] - int.range[,"start"]
        mean.int <- as.integer(mean(int.ran))
        adj.value <- 10^as.integer(nchar(as.character(mean.int))-3)
        adj.diff.value <- int.ran - as.integer(int.ran/adj.value)
        adj.mat <- lapply(1:num.rows,function(re.nums){
            each.mat <- rbind(re.T.ex[re.nums,])
            dw.nums <- c(1:num.rows)[1:num.rows < re.nums]
            te.st <- each.mat[,"start"] <= total.exs[,"start"]
            te.en <- each.mat[,"end"] >= total.exs[,"end"]
            over.ex <- rbind(total.exs[te.st & te.en,])
            over.ex <- over.ex[!is.element(paste(over.ex,collapse=""),
                paste(each.mat,collapse="")),]
            if (any(length(dw.nums))){
                sum.adj <- sum(adj.diff.value[dw.nums])
                if (any(length(over.ex))){
                    cbind(over.ex,sum.adj)
                }
                else    c(each.mat,sum.adj)
            }
            else{
                if (any(length(over.ex))){
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
    p.adj <- paste(adj.mat[,"start"],adj.mat[,"end"])
    p.tex <- paste(re.T.ex[,"start"],re.T.ex[,"end"])
    re.T.ex <- cbind(re.T.ex,adj=adj.mat[is.element(p.adj,p.tex),"adj"])
    total.adj <- merge(total.exs,adj.loci,by.x="start",by.y="ori.pos")
    total.adj <- merge(total.adj,adj.loci,by.x="end",by.y="ori.pos")
    total.adj <- unique(total.adj[,c("adj.pos.x","adj.pos.y")])
    adj.int.range <- IRanges(total.adj[,"adj.pos.x"],total.adj[,"adj.pos.y"])
    total.exs.adj.range <- GRanges(seqnames=Rle("*"),ranges=adj.int.range)
    dis.repre.mat <- disjoin(total.exs.adj.range)
    dis.repre.mat <- cbind(start=start(dis.repre.mat),end=end(dis.repre.mat))
    SNPre <- omics.test(testsQTLs,snplocus,"SNP")
    Mere <-omics.test(testMesQTLs,methyllocus,"Methylation")
    snp.granges <- SNPre[[1]]
    omic.mat.snp <- SNPre[[2]]
    me.granges <- Mere[[1]]
    omic.mat.me <- Mere[[2]]
    
    each.mat <- testASmodel
    i.des <- c("1st_des","2nd_des","Do_des","Up_des",
        "Short_des","Long_des","Neighbor_des","Retain_des")
    text.gplot<- rbind(each.mat[,!is.element(colnames(each.mat),i.des)])
    text.gplot <- rbind(text.gplot[,text.gplot != "NA"])
    txids <- txTable[txTable[,"GENEID"] == each.mat[,"EnsID"],"TXID"]
    u.exon.range <- unlist(exon.range)
    tx.exon.range <- u.exon.range[is.element(names(u.exon.range),txids),]
    tx.exon.locus <- paste(start(tx.exon.range),end(tx.exon.range),sep="-")
    names(tx.exon.locus) <- names(tx.exon.range)
    total.txs <- length(unique(names(tx.exon.locus)))
    e.mat.exs <- each.mat[,is.element(colnames(each.mat),i.des)]
    e.mat.exs <- strsplit(e.mat.exs,",")
    if (testType == "ASS"){
        sh.lo.de <- c(e.mat.exs$Short_des,e.mat.exs$Long_des)
        e.mat.exs <- list(sh.lo.de,e.mat.exs$Neighbor_des)
        names(e.mat.exs) <- c("Do_des","1st_des")
    }
    e.mat.exs <- e.mat.exs[e.mat.exs != "NA"]
    each.mat.dw <- NULL
    each.mat.1st <- NULL
    each.mat.2nd <- NULL
    each.mat.up <- NULL
    ex.nms <- names(e.mat.exs)
    if (any(seq_along(intersect(names(e.mat.exs),"Do_des")))){
        each.mat.dw <- e.mat.exs[[intersect(ex.nms,"Do_des")]]
    }
    if (any(seq_along(intersect(names(e.mat.exs),"1st_des")))){
        each.mat.1st <- e.mat.exs[[intersect(ex.nms,"1st_des")]]
    }
    if (any(seq_along(intersect(names(e.mat.exs),"2nd_des")))){
        each.mat.2nd <- e.mat.exs[[intersect(ex.nms,"2nd_des")]]
    }
    if (any(seq_along(intersect(names(e.mat.exs),"Up_des")))){
        each.mat.up <- e.mat.exs[[intersect(ex.nms,"Up_des")]]
    }
    IR.exs <- e.mat.exs$Retain_des
    T.nm <- c("Do_des","1st_des","Up_des","2nd_des")
    fi.num <- seq_len(length(each.mat.dw))
    se.num <- seq_len(length(each.mat.1st))
    th.num <- seq_len(length(each.mat.up))
    fo.num <- seq_len(length(each.mat.2nd))
    if (!any(fi.num))  fi.num <- 0
    if (!any(se.num))  se.num <- 0
    if (!any(th.num))  th.num <- 0
    if (!any(fo.num))  fo.num <- 0
    allcomb <- rbind(expand.grid(fi.num,se.num,th.num,fo.num))
    allcomb <- rbind(allcomb[,apply(allcomb,2,min) != 0])
    colnames(allcomb) <- T.nm[is.element(T.nm,names(e.mat.exs))]
    count.num <- 1
    final.re <- NULL
    each.re <- NULL
    final.re <- lapply(seq_len(nrow(allcomb)),function(ea.nm){
        each.re <- adjFun(ea.nm,allcomb)
        each.re
    })
    final.re <- do.call(rbind,final.re)
    p.f.re <- paste(final.re[,"Ex"],final.re[,"Txs"])
    pre.re <- tapply(p.f.re,final.re[,"Txs"],function(x){
        pre.each.mat <- do.call(rbind,strsplit(x," "))
        s.mat <- unlist(strsplit(pre.each.mat[,1],"~"))
        s.each.ex <- paste(unique(s.mat),collapse="~")
        c(s.each.ex,unique(pre.each.mat[,2]))
    })
    pre.re <- do.call(rbind,pre.re)
    colnames(pre.re) <- c("Ex","Txs")
    p.p.re <- paste(pre.re[,"Ex"],pre.re[,"Txs"])
    final.re <- tapply(p.p.re,pre.re[,"Ex"],function(x){
        pre.each.mat <- do.call(rbind,strsplit(x," "))
        s.each.ex <- unlist(strsplit(pre.each.mat[,1],"~"))
        s.each.ex <- unique(s.each.ex)
        cbind(s.each.ex,paste(pre.each.mat[,2],collapse=","))
    })
    final.re <- do.call(rbind,final.re)
    colnames(final.re) <- c("Ex","Txs")
    repre.ex <- cbind(dis.repre.mat,"RepresentativeExon")
    f.ex.mat <- do.call(rbind,strsplit(final.re[,"Ex"],"-"))
    final.re <- cbind(f.ex.mat,final.re[,-which(colnames(final.re) == "Ex")])
    final.re <- rbind(final.re,repre.ex)
    colnames(final.re) <- c("Exstart","Exend","Txs")
    p.f.re <- paste(final.re[,"Exstart"],final.re[,"Exend"],final.re[,"Txs"])
    final.re <- tapply(p.f.re,final.re[,"Txs"],function(x){
        pre.each.mat <- do.call(rbind,strsplit(x," "))
        unique(pre.each.mat)
    })
    final.re <- do.call(rbind,final.re)
    colnames(final.re) <- c("Exstart","Exend","Txs")
    u.final.txs <- unique(final.re[,"Txs"])
    p.x <- NULL
    txs <- NULL
    adj.txs.name <- lapply(u.final.txs,function(each.u){
        s.u.tx <- unlist(strsplit(each.u,","))
        N.s.tx <- length(s.u.tx)
        odd.even <- N.s.tx%%2
        if (N.s.tx > 2) {
            num.arr <- seq(2,N.s.tx,2)
            txs <- lapply(num.arr,function(x){
                p.x <- paste(s.u.tx[c(as.integer(x-1),x)],collapse=",")
                if (odd.even == 1)    p.x <- paste(p.x,",","\n",sep="")
                })
            if (odd.even == 1){
                txs <- paste(unlist(txs),s.u.tx[N.s.tx],sep="")
            }
            else if (odd.even == 0){
                txs <- paste(unlist(txs),collapse=",\n")
            }
            c(each.u,txs)
        }
    })
    
    adj.txs.name <- do.call(rbind,adj.txs.name)
    ov.tx.nm <- is.element(final.re[,"Txs"],adj.txs.name[,1])
    final.re[ov.tx.nm,"Txs"] <- adj.txs.name[,2]
    final.re <- data.frame(as.integer(final.re[,"Exstart"]),
        as.integer(final.re[,"Exend"]),final.re[,"Txs"])
    colnames(final.re) <- c("Exstart","Exend","Txs")
    min.ex <- as.integer(min(final.re[,"Exstart"])) - 200
    max.ex <- as.integer(max(final.re[,"Exend"])) + 200
    vline.mat <- cbind(unique(c(repre.ex[,c("start","end")])),"exons")
    vline.mat <- data.frame(as.integer(vline.mat[,1]),vline.mat[,2])
    colnames(vline.mat) <- c("pos","types")
    v.adj.ov <- is.element(adj.loci[,"adj.pos"],vline.mat[,"pos"])
    adj.loci <- rbind(adj.loci[v.adj.ov,])
    adj.loci <- data.frame(ori.pos=as.integer(adj.loci[,"ori.pos"]),
        adj.pos=adj.loci[,"adj.pos"])
    omics.names <- NULL
    if (length(snp.granges)){
        omics.names <- c("SNPid","SNP")
    }
    if (length(me.granges)){
        omics.names <- c("Meid","Methylation",omics.names)
    }
    ea.des <- is.element(colnames(each.mat),i.des)
    ori.exs.mat <- each.mat[,each.mat != "NA" & ea.des]
    ori.exs.mat <- unlist(strsplit(ori.exs.mat,","))
    s.exs <- do.call(rbind,strsplit(names(ori.exs.mat),"_"))[,1]
    names(ori.exs.mat) <- paste(s.exs,"Ex",sep="")
    ori.exs.mat <- cbind(ori.exs.mat,names(ori.exs.mat))
    colnames(ori.exs.mat) <- c("orExs","ExTypes")
    adj.exs.mat <- cbind(adj.mat[,s.e],adj.mat[,s.e] - adj.mat[,"adj"])
    colnames(adj.exs.mat) <- c("ori.start","ori.end","adj.start","adj.end")
    op <- paste(adj.exs.mat[,"ori.start"],adj.exs.mat[,"ori.end"],sep="-")
    ap <- paste(adj.exs.mat[,"adj.start"],adj.exs.mat[,"adj.end"],sep="-")
    adj.exs.mat <- cbind(op,ap)
    colnames(adj.exs.mat) <- c("orExs","adj.Exs")
    adj.exs.mat <- merge(ori.exs.mat,adj.exs.mat,by.x="orExs",by.y="orExs")
    mn.ran <- do.call(rbind,strsplit(as.matrix(adj.exs.mat[,"adj.Exs"]),"-"))
    me.r <- apply(mn.ran,1,function(x) as.integer(median(as.integer(x))))
    mn.ran <- cbind(as.matrix(adj.exs.mat[,"ExTypes"]),me.r)
    adj.exs.mat <- cbind(mn.ran,as.matrix(adj.exs.mat[,"adj.Exs"]))
    colnames(adj.exs.mat) <- c("Types","median.EX","Exs")
    p.final.re <- paste(final.re[,"Exstart"],final.re[,"Exend"],sep="-")
    text.exs <- cbind(p.final.re,as.matrix(final.re[,"Txs"]))
    colnames(text.exs) <- c("Exs","Txs")
    text.exs <- rbind(text.exs[text.exs[,"Txs"] != "RepresentativeExon",])
    each.exs.mat <- NULL
    final.text <- lapply(1:nrow(adj.exs.mat),function(each.nums){
        each.exs.mat <- adj.exs.mat[each.nums,]
        te.ea.ov <- is.element(text.exs[,"Exs"],each.exs.mat["Exs"])
        each.exs.mat <- cbind(rbind(text.exs[te.ea.ov,]),
        each.exs.mat["Types"],each.exs.mat["median.EX"])
        each.exs.mat
    })
    final.text <- do.call(rbind,final.text)
    colnames(final.text) <- c("Exs","Txs","Types","locus")
    final.text <- data.frame(final.text[,"Exs"],final.text[,"Txs"],
        final.text[,"Types"],as.integer(final.text[,"locus"]))
    colnames(final.text) <- c("Exs","Txs","Types","locus")
    total.txs <- unique(unlist(strsplit(as.matrix(final.re[,"Txs"]),",")))
    total.txs <- length(total.txs) - 1
    y.orders <- sort(unique(as.matrix(final.re[,"Txs"])),decreasing=TRUE)
    y.orders <- c("",omics.names,y.orders)
    int.lines <- rbind(cbind(min.ex,y.orders),cbind(max.ex,y.orders))
    int.lines <- int.lines[order(int.lines[,"y.orders"]),]
    int.lines <- data.frame(as.integer(int.lines[,"min.ex"]),
        int.lines[,"y.orders"])
    colnames(int.lines) <- c("v1","v2")
    omic.cn <- c("SNP","SNPid","Meid","Methylation","")
    int.lines <- int.lines[!is.element(int.lines[,"v2"],omic.cn),]
    ggplot.mat <- ggplot() + 
    geom_line(data=int.lines,aes(x=v1,y=v2),color="#363636")
    des.cn <- c("1st_des","2nd_des","Short_des","Long_des")
    p.adj <- paste(adj.mat[,"start"],adj.mat[,"end"],sep="-")
    ov.des <- intersect(colnames(each.mat),des.cn)
    adj.tar <- rbind(adj.mat[is.element(p.adj,each.mat[,ov.des]),])
    adj.tar.r <- rbind(adj.tar[,c("start","end")] - adj.tar[,"adj"])
    adj.tar.r <- paste(adj.tar.r[,"start"],adj.tar.r[,"end"],sep="-")
    T.l <- as.integer(total.txs)
    Col.mat <- colorRampPalette(c("red","grey"))(101)
    for (i in 1:nrow(final.re)){
        ef.re <- rbind(final.re[i,])
        p.e.re <- paste(ef.re[,"Exstart"],ef.re[,"Exend"],sep="-")
        test.tar <- is.element(p.e.re,adj.tar.r)
        ex.col <- "#999999"
        if (test.tar)    ex.col <- "#1d8b95"
        if (ef.re[,"Txs"] == "RepresentativeExon"){
            over.re <- final.re[final.re[,"Exstart"] <= ef.re[,"Exstart"]
                & final.re[,"Exend"] >= ef.re[,"Exend"],"Txs"]
            over.re <- unique(unlist(strsplit(as.matrix(over.re),",")))
            over.re <- length(which(over.re != "RepresentativeExon")) / T.l
            ex.col <- Col.mat[as.integer(100-over.re*100+1)]
        }
        exs <- data.frame(cbind(rbind(ef.re[,"Exstart"],ef.re[,"Exend"]),
            rbind(as.matrix(ef.re[,"Txs"]),as.matrix(ef.re[,"Txs"]))))
        colnames(exs) <- c("v1","v2")
        ggplot.mat <- ggplot.mat + 
            geom_line(data=exs,aes(x=as.integer(as.matrix(v1)),y=v2,group=1),
            size=15,color=ex.col)
    }
    spli.nm <- paste("Alternative Splicing Model for Index ",CalIndex,sep="")
    the.sp <- theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),axis.ticks = element_blank(),
        panel.grid.minor=element_blank(),panel.background=element_blank(),
        axis.line = element_line(colour = "white"))
    ggplot.mat <- ggplot.mat + scale_y_discrete(limits=y.orders) + 
        geom_vline(data=vline.mat,aes(xintercept=pos),color="lightgrey",
        linetype=2) + labs(x = "Genomic Coordinates", y ="", title=spli.nm) +
        geom_text(data=final.text,aes(x=locus,y=Txs,label=Types),
            check_overlap=TRUE,color="white",size=3) +
        geom_text(data=adj.loci,aes(x=adj.pos,y="",label=paste(" ",ori.pos,
            sep="")),size=2.5,angle=90,check_overlap=TRUE) + the.sp
    if (any(length(snp.granges))){
        si.SNP <- sQTLs[as.double(sQTLs[,"FdrByGeno"]) < 0.05,"SNP"]
        snp.granges[is.element(snp.granges[,"SNP"],si.SNP),"color"] <- "red"
        snp.granges <- data.frame(SNP=snp.granges[,"SNP"],
            adj=as.integer(snp.granges[,"adj"]),types=snp.granges[,"types"],
            color=snp.granges[,"color"])
        ggplot.mat <- ggplot.mat + 
            geom_point(data=snp.granges,aes(x=adj,y=types,shape=types),
            color=snp.granges[,"color"],pch=17,cex=2)+geom_text(data=
            omic.mat.snp,aes(x=adj.pos,y="SNPid",label=id),angle=90,size=2)
    }
    if (any(length(me.granges))){
        si.ME <- MesQTLs[as.double(MesQTLs[,"fdrByMet"]) < 0.05,"MeID"]
        me.granges[is.element(me.granges[,"Methyl"],si.ME),"color"] <- "red"
        me.granges <- data.frame(Methyl=me.granges[,"Methyl"],
            adj=as.integer(me.granges[,"adj"]),
            types=me.granges[,"types"],color=me.granges[,"color"])
        ggplot.mat <- ggplot.mat + 
            geom_point(data=me.granges,aes(x=adj,y=types),color=
            me.granges[,"color"],shape="|",cex=2)+geom_text(data=omic.mat.me,
            aes(x=adj.pos,y="Meid",label=id),angle=90,size=2)
    }
    ggplot.mat <- list(plot=ggplot.mat,text=gsub("-","-\n",text.gplot))
    mytheme <- ttheme_default(base_size = 8, base_colour = "black",
        colhead=list(fg_params = list(parse=TRUE)))
    mytheme <- ttheme_default(7)
    print ("Writing AS model figure into PDF format")
    pdf(file=paste(out.dir,"/",CalIndex,".pdf",sep=""),
        width = 2000, height = 2000,paper="a4r")
    text.p.box <- tableGrob(ggplot.mat$"text",NULL,theme=ttheme_default(7))
    grid.arrange(ggplot.mat$"plot",text.p.box,heights=c(1,0.3))
    if (any(length(testDiffgroups))){
        dnm <- paste("Differential PSI Levels between Conditions for Index ",
            CalIndex,sep="")
        print ("Writing Differential PSI level figure into PDF format")
        diff.plot <- BoxforGroup(testDiffgroups,testRatio)
        text.mat <- cbind(diff.plot$"text",Nofsam=Num.samples)
        text.plot.box <- tableGrob(text.mat,theme = mytheme)
        grid.arrange(diff.plot$"plot",text.plot.box,heights=c(1,0.3),
            top=textGrob(dnm,gp=gpar(fontsize=15,fontface = "bold")))
    }
    if (any(length(testsQTLs))){
        print ("Writing sQTLs figure into PDF format")
        s.sQTLs <- rbind(testsQTLs[as.double(testsQTLs[,"FdrByGeno"])<0.05,])
        sQTLs.plot <- BoxforsQTLs(s.sQTLs,snpdata)
        processing <- lapply(sQTLs.plot,function(each.plot){
            if (length(each.plot)){
                sq.nm <- paste("sQTLs Results of ",each.plot$"text"[,"SNP"],
                    " for Index ",CalIndex,sep="")
                text.mat <- cbind(each.plot$"text",Nofsam=Num.samples)
                text.plot.box <- tableGrob(text.mat,theme = mytheme)
                text.plot.box$widths <- unit(rep(1/ncol(text.plot.box), 
                    ncol(text.plot.box)), "npc")
                grid.arrange(each.plot$"plot",text.plot.box,heights=c(1,0.3),
                    top=textGrob(sq.nm,gp=gpar(fontsize=15,fontface="bold")))
            }
        })
    }
    if (any(length(testMesQTLs))){
        print ("Writing Me-sQTLs figure into PDF format")
        s.Me <- rbind(testMesQTLs[as.double(testMesQTLs[,"pByMet"]) < 0.05,])
        ME.plots <- BoxforMe(s.Me,testRatio)
        processing <- lapply(ME.plots,function(each.plot){
            if (length(each.plot)){
                me.nm <- paste("Me-sQTLs Results of ",
                    each.plot$"text"[,"MeID"]," for Index ",CalIndex,sep="")
                text.mat <- cbind(each.plot$"text",Nofsam=Num.samples)
                text.p.box <- tableGrob(text.mat,rows = NULL,theme = mytheme)
                text.p.box$widths <- unit(rep(1/ncol(text.p.box),
                    ncol(text.p.box)), "npc")
                grid.arrange(each.plot$"plot",text.p.box,heights=c(1,0.3),
                    top=textGrob(me.nm,gp=gpar(fontsize=15,fontface="bold")))
            }
        })
    }
    if (any(length(testClinical))){
        print ("Writing Clinical analysis figure into PDF format")
        Cli.plots <- BoxforClinical(CalIndex)
        text.mat <- cbind(Cli.plots$"text",Nofsam=Num.samples)
        text.plot.box <- tableGrob(text.mat,rows=NULL,theme = mytheme)
        text.plot.box$widths <- unit(rep(1/ncol(text.plot.box),
            ncol(text.plot.box)), "npc")
        grid.arrange(Cli.plots$"plot",text.plot.box,heights=c(1,0.3),
            top=textGrob(paste("Clinical Outcomes for Index ",
            CalIndex,sep=""),gp=gpar(fontsize=15,fontface = "bold")))
    }
    dev.off()
}
