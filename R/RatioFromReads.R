RatioFromReads <- function(ASdb=NULL,Total.bamfiles=NULL,readsInfo=
    c("paired","single"),readLen=NULL,inserSize=NULL,minr=3,CalIndex=NULL,
    Ncor=1,out.dir=NULL){
    splitSplice <- function(EX1,EX2){
        s.EX1 <- strsplit(unlist(strsplit(EX1,",")),"-")
        s.EX2 <- strsplit(unlist(strsplit(EX2,",")),"-")
        s.EX1 <- matrix(as.integer(do.call(rbind,s.EX1)),ncol=2)
        s.EX2 <- matrix(as.integer(do.call(rbind,s.EX2)),ncol=2)
        colnames(s.EX1) <- c("start","end")
        colnames(s.EX2) <- c("start","end")
        final.re <- unlist(lapply(s.EX1[,"end"],function(each.fi.s.EX){
            pre.re <- lapply(s.EX2[,"start"],function(each.se.s.EX){
                paste(sort(c(each.fi.s.EX,each.se.s.EX)),collapse="-")
            })
            do.call(rbind,pre.re)
        }))
        return (unique(final.re))
    }
    coorEX <- function(spl.re,g.Info,readLen,inseSize,min.r=5,AStype){
        Normalized.values <- function(exons.l){
            normal.max.in <- 0
            normal.max.skip <- 0
            normal.in <- 0
            normal.skip <- 0
            jun.in.normal <- 0
            jun.skip.normal <- 0
            if (AStype == "ES" | AStype == "MXE"){
                max.p.reads <- exons.l[1:3] + (read.l-1)
                fi.in.value <- abs(min((1+inse+read.l)-(exons.l[1]+1),0))
                fi.nex.value <- abs(min(exons.l[2]-(read.l + inse),0))
                se.in.value <- abs(min((1+inse+read.l)-(exons.l[2]+1),0))
                se.nex.value <- abs(min(exons.l[3]-(read.l + inse),0))
                fi.in.minus <- fi.in.value + fi.nex.value
                se.in.minus <- se.in.value + se.nex.value
                se.fi <- abs(min(abs(read.l+inse-exons.l[2])-(read.l-1),0))
                sk.in.value <-abs(min((1+inse+read.l)-(exons.l[1]+1),0))
                sk.nex.value <- abs(min(exons.l[3]-(read.l + inse),0))
                sk.minus.value <- sk.in.value + sk.nex.value
                normal.fi.in <- max((max.p.reads[1] - fi.in.minus),0)
                normal.se.in <- max((max.p.reads[2] - se.in.minus),0)
                normal.max.in <- max(normal.fi.in + normal.se.in - se.fi,0)
                normal.max.skip <- max(max.p.reads[1] - sk.minus.value,0)
                exon.r.a <- abs(min(1+inse-exons.l[1],0))
                exon.r.b <- abs(min(exons.l[2]-(read.l+inse),0))
                exon.r.c <- abs(min((1+inse)-(exons.l[2]+exons.l[3]),0))
                exon.r.d <- abs(min(read.l+inse-exons.l[2],0))
                exon.r.e <- abs(min((inse+1)-exons.l[3],0))
                jun.pair.a <- abs(min(read.l-1,exon.r.a))
                jun.pair.b <- abs(min(read.l-1,exons.l[1],exon.r.b))
                jun.pair.c <- abs(min(read.l-1,exon.r.c))
                jun.pair.d <- abs(min(read.l-1,exon.r.d))
                jun.pair.e <- abs(min(read.l-1,exon.r.e))
                jun.sk.pa <- sum(jun.pair.a,jun.pair.e)
                jun.in.pa <- sum(jun.sk.pa,jun.pair.b,jun.pair.c,jun.pair.d)
                normal.in <- max(normal.max.in - jun.in.pa,0)
                normal.skip <- max(normal.max.skip - jun.sk.pa,0)
                jun.in.normal <- max(2*(2*(read.l - 2*an.size + 1)) - 
                    (2*abs(min(exons.l[2] - (read.l - 2 + 1),0)) + 
                    abs(min(abs(read.l+inse-exons.l[2]) - (read.l-1),0))),0)
                jun.skip.normal <- max(2*(read.l - 2*an.size + 1),0)
            }
            else if (AStype == "ASS"){
                normal.in <- 2*exons.l[2]
                jun.in.normal <- max(2*(read.l - 2*an.size + 1),0)
                jun.skip.normal <- max(2*(read.l - 2*an.size + 1),0)
            }
            else if (AStype == "IR"){
                jun.in.normal <- max(2*(read.l - 2*an.size + 1),0)
                normal.skip <- 2*(exons.l[2] + (read.l-1))
            }
            return.mat <- c(normal.max.in,normal.max.skip,
                normal.in,normal.skip,jun.in.normal,jun.skip.normal)
            cn <- c("pair.in","pair.sk","pairwojun.in",
                "pairwojun.sk","jun.in","jun.sk")
            names(return.mat) <- cn
            return (return.mat)
        }
        g.1.info <- g.Info[[1]]
        g.2.info <- g.Info[[2]]
        paired.r <- spl.re$pairedInfo
        exon.r <- spl.re$exonInfo
        junction.r <- spl.re$junctionInfo
        inse <- inseSize
        total.exon.l <- 0
        if (sum(c(paired.r,exon.r,junction.r)) == 0){
            return ("NA")
        }
        total.exons <- do.call(rbind,strsplit(names(exon.r),"-"))
        total.exon.l <- apply(total.exons,1,function(x) diff(as.double(x)+1))
        g.1.p.r <- paired.r[is.element(names(paired.r),g.1.info$paired)]
        g.2.p.r <- paired.r[is.element(names(paired.r),g.2.info$paired)]
        g.1.e.r <- exon.r[is.element(names(exon.r),g.1.info$exon[2])]
        g.2.e.r <- 0
        j.nm <- names(junction.r)
        g.1.j.r <- junction.r[is.element(j.nm,g.1.info$junction)]
        g.2.j.r <- junction.r[is.element(j.nm,g.2.info$junction)]

        if (!length(g.1.p.r) & !length(g.2.p.r) & !length(g.1.e.r)){
            if (any(length(g.1.j.r)) & any(length(g.2.j.r))){
                g.1.p.r <- 0
                g.2.p.r <- 0
                g.1.e.r <- 0
            }
            else    return ("NA")
        }
        if (!any(length(g.1.e.r))){
            g.1.e.r <- 0
        }
        if (!any(length(g.1.j.r))) g.1.j.r <- 0
        if (!any(length(g.2.j.r))) g.2.j.r <- 0
        read.l <- as.integer(readLen)
        inse <- inseSize
        an.size=1
        normal.values <- Normalized.values(total.exon.l)
        if (AStype == "ESse") return ("NA")
        if (AStype == "ES" | AStype == "MXE"){
            if (AStype == "MXE"){
                fi.normal.values <- Normalized.values(total.exon.l[c(1,2,4)])
                se.normal.values <- Normalized.values(total.exon.l[c(1,3,4)])
                fi.normal.values[c("pair.sk","pairwojun.sk","jun.sk")] <- 
                    se.normal.values[c("pair.in","pairwojun.in","jun.in")]
                normal.values <- fi.normal.values
            }
            total.reads <- sum(g.1.p.r,g.1.j.r,g.2.p.r,g.2.j.r)
            g.1.nor.num <- sum(normal.values[c("pairwojun.in","jun.in")])
            g.2.nor.num <- sum(normal.values[c("pairwojun.sk","jun.sk")])
            group.1.read.count <- sum(g.1.p.r,g.1.j.r)
            group.2.read.count <- sum(g.2.p.r,g.2.j.r)
            read.num.test.1 <- sum(g.1.p.r) == 0 & sum(g.1.j.r) == 0
            read.num.test.2 <- sum(g.2.p.r) == 0 & sum(g.2.j.r) == 0
            if ((read.num.test.1) | (read.num.test.2) | total.reads < min.r){
                return ("NA")
                }
            if (sum(g.1.j.r) == 0 | sum(g.2.j.r) == 0){
                g.1.nor.num <- normal.values["pair.in"]
                g.2.nor.num <- normal.values["pair.sk"]
                group.1.read.count <- sum(g.1.p.r)
                group.2.read.count <- sum(g.2.p.r)
            }
            else if (sum(g.1.p.r) == 0 | sum(g.2.p.r) == 0){
                g.1.nor.num <- normal.values["jun.in"]
                g.2.nor.num <- normal.values["jun.sk"]
                group.1.read.count <- sum(g.1.j.r)
                group.2.read.count <- sum(g.2.j.r)
            }
        }
        else if (AStype == "IR"){
            g.1.nor.num <- normal.values["jun.in"]
            g.2.nor.num <- normal.values["pairwojun.sk"]
            group.1.read.count <- sum(g.1.j.r)
            group.2.read.count <- sum(g.1.e.r)
        }
        else if (AStype == "ASS"){
            g.1.nor.num <- normal.values["jun.in"]
            g.2.nor.num <- normal.values["jun.sk"]
            group.1.read.count <- sum(g.1.j.r)
            group.2.read.count <- sum(g.2.j.r)
            if(g.1.e.r > 5){
                g.1.nor.num <- sum(normal.values[c("pairwojun.in","jun.in")])
                group.1.read.count <- sum(g.1.e.r,g.1.j.r)
            }
        }
        ratio.g.1 <- group.1.read.count/g.1.nor.num
        ratio.g.2 <- group.2.read.count/g.2.nor.num
        group.total <- sum(ratio.g.1,ratio.g.2)
        group.1.2.ratio <- (group.2.read.count/g.2.nor.num)/group.total
        if (group.1.read.count == 0 | group.2.read.count == 0)  return ("NA")
        if (group.1.read.count < min.r & group.2.read.count < min.r){
            group.1.2.ratio <- "NA"
            }
        return (group.1.2.ratio)
    }
    splitEnv <- environment(splitSplice)
    coorEnv <- environment(coorEX)
    Each.Cal.ratio <- function(bamfiles=NULL,splicingInfo=NULL){
        ExReads <- function(t.ex,t.sp,g.e1,g.e2,g1,g2,s1,s2,ch,er,met,alt){
            coor.re <- NULL
            inse <- inserSize
            pre.bam.re <- lapply(seq_along(bamfiles),function(ebam){
                T.r <- SplicingReads(bamfiles[ebam],t.ex,t.sp,er,ch,met,inse)
                group.1.list <- list(g1,g.e1,s1)
                group.2.list <- list(g2,g.e2,s2)
                names(group.1.list) <- c("paired","exon","junction")
                names(group.2.list) <- c("paired","exon","junction")
                t.g.li <- list(group.1.list,group.2.list)
                names(t.g.li) <- c("Inclu","Skip")
                coor.re <- coorEnv$coorEX(T.r,t.g.li,readLen,inserSize,minr,alt)
                coor.re
            })
        pre.bam.re <- do.call(cbind,pre.bam.re)    
        pre.bam.re[is.na(pre.bam.re)] <- "NA"
        return (pre.bam.re)
        }
        if (!any(length(splicingInfo))) return (NULL)
        final.ES.result <- NULL
        final.ASS.result <- NULL
        final.IR.result <- NULL
        typenm <- names(splicingInfo)[lengths(splicingInfo) > 1]
        if (is.element("ES",typenm)){
            f.fi.result <- NULL
            f.se.result <- NULL
            f.MXE.result <- NULL
            each.re <- splicingInfo$"ES"
            ES.f.te <- each.re[,"2ndEX"] == "NA" & each.re[,"Types"] == "ES"
            ES.s.te <- each.re[,"2ndEX"] != "NA" & each.re[,"Types"] == "ES"
            MXE.te <- each.re[,"2ndEX"] != "NA" & each.re[,"Types"] == "MXE"
            ES.fi.result <- rbind(each.re[ES.f.te,])
            ES.se.result <- rbind(each.re[ES.s.te,])
            MXE.result <- rbind(each.re[MXE.te,])
            do.st <- do.call(rbind,strsplit(each.re[,"DownEX"],"-"))[,1]
            up.en <- do.call(rbind,strsplit(each.re[,"UpEX"],"-"))[,2]
            position.range <- rbind(unique(cbind(do.st,up.en)))
            colnames(position.range) <- c("start","end")
            if (any(seq_len(nrow(ES.fi.result)))){
                ES.fi.r <- foreach(ES.num=seq_len(nrow(ES.fi.result)),
                    .packages=called.packages,.combine=rbind) %dopar% {
                    ES.re <- rbind(ES.fi.result[ES.num,])
                    e.dw.st <- do.call(rbind,strsplit(ES.re[,"DownEX"],"-"))
                    e.up.en <- do.call(rbind,strsplit(ES.re[,"UpEX"],"-"))
                    each.ran <- rbind(unique(cbind(e.dw.st[,1],e.up.en[,2])))
                    fi.pr <- paste(ES.re[,c("DownEX","1stEX")],collapse="~")
                    se.pr <- paste(ES.re[,c("1stEX","UpEX")],collapse="~")
                    g.1.p <- rbind(fi.pr,se.pr)
                    g.2.p <- paste(ES.re[,c("DownEX","UpEX")],collapse="~")
                    fi.s1 <- splitEnv$splitSplice(ES.re[,"Do_des"],ES.re[,"1st_des"])
                    fi.s2 <- splitEnv$splitSplice(ES.re[,"1st_des"],ES.re[,"Up_des"])
                    g.1.s <- c(fi.s1,fi.s2)
                    g.2.s <- splitEnv$splitSplice(ES.re[,"Do_des"],ES.re[,"Up_des"])
                    t.ex <- ES.re[,c("DownEX","1stEX","UpEX")]
                    t.sp <- c(g.1.s,g.2.s)
                    e.chr <- unique(ES.re[,"Nchr"])
                    pr.re <- ExReads(t.ex,t.sp,t.ex,t.ex,g.1.p,g.2.p,
                        g.1.s,g.2.s,e.chr,each.ran,readsInfo,"ES")
                    pr.re
                }
                fi.re <- cbind(rbind(ES.fi.result[,c("Index","EnsID","Nchr",
                    "1stEX","2ndEX","DownEX","UpEX","Types")]),ES.fi.r)
                colnames(fi.re) <- c("Index","EnsID","Nchr","1stEX",
                    "2ndEX","DownEX","UpEX","Types",sample.names)
                f.fi.result <- fi.re
            }
            if (any(seq_len(nrow(ES.se.result)))){
                ES.se.r <- foreach(ES.num=seq_len(nrow(ES.se.result)),
                    .packages=called.packages,.combine=rbind) %dopar% {
                    ES.re <- rbind(ES.se.result[ES.num,])
                    do.st <- strsplit(ES.re[,"DownEX"],"-")
                    up.en <- strsplit(ES.re[,"UpEX"],"-")
                    do.st <- do.call(rbind,do.st)[,1]
                    up.en <- do.call(rbind,up.en)[,2]
                    each.ran <- rbind(unique(cbind(do.st,up.en)))
                    fi.pr <- paste(ES.re[,c("DownEX","1stEX")],collapse="~")
                    se.pr <- paste(ES.re[,c("1stEX","2ndEX")],collapse="~")
                    th.pr <- paste(ES.re[,c("2ndEX","UpEX")],collapse="~")
                    g.1.p <- rbind(fi.pr,se.pr,th.pr)
                    g.2.p <- paste(ES.re[,c("DownEX","UpEX")],collapse="~")
                    fi.s1 <- splitEnv$splitSplice(ES.re[,"Do_des"],ES.re[,"1st_des"])
                    fis2 <- splitEnv$splitSplice(ES.re[,"1st_des"],ES.re[,"2nd_des"])
                    fi.s3 <- splitEnv$splitSplice(ES.re[,"2nd_des"],ES.re[,"Up_des"])
                    g.1.s <- c(fi.s1,fis2,fi.s3)
                    g.2.s <- splitEnv$splitSplice(ES.re[,"Do_des"],ES.re[,"Up_des"])
                    t.sp <- c(g.1.s,g.2.s)
                    t.ex <- c(ES.re[,"DownEX"],ES.re[,"1stEX"],
                        ES.re[,"2ndEX"],ES.re[,"UpEX"])
                    e.chr <- unique(ES.re[,"Nchr"])
                    pr.re <- ExReads(t.ex,t.sp,t.ex,t.ex,g.1.p,g.2.p,
                        g.1.s,g.2.s,e.chr,each.ran,readsInfo,"ESse")
                    pr.re
                }
                se.re <- cbind(rbind(ES.se.result[,c("Index","EnsID","Nchr",
                    "1stEX","2ndEX","DownEX","UpEX","Types")]),ES.se.r)
                colnames(se.re) <- c("Index","EnsID","Nchr","1stEX",
                    "2ndEX","DownEX","UpEX","Types",sample.names)
                f.se.result <- se.re
            }
            if (any(seq_len(nrow(MXE.result)))){
                MXE.ratio <- foreach(ES.num=seq_len(nrow(MXE.result)),
                    .packages=called.packages,.combine=rbind) %dopar% {
                    ES.re <- rbind(MXE.result[ES.num,])
                    do.st <- do.call(rbind,strsplit(ES.re[,"DownEX"],"-"))
                    up.en <- do.call(rbind,strsplit(ES.re[,"UpEX"],"-"))
                    do.st <- do.st[,1]
                    up.en <- up.en[,2]
                    each.ran <- rbind(unique(cbind(do.st,up.en)))
                    fi.p.1 <- paste(ES.re[,c("DownEX","1stEX")],collapse="~")
                    fi.p.2 <- paste(ES.re[,c("DownEX","2ndEX")],collapse="~")
                    se.p.1 <- paste(ES.re[,c("1stEX","UpEX")],collapse="~")
                    se.p.2 <- paste(ES.re[,c("2ndEX","UpEX")],collapse="~")
                    g.1.p <- rbind(fi.p.1,se.p.1)
                    g.2.p <- rbind(fi.p.2,se.p.2)
                    fi.s1 <- splitEnv$splitSplice(ES.re[,"Do_des"],ES.re[,"1st_des"])
                    fi.s2 <- splitEnv$splitSplice(ES.re[,"1st_des"],ES.re[,"Up_des"])
                    se.s1 <- splitEnv$splitSplice(ES.re[,"Do_des"],ES.re[,"2nd_des"])
                    se.s2 <- splitEnv$splitSplice(ES.re[,"2nd_des"],ES.re[,"Up_des"])
                    g.1.s <- c(fi.s1,fi.s2)
                    g.2.s <- c(se.s1,se.s2)
                    t.sp <- c(g.1.s,g.2.s)
                    g.1.e <- c(ES.re[,"DownEX"],ES.re[,"1stEX"],
                        ES.re[,"UpEX"])
                    g.2.e <- c(ES.re[,"DownEX"],ES.re[,"2ndEX"],
                        ES.re[,"UpEX"])
                    t.ex <- ES.re[,c("DownEX","1stEX","2ndEX","UpEX")]
                    e.chr <- unique(ES.re[,"Nchr"])
                    pr.re <- ExReads(t.ex,t.sp,g.1.e,g.2.e,g.1.p,g.2.p,
                        g.1.s,g.2.s,e.chr,each.ran,readsInfo,"MXE")
                    pr.re
                }
                MXE.re <- cbind(rbind(MXE.result[,c("Index","EnsID","Nchr",
                    "1stEX","2ndEX","DownEX","UpEX","Types")]),MXE.ratio)
                colnames(MXE.re) <- c("Index","EnsID","Nchr","1stEX",
                    "2ndEX","DownEX","UpEX","Types",sample.names)
                f.MXE.result <- MXE.re
            }
            final.ES.result <- rbind(f.fi.result,f.se.result,f.MXE.result)
            if (!any(length(final.ES.result)))    final.ES.result <- NULL
            else    rownames(final.ES.result) <- 1:nrow(final.ES.result)
        }
        if (is.element("IR",typenm)){
            IR.re <- splicingInfo[["IR"]]
            do.st <- do.call(rbind,strsplit(IR.re[,"DownEX"],"-"))[,1]
            up.en <- do.call(rbind,strsplit(IR.re[,"UpEX"],"-"))[,2]
            position.range <- rbind(unique(cbind(do.st,up.en)))
            colnames(position.range) <- c("start","end")
            IR.ratio <- foreach(IR.num=seq_len(nrow(IR.re)),
                .packages=called.packages,.combine=rbind) %dopar% {
                ea.re <- rbind(IR.re[IR.num,])
                s.down <- strsplit(ea.re[,"DownEX"],"-")
                s.up <- strsplit(ea.re[,"UpEX"],"-")
                ex.sp <- paste(unlist(s.down)["DownEX2"],
                    unlist(s.up)["UpEX1"],sep="-")
                do.st <- do.call(rbind,s.down)[,1]
                up.en <- do.call(rbind,s.up)[,2]
                each.ran <- rbind(unique(cbind(do.st,up.en)))
                g.1.p <- paste(ea.re[,c("DownEX","UpEX")],collapse="~")
                g.1.p <- rbind(g.1.p)
                fi.pr.2 <- paste(c(ea.re[,"DownEX"],ex.sp),collapse="~")
                se.pr.2 <- paste(c(ex.sp,IR.re[,"UpEX"]),collapse="~")
                g.2.p <- rbind(fi.pair.2,se.pair.2)
                g.1.s <- c(splitEnv$splitSplice(ea.re[,"Do_des"],ea.re[,"Up_des"]))
                g.2.s <- "NA"
                t.ex <- c(ea.re[,"DownEX"],ex.sp,ea.re[,"UpEX"])
                e.chr <- unique(ea.re[,"Nchr"])
                pr.re <- ExReads(t.ex,g.1.s,t.ex,t.ex,g.1.p,g.2.p,
                        g.1.s,g.2.s,e.chr,each.ran,"exon","IR")
                pr.re
            }
            IR.result <- cbind(rbind(IR.result[,c("Index","EnsID","Nchr",
                "RetainEX","DownEX","UpEX","Types")]),IR.ratio)
            colnames(IR.result) <- c("Index","EnsID","Nchr","RetainEX",
                "DownEX","UpEX","Types",sample.names)
            final.IR.result <- IR.result
            if (!any(length(final.IR.result)))    final.IR.result <- NULL
            else    rownames(final.IR.result) <- 1:nrow(final.IR.result)
        }
        if (is.element("ASS",typenm)){
            sh.nm <- c("ShortNeighborEX","NeighborEX")
            lo.nm <- c("LongNeighborEX","NeighborEX")
            sh.des.nm <- c("ShortNeighbor_des","NeighborEX")
            lo.des.nm <- c("LongNeighbor_des","NeighborEX")
            final.A5SS.result <- NULL
            final.A3SS.result <- NULL
            ASS.result <- splicingInfo[["ASS"]]
            ASS.num <- grep("A[0-9]SS",ASS.result[,"Types"])
            ASS.re <- rbind(ASS.result[ASS.num,])
            ASS.ratio <- foreach(ASS.num=seq_len(nrow(ASS.re)),
                .packages=called.packages,.combine=rbind) %dopar% {
                ea.ty <- ASS.result[ASS.num,"Types"]
                ea.re <- rbind(ASS.result[ASS.num,])
                nei.info <- grep("NeighborEX",colnames(ea.re))
                ASS.nei <- unlist(strsplit(ea.re[,nei.info],"-"))
                sh.ex <- unlist(strsplit(ea.re[,"ShortEX"],"-"))
                lo.ex <- unlist(strsplit(ea.re[,"LongEX"],"-"))
                sh.ex.1 <- sh.ex["ShortEX1"]
                sh.ex.2 <- sh.ex["ShortEX2"]
                lo.ex.1 <- lo.ex["LongEX1"]
                lo.ex.2 <- lo.ex["LongEX2"]
                ASS.cn <- colnames(ea.re)
                lo.nei <- ea.re[,is.element(ASS.cn,lo.nm)]
                sh.nei <- ea.re[,is.element(ASS.cn,sh.nm)]
                lo.des <- ea.re[,is.element(ASS.cn,lo.des.nm)]
                sh.des <- ea.re[,is.element(ASS.cn,sh.des.nm)]
                if (ea.ty == "A5SS"){
                    min.r <- max(as.integer(sh.ex)) - 20
                    max.r <- min(as.integer(c(ASS.nei))) + 20
                    ASS.sp <- paste(sh.ex.2,lo.ex.2,sep="-")
                    names(ASS.sp) <- "Alt.ex"
                    g.2.p <- c(ea.re[,"ShortEX"],sh.nei)
                    g.2.p <- rbind(paste(g.2.p,collapse="~"))
                    g.1.s <- c(splitEnv$splitSplice(ea.re[,"Long_des"],lo.des))
                    g.2.s <- c(splitEnv$splitSplice(ea.re[,"Short_des"],sh.des))
                    g.1.e <- c(ASS.sp,sh.nei)
                    g.2.e <- c(ASS.sp,lo.nei)
                    t.ex <- c(ASS.sp,lo.nei)
                    }
                else if (ea.ty == "A3SS"){
                    max.r <- min(as.integer(sh.ex)) + 20
                    min.r <- max(as.integer(c(ASS.nei))) - 20
                    ASS.sp <- paste(lo.ex.1,sh.ex.1,sep="-")
                    names(ASS.sp) <- "Alt.ex"
                    g.2.p <- c(sh.nei,ea.re[,"ShortEX"])
                    g.2.p <- rbind(paste(g.2.p,collapse="~"))
                    g.1.s <- c(splitEnv$splitSplice(lo.des,ea.re[,"Long_des"]))
                    g.2.s <- c(splitEnv$splitSplice(sh.des,ea.re[,"Short_des"]))
                    g.1.e <- c(sh.nei,ASS.sp)
                    g.2.e <- c(lo.nei,ASS.sp)
                    t.ex <- c(lo.nei,ASS.sp)
                    }
                each.ran <- cbind(min.r,max.r)
                colnames(each.ran) <- c("start","end")
                g.1.p <- NULL
                t.sp <- c(g.1.s,g.2.s)
                e.chr <- unique(ea.re[,"Nchr"])
                pr.re <- ExReads(t.ex,g.1.s,g.1.e,g.2.e,g.1.p,g.2.p,
                        g.1.s,g.2.s,e.chr,each.ran,"exon","ASS")
                pr.re
            }
            te.cn <- "Index|EnsID|Nchr|ShortEX|LongEX|NeighborEX|Types"
            cn <- colnames(ASS.result)[grep(te.cn,colnames(ASS.result))]
            ASS.result <- cbind(rbind(ASS.result[,cn]),ASS.ratio)
            colnames(ASS.result) <- c(cn,sample.names)
            final.ASS.result <- ASS.result
            if (!any(length(final.ASS.result)))    final.ASS.result <- NULL
            else    rownames(final.ASS.result) <- 1:nrow(final.ASS.result)
        }
        final.re <- list(final.ES.result,final.ASS.result,final.IR.result)
        names(final.re) <- c("ES","ASS","IR")
        return (final.re)
    }
    ea.re <- NULL
    ES.num <- NULL
    IR.num <- NULL
    ASS.num <- NULL
    g.1.p <- NULL
    g.2.p <- NULL
    fi.pr <- NULL
    se.pr <- NULL
    pr.re <- NULL
    fi.pair.2 <- NULL
    se.pair.2 <- NULL
    final.re <- NULL
    registerDoParallel(cores=Ncor)
    called.packages <- c("GenomicRanges","GenomicFeatures")
    sample.files <- rbind(Total.bamfiles[,"path"])
    sample.names <- rbind(Total.bamfiles[,"names"])
    test.mat <- list(NULL,NULL,NULL)
    names(test.mat) <- c("ES","ASS","IR")
    T.spl <- ASdb@"SplicingModel"
    if (length(CalIndex)){
        if (ncol(T.spl$"ES") != 1){
            ES.re <- T.spl$"ES"[is.element(T.spl$"ES"[,"Index"],CalIndex),]
            test.mat$"ES" <- rbind(ES.re)
        }
        if (ncol(T.spl$"ASS") != 1){
            AS.re <- T.spl$"ASS"[is.element(T.spl$"ASS"[,"Index"],CalIndex),]
            test.mat$"ASS" <- rbind(AS.re)
        }
        if (ncol(T.spl$"IR") != 1){
            IR.re <- T.spl$"IR"[is.element(T.spl$"IR"[,"Index"],CalIndex),]
            test.mat$"IR" <- rbind(IR.re)
        }
    }
    else    test.mat <- T.spl
    final.ra <- Each.Cal.ratio(sample.files,test.mat)
    if (!any(length(final.ra$"ES")))    final.ra$"ES" <- as.matrix("NA")
    if (!any(length(final.ra$"ASS")))    final.ra$"ASS" <- as.matrix("NA")
    if (!any(length(final.ra$"IR")))    final.ra$"IR" <- as.matrix("NA")
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=final.ra,
        GroupDiff=ASdb@"GroupDiff",sQTLs=ASdb@"sQTLs",
        Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
    if (length(out.dir)){
        p.out <- paste(out.dir,"/AS_Ratio/",sep="")
        system(paste("mkdir -p ",p.out,sep=""))
        write.table(final.ra[["ES"]],
            paste(p.out,"ES_Ratio.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.ra[["ASS"]],
            paste(p.out,"ASS_Ratio.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.ra[["IR"]],
            paste(p.out,"IR_Ratio.txt",sep=""),sep='\t',quote=FALSE)
    }
    return(ASdb)
}

