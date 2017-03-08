RatioFromReads <- function(ASdb=NULL,Total.bamfiles=NULL,readsInfo=c("paired","single"),readLen=NULL,inserSize=NULL,minr=3,CalIndex=NULL,Ncor=1,out.dir=NULL){
    split.splice <- function(EX1,EX2){
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
    coorEX <- function(readsSplcing=NULL,groupInfo=NULL,readLen=NULL,inseSize=NULL,min.read=5,AStype=NULL,readsinfo=NULL){
        Normalized.values <- function(exons.l){
            normal.max.in <- 0
            normal.max.skip <- 0
            normal.in <- 0
            normal.skip <- 0
            jun.in.normal <- 0
            jun.skip.normal <- 0
            if (AStype == "ES" | AStype == "MXE"){
                max.paired.reads <- exons.l[1:3] + (read.l-1)
                fi.in.value <- abs(min((1+inse+read.l)-(exons.l[1]+1),0))
                fi.nex.value <- abs(min(exons.l[2]-(read.l + inse),0))
                se.in.value <- abs(min((1+inse+read.l)-(exons.l[2]+1),0))
                se.nex.value <- abs(min(exons.l[3]-(read.l + inse),0))
                fi.in.minus.value <- fi.in.value + fi.nex.value
                se.in.minus.value <- se.in.value + se.nex.value
                over.se.fi <- abs(min(abs(read.l+inse-exons.l[2])-(read.l-1),0))
                sk.in.value <-abs(min((1+inse+read.l)-(exons.l[1]+1),0))
                sk.nex.value <- abs(min(exons.l[3]-(read.l + inse),0))
                sk.minus.value <- sk.in.value + sk.nex.value
                normal.fi.in <- max((max.paired.reads[1] - fi.in.minus.value),0)
                normal.se.in <- max((max.paired.reads[2] - se.in.minus.value),0)
                normal.max.in <- max(normal.fi.in + normal.se.in - over.se.fi,0)
                normal.max.skip <- max(max.paired.reads[1] - sk.minus.value,0)
                jun.pair.a <- abs(min(read.l-1,abs(min(1+inse-exons.l[1],0))))
                jun.pair.b <- abs(min(read.l-1,exons.l[1],abs(min(exons.l[2]-(read.l+inse),0))))
                jun.pair.c <- abs(min(read.l-1,abs(min((1+inse)-(exons.l[2]+exons.l[3]),0))))
                jun.pair.d <- abs(min(read.l-1,abs(min(read.l+inse-exons.l[2],0))))
                jun.pair.e <- abs(min(read.l-1,abs(min((inse+1)-exons.l[3],0))))
                jun.in.pair <- sum(jun.pair.a,jun.pair.b,jun.pair.c,jun.pair.d,jun.pair.e)
                jun.skip.pair <- sum(jun.pair.a,jun.pair.e)
                normal.in <- max(normal.max.in - jun.in.pair,0)
                normal.skip <- max(normal.max.skip - jun.skip.pair,0)
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
            return.mat <- c(normal.max.in,normal.max.skip,normal.in,normal.skip,jun.in.normal,jun.skip.normal)
            names(return.mat) <- c("pair.in","pair.sk","pairwojun.in","pairwojun.sk","jun.in","jun.sk")
            return (return.mat)
        }
        group.1.info <- groupInfo[[1]]
        group.2.info <- groupInfo[[2]]
        paired.r <- readsSplcing$pairedInfo
        exon.r <- readsSplcing$exonInfo
        junction.r <- readsSplcing$junctionInfo
        inse <- inseSize
        total.exon.l <- 0
        if (sum(paired.r) == 0 & sum(exon.r) == 0 & sum(junction.r) == 0) return ("NA")
        total.exons <- do.call(rbind,strsplit(names(exon.r),"-"))
        total.exon.l <- apply(total.exons,1,function(x) diff(as.double(x)+1))
        group.1.p.r <- paired.r[is.element(names(paired.r),group.1.info$paired)]
        group.2.p.r <- paired.r[is.element(names(paired.r),group.2.info$paired)]
        group.1.e.r <- exon.r[is.element(names(exon.r),group.1.info$exon[2])]
        group.2.e.r <- 0
        group.1.j.r <- junction.r[is.element(names(junction.r),group.1.info$junction)]
        group.2.j.r <- junction.r[is.element(names(junction.r),group.2.info$junction)]
        if (length(group.1.p.r) == 0 & length(group.2.p.r) == 0 & length(group.1.e.r) == 0){
            if (length(group.1.j.r) != 0 & length(group.2.j.r) != 0){
                group.1.p.r <- 0
                group.2.p.r <- 0
                group.1.e.r <- 0
            }
            else    return ("NA")
        }
        if (length(group.1.e.r) == 0){
            group.1.e.r <- 0
        }
        if (length(group.1.j.r) == 0) group.1.j.r <- 0
        if (length(group.2.j.r) == 0) group.2.j.r <- 0
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
            total.reads <- sum(group.1.p.r,group.1.j.r,group.2.p.r,group.2.j.r)
            group.1.nor.num <- sum(normal.values[c("pairwojun.in","jun.in")])
            group.2.nor.num <- sum(normal.values[c("pairwojun.sk","jun.sk")])
            group.1.read.count <- sum(group.1.p.r,group.1.j.r)
            group.2.read.count <- sum(group.2.p.r,group.2.j.r)
            read.num.test.1 <- sum(group.1.p.r) == 0 & sum(group.1.j.r) == 0
            read.num.test.2 <- sum(group.2.p.r) == 0 & sum(group.2.j.r) == 0
            if ((read.num.test.1) | (read.num.test.2) | total.reads < min.read)    return ("NA")
            if (sum(group.1.j.r) == 0 | sum(group.2.j.r) == 0){
                group.1.nor.num <- normal.values["pair.in"]
                group.2.nor.num <- normal.values["pair.sk"]
                group.1.read.count <- sum(group.1.p.r)
                group.2.read.count <- sum(group.2.p.r)
            }
            else if (sum(group.1.p.r) == 0 | sum(group.2.p.r) == 0){
                group.1.nor.num <- normal.values["jun.in"]
                group.2.nor.num <- normal.values["jun.sk"]
                group.1.read.count <- sum(group.1.j.r)
                group.2.read.count <- sum(group.2.j.r)
            }
        }
        else if (AStype == "IR"){
            group.1.nor.num <- normal.values["jun.in"]
            group.2.nor.num <- normal.values["pairwojun.sk"]
            group.1.read.count <- sum(group.1.j.r)
            group.2.read.count <- sum(group.1.e.r)
        }
        else if (AStype == "ASS"){
            group.1.nor.num <- normal.values["jun.in"]
            group.2.nor.num <- normal.values["jun.sk"]
            group.1.read.count <- sum(group.1.j.r)
            group.2.read.count <- sum(group.2.j.r)
            if(group.1.e.r > 5){
                group.1.nor.num <- sum(normal.values[c("pairwojun.in","jun.in")])
                group.1.read.count <- sum(group.1.e.r,group.1.j.r)
            }
        }
        group.total <- sum(group.1.read.count/group.1.nor.num,group.2.read.count/group.2.nor.num)
        group.1.2.ratio <- (group.2.read.count/group.2.nor.num)/group.total
        if (group.1.read.count == 0 | group.2.read.count == 0)    return ("NA")
        if (group.1.read.count < min.read & group.2.read.count < min.read)    group.1.2.ratio <- "NA"
        return (group.1.2.ratio)
    }
    Each.Cal.ratio <- function(bamfiles=NULL,splicingInfo=NULL){
        if (length(splicingInfo) == 0 | length(splicingInfo) == 0) return (NULL)
        final.ES.result <- NULL
        final.ASS.result <- NULL
        final.IR.result <- NULL
        typenames.t <- unlist(lapply(splicingInfo,function(x) length(x[x[,1] != "NA",1])))
        typenames <- names(typenames.t[which(typenames.t > 0)])
        if (is.element("ES",typenames)){
            final.ES.fi.result <- NULL
            final.ES.se.result <- NULL
            final.MXE.result <- NULL
            each.splice.result <- splicingInfo[["ES"]]
            ES.fi.test <- each.splice.result[,"2ndEX"] == "NA" & each.splice.result[,"Types"] == "ES"
            ES.se.test <- each.splice.result[,"2ndEX"] != "NA" & each.splice.result[,"Types"] == "ES"
            MXE.test <- each.splice.result[,"2ndEX"] != "NA" & each.splice.result[,"Types"] == "MXE"
            ES.fi.result <- rbind(each.splice.result[ES.fi.test,])
            ES.se.result <- rbind(each.splice.result[ES.se.test,])
            MXE.result <- rbind(each.splice.result[MXE.test,])
            do.st <- do.call(rbind,strsplit(each.splice.result[,"DownEX"],"-"))[,1]
            up.en <- do.call(rbind,strsplit(each.splice.result[,"UpEX"],"-"))[,2]
            position.range <- rbind(unique(cbind(do.st,up.en)))
            colnames(position.range) <- c("start","end")
            if (nrow(ES.fi.result) != 0){
                ES.fi.ratio <- foreach(ES.num=1:nrow(ES.fi.result),.packages=called.packages,.combine=rbind) %dopar% {
                    each.ES.re <- rbind(ES.fi.result[ES.num,])
                    es.dw.st <- do.call(rbind,strsplit(each.ES.re[,"DownEX"],"-"))[,1]
                    es.up.en <- do.call(rbind,strsplit(each.ES.re[,"UpEX"],"-"))[,2]
                    each.ranges <- rbind(unique(cbind(es.dw.st,es.up.en)))
                    pre.bam.re <- lapply(1:length(bamfiles),function(each.bam){
                        fi.pair <- paste(c(each.ES.re[,"DownEX"],each.ES.re[,"1stEX"]),collapse="~")
                        se.pair <- paste(c(each.ES.re[,"1stEX"],each.ES.re[,"UpEX"]),collapse="~")
                        tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.ES.re[,"Nchr"]),each.ranges,0)
                        group.1.pair <- rbind(fi.pair,se.pair)
                        group.2.pair <- paste(c(each.ES.re[,"DownEX"],each.ES.re[,"UpEX"]),collapse="~")
                        fi.sp.1 <- split.splice(each.ES.re[,"Do_des"],each.ES.re[,"1st_des"])
                        fi.sp.2 <- split.splice(each.ES.re[,"1st_des"],each.ES.re[,"Up_des"])
                        group.1.spl <- c(fi.sp.1,fi.sp.2)
                        group.2.spl <- split.splice(each.ES.re[,"Do_des"],each.ES.re[,"Up_des"])
                        total.reads <- SplicingReads(tested.reads,each.ES.re[,c("DownEX","1stEX","UpEX")],c(group.1.spl,group.2.spl),readsInfo,inserSize)
                        group.exon <- c(each.ES.re[,"DownEX"],each.ES.re[,"1stEX"],each.ES.re[,"UpEX"])
                        group.1.list <- list(group.1.pair,group.exon,group.1.spl)
                        group.2.list <- list(group.2.pair,group.exon,group.2.spl)
                        names(group.1.list) <- c("paired","exon","junction")
                        names(group.2.list) <- c("paired","exon","junction")
                        total.group.list <- list(group.1.list,group.2.list)
                        names(total.group.list) <- c("Inclu","Skip")
                        coorEX(total.reads,total.group.list,readLen,inserSize,minr,"ES")
                    })
                    pre.bam.re <- do.call(cbind,pre.bam.re)
                    pre.bam.re[is.na(pre.bam.re)] <- "NA"
                    pre.bam.re
                }
                ES.fi.result <- cbind(rbind(ES.fi.result[,c("Index","EnsID","Nchr","1stEX",
                    "2ndEX","DownEX","UpEX","Types")]),ES.fi.ratio)
                colnames(ES.fi.result) <- c("Index","EnsID","Nchr","1stEX","2ndEX",
                    "DownEX","UpEX","Types",sample.names)
                final.ES.fi.result <- ES.fi.result
            }
            if (nrow(ES.se.result) != 0){
                ES.se.ratio <- foreach(ES.num=1:nrow(ES.se.result),.packages=called.packages,.combine=rbind) %dopar% {
                    each.ES.re <- rbind(ES.se.result[ES.num,])
                    do.st <- do.call(rbind,strsplit(each.ES.re[,"DownEX"],"-"))[,1]
                    up.en <- do.call(rbind,strsplit(each.ES.re[,"UpEX"],"-"))[,2]
                    each.ranges <- rbind(unique(cbind(do.st,up.en)))
                    pre.bam.re <- lapply(1:length(bamfiles),function(each.bam){
                        tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.ES.re[,"Nchr"]),each.ranges,0)
                        group.1.pair <- rbind(paste(c(each.ES.re[,"DownEX"],each.ES.re[,"1stEX"]),collapse="~"),
                            paste(c(each.ES.re[,"1stEX"],each.ES.re[,"2ndEX"]),collapse="~"),
                            paste(c(each.ES.re[,"2ndEX"],each.ES.re[,"UpEX"]),collapse="~"))
                        group.2.pair <- paste(c(each.ES.re[,"DownEX"],each.ES.re[,"UpEX"]),collapse="~")
                        group.1.spl <- c(split.splice(each.ES.re[,"Do_des"],each.ES.re[,"1st_des"]),
                            split.splice(each.ES.re[,"1st_des"],each.ES.re[,"2nd_des"]),
                            split.splice(each.ES.re[,"2nd_des"],each.ES.re[,"Up_des"]))
                        group.2.spl <- split.splice(each.ES.re[,"Do_des"],each.ES.re[,"Up_des"])
                        total.reads <- SplicingReads(tested.reads,each.ES.re[,c("DownEX","1stEX","2ndEX","UpEX")],c(group.1.spl,group.2.spl),readsInfo,inserSize)
                        group.exon <- c(each.ES.re[,"DownEX"],each.ES.re[,"1stEX"],
                            each.ES.re[,"2ndEX"],each.ES.re[,"UpEX"])
                        group.1.list <- list(group.1.pair,group.exon,group.1.spl)
                        group.2.list <- list(group.2.pair,group.exon,group.2.spl)
                        names(group.1.list) <- c("paired","exon","junction")
                        names(group.2.list) <- c("paired","exon","junction")
                        total.group.list <- list(group.1.list,group.2.list)
                        names(total.group.list) <- c("Inclu","Skip")
                        coorEX(total.reads,total.group.list,readLen,inserSize,minr,"ESse")
                    })
                    pre.bam.re <- do.call(cbind,pre.bam.re)
                    pre.bam.re[is.na(pre.bam.re)] <- "NA"
                    pre.bam.re
                }
                ES.se.result <- cbind(rbind(ES.se.result[,c("Index","EnsID","Nchr","1stEX",
                    "2ndEX","DownEX","UpEX","Types")]),ES.se.ratio)
                colnames(ES.se.result) <- c("Index","EnsID","Nchr","1stEX","2ndEX",
                    "DownEX","UpEX","Types",sample.names)
                final.ES.se.result <- ES.se.result
            }
            if (nrow(MXE.result) != 0){
                MXE.ratio <- foreach(ES.num=1:nrow(MXE.result),.packages=called.packages,.combine=rbind) %dopar% {
                    each.ES.re <- rbind(MXE.result[ES.num,])
                    do.st <- do.call(rbind,strsplit(each.ES.re[,"DownEX"],"-"))[,1]
                    up.en <- do.call(rbind,strsplit(each.ES.re[,"UpEX"],"-"))[,2]
                    each.ranges <- rbind(unique(cbind(do.st,up.en)))
                    pre.bam.re <- lapply(1:length(bamfiles),function(each.bam){
                        tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.ES.re[,"Nchr"]),each.ranges,0)
                        fi.pair.1 <- paste(c(each.ES.re[,"DownEX"],each.ES.re[,"1stEX"]),collapse="~")
                        se.pair.1 <- paste(c(each.ES.re[,"1stEX"],each.ES.re[,"UpEX"]),collapse="~")
                        fi.pair.2 <- paste(c(each.ES.re[,"DownEX"],each.ES.re[,"2ndEX"]),collapse="~")
                        se.pair.2 <- paste(c(each.ES.re[,"2ndEX"],each.ES.re[,"UpEX"]),collapse="~")
                        group.1.pair <- rbind(fi.pair,se.pair)
                        group.2.pair <- rbind(fi.pair.2,se.pair.2)
                        fi.sp.1 <- split.splice(each.ES.re[,"Do_des"],each.ES.re[,"1st_des"])
                        fi.sp.2 <- split.splice(each.ES.re[,"1st_des"],each.ES.re[,"Up_des"])
                        se.sp.1 <- split.splice(each.ES.re[,"Do_des"],each.ES.re[,"2nd_des"])
                        se.sp.2 <- split.splice(each.ES.re[,"2nd_des"],each.ES.re[,"Up_des"])
                        group.1.spl <- c(fi.sp.1,fi.sp.2)
                        group.2.spl <- c(se.sp.1,se.sp.2)
                        total.reads <- SplicingReads(tested.reads,each.ES.re[,c("DownEX","1stEX","2ndEX","UpEX")],c(group.1.spl,group.2.spl),readsInfo,inserSize)
                        group.1.exon <- c(each.ES.re[,"DownEX"],each.ES.re[,"1stEX"],each.ES.re[,"UpEX"])
                        group.2.exon <- c(each.ES.re[,"DownEX"],each.ES.re[,"2ndEX"],each.ES.re[,"UpEX"])
                        group.1.list <- list(group.1.pair,group.1.exon,group.1.spl)
                        group.2.list <- list(group.2.pair,group.2.exon,group.2.spl)
                        names(group.1.list) <- c("paired","exon","junction")
                        names(group.2.list) <- c("paired","exon","junction")
                        total.group.list <- list(group.2.list,group.1.list)
                        names(total.group.list) <- c("Inclu","Skip")
                        coorEX(total.reads,total.group.list,readLen,inserSize,minr,"MXE")
                    })
                    pre.bam.re <- do.call(cbind,pre.bam.re)
                    pre.bam.re[is.na(pre.bam.re)] <- "NA"
                    pre.bam.re
                }
                MXE.result <- cbind(rbind(MXE.result[,c("Index","EnsID","Nchr","1stEX",
                    "2ndEX","DownEX","UpEX","Types")]),MXE.ratio)
                colnames(MXE.result) <- c("Index","EnsID","Nchr","1stEX","2ndEX",
                    "DownEX","UpEX","Types",sample.names)
                final.MXE.result <- MXE.result
            }
            final.ES.result <- rbind(final.ES.fi.result,final.ES.se.result,final.MXE.result)
            if (length(final.ES.result) == 0)    final.ES.result <- NULL
            else    rownames(final.ES.result) <- 1:nrow(final.ES.result)
        }
        if (is.element("IR",typenames)){
            IR.result <- splicingInfo[["IR"]]
            each.splice.result <- IR.result
            do.st <- do.call(rbind,strsplit(each.splice.result[,"DownEX"],"-"))[,1]
            up.en <- do.call(rbind,strsplit(each.splice.result[,"UpEX"],"-"))[,2]
            position.range <- rbind(unique(cbind(do.st,up.en)))
            colnames(position.range) <- c("start","end")
            IR.ratio <- foreach(IR.num=1:nrow(each.splice.result),.packages=called.packages,.combine=rbind) %dopar% {
                each.IR.re <- rbind(each.splice.result[IR.num,])
                IR.do.st <- do.call(rbind,strsplit(each.IR.re[,"DownEX"],"-"))[,1]
                IR.up.en <- do.call(rbind,strsplit(each.IR.re[,"UpEX"],"-"))[,2]
                each.ranges <- rbind(unique(cbind(IR.do.st,IR.up.en)))
                pre.bam.re <- lapply(1:length(bamfiles),function(each.bam){
                    tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.IR.re[,"Nchr"]),each.ranges,0)
                    IR.spli.ex <- paste(unlist(strsplit(each.IR.re[,"DownEX"],"-"))["DownEX2"],
                        unlist(strsplit(each.IR.re[,"UpEX"],"-"))["UpEX1"],sep="-")
                    reads.num.test <- length(grep("N",tested.reads[,"cigar"])) >    5
                    IR.test.do <- unlist(strsplit(each.IR.re[,"DownEX"],"-"))["DownEX2"]
                    IR.test.up <- unlist(strsplit(each.IR.re[,"UpEX"],"-"))["UpEX1"]
                    if ((reads.num.test)| unlist(IR.test.do < IR.test.up)){
                        group.1.pair <- rbind(paste(c(each.IR.re[,"DownEX"],each.IR.re[,"UpEX"]),collapse="~"))
                        fi.pair.2 <- paste(c(each.IR.re[,"DownEX"],IR.spli.ex),collapse="~")
                        se.pair.2 <- paste(c(IR.spli.ex,each.IR.re[,"UpEX"]),collapse="~")
                        group.2.pair <- rbind(fi.pair.2,se.pair.2)
                        group.1.spl <- c(split.splice(each.IR.re[,"Do_des"],each.IR.re[,"Up_des"]))
                        group.2.spl <- "NA"
                        total.reads <- SplicingReads(tested.reads,c(each.IR.re[,"DownEX"],IR.spli.ex,each.IR.re[,"UpEX"]),group.1.spl,"exon",inserSize)
                        group.1.exon <- c(each.IR.re[,"DownEX"],IR.spli.ex,each.IR.re[,"UpEX"])
                        group.2.exon <- c(each.IR.re[,"DownEX"],IR.spli.ex,each.IR.re[,"UpEX"])
                        group.1.list <- list(group.1.pair,group.1.exon,group.1.spl)
                        group.2.list <- list(group.2.pair,group.2.exon,group.2.spl)
                        names(group.1.list) <- c("paired","exon","junction")
                        names(group.2.list) <- c("paired","exon","junction")
                        total.group.list <- list(group.1.list,group.2.list)
                        names(total.group.list) <- c("Inclu","Skip")
                        coorEX(total.reads,total.group.list,readLen,inserSize,minr,"IR")
                    }
                    else    "NA"
                })
                pre.bam.re <- do.call(cbind,pre.bam.re)
                pre.bam.re[is.na(pre.bam.re)] <- "NA"
                pre.bam.re
            }
            IR.result <- cbind(rbind(IR.result[,c("Index","EnsID","Nchr","RetainEX",
                "DownEX","UpEX","Types")]),IR.ratio)
            colnames(IR.result) <- c("Index","EnsID","Nchr","RetainEX",
                "DownEX","UpEX","Types",sample.names)
            final.IR.result <- IR.result
            if (length(final.IR.result) == 0)    final.IR.result <- NULL
            else    rownames(final.IR.result) <- 1:nrow(final.IR.result)
        }
        if (is.element("ASS",typenames)){
            final.A5SS.result <- NULL
            final.A3SS.result <- NULL
            ASS.result <- splicingInfo[["ASS"]]
            A3SS.result <- rbind(ASS.result[ASS.result[,"Types"] == "A3SS",])
            A5SS.result <- rbind(ASS.result[ASS.result[,"Types"] == "A5SS",])
            if (nrow(A5SS.result) != 0){
                position.range <- lapply(1:nrow(A5SS.result),function(ASS.nums){
                    short.ex.st <- unlist(strsplit(A5SS.result[ASS.nums,"ShortEX"],"-"))
                    long.ex.st <- unlist(strsplit(A5SS.result[ASS.nums,"LongEX"],"-"))
                    nei.ex.st <- unlist(strsplit(A5SS.result[ASS.nums,grep("NeighborEX",colnames(A5SS.result))],"-"))
                    min.r <- min(c(short.ex.st,long.ex.st))
                    max.r <- max(c(nei.ex.st))
                    cbind(min.r,max.r)
                })
                position.range <- unique(do.call(rbind,position.range))
                colnames(position.range) <- c("start","end")
                A5SS.ratio <- foreach(A5SS.num=1:nrow(A5SS.result),.packages=called.packages,.combine=rbind) %dopar% {
                    each.A5SS.re <- rbind(A5SS.result[A5SS.num,])
                    A5SS.nei.ex.st <- unlist(strsplit(each.A5SS.re[,grep("NeighborEX",colnames(each.A5SS.re))],"-"))
                    min.r <- max(as.integer(unlist(strsplit(each.A5SS.re[,"ShortEX"],"-")))) - 20
                    max.r <- min(as.integer(c(A5SS.nei.ex.st))) + 20
                    each.ranges <- cbind(min.r,max.r)
                    colnames(each.ranges) <- c("start","end")
                    pre.bam.re <- lapply(1:length(bamfiles),function(each.bam){
                        tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.A5SS.re[,"Nchr"]),each.ranges,0)
                        if (length(grep("N",tested.reads[,"cigar"])) >    5){
                            ASS5.spli.ex <- paste(unlist(strsplit(each.A5SS.re[,"ShortEX"],"-"))["ShortEX2"],unlist(strsplit(each.A5SS.re[,"LongEX"],"-"))["LongEX2"],sep="-")
                            names(ASS5.spli.ex) <- "Alt.ex"
                            A5SS.cn <- colnames(each.A5SS.re)
                            A5SS.long.nei <- each.A5SS.re[,is.element(A5SS.cn,c("LongNeighborEX","NeighborEX"))]
                            A5SS.short.nei <- each.A5SS.re[,is.element(A5SS.cn,c("ShortNeighborEX","NeighborEX"))]
                            A5SS.long.nei.des <- each.A5SS.re[,is.element(A5SS.cn,c("LongNeighbor_des","Neighbor_des"))]
                            A5SS.short.nei.des <- each.A5SS.re[,is.element(A5SS.cn,c("ShortNeighbor_des","Neighbor_des"))]
                            group.2.pair <- rbind(paste(c(each.A5SS.re[,"ShortEX"],A5SS.short.nei),collapse="~"))
                            group.1.spl <- c(split.splice(each.A5SS.re[,"Long_des"],A5SS.long.nei.des))
                            group.2.spl <- c(split.splice(each.A5SS.re[,"Short_des"],A5SS.short.nei.des))
                            total.reads <- SplicingReads(tested.reads,c(ASS5.spli.ex,A5SS.long.nei),c(group.1.spl,group.2.spl),"exon",inserSize)
                            group.1.exon <- c(each.A5SS.re[,"LongEX"],ASS5.spli.ex,A5SS.short.nei)
                            group.2.exon <- c(each.A5SS.re[,"ShortEX"],ASS5.spli.ex,A5SS.long.nei)
                            group.1.list <- list(group.1.pair,group.1.exon,group.1.spl)
                            group.2.list <- list(group.2.pair,group.2.exon,group.2.spl)
                            names(group.1.list) <- c("paired","exon","junction")
                            names(group.2.list) <- c("paired","exon","junction")
                            total.group.list <- list(group.1.list,group.2.list)
                            names(total.group.list) <- c("Inclu","Skip")
                            coorEX(total.reads,total.group.list,readLen,inserSize,minr,"ASS")
                        }
                        else "NA"
                    })
                    pre.bam.re <- do.call(cbind,pre.bam.re)
                    pre.bam.re[is.na(pre.bam.re)] <- "NA"
                    pre.bam.re
                }
                cn <- colnames(A5SS.result)[grep("Index|EnsID|Nchr|ShortEX|LongEX|NeighborEX|Types",colnames(A3SS.result))]
                A5SS.result <- cbind(rbind(A5SS.result[,cn]),A5SS.ratio)
                colnames(A5SS.result) <- c(cn,sample.names)
                final.A5SS.result <- A5SS.result
            }
            if (nrow(A3SS.result) != 0){
                position.range <- lapply(1:nrow(A3SS.result),function(ASS.nums){
                    short.ex.st <- unlist(strsplit(A3SS.result[ASS.nums,"ShortEX"],"-"))
                    long.ex.st <- unlist(strsplit(A3SS.result[ASS.nums,"LongEX"],"-"))
                    nei.ex.st <- unlist(strsplit(A3SS.result[ASS.nums,grep("NeighborEX",colnames(A3SS.result))],"-"))
                    max.r <- max(c(short.ex.st,long.ex.st))
                    min.r <- min(c(nei.ex.st))
                    cbind(min.r,max.r)
                })
                position.range <- unique(do.call(rbind,position.range))
                colnames(position.range) <- c("start","end")
                A3SS.ratio <-    foreach(A3SS.num=1:nrow(A3SS.result),.packages=called.packages,.combine=rbind) %dopar% {
                    each.A3SS.re <- rbind(A3SS.result[A3SS.num,])
                    A3SS.nei.ex.st <- unlist(strsplit(each.A3SS.re[,grep("NeighborEX",colnames(each.A3SS.re))],"-"))
                    max.r <- min(as.integer(unlist(strsplit(each.A3SS.re[,"ShortEX"],"-")))) + 20
                    min.r <- max(as.integer(c(A3SS.nei.ex.st))) - 20
                    each.ranges <- cbind(min.r,max.r)
                    colnames(each.ranges) <- c("start","end")
                    pre.bam.re <- lapply(1:length(bamfiles),function(each.bam){
                        tested.reads <- ReadBamToPosition(bamfiles[each.bam],unique(each.A3SS.re[,"Nchr"]),each.ranges,0)
                        if (length(grep("N",tested.reads[,"cigar"])) >    5){
                            ASS3.spli.ex <- paste(unlist(strsplit(each.A3SS.re[,"LongEX"],"-"))["LongEX1"],unlist(strsplit(each.A3SS.re[,"ShortEX"],"-"))["ShortEX1"],sep="-")
                            names(ASS3.spli.ex) <- "Alt.ex"
                            A3SS.cn <- colnames(each.A5SS.re)
                            A3SS.long.nei <- each.A3SS.re[,is.element(A3SS.cn,c("LongNeighborEX","NeighborEX"))]
                            A3SS.short.nei <- each.A3SS.re[,is.element(A3SS.cn,c("ShortNeighborEX","NeighborEX"))]
                            A3SS.long.nei.des <- each.A3SS.re[,is.element(A3SS.cn,c("LongNeighbor_des","Neighbor_des"))]
                            A3SS.short.nei.des <- each.A3SS.re[,is.element(A3SS.cn,c("ShortNeighbor_des","Neighbor_des"))]
                            group.2.pair <- rbind(paste(c(A3SS.short.nei,each.A3SS.re[,"ShortEX"]),collapse="~"))
                            group.1.spl <- c(split.splice(A3SS.long.nei.des,each.A3SS.re[,"Long_des"]))
                            group.2.spl <- c(split.splice(A3SS.short.nei.des,each.A3SS.re[,"Short_des"]))
                            total.reads <- SplicingReads(tested.reads,c(A3SS.cn,ASS3.spli.ex),c(group.1.spl,group.2.spl),"exon",inserSize)
                            group.1.exon <- c(A3SS.short.nei,each.A3SS.re[,"LongEX"],ASS3.spli.ex)
                            group.2.exon <- c(A3SS.long.nei,each.A3SS.re[,"ShortEX"],ASS3.spli.ex)
                            group.1.exon <- c(A3SS.short.nei,ASS3.spli.ex)
                            group.2.exon <- c(A3SS.long.nei,ASS3.spli.ex)
                            group.1.list <- list(group.1.pair,group.1.exon,group.1.spl)
                            group.2.list <- list(group.2.pair,group.2.exon,group.2.spl)
                            names(group.1.list) <- c("paired","exon","junction")
                            names(group.2.list) <- c("paired","exon","junction")
                            total.group.list <- list(group.1.list,group.2.list)
                            names(total.group.list) <- c("Inclu","Skip")
                            coorEX(total.reads,total.group.list,readLen,inserSize,minr,"ASS")
                        }
                        else "NA"
                    })
                    pre.bam.re <- do.call(cbind,pre.bam.re)
                    pre.bam.re[is.na(pre.bam.re)] <- "NA"
                    pre.bam.re
                }
                cn <- colnames(A3SS.result)[grep("Index|EnsID|Nchr|ShortEX|LongEX|NeighborEX|Types",colnames(A3SS.result))]
                A3SS.result <- cbind(rbind(A3SS.result[,cn]),A3SS.ratio)
                colnames(A3SS.result) <- c(cn,sample.names)
                final.A3SS.result <- A3SS.result
            }
            final.ASS.result <- rbind(final.A5SS.result,final.A3SS.result)
            if (length(final.ASS.result) == 0)    final.ASS.result <- NULL
            else    rownames(final.ASS.result) <- 1:nrow(final.ASS.result)
        }
        final.result <- list(final.ES.result,final.ASS.result,final.IR.result)
        names(final.result) <- c("ES","ASS","IR")
        return (final.result)
    }
    each.ES.re <- NULL
    ES.num <- NULL
    IR.num <- NULL
    A5SS.num <- NULL
    A3SS.num <- NULL
    group.1.pair <- NULL
    group.2.pair <- NULL
    fi.pair <- NULL
    se.pair <- NULL
    registerDoParallel(cores=Ncor)
    called.packages <- c("GenomicRanges","GenomicFeatures")
    sample.files <- rbind(Total.bamfiles[,"path"])
    sample.names <- rbind(Total.bamfiles[,"names"])
    Total.splicingInfo <- ASdb@"SplicingModel"
    if (length(CalIndex) != 0){
        if (ncol(Total.splicingInfo$"ES") != 1){
            each.result <- Total.splicingInfo$"ES"[is.element(Total.splicingInfo$"ES"[,"Index"],CalIndex),]
            Total.splicingInfo$"ES" <- rbind(each.result)
        }
        if (ncol(Total.splicingInfo$"ASS") != 1){
            each.result <- Total.splicingInfo$"ASS"[is.element(Total.splicingInfo$"ASS"[,"Index"],CalIndex),]
            Total.splicingInfo$"ASS" <- rbind(each.result)
        }
        if (ncol(Total.splicingInfo$"IR") != 1){
            each.result <- Total.splicingInfo$"IR"[is.element(Total.splicingInfo$"IR"[,"Index"],CalIndex),]
            Total.splicingInfo$"IR" <- rbind(each.result)
        }
    }
    final.total.ratio <- Each.Cal.ratio(sample.files,Total.splicingInfo)
    if (length(final.total.ratio$"ES") == 0)    final.total.ratio$"ES" <- as.matrix("NA")
    if (length(final.total.ratio$"ASS") == 0)    final.total.ratio$"ASS" <- as.matrix("NA")
    if (length(final.total.ratio$"IR") == 0)    final.total.ratio$"IR" <- as.matrix("NA")
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=final.total.ratio,
        GroupDiff=ASdb@"GroupDiff",sQTLs=ASdb@"sQTLs",Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
    if (length(out.dir) != 0){
        system(paste("mkdir -p ",out.dir,"/AS_Ratio",sep=""))
        write.table(final.total.ratio[["ES"]],paste(out.dir,"/AS_Ratio/ES_Ratio.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.total.ratio[["ASS"]],paste(out.dir,"/AS_Ratio/ASS_Ratio.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.total.ratio[["IR"]],paste(out.dir,"/AS_Ratio/IR_Ratio.txt",sep=""),sep='\t',quote=FALSE)
    }
    return(ASdb)
}

