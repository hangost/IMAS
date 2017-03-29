ExonsCluster <- function(ASdb,GTFdb,Ncor=1){
    sortEX <- function(test.range,outtype="mat"){
        over.ranges <- findOverlaps(test.range,test.range,select="all")
        over.ranges <- unique(as.matrix(over.ranges))
        of <- NULL
        sort.firstEX <- tapply(over.ranges[,2],over.ranges[,1],function(of){
            paste(of,collapse=",")
        })
        u.sort.firstEX <- unique(sort.firstEX)
        pre.result <- NULL
        final.sorted.EX <- lapply(u.sort.firstEX,function(u.fe){
            over.num <- as.double(unlist(strsplit(u.fe,",")))
            if (outtype == "ranges"){
                pre.result <- test.range[over.num]
                pre.result
            }
            else {
                start.ran <- as.matrix(start(test.range[over.num]))
                end.ran <- as.matrix(end(test.range[over.num]))
                pre.result <- cbind(start.ran,end.ran)
                colnames(pre.result) <- c("start","end")
                pre.result
            }
        })
        return(final.sorted.EX)
    }
    p.st <- function(pasted.string = NULL,pas.string = ","){
        return (paste(sort(unique(pasted.string)),collapse=pas.string))
    }
    merge.mat <- function(mat1,mat2,mat3,AtTypes,mat4=NULL){
        merge.mat <- NULL
        if (AtTypes == "A5SS" | AtTypes == "A3SS"){
            p.mat1 <- paste(mat1[,"start"],mat1[,"end"],sep="-")
            p.mat2 <- paste(mat2[,"start"],mat2[,"end"],sep="-")
            p.mat3 <- paste(mat3[,"start"],mat3[,"end"],sep="-")
            dse.mat1 <- unlist(strsplit(mat1[,"des"],","))
            dse.mat2 <- unlist(strsplit(mat2[,"des"],","))
            dse.mat3 <- unlist(strsplit(mat3[,"des"],","))
            dse.mat1 <- do.call(rbind,strsplit(dse.mat1,"-"))
            dse.mat2 <- do.call(rbind,strsplit(dse.mat2,"-"))
            dse.mat3 <- do.call(rbind,strsplit(dse.mat3,"-"))
            colnames(dse.mat1) <- c("start","end")
            colnames(dse.mat2) <- c("start","end")
            colnames(dse.mat3) <- c("start","end")
            if (AtTypes == "A5SS") mat.stan <- c("end","end","start")
            if (AtTypes == "A3SS") mat.stan <- c("start","start","end")
            u.dse.mat1 <- unique(dse.mat1[,mat.stan[1]])
            u.dse.mat2 <- unique(dse.mat2[,mat.stan[2]])
            u.dse.mat3 <- sort(as.double(unique(dse.mat3[,mat.stan[3]])))
            total.des.mat <- lapply(u.dse.mat1,function(e.d.mat1){
                only.des2 <- u.dse.mat2[u.dse.mat2!=e.d.mat1]
                only.des2 <- paste(sort(only.des2),collapse=",")
                only.des3 <- paste(sort(u.dse.mat3),collapse=",")
                p.u.dse.mat3 <- paste(u.dse.mat3,collapse=",")
                if (AtTypes == "A5SS"){
                    p.des.2.3 <- paste(only.des2,only.des3,sep="-")
                    merge.des.mat <- paste(e.d.mat1,p.u.dse.mat3,sep="-")
                    merge.des.mat <- c(merge.des.mat,p.des.2.3)
                }
                if (AtTypes == "A3SS"){
                    p.des.3.2 <- paste(only.des3,only.des2,sep="-")
                    merge.des.mat <- paste(p.u.dse.mat3,e.d.mat1,sep="-")
                    merge.des.mat <- c(merge.des.mat,p.des.3.2)
                }
                cn <- c("start","end","des")
                ex1 <- rbind(mat1[mat1[,mat.stan[1]]==e.d.mat1,cn])
                ex2 <- rbind(mat2[mat2[,mat.stan[2]]!=e.d.mat1,cn])
                ex3 <- rbind(mat3[,c("start","end")])
                mat1.ex.des <- sort(ex1[,"des"])
                mat2.ex.des <- sort(ex2[,"des"])
                mat3.ex.des <- sort(paste(ex3[,"start"],ex3[,"end"],sep="-"))
                mat1.p <- paste(min(ex1[,"start"]),max(ex1[,"end"]),sep="-")
                mat2.p <- paste(min(ex2[,"start"]),max(ex2[,"end"]),sep="-")
                mat3.p <- paste(min(ex3[,"start"]),max(ex3[,"end"]),sep="-")
                only.des1 <- paste(mat1.ex.des,collapse=",")
                only.des2 <- paste(mat2.ex.des,collapse=",")
                only.des3 <- paste(mat3.ex.des,collapse=",")
                merge.mat <- c(mat1.p,mat2.p,mat3.p,
                    only.des1,only.des2,only.des3,merge.des.mat)
            })
            final.mat <- do.call(rbind,total.des.mat)
        }
        if (AtTypes == "ES" | AtTypes == "MXE"){
            p.mat1 <- paste(mat1[,"start"],mat1[,"end"],sep="-")
            p.mat2 <- paste(mat2[,"start"],mat2[,"end"],sep="-")
            p.mat3 <- paste(mat3[,"start"],mat3[,"end"],sep="-")
            dse.mat1 <- unlist(strsplit(mat1[,"des"],","))
            dse.mat2 <- unlist(strsplit(mat2[,"des"],","))
            dse.mat3 <- unlist(strsplit(mat3[,"des"],","))
            dse.mat1 <- do.call(rbind,strsplit(dse.mat1,"-"))
            dse.mat2 <- do.call(rbind,strsplit(dse.mat2,"-"))
            dse.mat3 <- do.call(rbind,strsplit(dse.mat3,"-"))
            colnames(dse.mat1) <- c("start","end")
            colnames(dse.mat2) <- c("start","end")
            colnames(dse.mat3) <- c("start","end")
            p.dse1.s <- p.st(dse.mat1[,"start"])
            p.ds1.e <- p.st(dse.mat1[,"end"])
            p.dse2.s <- p.st(dse.mat2[,"start"])
            p.ds2.e <- p.st(dse.mat2[,"end"])
            p.dse3.s <- p.st(dse.mat3[,"start"])
            p.ds3.e <- p.st(dse.mat3[,"end"])
            if (!any(seq_along(mat4))){
                final.des.mat <- paste(p.ds2.e,p.dse1.s,sep="-")
                final.des.mat <- c(final.des.mat,"NA")
                final.mat <- c(final.des.mat,paste(p.ds1.e,p.dse3.s,sep="-"))
                final.mat <- c(p.mat1,"NA",p.mat2,p.mat3,mat1[,"des"],
                    "NA",mat2[,"des"],mat3[,"des"],final.mat,"ES")
            }
            if(any(seq_along(mat4))){
                p.mat4 <- paste(mat4[,"start"],mat4[,"end"],sep="-")
                dse.mat4 <- unlist(strsplit(mat4[,"des"],","))
                dse.mat4 <- do.call(rbind,strsplit(dse.mat4,"-"))
                colnames(dse.mat4) <- c("start","end")
                p.dse4.s <- p.st(dse.mat4[,"start"])
                p.dse4.e <- p.st(dse.mat4[,"end"])
                if(AtTypes == "ES"){
                    f.mat <- paste(p.ds3.e,p.dse1.s,sep="-")
                    f.mat <- c(f.mat,paste(p.ds1.e,p.dse2.s,sep="-"))
                    final.des.mat <- c(f.mat,paste(p.ds2.e,p.dse4.s,sep="-"))
                }
                if(AtTypes == "MXE"){
                    fir.spli <- paste(p.ds3.e,p.dse1.s,sep="-")
                    sec.spli <- paste(p.ds1.e,p.dse4.s,sep="-")
                    thi.spli <- paste(p.ds3.e,p.dse2.s,sep="-")
                    for.spli <- paste(p.ds2.e,p.dse4.s,sep="-")
                    p.thi.for <- paste(thi.spli,"|",for.spli,sep="")
                    f.mat <- paste(fir.spli,"|",sec.spli,sep="")
                    f.mat <- c(f.mat,"NA")
                    final.des.mat <- c(f.mat,p.thi.for)
                }
                p.mat <- c(p.mat1,p.mat2,p.mat3,p.mat4)
                des.mat <- c(mat1[,"des"],mat2[,"des"],
                    mat3[,"des"],mat4[,"des"])
                final.mat <- c(p.mat,des.mat,final.des.mat,AtTypes)
            }
        }
        return (final.mat)
    }
    ASS.Alt.result <- function(altSplice){
        mer.Ex <- function(Ex.mat,std.lo,Alt.type){
            test.re <- NULL
            merged.mat <- tapply(Ex.mat,std.lo,function(total.sem){
                sem <- do.call(rbind,strsplit(total.sem,"-"))
                colnames(sem) <- c("start","end")
                merge.des <- paste(sort(total.sem),collapse=",")
                order.lo <- order(as.double(sem[,"start"]))
                test.re <- cbind(rbind(sem[order.lo[1],]),rbind(merge.des))
                test.re
            })
            return (do.call(rbind,merged.mat))
        }
        ASS.merge.f <- function(ASS.result,Alt.type){
            neig.nm <- grep("Neighbor_des",colnames(ASS.result))
            neig.nm <- colnames(ASS.result)[neig.nm]
            s.Short.EX <- strsplit(ASS.result[,"ShortEX"],"-")
            s.long.EX <- strsplit(ASS.result[,"LongEX"],"-")
            Short.EX <- unique(do.call(rbind,s.Short.EX))
            long.EX <- unique(do.call(rbind,s.long.EX))
            colnames(Short.EX) <- c("start","end")
            colnames(long.EX) <- c("start","end")
            total.EX <- rbind(Short.EX,long.EX)
            t.st <- as.integer(total.EX[,"start"])
            t.en <- as.integer(total.EX[,"end"])
            each.ranges <- IRanges(start=t.st,end=t.en)
            total.EX.ranges <- GRanges(seqnames="*",ranges=each.ranges)
            Short.EX <- unlist(strsplit(ASS.result[,"Short_des"],","))
            Short.EX <- strsplit(Short.EX,"-")
            Short.EX <- unique(do.call(rbind,Short.EX))
            long.EX <- unlist(strsplit(ASS.result[,"Long_des"],","))
            long.EX <- strsplit(long.EX,"-")
            long.EX <- unique(do.call(rbind,long.EX))
            neighbor.EX <- unlist(strsplit(ASS.result[,neig.nm],","))
            neighbor.EX <- strsplit(neighbor.EX,"-")
            neighbor.EX <- unique(do.call(rbind,neighbor.EX))
            colnames(Short.EX) <- c("start","end")
            colnames(long.EX) <- c("start","end")
            colnames(neighbor.EX) <- c("start","end")
            s.pa.ex.mat <- paste(Short.EX[,"start"],Short.EX[,"end"],sep="-")
            l.pa.ex.mat <- paste(long.EX[,"start"],long.EX[,"end"],sep="-")
            if (Alt.type == "A5SS"){
                A5.num <- altSplice[,"Types"] == "A5SS"
                A5.re <- rbind(altSplice[A5.num,])
                short.m <- mer.Ex(s.pa.ex.mat,Short.EX[,"end"],"A5SS")
                long.m <- mer.Ex(l.pa.ex.mat,long.EX[,"end"],"A5SS")
            }
            else if(Alt.type == "A3SS"){
                A3.num <- altSplice[,"Types"] == "A3SS"
                A3.re <- rbind(altSplice[A3.num,])
                short.m <- mer.Ex(s.pa.ex.mat,Short.EX[,"start"],"A3SS")
                long.m <- mer.Ex(l.pa.ex.mat,long.EX[,"start"],"A3SS")
            }
            p.n.E <- paste(neighbor.EX[,"start"],neighbor.EX[,"end"],sep="-")
            neighbor.ex.merge <- cbind(neighbor.EX,p.n.E)
            colnames(neighbor.ex.merge) <- c("start","end","des")
            colnames(short.m) <- c("start","end","des")
            colnames(long.m) <- c("start","end","des")
            ASS.final.result <- cbind(merge.mat(short.m,
                long.m,neighbor.ex.merge,Alt.type),Alt.type)
            cn <- c("ShortEX","LongEX","NeighborEX","Short_des","Long_des",
                "Neighbor_des","splicing in 1EX","splicing in 2EX","Types")
            colnames(ASS.final.result) <- cn
            return (ASS.final.result)
        }
        ASS.final.result <- NULL
        merged.mat <- NULL
        test.pos.re <- NULL
        A5.num <- altSplice[,"Types"] == "A5SS"
        A3.num <- altSplice[,"Types"] == "A3SS"
        if (is.element("TRUE",A5.num)){
            A5SS.f.result <- ASS.merge.f(rbind(altSplice[A5.num,]),"A5SS")
        }
        if (is.element("TRUE",A3.num)){
            A3SS.f.result <- ASS.merge.f(rbind(altSplice[A3.num,]),"A3SS")
        }
        final.result <- rbind(A5SS.f.result,A3SS.f.result)
        return (final.result)
    }
    
    ES.Alt.result <- function(altSplice){
        ES.merge.f <- function(ES.result,Alt.type){
            firstEX <- strsplit(ES.result[,"1stEX"],"-")
            firstEX <- unique(do.call(rbind,firstEX))
            colnames(firstEX) <- c("start","end")
            fi.s <- as.integer(firstEX[,"start"])
            fi.e <- as.integer(firstEX[,"end"])
            each.ranges <- IRanges(start=fi.s,end=fi.e)
            firstEX.range <- GRanges(seqnames="*",ranges=each.ranges)
            merged.final.result <- NULL
            sorted.first.EX <- sortEX(firstEX.range,"ranges")
            se.test <- any(which(ES.result[,"2ndEX"] != "NA"))
            final.result <- NULL
            cn <- c("1stEX","2ndEX","DownEX","UpEX","1st_des","2nd_des"
                ,"Do_des","Up_des","1stSpl","2ndSpl","3rdSpl","Types")
            m.final.result <- lapply(sorted.first.EX,function(sfe){
                ea.fi.EX <- cbind(start(sfe),end(sfe))
                colnames(ea.fi.EX) <- c("start","end")
                p.fi.ex <- paste(ea.fi.EX[,"start"],ea.fi.EX[,"end"],sep="-")
                over.ex <- is.element(ES.result[,"1stEX"],p.fi.ex)
                sub.ES <- rbind(ES.result[over.ex,])
                DoEX <- unique(sub.ES[,"Do_des"])
                DoEX <- unlist(strsplit(DoEX,","))
                DoEX <- unique(do.call(rbind,strsplit(DoEX,"-")))
                UpEX <- unique(sub.ES[,"Up_des"])
                UpEX <- unlist(strsplit(UpEX,","))
                UpEX <- unique(do.call(rbind,strsplit(UpEX,"-")))
                colnames(DoEX) <- c("start","end")
                colnames(UpEX) <- c("start","end")
                Do.des <- p.st(paste(DoEX[,"start"],DoEX[,"end"],sep="-"))
                up.des <- p.st(paste(UpEX[,"start"],UpEX[,"end"],sep="-"))
                tar.des <- p.st(p.fi.ex)
                merged.Do <- cbind(rbind(c(min(as.double(DoEX[,"start"])),
                    max(as.double(DoEX[,"end"])))),rbind(Do.des))
                merged.Up <- cbind(rbind(c(min(as.double(UpEX[,"start"])),
                    max(as.double(UpEX[,"end"])))),rbind(up.des))
                mer.tar <- cbind(rbind(c(min(as.double(ea.fi.EX[,"start"])),
                    max(as.double(ea.fi.EX[,"end"])))),rbind(tar.des))
                colnames(mer.tar) <- c("start","end","des")
                colnames(merged.Up) <- c("start","end","des")
                colnames(merged.Do) <- c("start","end","des")
                if (Alt.type == "ES" & !se.test){
                    final.result <- c(merge.mat(mer.tar,
                        merged.Do,merged.Up,Alt.type))
                    names(final.result) <- cn
                    final.result
                }
                else if ((Alt.type == "MXE" | Alt.type == "ES") & se.test){
                    secondEX <- strsplit(sub.ES[,"2ndEX"],"-")
                    secondEX <- unique(do.call(rbind,secondEX))
                    colnames(secondEX) <- c("start","end")
                    se.s <- as.integer(secondEX[,"start"])
                    se.e <- as.integer(secondEX[,"end"])
                    se.ea.ran <- IRanges(start=se.s,end=se.e)
                    secondEX.range <- GRanges(seqnames="*",ranges=se.ea.ran)
                    sorted.secondEX <- sortEX(secondEX.range,"ranges")
                    if (length(sorted.secondEX) > 1)  return (NULL)
                    sorted.secondEX <- sorted.secondEX[[1]]
                    ov.ex.ran <- findOverlaps(sfe,sorted.secondEX)
                    over.nums <- unique(as.matrix(ov.ex.ran)[,"subjectHits"])
                    if (any(length(over.nums))){
                        secondEX <- rbind(secondEX[-over.nums,])
                    }
                    if (!any(length(secondEX)))  return (NULL)
                    pse <- paste(secondEX[,"start"],secondEX[,"end"],sep="-")
                    second.des <- p.st(pse)
                    m.se <- cbind(rbind(c(min(as.double(secondEX[,"start"])),
                        max(as.double(secondEX[,"end"])))),rbind(second.des))
                    colnames(m.se) <- c("start","end","des")
                    final.result <- c(merge.mat(mer.tar,
                        m.se,merged.Do,Alt.type,merged.Up))
                    names(final.result) <- cn
                    final.result
                }
            })
            if (any(length(m.final.result))){
                m.final.result <- do.call(rbind,m.final.result)
                }
            return (m.final.result)
        }
        fi.ES <- rbind(altSplice[altSplice[,"2ndEX"] == "NA",])
        se.ES <- rbind(altSplice[altSplice[,"2ndEX"] != "NA",])
        se.ES <- rbind(altSplice[altSplice[,"Types"] != "ES",])
        mxe.ES <- rbind(altSplice[altSplice[,"Types"] == "MXE",])
        fi.ES.re <- NULL
        se.ES.re <- NULL
        mxe.ES.re <- NULL
        if (any(length(fi.ES))){
            fi.ES.re <- ES.merge.f(fi.ES,"ES")
        }
        if (any(length(se.ES))){
            se.ES.re <- ES.merge.f(se.ES,"ES")
        }
        if (any(length(mxe.ES))){
            mxe.ES.re <- ES.merge.f(mxe.ES,"MXE")
        }
        final.result <- rbind(fi.ES.re,se.ES.re,mxe.ES.re)
        return (final.result)
    }
    IR.Alt.result <- function(altSplice){
        rm.alt.result <- NULL
        IR.result <- altSplice
        Do.test.ex <- altSplice[,"DownEX"]
        Up.test.ex <- altSplice[,"UpEX"]
        Do.test.ex <- do.call(rbind,strsplit(Do.test.ex,"-"))
        Up.test.ex <- do.call(rbind,strsplit(Up.test.ex,"-"))
        colnames(Do.test.ex) <- c("start","end")
        colnames(Up.test.ex) <- c("start","end")
        firstEX <- cbind(Do.test.ex[,"end"],Up.test.ex[,"start"])
        colnames(firstEX) <- c("start","end")
        fi.st <- as.integer(firstEX[,"start"])
        fi.en <- as.integer(firstEX[,"end"])
        each.ranges <- IRanges(start=fi.st,end=fi.en)
        firstEX.range <- GRanges(seqnames="*",ranges=each.ranges)
        Do.test.ex <- IR.result[,"DownEX"]
        Up.test.ex <- IR.result[,"UpEX"]
        Do.test.ex <- do.call(rbind,strsplit(Do.test.ex,"-"))
        Up.test.ex <- do.call(rbind,strsplit(Up.test.ex,"-"))
        colnames(Do.test.ex) <- c("start","end")
        colnames(Up.test.ex) <- c("start","end")
        firstEX <- cbind(Do.test.ex[,"end"],Up.test.ex[,"start"])
        colnames(firstEX) <- c("start","end")
        fi.st <- as.integer(firstEX[,"start"])
        fi.en <- as.integer(firstEX[,"end"])
        each.ranges <- IRanges(start=fi.st,end=fi.en)
        firstEX.range <- GRanges(seqnames="*",ranges=each.ranges)
        sorted.first.EX <- sortEX(firstEX.range,"mat")
        merged.result <- NULL
        final.result <- lapply(sorted.first.EX,function(sfe){
            IR.dw <- do.call(rbind,strsplit(IR.result[,"DownEX"],"-"))[,2]
            IR.up <- do.call(rbind,strsplit(IR.result[,"UpEX"],"-"))[,1]
            each.do.up.ex <- paste(IR.dw,IR.up,sep="-")
            sf.ex <- paste(rbind(sfe)[,"start"],rbind(sfe)[,"end"],sep="-")
            each.result <- rbind(IR.result[is.element(each.do.up.ex,sf.ex),])
            re.ex <- unlist(strsplit(each.result[,"Retain_des"],","))
            do.ex <- unlist(strsplit(each.result[,"Do_des"],","))
            up.ex <- unlist(strsplit(each.result[,"Up_des"],","))
            re.ex <- do.call(rbind,strsplit(sort(re.ex),"-"))
            do.ex <- do.call(rbind,strsplit(sort(do.ex),"-"))
            up.ex <- do.call(rbind,strsplit(sort(up.ex),"-"))
            colnames(re.ex) <- c("start","end")
            colnames(do.ex) <- c("start","end")
            colnames(up.ex) <- c("start","end")
            p.re.ex <- paste(min(re.ex[,"start"]),max(re.ex[,"end"]),sep="-")
            p.do.ex <- paste(min(do.ex[,"start"]),max(do.ex[,"end"]),sep="-")
            p.up.ex <- paste(min(up.ex[,"start"]),max(up.ex[,"end"]),sep="-")
            add.re.ex <- unique(paste(do.ex[,"start"],up.ex[,"end"],sep="-"))
            re.des <- unlist(strsplit(each.result[,"Retain_des"],","))
            re.des <- p.st(c(re.des,add.re.ex))
            do.des <- p.st(c(unlist(strsplit(each.result[,"Do_des"],","))))
            up.des <- p.st(c(unlist(strsplit(each.result[,"Up_des"],","))))
            do.re.ex <- c(do.ex[,"start"],re.ex[,"start"])
            up.re.ex <- c(up.ex[,"end"],re.ex[,"end"])
            out.sd <- unique(paste(p.st(do.re.ex),p.st(up.re.ex),sep="-"))
            in.sd <- paste(p.st(do.ex[,"end"]),p.st(up.ex[,"start"]),sep="-")
            in.sd <- unique(in.sd)
            merged.result <- rbind(c(p.re.ex,p.do.ex,p.up.ex,re.des,do.des,
                up.des,out.sd,in.sd,"IR"))
        })
        final.result <- do.call(rbind,final.result)
        colnames(final.result) <- c("RetainEX","DownEX","UpEX","Retain_des",
            "Do_des","Up_des","Outter_splice","Inner_splice","Types")
        return (final.result)
    }
    for.f <- function(Alt.mat,alt.gene,alt.type){
        each.result <- NULL
        final.result <- foreach(alt.num=seq_along(alt.gene),
            .packages=called.packages,.combine=rbind) %dopar% {
            over.num <- Alt.mat[,"EnsID"] == alt.gene[alt.num]
            Alt.mat <- rbind(Alt.mat[over.num,])
            if (alt.type == "ES"){
                each.result <- rbind(ES.Alt.result(Alt.mat))
            }
            else if (alt.type == "ASS"){
                each.result <- rbind(ASS.Alt.result(Alt.mat))
            }
            else if (alt.type == "IR"){
                each.result <- rbind(IR.Alt.result(Alt.mat))
            }
            spl.nums <- grep("Spl|spl",colnames(each.result))
            e.chr <- unique(Alt.mat[,"Nchr"])
            each.result <- rbind(each.result[,-spl.nums])
            each.result <- cbind(EnsID=alt.gene[alt.num],Nchr=e.chr,
                Strand=unique(Alt.mat[,"Strand"]),each.result)
            each.result
        }
    return (final.result)
    }
    GTFdb <- chrseparate(GTFdb,1:22)
    registerDoParallel(cores=Ncor)
    trans.intron.range <- intronsByTranscript(GTFdb)
    tx.cns <- c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND")
    txTable <- try(select(GTFdb, keys=names(trans.intron.range),
        columns=tx.cns, keytype="TXID"),silent=TRUE)
    txTable <- gsub(" ","",as.matrix(txTable))
    trans.intron.range <- unlist(trans.intron.range)
    called.packages <- c("GenomicRanges","GenomicFeatures")
    Alt.splice.result <- ASdb@SplicingModel
    out.sd <- NULL
    in.sd <- NULL
    alt.num <- NULL
    ES.gene <- NULL
    ASS.gene <- NULL
    IR.gene <- NULL
    ES.num <- NULL
    ASS.num <- NULL
    IR.num <- NULL
    final.result <- NULL
    if (ncol(Alt.splice.result$ES) != 1){
        ES.gene <- unique(Alt.splice.result$ES[,"EnsID"])
    }
    if (ncol(Alt.splice.result$ASS) != 1){
        ASS.gene <- unique(Alt.splice.result$ASS[,"EnsID"])
    }
    if (ncol(Alt.splice.result$IR) != 1){
        IR.gene <- unique(Alt.splice.result$IR[,"EnsID"])
    }
    rm.num <- NULL
    final.ES.result <- NULL
    final.ASS.result <- NULL
    final.IR.result <- NULL
    total.list <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
    names(total.list) <- c("ES","ASS","IR")
    if(any(seq_along(ES.gene))){
        each.mat <- Alt.splice.result$ES
        each.re <- for.f(each.mat,ES.gene,"ES")
        if (any(length(each.re))){
            rownames(each.re) <- seq_len(nrow(each.re))
            ES.nms <- paste("ES",seq_len(nrow(each.re)),sep="")
            each.re <- cbind(Index=ES.nms,each.re)
            total.list$"ES" <- each.re
        }
    }
    if(any(seq_along(ASS.gene))){
        each.mat <- Alt.splice.result$ASS
        each.re <- for.f(each.mat,ASS.gene,"ASS")
        if (any(length(each.re))){
            rownames(each.re) <- seq_len(nrow(each.re))
            ASS.nms <- paste("ASS",seq_len(nrow(each.re)),sep="")
            each.re <- cbind(Index=ASS.nms,each.re)
            total.list$"ASS" <- each.re
        }
    }
    if(any(seq_along(IR.gene))){
        each.mat <- Alt.splice.result$IR
        each.re <- for.f(each.mat,IR.gene,"IR")
        if (any(length(each.re))){
            rownames(each.re) <- seq_len(nrow(each.re))
            IR.nms <- paste("IR",seq_len(nrow(each.re)),sep="")
            each.re <- cbind(Index=IR.nms,each.re)
            total.list$"IR" <- each.re
        }
    }
    ASdb <- new("ASdb",SplicingModel=total.list,Ratio=ASdb@"Ratio",
        GroupDiff=ASdb@"GroupDiff",sQTLs=ASdb@"sQTLs",
        Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
    return (ASdb)
}



