ExonsCluster <- function(ASdb,GTFdb,Ncor=1){
    sortEX <- function(test.range,outtype="mat"){
        over.ranges <- unique(as.matrix(findOverlaps(test.range,test.range,select="all")))
        sort.firstEX <- tapply(over.ranges[,2],over.ranges[,1],function(of){
            paste(of,collapse=",")
        })
        u.sort.firstEX <- unique(sort.firstEX)
        final.sorted.EX <- lapply(u.sort.firstEX,function(u.fe){
            over.num <- as.double(unlist(strsplit(u.fe,",")))
            if (outtype == "ranges"){
                test.range[over.num]
            }
            else {
                pre.result <- cbind(as.matrix(start(test.range[over.num])),as.matrix(end(test.range[over.num])))
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
        if (AtTypes == "A5SS" | AtTypes == "A3SS"){
            final.mat <- NULL
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
            total.des.mat <- lapply(u.dse.mat1,function(each.dse.mat1){
                if (AtTypes == "A5SS"){
                    merge.des.mat <- paste(each.dse.mat1,paste(u.dse.mat3,collapse=","),sep="-")
                    only.des2 <- paste(sort(u.dse.mat2[u.dse.mat2!=each.dse.mat1]),collapse=",")
                    only.des3 <- paste(sort(u.dse.mat3),collapse=",")
                    merge.des.mat <- c(merge.des.mat,paste(only.des2,only.des3,sep="-"))
                }
                if (AtTypes == "A3SS"){
                    merge.des.mat <- paste(paste(u.dse.mat3,collapse=","),each.dse.mat1,sep="-")
                    only.des2 <- paste(sort(u.dse.mat2[u.dse.mat2!=each.dse.mat1]),collapse=",")
                    only.des3 <- paste(sort(u.dse.mat3),collapse=",")
                    merge.des.mat <- c(merge.des.mat,paste(only.des3,only.des2,sep="-"))
                }
                mat1.ex.mat <- rbind(mat1[mat1[,mat.stan[1]]==each.dse.mat1,c("start","end","des")])
                mat2.ex.mat <- rbind(mat2[mat2[,mat.stan[2]]!=each.dse.mat1,c("start","end","des")])
                mat3.ex.mat <- rbind(mat3[,c("start","end")])
                mat1.ex.des <- sort(mat1.ex.mat[,"des"])
                mat2.ex.des <- sort(mat2.ex.mat[,"des"])
                mat3.ex.des <- sort(paste(mat3.ex.mat[,"start"],mat3.ex.mat[,"end"],sep="-"))
                mat1.ex.pa <- paste(min(mat1.ex.mat[,"start"]),max(mat1.ex.mat[,"end"]),sep="-")
                mat2.ex.pa <- paste(min(mat2.ex.mat[,"start"]),max(mat2.ex.mat[,"end"]),sep="-")
                mat3.ex.pa <- paste(min(mat3.ex.mat[,"start"]),max(mat3.ex.mat[,"end"]),sep="-")
                only.des1 <- paste(mat1.ex.des,collapse=",")
                only.des2 <- paste(mat2.ex.des,collapse=",")
                only.des3 <- paste(mat3.ex.des,collapse=",")
                merge.mat <- c(mat1.ex.pa,mat2.ex.pa,mat3.ex.pa,only.des1,only.des2,only.des3,merge.des.mat)
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
            if (length(mat4) == 0){
                final.des.mat <- paste(p.st(dse.mat2[,"end"]),p.st(dse.mat1[,"start"]),sep="-")
                final.des.mat <- c(final.des.mat,"NA")
                final.des.mat <- c(final.des.mat,paste(p.st(dse.mat1[,"end"]),p.st(dse.mat3[,"start"]),sep="-"))
                final.mat <- c(p.mat1,"NA",p.mat2,p.mat3,mat1[,"des"],"NA",mat2[,"des"],mat3[,"des"],final.des.mat,"ES")
            }
            if(length(mat4) != 0){
                p.mat4 <- paste(mat4[,"start"],mat4[,"end"],sep="-")
                dse.mat4 <- unlist(strsplit(mat4[,"des"],","))
                dse.mat4 <- do.call(rbind,strsplit(dse.mat4,"-"))
                colnames(dse.mat4) <- c("start","end")
                if(AtTypes == "ES"){
                    final.des.mat <- paste(p.st(dse.mat3[,"end"]),p.st(dse.mat1[,"start"]),sep="-")
                    final.des.mat <- c(final.des.mat,paste(p.st(dse.mat1[,"end"]),p.st(dse.mat2[,"start"]),sep="-"))
                    final.des.mat <- c(final.des.mat,paste(p.st(dse.mat2[,"end"]),p.st(dse.mat4[,"start"]),sep="-"))
                }
                if(AtTypes == "MXE"){
                    fir.spli <- paste(p.st(dse.mat3[,"end"]),p.st(dse.mat1[,"start"]),sep="-")
                    sec.spli <- paste(p.st(dse.mat1[,"end"]),p.st(dse.mat4[,"start"]),sep="-")
                    thi.spli <- paste(p.st(dse.mat3[,"end"]),p.st(dse.mat2[,"start"]),sep="-")
                    for.spli <- paste(p.st(dse.mat2[,"end"]),p.st(dse.mat4[,"start"]),sep="-")
                    final.des.mat <- paste(fir.spli,"|",sec.spli,sep="")
                    final.des.mat <- c(final.des.mat,"NA")
                    final.des.mat <- c(final.des.mat,paste(thi.spli,"|",for.spli,sep=""))
                }
                p.mat <- c(p.mat1,p.mat2,p.mat3,p.mat4)
                des.mat <- c(mat1[,"des"],mat2[,"des"],mat3[,"des"],mat4[,"des"])
                final.mat <- c(p.mat,des.mat,final.des.mat,AtTypes)
            }
        }
        return (final.mat)
    }
    
    ASS.Alt.result <- function(altSplice,total.intron.mat){
        neig.nm <- c("ShortNeighbor_des","LongNeighborEX")
        ASS.type <- unique(altSplice[,"Types"])
        rm.alt.result <- NULL
        A5SS.final.result <- NULL
        A3SS.final.result <- NULL
        if(is.element("A5SS",ASS.type)){
            A5SS.result <- rbind(altSplice[altSplice[,"Types"] == "A5SS",])
            Short.EX <- unique(do.call(rbind,strsplit(A5SS.result[,"ShortEX"],"-")))
            long.EX <- unique(do.call(rbind,strsplit(A5SS.result[,"LongEX"],"-")))
            colnames(Short.EX) <- c("start","end")
            colnames(long.EX) <- c("start","end")
            total.EX <- rbind(Short.EX,long.EX)
            each.ranges <- IRanges(start=as.integer(total.EX[,"start"]),end=as.integer(total.EX[,"end"]))
            total.EX.ranges <- GRanges(seqnames="*",ranges=each.ranges)
            Short.EX <- strsplit(unlist(strsplit(A5SS.result[,"Short_des"],",")),"-")
            Short.EX <- unique(do.call(rbind,Short.EX))
            long.EX <- strsplit(unlist(strsplit(A5SS.result[,"Long_des"],",")),"-")
            long.EX <- unique(do.call(rbind,long.EX))
            neighbor.EX <- strsplit(unlist(strsplit(A5SS.result[,neig.nm],",")),"-")
            neighbor.EX <- unique(do.call(rbind,neighbor.EX))
            colnames(Short.EX) <- c("start","end")
            colnames(long.EX) <- c("start","end")
            colnames(neighbor.EX) <- c("start","end")
            pa.ex.mat <- paste(Short.EX[,"start"],Short.EX[,"end"],sep="-")
            short.ex.merge <- tapply(pa.ex.mat,Short.EX[,"end"],function(total.sem){
                sem <- do.call(rbind,strsplit(total.sem,"-"))
                colnames(sem) <- c("start","end")
                merge.descrip <- paste(sort(total.sem),collapse=",")
                cbind(rbind(sem[order(as.double(sem[,"start"]),decreasing=FALSE)[1],]),rbind(merge.descrip))
            })
            pa.ex.mat <- paste(long.EX[,"start"],long.EX[,"end"],sep="-")
            long.ex.merge <- tapply(pa.ex.mat,long.EX[,"end"],function(total.sem){
                sem <- do.call(rbind,strsplit(total.sem,"-"))
                colnames(sem) <- c("start","end")
                merge.descrip <- paste(sort(total.sem),collapse=",")
                cbind(rbind(sem[order(as.double(sem[,"start"]),decreasing=FALSE)[1],]),rbind(merge.descrip))
            })
            short.ex.merge <- do.call(rbind,short.ex.merge)
            long.ex.merge <- do.call(rbind,long.ex.merge)
            short.possible.result <- lapply(short.ex.merge[,"end"],function(s.e){
                short.possible <- cbind(s.e,neighbor.EX[,"start"],"possible")
                colnames(short.possible) <- c("start","end","status")
                short.test.int <- paste(short.possible[,"start"],short.possible[,"end"],sep="-")
                short.possible[is.element(short.test.int,total.intron.mat),"status"] <- "exist"
                paste(short.possible[,"start"],"-",short.possible[,"end"],sep="")
            })
            long.possible.result <- lapply(long.ex.merge[,"end"],function(s.e){
                long.possible <- cbind(s.e,neighbor.EX[,"start"],"possible")
                colnames(long.possible) <- c("start","end","status")
                long.test.int <- paste(long.possible[,"start"],long.possible[,"end"],sep="-")
                long.possible[is.element(long.test.int,total.intron.mat),"status"] <- "exist"
                paste(long.possible[,"start"],"-",long.possible[,"end"],sep="")
            })
            short.possible.result <- do.call(rbind,short.possible.result)
            long.possible.result <- do.call(rbind,long.possible.result)
            short.long.des <- c(p.st(short.possible.result),p.st(long.possible.result))
            neighbor.ex.merge <- cbind(neighbor.EX,paste(neighbor.EX[,"start"],neighbor.EX[,"end"],sep="-"))
            colnames(neighbor.ex.merge) <- c("start","end","des")#,"des.1","des.2")
            colnames(short.ex.merge) <- c("start","end","des")
            colnames(long.ex.merge) <- c("start","end","des")
            A5SS.final.result <- cbind(merge.mat(short.ex.merge,long.ex.merge,neighbor.ex.merge,"A5SS"),"A5SS")
            colnames(A5SS.final.result) <- c("ShortEX","LongEX","NeighborEX","Short_des","Long_des",
            "Neighbor_des","splicing in 1stEX ","splicing in 2ndEX","Types")
        }
        
        
        if(is.element("A3SS",ASS.type)){
            A3SS.result <- rbind(altSplice[altSplice[,"Types"] == "A3SS",])
            Short.EX <- unique(do.call(rbind,strsplit(A3SS.result[,"ShortEX"],"-")))
            long.EX <- unique(do.call(rbind,strsplit(A3SS.result[,"LongEX"],"-")))
            colnames(Short.EX) <- c("start","end")
            colnames(long.EX) <- c("start","end")
            total.EX <- rbind(Short.EX,long.EX)
            each.ranges <- IRanges(start=as.integer(total.EX[,"start"]),end=as.integer(total.EX[,"end"]))
            total.EX.ranges <- GRanges(seqnames="*",ranges=each.ranges)
            Short.EX <- unique(do.call(rbind,strsplit(unlist(strsplit(A3SS.result[,"Short_des"],",")),"-")))
            long.EX <- unique(do.call(rbind,strsplit(unlist(strsplit(A3SS.result[,"Long_des"],",")),"-")))
            neighbor.EX <- unique(do.call(rbind,strsplit(unlist(strsplit(A3SS.result[,neig.nm],",")),"-")))
            colnames(Short.EX) <- c("start","end")
            colnames(long.EX) <- c("start","end")
            colnames(neighbor.EX) <- c("start","end")
            pa.ex.mat <- paste(Short.EX[,"start"],Short.EX[,"end"],sep="-")
            short.ex.merge <- tapply(pa.ex.mat,Short.EX[,"start"],function(total.sem){
                sem <- do.call(rbind,strsplit(total.sem,"-"))
                colnames(sem) <- c("start","end")
                merge.descrip <- paste(sort(total.sem),collapse=",")
                cbind(rbind(sem[order(as.double(sem[,"end"]),decreasing=TRUE)[1],]),rbind(merge.descrip))
            })
            pa.ex.mat <- paste(long.EX[,"start"],long.EX[,"end"],sep="-")
            long.ex.merge <- tapply(pa.ex.mat,long.EX[,"start"],function(total.sem){
                sem <- do.call(rbind,strsplit(total.sem,"-"))
                colnames(sem) <- c("start","end")
                merge.descrip <- paste(sort(total.sem),collapse=",")
                cbind(rbind(sem[order(as.double(sem[,"end"]),decreasing=TRUE)[1],]),rbind(merge.descrip))
            })
            short.ex.merge <- do.call(rbind,short.ex.merge)
            long.ex.merge <- do.call(rbind,long.ex.merge)
            short.possible.result <- lapply(short.ex.merge[,"start"],function(s.e){
                short.possible <- cbind(neighbor.EX[,"end"],s.e,"possible")
                colnames(short.possible) <- c("start","end","status")
                short.test.int <- paste(short.possible[,"start"],short.possible[,"end"],sep="-")
                short.possible[is.element(short.test.int,total.intron.mat),"status"] <- "exist"
                paste(short.possible[,"start"],"-",short.possible[,"end"],sep="")
            })
            long.possible.result <- lapply(long.ex.merge[,"start"],function(s.e){
                long.possible <- cbind(neighbor.EX[,"end"],s.e,"possible")
                colnames(long.possible) <- c("start","end","status")
                long.test.int <- paste(long.possible[,"start"],long.possible[,"end"],sep="-")
                long.possible[is.element(long.test.int,total.intron.mat),"status"] <- "exist"
                paste(long.possible[,"start"],"-",long.possible[,"end"],sep="")
            })
            short.possible.result <- do.call(rbind,short.possible.result)
            long.possible.result <- do.call(rbind,long.possible.result)
            short.long.des <- c(paste(short.possible.result,collapse=","),paste(long.possible.result,collapse=","))
            neighbor.ex.merge <- cbind(neighbor.EX,paste(neighbor.EX[,"start"],neighbor.EX[,"end"],sep="-"))
            colnames(neighbor.ex.merge) <- c("start","end","des")#,"des.1","des.2")
            colnames(short.ex.merge) <- c("start","end","des")
            colnames(long.ex.merge) <- c("start","end","des")
            A3SS.final.result <- cbind(merge.mat(short.ex.merge,long.ex.merge,neighbor.ex.merge,"A3SS"),"A3SS")
            colnames(A3SS.final.result) <- c("ShortEX","LongEX","NeighborEX","Short_des","Long_des",
            "Neighbor_des","splicing in 1stEX ","splicing in 2ndEX","Types")
        }
        final.result <- rbind(A5SS.final.result,A3SS.final.result)
        return (final.result)
    }
    ES.Alt.result <- function(altSplice){
        rm.alt.result <- NULL
        #ES.result <- altSplice[altSplice[,"2ndEX"] == "NA",]
        ES.result <- altSplice
        firstEX <- altSplice[,"1stEX"]
        firstEX <- unique(do.call(rbind,strsplit(firstEX,"-")))
        colnames(firstEX) <- c("start","end")
        each.ranges <- IRanges(start=as.integer(firstEX[,"start"]),end=as.integer(firstEX[,"end"]))
        firstEX.range <- GRanges(seqnames="*",ranges=each.ranges)
        first.ES.result <- rbind(ES.result[ES.result[,"2ndEX"] == "NA",])
        firstEX <- first.ES.result[,"1stEX"]
        firstEX <- unique(do.call(rbind,strsplit(firstEX,"-")))
        merged.final.result <- NULL
        if (length(firstEX) != 0){ # replace
            colnames(firstEX) <- c("start","end")
            each.ranges <- IRanges(start=as.integer(firstEX[,"start"]),end=as.integer(firstEX[,"end"]))
            firstEX.range <- GRanges(seqnames="*",ranges=each.ranges)
            sorted.first.EX <- sortEX(firstEX.range,"mat")
            merged.final.result <- lapply(sorted.first.EX,function(sfe){
                rownames(sfe) <- 1:nrow(sfe)
                each.first.EX <- sfe
                over.ex <- is.element(first.ES.result[,"1stEX"],paste(each.first.EX[,"start"],each.first.EX[,"end"],sep="-"))
                Do.Up.EX <- rbind(first.ES.result[over.ex,c("Do_des","Up_des")])
                DoEX <- unique(Do.Up.EX[,"Do_des"])
                DoEX <- unlist(strsplit(DoEX,","))
                DoEX <- unique(do.call(rbind,strsplit(DoEX,"-")))
                UpEX <- unique(Do.Up.EX[,"Up_des"])
                UpEX <- unlist(strsplit(UpEX,","))
                UpEX <- unique(do.call(rbind,strsplit(UpEX,"-")))
                colnames(DoEX) <- c("start","end")
                colnames(UpEX) <- c("start","end")
                Do.des <- p.st(paste(DoEX[,"start"],DoEX[,"end"],sep="-"))
                up.des <- p.st(paste(UpEX[,"start"],UpEX[,"end"],sep="-"))
                tar.des <- p.st(paste(each.first.EX[,"start"],each.first.EX[,"end"],sep="-"))
                merged.Do <- cbind(rbind(c(min(as.double(DoEX[,"start"])),
                    max(as.double(DoEX[,"end"])))),rbind(Do.des))
                merged.Up <- cbind(rbind(c(min(as.double(UpEX[,"start"])),
                    max(as.double(UpEX[,"end"])))),rbind(up.des))
                merged.tar <- cbind(rbind(c(min(as.double(each.first.EX[,"start"])),
                    max(as.double(each.first.EX[,"end"])))),rbind(tar.des))
                colnames(merged.tar) <- c("start","end","des")
                colnames(merged.Up) <- c("start","end","des")
                colnames(merged.Do) <- c("start","end","des")
                final.result <- c(merge.mat(merged.tar,merged.Do,merged.Up,"ES"))
                names(final.result) <- c("1stEX","2ndEX","DownEX","UpEX","1st_des",
                    "2nd_des","Do_des","Up_des","1stSpl","2ndSpl","3rdSpl","Types")
                final.result
            })
            merged.final.result <- do.call(rbind,merged.final.result)
        }
        sec.test <- ES.result[,"2ndEX"] != "NA" & ES.result[,"Types"] == "ES"
        MXE.test <- ES.result[,"2ndEX"] != "NA" & ES.result[,"Types"] == "MXE"
        sec.ES.result <- rbind(ES.result[sec.test,])
        MXE.ES.result <- rbind(ES.result[MXE.test,])
        SE.MXE <- unique(ES.result[ES.result[,"2ndEX"] != "NA","Types"])
        if (length(SE.MXE) == 0) return (list(merged.final.result,rm.alt.result))
        sec.merged.final.result <- lapply(SE.MXE,function(se.mxe.type){
            sec.test <- ES.result[,"2ndEX"] != "NA" & ES.result[,"Types"] == se.mxe.type
            sec.ES.result <- rbind(ES.result[sec.test,])
            firstEX <- sec.ES.result[,"1stEX"]
            secondEX <- sec.ES.result[,"2ndEX"]
            firstEX <- unique(do.call(rbind,strsplit(firstEX,"-")))
            secondEX <- unique(do.call(rbind,strsplit(secondEX,"-")))
            colnames(firstEX) <- c("start","end")
            colnames(secondEX) <- c("start","end")
            fi.each.range <- IRanges(start=as.integer(firstEX[,"start"]),end=as.integer(firstEX[,"end"]))
            se.each.range <- IRanges(start=as.integer(secondEX[,"start"]),end=as.integer(secondEX[,"end"]))
            firstEX.range <- GRanges(seqnames="*",ranges=fi.each.range)
            secondEX.range <- GRanges(seqnames="*",ranges=se.each.range)
            firstEX.range <- sortEX(firstEX.range,"ranges")
            sorted.secondEX <- do.call(rbind,(sortEX(secondEX.range)))
            each.ranges <- IRanges(start=as.integer(sorted.secondEX[,"start"]),end=as.integer(sorted.secondEX[,"end"]))
            sorted.secondEX.range <- GRanges(seqnames="*",ranges=each.ranges)
            pre.merged.final.result <- lapply(firstEX.range,function(fi.es.ex){
                over.ex.ranges <- findOverlaps(fi.es.ex,sorted.secondEX.range,select="all")
                over.ex.ranges <- unique(as.matrix(over.ex.ranges)[,"subjectHits"])
                sub.sorted.secondEX <- sorted.secondEX
                sub.firstEX <- cbind(start(fi.es.ex),end(fi.es.ex))
                colnames(sub.firstEX) <- c("start","end")
                if(length(over.ex.ranges) != 0){
                    sub.sorted.secondEX <- rbind(sorted.secondEX[-over.ex.ranges,])
                }
                sub.fi.ex <- paste(sub.firstEX[,"start"],sub.firstEX[,"end"],sep="-")
                sub.sec.ex <- paste(sub.sorted.secondEX[,"start"],sub.sorted.secondEX[,"end"],sep="-")
                Do.Up.fir.EX <- rbind(sec.ES.result[is.element(sec.ES.result[,"1stEX"],sub.fi.ex),])
                sub.sorted.secondEX <- sub.sorted.secondEX[is.element(sub.sec.ex,Do.Up.fir.EX[,"2ndEX"]),]
                sub.sorted.secondEX <- rbind(sub.sorted.secondEX)
                each.merged.final.result <- lapply(sub.sec.ex,function(se.es.ex){
                    sub.secondEX <- rbind(unlist(strsplit(se.es.ex,"-")))
                    colnames(sub.secondEX) <- c("start","end")
                    total.sub.sec.ex <- paste(sub.secondEX[,"start"],sub.secondEX[,"end"],sep="-")
                    Do.Up.EX <- rbind(Do.Up.fir.EX[is.element(Do.Up.fir.EX[,"2ndEX"],total.sub.sec.ex),])
                    DoEX <- unique(Do.Up.EX[,"Do_des"])
                    DoEX <- unlist(strsplit(DoEX,","))
                    DoEX <- unique(do.call(rbind,strsplit(DoEX,"-")))
                    UpEX <- unique(Do.Up.EX[,"Up_des"])
                    UpEX <- unlist(strsplit(UpEX,","))
                    UpEX <- unique(do.call(rbind,strsplit(UpEX,"-")))
                    colnames(DoEX) <- c("start","end")
                    colnames(UpEX) <- c("start","end")
                    colnames(sub.firstEX) <- c("start","end")
                    colnames(sub.firstEX) <- c("start","end")
                    Do.des <- paste(sort(unique(paste(DoEX[,"start"],DoEX[,"end"],sep="-"))),collapse=",")
                    up.des <- paste(sort(unique(paste(UpEX[,"start"],UpEX[,"end"],sep="-"))),collapse=",")
                    first.des <- paste(sort(unique(paste(sub.firstEX[,"start"],sub.firstEX[,"end"],sep="-"))),collapse=",")
                    second.des <- p.st(total.sub.sec.ex)
                    merged.Do <- cbind(rbind(c(min(as.double(DoEX[,"start"])),max(as.double(DoEX[,"end"])))),rbind(Do.des))
                    merged.Up <- cbind(rbind(c(min(as.double(UpEX[,"start"])),max(as.double(UpEX[,"end"])))),rbind(up.des))
                    merged.fir <- cbind(min(sub.firstEX[,"start"]),max(sub.firstEX[,"end"]),rbind(first.des))
                    merged.sec <- cbind(min(sub.secondEX[,"start"]),max(sub.secondEX[,"end"]),rbind(second.des))
                    colnames(merged.fir) <- c("start","end","des")
                    colnames(merged.sec) <- c("start","end","des")
                    colnames(merged.Up) <- c("start","end","des")
                    colnames(merged.Do) <- c("start","end","des")
                    final.result <- c(merge.mat(merged.fir,merged.sec,merged.Do,se.mxe.type,merged.Up))
                    names(final.result) <- rbind(c("1stEX","2ndEX","DownEX","UpEX","1st_des",
                    "2nd_des","Do_des","Up_des","1stSpl","2ndSpl","3rdSpl","Types"))
                    final.result
                })
                do.call(rbind,each.merged.final.result)
            })
            pre.merged.final.result <- do.call(rbind,pre.merged.final.result)
            pre.merged.final.result
        })
        sec.merged.final.result <- do.call(rbind,sec.merged.final.result)
        merged.final.result <- rbind(merged.final.result,sec.merged.final.result)
        return (merged.final.result)
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
        each.ranges <- IRanges(start=as.integer(firstEX[,"start"]),end=as.integer(firstEX[,"end"]))
        firstEX.range <- GRanges(seqnames="*",ranges=each.ranges)
        Do.test.ex <- IR.result[,"DownEX"]
        Up.test.ex <- IR.result[,"UpEX"]
        Do.test.ex <- do.call(rbind,strsplit(Do.test.ex,"-"))
        Up.test.ex <- do.call(rbind,strsplit(Up.test.ex,"-"))
        colnames(Do.test.ex) <- c("start","end")
        colnames(Up.test.ex) <- c("start","end")
        firstEX <- cbind(Do.test.ex[,"end"],Up.test.ex[,"start"])
        colnames(firstEX) <- c("start","end")
        firstEX.range <- GRanges(seqnames="*",ranges=IRanges(start=as.integer(firstEX[,"start"]),end=as.integer(firstEX[,"end"])))
        sorted.first.EX <- sortEX(firstEX.range,"mat")
        merged.final.result <- lapply(sorted.first.EX,function(sfe){
            IR.dw <- do.call(rbind,strsplit(IR.result[,"DownEX"],"-"))[,2]
            IR.up <- do.call(rbind,strsplit(IR.result[,"UpEX"],"-"))[,1]
            each.do.up.ex <- paste(IR.dw,IR.up,sep="-")
            sfe.ex <- paste(rbind(sfe)[,"start"],rbind(sfe)[,"end"],sep="-")
            each.result <- rbind(IR.result[is.element(each.do.up.ex,sfe.ex),])
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
            re.des <- p.st(c(unlist(strsplit(each.result[,"Retain_des"],",")),add.re.ex))
            do.des <- p.st(c(unlist(strsplit(each.result[,"Do_des"],","))))
            up.des <- p.st(c(unlist(strsplit(each.result[,"Up_des"],","))))
            out.splicing.des <- unique(paste(paste(sort(unique(c(do.ex[,"start"],re.ex[,"start"]))),collapse=","),
            paste(sort(unique(c(up.ex[,"end"],re.ex[,"end"]))),collapse=","),sep="-"))
            in.splicing.des <- unique(paste(paste(sort(unique(do.ex[,"end"])),collapse=","),
            paste(sort(unique(up.ex[,"start"])),collapse=","),sep="-"))
            merged.result <- rbind(c(p.re.ex,p.do.ex,p.up.ex,re.des,do.des,
            up.des,out.splicing.des,in.splicing.des,"IR"))
        })
        merged.final.result <- do.call(rbind,merged.final.result)
        colnames(merged.final.result) <- c("RetainEX","DownEX","UpEX","Retain_des",
            "Do_des","Up_des","Outter_splice","Inner_splice","Types")
        return (merged.final.result)
    }
    GTFdb <- chrseparate(GTFdb,1:22)
    registerDoParallel(cores=Ncor)
    trans.intron.range <- intronsByTranscript(GTFdb)
    tx.cns <- c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND")
    txTable <- try(select(GTFdb, keys=names(trans.intron.range), columns=tx.cns, keytype="TXID"),silent=TRUE)
    txTable <- gsub(" ","",as.matrix(txTable))
    trans.intron.range <- unlist(trans.intron.range)
    called.packages <- c("GenomicRanges","GenomicFeatures")
    Alt.splice.result <- ASdb@SplicingModel
    ES.gene <- NULL
    ASS.gene <- NULL
    IR.gene <- NULL
    ES.num <- NULL
    ASS.num <- NULL
    IR.num <- NULL
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
    if(length(ES.gene != 0)){
        final.ES.result <- foreach(ES.num=1:length(ES.gene),.packages=called.packages,.combine=rbind) %dopar% { 
            ES.mat <- rbind(Alt.splice.result$ES[Alt.splice.result$ES[,"EnsID"] == ES.gene[ES.num],])
            total.each.result <- rbind(ES.Alt.result(ES.mat))
            total.each.result <- rbind(total.each.result[,-grep("Spl",colnames(total.each.result))])
            total.each.result <- cbind(EnsID=ES.gene[ES.num],Nchr=unique(ES.mat[,"Nchr"]),Strand=unique(ES.mat[,"Strand"]),total.each.result)
        }
    }
    if(length(ASS.gene != 0)){
        final.ASS.result <- foreach(ASS.num=1:length(ASS.gene),.packages=called.packages,.combine=rbind) %dopar% {
            TX.nums <- txTable[txTable[,"GENEID"] == ASS.gene[ASS.num],"TXID"]
            intron.mat <- trans.intron.range[is.element(names(trans.intron.range),TX.nums),]
            intron.mat <- paste(start(intron.mat)-1,end(intron.mat)+1,sep="-")
            ASS.mat <- rbind(Alt.splice.result$ASS[Alt.splice.result$ASS[,"EnsID"] == ASS.gene[ASS.num],])
            total.each.result <- rbind(ASS.Alt.result(ASS.mat,intron.mat))
            total.each.result <- rbind(total.each.result[,-grep("spl",colnames(total.each.result))])
            total.each.result <- cbind(EnsID=ASS.gene[ASS.num],Nchr=unique(ASS.mat[,"Nchr"]),Strand=unique(ASS.mat[,"Strand"]),total.each.result)
        }
    }
    if(length(IR.gene != 0)){
        final.IR.result <- foreach(IR.num=1:length(IR.gene),.packages=called.packages,.combine=rbind) %dopar% {
            IR.mat <- rbind(Alt.splice.result$IR[Alt.splice.result$IR[,"EnsID"] == IR.gene[IR.num],])
            total.each.result <- rbind(IR.Alt.result(IR.mat))
            total.each.result <- rbind(total.each.result[,-grep("Spl",colnames(total.each.result))])
            total.each.result <- cbind(EnsID=IR.gene[IR.num],Nchr=unique(IR.mat[,"Nchr"]),Strand=unique(IR.mat[,"Strand"]),total.each.result)
        }
    }
    if (length(final.ES.result) != 0){
        final.ES.result <- cbind(Index=paste("ES",1:nrow(final.ES.result),sep=""),final.ES.result)
        rownames(final.ES.result) <- 1:nrow(final.ES.result)
        total.list$"ES" <- final.ES.result
    }
    if (length(final.ASS.result) != 0){
        final.ASS.result <- cbind(Index=paste("ASS",1:nrow(final.ASS.result),sep=""),final.ASS.result)
        total.list$"ASS" <- final.ASS.result
    }
    if (length(final.IR.result) != 0){
        final.IR.result <- cbind(Index=paste("IR",1:nrow(final.IR.result),sep=""),final.IR.result)
        total.list$"IR" <- final.IR.result
    }
    ASdb <- new("ASdb",SplicingModel=total.list,Ratio=ASdb@"Ratio",GroupDiff=ASdb@"GroupDiff",
        sQTLs=ASdb@"sQTLs",Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
    return (ASdb)
}
