SplicingReads <- function(samfile=NULL,test.exon=NULL,spli.jun=NULL,readsinfo="paired",inse=40){
    ReadCover <- function(paired.read,met="paired"){
        read.count <- nrow(paired.read)
        startnum <- NULL                                
        pre.final.g.range <- lapply(as.integer(c(1:read.count)),function(read.num){
            final.g.mat <- NULL
            meta.info <- NULL
            each.read <- rbind(paired.read[read.num,])
            each.name <- each.read[,"qname"]
            cigar.info <- cigarParsing(each.read)
            each.g.range <- NULL
            if (length(cigar.info) != 0){
                N.final <- NULL
                chrnum <- as.matrix(each.read[,"rname"])
                startnum <- as.matrix(each.read[,"pos"])
                meta.count <- 0
                cigar.count <- 1:nrow(cigar.info)
                final.g.mat <- do.call(rbind,lapply(cigar.count,function(j){
                    g.mat <- NULL
                    name.cigar.info <- rownames(cigar.info)
                    pre.info <- cigar.count[cigar.count < j]
                    t.ci <- c("=","M","X","D","N")
                    testable.cigar <- pre.info[is.element(rownames(cigar.info)[pre.info],t.ci)]
                    if (length(testable.cigar) != 0){
                        startnum <- as.integer(startnum) + sum(as.integer(cigar.info[testable.cigar]))
                    }
                    if (met != "N"){
                        if (is.element(name.cigar.info[j],c("=","M"))){
                            g.mat <- c(startnum,as.integer(startnum)+as.integer(cigar.info[j,1])-1)
                            g.mat
                        }
                    }
                    else if(met == "N"){
                        if (is.element(name.cigar.info[j],c("N"))){
                            st.n <- as.integer(startnum)-1
                            en.n <- as.integer(startnum)+as.integer(cigar.info[j,1])
                            N.final <- rbind(N.final,c(st.n,en.n))
                            N.final
                        }
                    }
                }))
                if (met == "over.rm"){
                    each.diff <- apply(final.g.mat,1,function(x) diff(as.integer(x))) + 1
                    final.g.mat <- final.g.mat[each.diff / sum(each.diff) > 0.1]
                    final.g.mat <- matrix(final.g.mat,ncol=2,byrow=FALSE)
                    ch.cig <- NULL
                    if (nrow(final.g.mat) >1){
                        ch.cig <- unlist(lapply(1:as.integer(nrow(final.g.mat)-1),function(g.num){
                            di.v <- diff(c(max(as.integer(final.g.mat[g.num,])),min(as.integer(final.g.mat[g.num+1,]))))
                            paste(diff(as.integer(final.g.mat[g.num,]))+1,"M",di.v-1,"N",sep="")
                        }))
                        ch.cig <- paste(ch.cig,collapse="")
                    }
                    ch.cig <- paste(ch.cig,diff(as.integer(final.g.mat[nrow(final.g.mat),]))+1,"M",sep="")
                    final.g.mat <- c(min(final.g.mat),max(final.g.mat),ch.cig)
                    final.g.mat <- cbind(c(read.num),c(chrnum),final.g.mat[1],final.g.mat[2],
                        each.name,each.read[,"mapq"],final.g.mat[3])
                    each.g.range <- final.g.mat
                    colnames(each.g.range) <- c("Num","chr","start","end","readsinfo","mapq","cigar")
                }
                if (length(final.g.mat) == 0){
                    final.g.mat <- cbind(as.matrix(startnum),as.matrix(startnum))
                }
                if (met != "over.rm"){
                    each.g.range <- cbind(c(read.num),c(chrnum),final.g.mat[,1],
                        final.g.mat[,2],each.name,each.read[,"mapq"])
                }
                if (met == "N"){
                    each.g.range <- cbind(each.name,final.g.mat)
                    colnames(each.g.range) <- c("reads","start","end")
                }
            }
            each.g.range})
        final.g.range <- do.call(rbind,pre.final.g.range)
        if (met == "N" | met == "over.rm"){
            return (final.g.range)
        }
        colnames(final.g.range) <- c("Num","chr","start","end","readsinfo","mapq")
        return (final.g.range)
    }
    final.exon.result <- NULL
    if (length(colnames(test.exon)) == 0 & length(test.exon) != 0){
        test.exon <- do.call(rbind,strsplit(test.exon,"-"))
        colnames(test.exon) <- c("start","end")
        te.en <- test.exon[is.element(rownames(test.exon),c("DownEX","1stEX","2ndEX")),"end"]
        te.st <- test.exon[is.element(rownames(test.exon),c("1stEX","2ndEX","UpEX")),"start"]
        spli.sites <- cbind(te.en,te.st)
        colnames(spli.sites) <- colnames(test.exon)
        final.exon.result <- c(1:nrow(test.exon))
        names(final.exon.result) <- paste(test.exon[,"start"],test.exon[,"end"],sep="-")
    }
    if (length(colnames(spli.jun)) == 0 & length(spli.jun) != 0){
        total.spli.sites <- unique(do.call(rbind,strsplit(spli.jun,"-")))
        colnames(total.spli.sites) <- c("start","end")
    }        
    return.mat <- list(0,0,0)
    names(return.mat) <- c("pairedInfo","exonInfo","junctionInfo")
    if (length(samfile) == 0) return (return.mat)
    samfile <- unique(rbind(samfile))
    if (nrow(samfile) < 10) return (return.mat)
    Am.result <- NULL
    t.sam.reads <- table(paste(samfile[,"qname"],samfile[,"seq"]))
    sam.q.n <- paste(samfile[,"qname"],samfile[,"seq"])
    u.reads <- names(which(t.sam.reads < 2))
    samfile <- unique(rbind(samfile[is.element(sam.q.n,u.reads),]))
    if (nrow(samfile) < 10) return (return.mat)
    N.num <- grep("N",samfile[,"cigar"])
    if (readsinfo == "single" | readsinfo == "paired" | readsinfo == "exon"){
        samfile2 <- samfile
        not.N.reads <- samfile
        if (length(N.num) != 0){
            not.N.reads <- unique(rbind(samfile[-N.num,]))
        }
    } 
    len.test1 <- (nrow(not.N.reads) < 15 & length(N.num) < 5)
    len.test2 <- (nrow(not.N.reads) < 5 & length(N.num) < 15)
    if (len.test1 | len.test2) return (return.mat)
    pre.cover.result <- NULL
    final.cover.result <- NULL
    paired.result <- NULL
    used.reads.nm <- NULL
    if (readsinfo != "single"){
        not.N.reads -> test.reads
        if (length(test.reads) != 0){
            u.reads.name <- table(as.matrix(paste(test.reads[,"qname"],test.reads[,"pos"],sep="_-_")))
            s.u.reads.name <- do.call(rbind,strsplit(names(u.reads.name)[u.reads.name == 1],"_-_"))[,1]
            test.reads <- rbind(test.reads[is.element(test.reads[,"qname"],s.u.reads.name),])
            test.reads.mat <- ReadCover(test.reads)
            te.chr.names <- as.character(test.reads.mat[1,"chr"])
            pa.chr.names <- test.reads.mat[,"chr"]
            te.iranges <- IRanges(start=as.integer(test.exon[,"start"]),end=as.integer(test.exon[,"end"]))
            pa.iranges <- IRanges(start=as.integer(test.reads.mat[,"start"]),end=as.integer(test.reads.mat[,"end"]))
            test.exon.range <- GRanges(seqnames=te.chr.names,ranges=te.iranges)
            paired.read.range <- GRanges(seqnames=pa.chr.names,ranges=pa.iranges,reads=test.reads.mat[,"readsinfo"])
            exon.cover.range <- as.matrix(findOverlaps(test.exon.range,paired.read.range))
            meta.reads <- test.reads.mat[,"readsinfo"]
            if ((readsinfo == "paired") & length(exon.cover.range) != 0){
                paired.result <- tapply(exon.cover.range[,2],exon.cover.range[,1],function(x){
                    meta.reads[x]
                })
                names(paired.result) <- unique(exon.cover.range[,"queryHits"])
                covered.ex <- length(paired.result)
                total.names.re <- names(paired.result)
                rm.inter.reads <- NULL
                if (length(total.names.re) > 1){
                    all.comb <- combn(1:length(total.names.re),2)
                    each.pre.re <- lapply(1:ncol(all.comb),function(all.x.num){
                        pa.1.re <- paired.result[[as.character(all.comb[1,all.x.num])]]
                        pa.2.re <- paired.result[[as.character(all.comb[2,all.x.num])]]
                        intersect(pa.1.re,pa.2.re)
                    })
                    names(each.pre.re) <- apply(all.comb,2,function(x) paste(x,collapse="-"))
                    pre.cover.result <- lapply(total.names.re,function(each.pa.re){
                        next.test.na <- total.names.re[as.integer(total.names.re) > as.integer(each.pa.re)]
                        if (length(next.test.na) != 0){
                            each.paired.length <- lapply(next.test.na,function(next.pa.re){
                                pre.each.comb <- NULL
                                if (each.pa.re != "1"){
                                    pre.each.comb <- combn(1:as.integer(each.pa.re),2)
                                    pre.each.comb <- apply(pre.each.comb,2,function(x) paste(x,collapse="-"))
                                }
                                each.comb <- combn(each.pa.re:next.pa.re,2)
                                each.comb <- c(pre.each.comb,apply(each.comb,2,function(x) paste(x,collapse="-")))
                                ea.over <- is.element(each.comb,paste(each.pa.re,next.pa.re,sep="-"))
                                rm.over <- !is.element(each.comb,paste(each.pa.re,next.pa.re,sep="-"))
                                ea.ne.re <- unlist(each.pre.re[each.comb[ea.over]])
                                rm.ne.re <- unlist(each.pre.re[each.comb[rm.over]])
                                final.ne.re <- !is.element(ea.ne.re,rm.ne.re)
                                if (length(final.ne.re) != 0)        ea.ne.re[final.ne.re] 
                                else if (length(ea.ne.re) == 0) NULL         
                                else        final.ne.re
                            })
                            pa.each.pa.re <- paste(test.exon[as.integer(each.pa.re),],collapse="-")
                            ne.test.ex <- rbind(test.exon[as.integer(next.test.na),])
                            pa.next.pa.re <- apply(ne.test.ex,1,function(x) paste(x,collapse="-"))
                            names(each.paired.length) <- paste(pa.each.pa.re,pa.next.pa.re,sep="~")
                            each.paired.length
                        }
                        else NULL
                    })
                }
            }
            else if ((readsinfo == "exon") & length(exon.cover.range) != 0){
                tar.reads.count <- length(unique(test.reads.mat[exon.cover.range[exon.cover.range[,1] == 2,2],"readsinfo"]))
            }
        }
        number.ex <- 1:nrow(test.exon)
        total.ex.spli <- unlist(lapply(number.ex,function(num.ex){
            next.ex.num <- number.ex[number.ex > num.ex]
            if (length(next.ex.num) != 0){
                lapply(next.ex.num,function(n){
                    paste(paste(test.exon[num.ex,],collapse="-"),paste(test.exon[n,],collapse="-"),sep="~")
                })
            }
        }))
        used.reads.nm <- unlist(pre.cover.result)
        pair.length <- unlist(lapply(pre.cover.result,function(each.cover){
            lapply(each.cover,function(se.each.cover){
                length(se.each.cover)
            })
        }))
        pre.cover.result <- rbind(pair.length)
        not.names <- total.ex.spli[!is.element(total.ex.spli,colnames(pre.cover.result))]
        final.cover.result <- rbind(rep(0,length(not.names)))
        colnames(final.cover.result) <- not.names
        final.cover.result <- cbind(final.cover.result,pre.cover.result)
        final.cover.result <- rbind(final.cover.result[,order(colnames(final.cover.result))])
        exon.count <- 1:nrow(test.exon)
        test.exon.count <- intersect(names(paired.result),exon.count)
        pre.exon.result <- unlist(lapply(test.exon.count,function(ex.num){
            p.test.ex <- paste(test.exon[ex.num,],collapse="-")
            sum.cover <- sum(final.cover.result[,grep(p.test.ex,colnames(final.cover.result))])
            length(paired.result[[as.character(ex.num)]]) - sum.cover
        }))
        final.exon.result <- rep(0,nrow(test.exon))
        final.exon.result[test.exon.count] <- pre.exon.result
        final.exon.result[is.na(final.exon.result)] <- 0
        names(final.exon.result) <- paste(test.exon[,"start"],test.exon[,"end"],sep="-")
        if (readsinfo == "single" | readsinfo == "exon"){
            final.exon.result[2] <- tar.reads.count
        }
        final.N.result <- rep(0,nrow(test.exon))
    }
    final.N.result <- NULL
    if (length(N.num) != 0){
        N.reads <- unique(rbind(samfile[N.num,]))
        if (nrow(N.reads) != 0){
            u.reads.name <- table(as.matrix(paste(N.reads[,"qname"],N.reads[,"pos"],sep="_-_")))
            N.ex.p <- paste(N.reads[,"qname"],N.reads[,"pos"],sep="_-_")
            u.reads.name <- names(u.reads.name)[u.reads.name == 1]
            test.reads <- rbind(N.reads[is.element(N.ex.p,u.reads.name),])
            if (nrow(test.reads) != 0){
                final.N.result <- ReadCover(test.reads,"N")
                final.N.result <- table(paste(final.N.result[,"start"],final.N.result[,"end"],sep="-"))
            }
        }
    }
    final.mat <- list(final.cover.result[1,],final.exon.result,final.N.result,Am.result)
    names(final.mat) <- c("pairedInfo","exonInfo","junctionInfo","trimming")
    return (final.mat)
}

