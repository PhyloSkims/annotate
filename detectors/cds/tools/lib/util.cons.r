#
# R consensus utilities
#

#
# compute consensus
#

cons.build <- function(seqs, backcount=1) {
  xx <- do.call(rbind, sapply(seqs, strsplit, "", USE.NAMES=F))
  lv <- c("a", "c", "g", "t", ".", "-")
  mx <- apply(xx, 2, function(x) table(factor(x, levels=lv)))[1:4,]
  cx <- colSums(mx)
  mx <- mx + backcount
  mx[,cx==0] <- 0
  apply(mx, 2, function(x) x / sum(x))
}

#
# score consensus
#
  
cons.score <- function(cons, seq, pos=1:ncol(cons)) {
   seq <- strsplit(seq, "")[[1]]
   if (length(seq) != ncol(cons)) {
     warning("incompatible seq and cons size")
     return(NA)
   }
   ppx <- sapply(pos, function(i) cons[seq[i],i])
   sum(log10(ppx+1e-6))
}

#
# logratio to uniform model score
#

cons.logratio <- function(cons, seq, m0=NULL, pos=1:ncol(cons)) {
  if (is.null(m0)) {
    m0 <- matrix(rep(0.25, 4), nrow=4, ncol=ncol(cons))
    rownames(m0) <- c('a', 'c', 'g', 't')
  }
  
  sc <- cons.score(cons, seq, pos=pos)
  sc0 <- cons.score(m0, seq, pos=pos)
  
  2 * (log(10^sc, base=2) - log(10^sc0, base=2))
}

#
# shuffle sequence
#

seq.shuf <- function(seq) {
  paste0(sample(strsplit(seq, "")[[1]], nchar(seq), replace=F), collapse="")
}

#
# compute confusion matrix between actual and shuffled sequences
#

cons.confusion <- function(cons, seq, m0=NULL, pos=1:ncol(cons), thresh=0) {
  som <- function(x) sum(x, na.rm=T)

  res <- list()
  res$l2scor <- l2scor <- sapply(seq, function(s) cons.logratio(cons, s, m0=m0, pos=pos))

  seq <- sapply(seq, seq.shuf)
  res$r2scor <- r2scor <- sapply(seq, function(s) cons.logratio(cons, s, m0=m0, pos=pos))

  res$conf <- conf <- matrix(c(som(l2scor >= thresh), som(l2scor < thresh),
                               som(r2scor >= thresh), som(r2scor < thresh)),
                             nrow=2, byrow=T)

  res$acc <- sum(diag(conf)) / sum(conf)
  res$sen <- conf[1,1] / sum(conf[1,])
  res$sel <- conf[1,1] / sum(conf[,1])

  res  
}

#
# plot consensus
#


cons.plot <- function(cons, main="consensus") {
  cols <- c("blue", "orange", "red", "green")
  bp <- barplot(cons, col=cols, ylim=c(0,1), main=main)
  plx <- apply(cons, 2, function(col) -sum(col * log(col+1e-6, base=4)))
  lines(bp, plx, type="b", pch=19)
  legend(0, 1.1, c("a","c","g","t"), fill=cols, horiz=T, xpd=NA, bty="n")
  legend(20, 1.1, "entropy", pch=19, horiz=T, xpd=NA, bty="n")
  invisible()
}

#
# plot confusion scores histograms
#

cons.histconf <- function(conf, main="junction logr score") {
  lrh <- hist(c(conf$l2scor, conf$r2scor), br=50, plot=F)
  lh <- hist(conf$l2scor, br=lrh$breaks, plot=F)
  rh <- hist(conf$r2scor, br=lrh$breaks, plot=F)
  xx <- rbind(lh$counts, rh$counts) / sum(lh$counts)
  colnames(xx) <- lrh$mids
  barplot(xx, col=c(3,2), beside=T, main=main)
  legend(0, 0.1, c("true", "shuffled"), fill=c(3,2), horiz=F, xpd=NA)
  invisible()
}
