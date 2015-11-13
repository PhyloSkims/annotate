#!/usr/bin/env Rscript
#
# plot summary graphics of comparisons
#
#

LIBDIR <- Sys.getenv("LIB_DIR")
if (LIBDIR == "") LIBDIR = "."

source(paste0(LIBDIR, "/util.plot.r"))

COLORS <- 2:10

#

OUT.DEV <- TRUE
OUT.TYPE <- "pdf"
OUT.FILE <- "compare"
if (OUT.DEV) uplot.init.dev(OUT.FILE, OUT.TYPE)

#

tab <- read.table("compare.txt", header=T, comment.char="", stringsAsFactors=F)

#

par(xpd=NA)

#

sel  <- c("cor", "alcor", "acc", "wrong", "over", "misstot")
tab$ptot <- rowSums(tab[,sel])
for (s in sel)
  tab[,paste0("p", s)] <- tab[,s] * 100 / tab$ptot

colors <- head(COLORS, length(sel))

cols <- paste0("p", sel)
ord  <- order(tab$pcor+tab$palcor+tab$pacc, decreasing=T)

barplot(t(tab[ord,cols]), names.arg=tab$X.org[ord], 
        ylim=c(0,100), col=colors, las=2, cex.names=0.5)

legend(0, 110, sel, fill=colors, cex=0.7, horiz=T)

#

sel  <- c("cor", "alcor", "acc", "wrong", "over", "misschlo")
tab$rtot <- rowSums(tab[,sel])
for (s in sel)
  tab[,paste0("r", s)] <- tab[,s] * 100 / tab$rtot

colors <- head(COLORS, length(sel))

cols <- paste0("r", sel)
ord  <- order(tab$rcor+tab$ralcor+tab$racc, decreasing=T)

barplot(t(tab[ord,cols]), names.arg=tab$X.org[ord], 
        ylim=c(0,100), col=colors, las=2, cex.names=0.5)

legend(0, 110, sel, fill=colors, cex=0.7, horiz=T)

#

if (OUT.DEV) {
  cat("# plot file:", paste0(OUT.FILE, ".", OUT.TYPE), "\n")
  invisible(dev.off())
}

quit(save='no')

