#!/usr/bin/env Rscript
#
# plots models previously computed by make.models.r
#
# source("plot.models.r")
#

require(vcd)
require(plotrix)

LIBDIR <- Sys.getenv("LIB_DIR")
if (LIBDIR == "") LIBDIR = "."

source(paste0(LIBDIR, "/util.base.r"))
source(paste0(LIBDIR, "/util.plot.r"))
source(paste0(LIBDIR, "/util.cons.r"))
source(paste0(LIBDIR, "/util.grid.r"))

# -------------------------------
# setup
# -------------------------------

OUT.DEV <- TRUE
OUT.TYPE <- "pdf"
OUT.FILE <- "models"

if (OUT.DEV) uplot.init.dev(OUT.FILE, OUT.TYPE)

# -------------------------------
# Load data
# -------------------------------

notify("loading DB data")
load("db.data.Rdata")

params  <- DB$params
chromo  <- DB$chromo
cds.lst <- DB$cds.lst
intron  <- DB$intron
cons    <- DB$cons

# -------------------------------
# Genomes infos
# -------------------------------

grd.titlepage("Species")
grd.textpage(lineno=1, "# org: ", nrow(chromo))

#
# general stats
#

grd.hist(chromo, "len", main="Histogram of chromosome length")
grd.hist(chromo, "gc", pos.quant=c(0.75, 0.6), main="Histogram of chromosome GC")
grd.hist(chromo, "nbCds")
grd.fplot(chromo, "len", "nbCds")

#
# nb cds no introns
#

chromo$nbCds_Mono <- chromo$nbCds_int0
chromo$nbCds_Poly <- chromo$nbCds_int1 + chromo$nbCds_intsup1
chromo$percentPoly <- round(chromo$nbCds_Poly * 100 / (chromo$nbCds_Poly + chromo$nbCds_Mono))

grd.hist(chromo, "nbCds_Mono", main="Histogram of monoexonic Cds")
grd.hist(chromo, "nbCds_Poly", main="Histogram of polyexonic Cds")
grd.hist(chromo, "percentPoly", pos.sum=c(0.23,0.6), main="Histogram of % polyexonic")

grd.fplot(chromo, "nbCds", "nbCds_Mono", TRUE, ablin=list(a=0, b=1, col=3))
grd.fplot(chromo, "nbCds", "nbCds_Poly")

#
# cds size
#

grd.hist(chromo, "meanCdsSize", pos.quant=NULL, main="Histogram of Cds size")

# -------------------------------
# CDS
# -------------------------------

cds.all <- do.call(rbind, cds.lst)

grd.titlepage("CDS")

grd.textpage(lineno=1, "# core cds core cutoff: ", params$CORE_NCDS_CUTOFF)
grd.textpage(lineno=2, "# core cds shell cutoff: ", params$SHEL_NCDS_CUTOFF)

grd.textpage(lineno=4, "# total cds: ", nrow(cds.all))
grd.textpage(lineno=5, "# core cds: ", nrow(cds.lst[["core"]]))
grd.textpage(lineno=6, "# shell cds: ", nrow(cds.lst[["shell"]]))
grd.textpage(lineno=7, "# dust cds: ", nrow(cds.lst[["dust"]]))

grd.textpage(lineno=9,  "# total org: ", length(unique(cds.all$X.locus)))
grd.textpage(lineno=10, "# core org: ",  length(unique(cds.lst[["core"]]$X.locus)))
grd.textpage(lineno=11, "# shell org: ", length(unique(cds.lst[["shell"]]$X.locus)))
grd.textpage(lineno=12, "# dust org: ",  length(unique(cds.lst[["dust"]]$X.locus)))


grd.textpage(lineno=14, "# total families: ", length(unique(cds.all$genefam)))
grd.textpage(lineno=15, "# core families: ",  length(unique(cds.lst[["core"]]$genefam)))
grd.textpage(lineno=16, "# shell families: ", length(unique(cds.lst[["shell"]]$genefam)))
grd.textpage(lineno=17, "# dust families: ",  length(unique(cds.lst[["dust"]]$genefam)))

uplot.setup(mfrow=c(2,2), xpd=NA)

x <- sapply(cds.lst, nrow)
uplot.pie(x, main="CDS", text.r=1.1, col=c(3,2,4))

x <- sapply(cds.lst, function(cds) length(unique(cds$X.locus)))
uplot.pie(x, main="ORG", text.r=1.1, col=c(3,2,4))

x <- sapply(cds.lst, function(cds) length(unique(cds$genefam)))
uplot.pie(x, main="FAM", text.r=1.1, col=c(3,2,4))

uplot.setup(xpd=F)

#
# plot genes cutoff
#

cds.all <- do.call(rbind, cds.lst)
cds.byfam <- split(cds.all, cds.all$genefam)

tab <- sort(sapply(cds.byfam, nrow), decreasing=T)
cols <- rep("red", length(tab))
cols[tab >= params$SHEL_NCDS_CUTOFF] <- "blue"
cols[tab >= params$CORE_NCDS_CUTOFF] <- "green"
barplot(tab, col=cols, border=NA, main="# genes")

cols <- cols[tab >= 50]
tab <- tab[tab >= 50]
barplot(tab, col=cols, border=NA, las=2, cex.names=0.5, main="# genes in core")
abline(h=params$CORE_NCDS_CUTOFF, col=1)
text(50, 200, "CORE_NCDS_CUTOFF", pos=3)

#
# cds length for core
#

invisible(sapply(c("core", "shell", "dust"), function(what) {

  cds <- cds.lst[[what]]

  x <- split(cds$length, cds$genefam)
  x <- x[order(sapply(x, mean))]
  
  uplot.setup(mfrow=c(2,1))
  
  boxplot(x, pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5, outcex = 0.1),
          las=2, cex.axis=0.5, main=paste0(what, " genes - length distribution"))
  
  boxplot(x, pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5, outcex = 0.1), ylim=c(0,2000),
          las=2, cex.axis=0.5, main=paste0(what, " genes - length distribution zoom"))
  
  uplot.setup()
}))

# -------------------------------
# starts & stops
# -------------------------------

cds <- cds.lst[["core"]]

grd.titlepage("Starts and Stops")

tab <- sort(table(cds$start), dec=T)
tab <- tab[tab >= 100]
tab <- tab / sum(tab) * 100

barplot(tab, log="y", las=1, ylim=c(0.1, 100), main="start frequencies (%)")
text(0.7, 50, round(tab[1], 2))
text(1.9, 2.3, round(tab[2], 2))
text(3.1, 2.3, round(tab[3], 2))

#
# start by org and gc
#

x <- split(cds$start, cds$X.locus)

fatg <- sapply(x, function(x) table(x)["atg"]/length(x)*100)
names(fatg) <- names(x)
chromo$fatg <- round(fatg[chromo$X.locus], 2)

fgtg <- sapply(x, function(x) table(x)["gtg"]/length(x)*100)
names(fgtg) <- names(x)
chromo$fgtg <- round(fgtg[chromo$X.locus], 2)

facg <- sapply(x, function(x) table(x)["acg"]/length(x)*100)
names(facg) <- names(x)
chromo$facg <- round(facg[chromo$X.locus], 2)

grd.hist(chromo, "fatg", pos.quant=c(0.5, 0.6), main="Histogram of atg freq. by org")
grd.hist(chromo, "fgtg", main="Histogram of gtg freq. by org")
grd.hist(chromo, "facg", pos.sum=c(0.3, 0.6), main="Histogram of acg freq. by org")

grd.fplot(chromo, "gc", "fatg", main="atg freq. by org GC", pos=c(0.2, 0.3))
grd.fplot(chromo, "gc", "fgtg", main="gtg freq. by org GC")
grd.fplot(chromo, "gc", "facg", main="acg freq. by org GC")

ter <- cbind(fatg, fgtg, facg)
colnames(ter) <- c("ATG", "GTG", "ACG")
igc <- cut(chromo$gc, breaks=quantile(chromo$gc, seq(0, 1, 0.1)), include.lowest=T, labels=1:10)
cols <- rainbow(10)[igc]
ternaryplot(ter, col=cols, cex=0.2, main="Start by org", labels="outside")

#
# start by common genes
#

x <- split(cds$start, cds$genefam)

fatg <- sapply(x, function(x) table(x)["atg"]/length(x)*100)
names(fatg) <- names(x)

fgtg <- sapply(x, function(x) table(x)["gtg"]/length(x)*100)
names(fgtg) <- names(x)

facg <- sapply(x, function(x) table(x)["acg"]/length(x)*100)
names(facg) <- names(x)

barplot(sort(fatg)[1:10], las=2, main="atg freq. by gene")
barplot(sort(fgtg, dec=T)[1:10], las=2, main="gtg freq. by gene")
barplot(sort(facg, dec=T)[1:10], las=2, main="acg freq. by gene")

ter <- cbind(fatg, fgtg, facg)
colnames(ter) <- c("ATG", "GTG", "ACG")
ternaryplot(ter, col=1, cex=0.5, id=rownames(ter), main="Starts by genes", labels="outside")

# -------------------------------
# stops
# -------------------------------

#
# stop by org and gc
#

tab <- sort(table(cds$stop), dec=T)
tab <- tab[tab >= 100]
tab <- tab / sum(tab) * 100

barplot(tab, log="y", las=1, ylim=c(0.1, 100), main="stop frequencies (%)")
text(0.7, 80, round(tab[1], 2))
text(1.9, 30, round(tab[2], 2))
text(3.1, 25, round(tab[3], 2))

x <- split(cds$stop, cds$X.locus)

ftaa <- sapply(x, function(x) table(x)["taa"]/length(x)*100)
names(ftaa) <- names(x)
chromo$ftaa <- round(ftaa[chromo$X.locus], 2)

ftag <- sapply(x, function(x) table(x)["tag"]/length(x)*100)
names(ftag) <- names(x)
chromo$ftag <- round(ftag[chromo$X.locus], 2)

ftga <- sapply(x, function(x) table(x)["tga"]/length(x)*100)
names(ftga) <- names(x)
chromo$ftga <- round(ftga[chromo$X.locus], 2)

grd.hist(chromo, "ftaa", pos.quant=c(0.7, 0.6), main="Histogram of taa freq. by org")
grd.hist(chromo, "ftag", pos.quant=c(0.8, 0.6), main="Histogram of tag freq. by org")
grd.hist(chromo, "ftga", pos.quant=c(0.8, 0.6), main="Histogram of tga freq. by org")

grd.fplot(chromo, "gc", "ftaa", main="taa freq. by org GC", pos=c(0.2, 0.3))
grd.fplot(chromo, "gc", "ftag", main="tag freq. by org GC")
grd.fplot(chromo, "gc", "ftga", main="tga freq. by org GC")

ter <- cbind(ftaa, ftag, ftga)
colnames(ter) <- c("TAA", "TAG", "TGA")
igc <- cut(chromo$gc, breaks=quantile(chromo$gc, seq(0, 1, 0.1)), include.lowest=T, labels=1:10)
cols <- rainbow(10)[igc]
ternaryplot(ter, col=cols, cex=0.2, main="Stops by org", labels="outside")

#
# stop by common genes
#

x <- split(cds$stop, cds$genefam)

ftaa <- sapply(x, function(x) table(x)["taa"]/length(x)*100)
names(ftaa) <- names(x)

ftag <- sapply(x, function(x) table(x)["tag"]/length(x)*100)
names(ftag) <- names(x)

ftga <- sapply(x, function(x) table(x)["tga"]/length(x)*100)
names(ftga) <- names(x)

barplot(sort(ftaa), las=2, cex.names=0.5, ylim=c(0,100), main="taa freq. by gene")
barplot(sort(ftag), las=2, cex.names=0.5, ylim=c(0,100), main="tag freq. by gene")
barplot(sort(ftga), las=2, cex.names=0.5, ylim=c(0,100), main="tga freq. by gene")

ter <- cbind(ftaa, ftag, ftga)
colnames(ter) <- c("TAA", "TAG", "TGA")
ternaryplot(ter, col=1, cex=0.3, id=rownames(ter), main="Stops by genes", labels="outside")

# -------------------------------
# splice junctions
# -------------------------------

grd.titlepage("Splice Junctions")

grd.textpage(lineno=1, "# intron in core: ", nrow(intron))

#
# intron size
#

intron$size <- intron$to - intron$from + 1

grd.hist(intron, "size", pos.quant=NULL, main="Histogram of intron size", br=1000, xlim=c(0,2000))

#
# nb intron / gene
#

x <- split(intron, intron$genefam)
x <- x[order(sapply(x, function(x) mean(x$intron_nb)), decreasing=T)]

nmax <- max(intron$intron_nb)
lintron <- lapply(x, function(x) x$intron_nb)
mintron <- sapply(lintron, function(x) table(factor(x, levels=1:nmax)))

lintron0 <- table(cds[cds$nexon == 1,"genefam"])[names(lintron)]
mintron <- rbind("0"=lintron0, mintron)
mintron <- t(t(mintron)/colSums(mintron))

mintron[mintron==0] <- NA

nn <- nrow(mintron)
xx <- mintron[nn:1,]
ll <- lapply(1:nn, function(i) xx[i,])
mintron <- mintron[,do.call(order, c(ll, decreasing=T))]

battleship.plot(mintron, maxxspan=0.3, maxyspan=0.3, 
                cex.labels=0.7,
                main="% intron per polyexonic gene")

#
# acceptors / donors
#

tab <- sort(table(intron$acc), dec=T)
tab <- tab[tab >= 100]
tab <- tab / sum(tab) * 100

barplot(tab, log="y", las=1, ylim=c(0.1, 100), main="acceptor frequencies (%)")
text(0.7, 50, round(tab[1], 2))
text(1.9, 3, round(tab[2], 2))
text(3.1, 2.3, round(tab[3], 2))

tab <- sort(table(intron$don), dec=T)
tab <- tab[tab >= 100]
tab <- tab / sum(tab) * 100

barplot(tab, log="y", las=1, ylim=c(0.1, 100), main="donor frequencies (%)")
text(0.7, 50, round(tab[1], 1))
text(1.9, 40, round(tab[2], 1))
text(3.1, 15, round(tab[3], 1))
text(4.3, 12, round(tab[4], 1))

#
# consensus all
#

cons$all <- cons.build(intron$acceptor.donor)

invisible(sapply(rev(names(cons)), function(what) {
  cons.plot(cons[[what]], paste0("consensus ", what))
}))

#
# default consensus score by consensus length
#

cons.def <- cons[["default"]]
cons.def <- cons.def[,! is.nan(colSums(cons.def))]
seq.def <- sapply(intron$acceptor.donor, function(s) gsub("[^acgt]", "", s))

epx <- apply(cons.def, 2, function(col) -sum(col * log(col, base=4)))
opx <- order(epx)

conf.def <- lapply(seq(2, length(opx), by=2), function(n) {
         pos <- head(opx, n)
         notify(n, "/", length(opx))
         cons.confusion(cons.def, seq.def, thresh=0, pos=pos)
})

acc <- sapply(conf.def, function(x) x$acc)
sen <- sapply(conf.def, function(x) x$sen)
sel <- sapply(conf.def, function(x) x$sel)

plot(sel, ylim=c(0.7, 1), pch=1, type="b", main="accuracy by nb consensus positions", ylab="")
lines(sen, type="b", pch=2)
lines(acc, type="b", pch=3)
legend(1, 0.95, c("sensit.", "select.", "accur."), pch=1:3, horiz=T, bty="n")

#
# default consensus score by genes
#

conf.def <- cons.confusion(cons.def, seq.def, thresh=0)
cons.histconf(conf.def)

sfam <- split(conf.def$l2scor, intron$genefam)
sfam <- sfam[order(sapply(sfam, median))]

boxplot(sfam, pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5, outcex = 0.1),
        las=2, cex.axis=0.7, main="default junction logr score by genes")
abline(h=0)

#
# end
#

if (OUT.DEV) {
  cat("+ plot file:", paste0(OUT.FILE, ".", OUT.TYPE), "\n")
  invisible(dev.off())
}

quit(save='no')
