#!/usr/bin/env Rscript
#
# compute start, stop, splice-junctions models for core DB
#
# source("make.models.r")
#

LIB_DIR <- Sys.getenv("LIB_DIR")
if (LIB_DIR == "") LIB_DIR = "."

source(paste0(LIB_DIR, "/util.base.r"))
source(paste0(LIB_DIR, "/util.cons.r"))
source(paste0(LIB_DIR, "/util.modelio.r"))

# -------------------------------
# parameters
# -------------------------------

# core cutoffs

source("db.models.params.txt")

# -------------------------------
# genome infos
# -------------------------------

notify("loading info table")
chromo <- read.table("db.info.txt", com="", head=T, stringsAsFactors=F)

# -------------------------------
# CDS
# -------------------------------

notify("loading cds table")
cds <- read.table("db.cds.txt", com="", header=T, stringsAsFactors=F)

cds$start <- as.factor(cds$start)
cds$stop <- as.factor(cds$stop)

cds <- cds[cds$status=="Ok",]
cds <- cds[cds$genefam!="none",]

cds$categ <- "dust"

x <- sort(table(cds$genefam), dec=T)
ok <- names(x[x >= CORE_NCDS_CUTOFF])

cds$categ[cds$genefam %in% ok] <- "core"

x <- x[! names(x) %in% ok]
ok <- names(x[x >= SHEL_NCDS_CUTOFF])

cds$categ[cds$genefam %in% ok] <- "shell"

#

cds.ori <- cds

cds.lst <- split(cds.ori, cds.ori$categ)

#
# write out families 
#

# patterns & names

invisible(lapply(cds.lst, function(cds) {

  x <- sort(table(cds$genefam), decreasing=T)
  tab <- paste0("^", names(x), "$")
  names(tab) <- names(x)

  y <- sapply(split(cds$gene, cds$genefam), function(g) {
    head(names(sort(table(g), decreasing=T)), 1)
  })
  
  tab <- cbind(tab, y[names(x)])

  categ <- unique(cds$categ)
  f <- paste0("db.", categ, ".pat.txt")
  notify("writing patterns for", categ, ":", f)
  write.table(tab, file=f, quote=F, col.names=F, row.names=T)
}))

# -------------------------------
# Start models (core only)
# -------------------------------

if (! "core" %in% names(cds.lst)) {
  notify("*** no gene found in core")
  notify("*** please change parameters")
  quit(save='no', status=1)
}

cds <- cds.lst[["core"]]

#
# start by genes
#

tab <- split(cds$start, cds$genefam)

fatg <- sapply(tab, function(x) table(x)["atg"]/length(x)*100)
names(fatg) <- names(tab)

start.dft <- names(which(fatg >= CORE_START_ATG_CUTOFF))
start.spc <- names(which(fatg < CORE_START_ATG_CUTOFF))

tab <- cds[cds$genefam %in% start.dft,]
tab <- table(tab$start)

# default model

x <- sort(tab[tab>=CORE_START_DFT_CUTOFF], decreasing=T)
write.model.start(x, "default")

# gene specific models

invisible(sapply(start.spc, function(g) {
  x <- cds[cds$genefam == g,]
  tx <- table(x$start)
  tx <- sort(tx[tx>=CORE_START_OTH_CUTOFF], decreasing=T)
  write.model.start(tx, g)
}))

# -------------------------------
# Stop models (core only)
# -------------------------------

# write default stop model

tab <- table(cds$stop)
x <- sort(tab[tab>=CORE_STOP_CUTOFF], decreasing=T)
write.model.stop(x, "default")

# -------------------------------
# splice junctions
# -------------------------------

notify("loading intron table")
intron <- read.table("db.intron.txt", com="", header=T, stringsAsFactors=F)

# remove invalid sequences

intron$seq <- gsub("\\.|-", "", intron$acceptor.donor)

lseq <- nchar(gsub("[^acgt]", "", intron$seq))

intron <- intron[lseq == 20,]

# remove genes out of core

intron <- intron[intron$genefam %in% cds$genefam,]

# acceptors / donors

intron$acc <- substr(intron$seq, 5, 6)
intron$don <- substr(intron$seq, 15, 16)

# consensus

cons.px <- cons.build(intron$acceptor.donor)
cons.px <- cons.px[,! is.nan(colSums(cons.px))]

seq.px <- sapply(intron$acceptor.donor, function(s) gsub("[^acgt]", "", s))

conf.px <- cons.confusion(cons.px, seq.px)

sfam <- split(conf.px$l2scor, intron$genefam)
sfam <- sfam[order(sapply(sfam, median))]

# extract splice exceptions

name.bad  <- names(which(sapply(sfam, median) < 0))
name.spc  <- names(which(sapply(sfam[name.bad], length) >= CORE_SPLICE_CUTOFF))
name.ok   <- setdiff(unique(intron$genefam), name.bad)
name.bad  <- setdiff(name.bad, name.spc)
name.list <- c(sapply(name.spc, function(x) x), list(default=name.ok))

cons <- lapply(name.list, function(x) cons.build(intron[intron$genefam %in% x, "acceptor.donor"]))

# write junction models

invisible(sapply(names(cons), function(n) write.model.splice3(cons[[n]], n)))
invisible(sapply(names(cons), function(n) write.model.splice5(cons[[n]], n)))

# use uniform model for bad guys

invisible(sapply(name.bad, function(n) write.unif.splice(3, n)))
invisible(sapply(name.bad, function(n) write.unif.splice(5, n)))
  
invisible(write.unif.splice('', "none"))

# -------------------------------
# keep data for plotting
# -------------------------------

DB <- list()

params <- list()

params$CORE_NCDS_CUTOFF      <- CORE_NCDS_CUTOFF
params$CORE_START_ATG_CUTOFF <- CORE_START_ATG_CUTOFF
params$CORE_START_DFT_CUTOFF <- CORE_START_DFT_CUTOFF
params$CORE_START_OTH_CUTOFF <- CORE_START_OTH_CUTOFF
params$CORE_STOP_CUTOFF      <- CORE_STOP_CUTOFF
params$CORE_SPLICE_CUTOFF    <- CORE_SPLICE_CUTOFF

params$SHEL_NCDS_CUTOFF      <- SHEL_NCDS_CUTOFF

DB$params  <- params
DB$chromo  <- chromo
DB$cds.lst <- cds.lst
DB$intron  <- intron
DB$cons    <- cons

notify("saving db.data.Rdata")
save(DB, file="db.data.Rdata")

# -------------------------------
# end
# -------------------------------

quit(save='no')
