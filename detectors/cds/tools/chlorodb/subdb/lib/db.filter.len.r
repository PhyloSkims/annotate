#!/usr/bin/env Rscript
#

args <- commandArgs(T)
path  <- if(length(args) > 0) args[1] else 'len.txt'
delta <- if(length(args) > 1) args[2] else 0.5

tab <- read.table(path, header=T)

lmed <- median(tab$len)

dlen <- lmed * as.numeric(delta)

tab$ok <- (abs(tab$len-lmed)/lmed) <= delta

write.table(tab, quote=F)

quit(save='no')

