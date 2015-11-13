#!/usr/bin/env Rscript
#

require(igraph, warn.conflicts=F)

args <- commandArgs(T)
path <- if(length(args) > 0) args[1] else 'graph.dl'

g <- read.graph(path, format='dl')

cc <- clusters(g)

res <- cbind(V(g)$name, membership(cc))

write.table(res, quote=FALSE, row.names=FALSE, col.names=FALSE)

quit(save="no")

