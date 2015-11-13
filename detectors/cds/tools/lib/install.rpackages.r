#!/usr/bin/env Rscript
#
# check and install required packages
#

out <- function(...) {
  cat(paste0('+ ', ..., '\n'), file=stderr())
}

installed <- function(package) {
  package %in% rownames(installed.packages())
}

check <- function(package, repos="http://cran.univ-lyon1.fr") {
  if (installed(package)) {
    out("R package ", package, " installed")
  } else {
    out("Installing R package ", package, " from ", repos)
    install.packages(package, repos=repos)
  }
  invisible(installed(package))
}

check("grid")
check("gridExtra")
check("vcd")
check("plotrix")
check("igraph")

quit(save='no', status=0)

