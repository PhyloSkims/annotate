#
# R models I/O utilities
#

#
# write start model
#

write.model.start <- function(frq, what) {
  dir.create("models", showWarnings=F)
  fil <- paste0("models/start.", what, ".frq")
  notify("writing start model:", fil)
  cat("# start model :", what, "\n", file=fil)
  for (x in names(frq))
    cat(x, frq[x]/sum(frq), frq[x], "\n", file=fil, append=T)
  invisible(fil)
}

#
# write stop model
#

write.model.stop <- function(frq, what) {
  dir.create("models", showWarnings=F)
  fil <- paste0("models/stop.", what, ".frq")
  notify("writing stop model:", fil)
  cat("# stop model :", what, "(freq. ignored)\n", file=fil)
  for (x in names(frq))
    cat(x, frq[x]/sum(frq), frq[x], "\n", file=fil, append=T)
  invisible(fil)
}

#
# write splice3 model
# [FIXME] positions are hard-coded
#

write.model.splice3 <- function(cons, what) {
  dir.create("models", showWarnings=F)
  fil <- paste0("models/splice3.", what, ".frq")
  notify("writing splice3 model:", fil)
  .catcons <- function(i) {
     cat(round(cons[c("a","c","g","t"), i]*100, 0), "\n",
         file=fil, append=T)
  }
  cat("# 3' splice model :", what, "\n", file=fil)
  cat("# A C G T\n", file=fil, append=T)
  sapply(seq.int(1, 4), .catcons)
  cat("splice\n", file=fil, append=T)
  sapply(seq.int(6, 11), .catcons)
  invisible(fil)
}

#
# write splice5 model
# [FIXME] positions are hard-coded
#

write.model.splice5 <- function(cons, what) {
  dir.create("models", showWarnings=F)
  fil <- paste0("models/splice5.", what, ".frq")
  notify("writing splice5 model:", fil)
  .catcons <- function(i) {
    cat(round(cons[c("a","c","g","t"), i]*100, 0), "\n",
        file=fil, append=T)
  }
  cat("# 5' splice model :", what, "\n", file=fil)
  cat("# A C G T\n", file=fil, append=T)
  sapply(seq.int(13, 18), .catcons)
  cat("splice\n", file=fil, append=T)
  sapply(seq.int(20, 23), .catcons)
  invisible(fil)
}

#
# write splice3/5 uniform model
#

write.unif.splice <- function(pos, what) {
  dir.create("models", showWarnings=F)
  fil <- paste0("models/splice", pos, ".", what, ".frq")
  notify("writing uniform splice", pos, "model:", fil)
  cat("# 3'/5' splice null model", file=fil)
  cat("# A C G T\n", file=fil, append=T)
  cat("25 25 25 25\n", file=fil, append=T)
  cat("splice\n", file=fil, append=T)
  cat("25 25 25 25\n", file=fil, append=T)
  invisible(fil)
}

