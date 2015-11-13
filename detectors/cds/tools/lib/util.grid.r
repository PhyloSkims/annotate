#
# R misc grid plotting
#

require(grid)
require(gridExtra)

#
# get line height
#

grd.lineheight <- function(s="X") {
  convertHeight(unit(1,"strheight", s), "native", valueOnly=T)
}

#
# quantile table
#

grd.qtab <- function(df, what, cols, n=5) {
  df <- df[order(df[,what], decreasing=T),cols]
  sep <- head(df,1)
  sep[] <- "-"
  rbind(head(df, n), sep, tail(df, n))
}

#
# histogram with tables
#

grd.hist <- function(df, what, cols = c(1,2, which(colnames(df) == what)),
                      breaks=50, pos.sum=c(0.2,0.6), pos.quant=c(0.7,0.6), cex=0.7,
                      main=paste0("Histogram of ", what), ...) {
  hist(df[,what], breaks=breaks, xlab=what, main=main, ...)
  if (! is.null(pos.sum)) {
    pushViewport(viewport(pos.sum[1], pos.sum[2], gp=gpar(cex=cex)))
    grid.table(x<-summary(df[,what]), rows=names(x))
    popViewport()
  }
  if (! is.null(pos.quant)) {
    pushViewport(viewport(pos.quant[1], pos.quant[2], gp=gpar(cex=cex)))
    grid.table(grd.qtab(df, what, cols), rows=NULL)
    popViewport()
  }
  invisible()  
}

#
# plot with fit
#

grd.fplot <- function(df, what.x, what.y, linfit=T, pos=c(0.2, 0.8), ablin=NULL, ...) {
  plot(df[,what.x], df[,what.y], xlab=what.x, ylab=what.y, ...)
  if (linfit) {
    fit <- lm(df[,what.y] ~ df[,what.x])
    abline(fit, col=2)
    pushViewport(viewport(gp=gpar(col=2)))
    a <- sprintf("%.2e", coef(fit)[2])
    b <- sprintf("%.2e", coef(fit)[1])
    grid.text(paste0(what.y, " = ", a, " * ", what.x, " + ", b),
              pos[1], pos[2], just="left")
    pos[2] = pos[2] - 2 * grd.lineheight()
    grid.text(paste0("R2=", round(summary(fit)$r.squared, 3)),
              pos[1], pos[2], just="left")
    popViewport()
  }
  if (! is.null(ablin))
    do.call(abline, ablin)
  invisible()  
}

#
# write text
#

grd.textpage <- function(..., lineno=0, left=0.1, top=0.9, cex=1, fact=1.4) {
  txt <- do.call(paste0, list(...))
  pushViewport(viewport(gp=gpar(cex=cex)))
  grid.text(txt, left, top-lineno*grd.lineheight()*fact, just="left")
  popViewport()
  invisible(txt)
}

#
# title page
#

grd.titlepage <- function(title, x=0.5, y=0.7, cex=3, ...) {
  notify("processing", title)
  grid.newpage()
  grid.text(title, x, y, gp=gpar(cex=cex), ...)
  invisible()  
}
