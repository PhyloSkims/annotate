#
# R plot utilities
#

#
# setup graphic device
# tiff: high resolution 600 dpi
# pdf
#

uplot.init.dev <- function(fname, type="pdf", width=7, height=7, resol=600, ...) {
  fname <- paste0(fname, ".", type)
  res <- NULL
  if (type == "tiff") {
    res <- tiff(fname, width=width, height=height, units="in", res=resol, ...)
  }
  if (type == "pdf") {
    res <- pdf(fname, width=width, height=height, ...)
  }
  invisible(res)
}

#
# convert pdf to tiff using ghostscript
#

uplot.convert2tiff <- function(fname, resol=600) {
  infile <- paste0(fname, ".pdf")
  oufile <- paste0(fname, ".tif")
  cmd <- paste0("echo quit | gs -r", resol, "-dBATCH -dNOPAUSE -sDEVICE=tiff12nc -sCompression=lzw -sOutputFile=", oufile, " ", infile)
  system(cmd)
}

#
# default plot setup
#

uplot.setup <- function(mfrow=c(1,1),
                       las=1,
                       mgp=c(2, 0.7, 0),
                       oma=c(0, 0, 0, 0),
                       mar=c(4, 3, 3, 2),
                       cex.main=1,
                       font.main=1,
                       family='Helvetica', ...) {
  par(mfrow=mfrow, las=las, mgp=mgp, oma=oma, mar=mar, cex.main=cex.main, font.main=font.main, family=family, ...)
}

#
# pie plot
#

uplot.pie <- function(tab, main="", labels=c("name", "val", "per"), text.r=0.5, text.col="black", text.cex=1, main.pos=c(0,0), main.col="black", ...) {
  pie(tab, edges=2000, main="", labels="", ...)
  text(main.pos[1], main.pos[2], main, cex=1.5, col=main.col)
  prop <- tab/sum(tab)
  theta <- 2*pi * (cumsum(prop) - prop/2)
  lab <- list(name=names(tab), val=tab, per=sprintf("%d%%", round(prop*100)))
  lab <- apply(data.frame(lab[labels]), 1, paste, collapse="\n") 
  if (length(lab) > 0)
    text(text.r*cos(theta), text.r*sin(theta), lab, cex=text.cex, col=text.col)
  invisible(NULL)
}

#
# plot utility : color representation of a table
#

uplot.table <- function(tab, col=heat.colors(100), with.lines=TRUE) {
  image(as.matrix(tab), xaxt="n", yaxt="n", col=col)
  nli <- nrow(tab)
  nco <- ncol(tab)
  dx  <- 0.5 / (nli-1)
  dy  <- 0.5 / (nco-1)
  xf <- (seq_len(nli)-1)/(nli-1) - dx
  yf <- (seq_len(nco)-1)/(nco-1) - dy
  if (with.lines) {
    segments(xf, -dy, xf, 1+dy)
    segments(-dx, yf, 1+dx, yf)
  }
  Axis(c(0,1), at=xf+dx, side=1, labels=rownames(tab), las=2, cex.axis=0.5, padj=0)
  Axis(c(0,1), at=yf+dy, side=2, labels=colnames(tab), las=2, cex.axis=0.5, padj=0)
  invisible(NULL)
}

#
# plot utility : identify points within user's rectangle
#

rect.identify <- function(data) {
  if (is.null(dim(data))) data <- cbind(seq_along(data), data)
  xy <- locator(n=2, type='n')
  r <- matrix(c(range(xy$x), range(xy$y)), ncol=2, byrow=TRUE)
  rect(r[1,1], r[2,1], r[1,2], r[2,2], border='red')
  .in.range <- function(p, r) {
    .in.int <- function(i) {
      (p[i] >= r[i,1]) && (p[i] <= r[i,2])
    }
    .in.int(1) && .in.int(2)
  }
  isel <- which(apply(data, 1, .in.range, r))
  points(data[isel,], col='red', pch=19)
  isel
}

