#
# get feature info from embl
#

# @include libgbk.awk

function GC(s, _local_, i, len) {
  s = toupper(s)
  len = length(s)
  gsub("G|C", "", s)
  return ((len - length(s)) * 100 / (len ? len : 1))
}

#
# rules
#

BEGIN {
  print "#locus orga len oklen gc nbCds nbCds_int0 nbCds_int1 nbCds_intsup1 perCds_noex meanCdsSize nbtRNA nbrRNA nboRNA"
}

/^ID/ {
  locus = $2
  gsub(";", "", locus)
  next
}

/^OS/ {
  orga = substr($0, 6)
  gsub(" ", "_", orga)
  next
}

/^FT   source/ {
  GetLoc($3, loc);
  len = loc[2];
  next
}

/^FT   CDS/ {
  meanCds = meanCds * nbCds + LenLocation($3)
  nbCds++
  meanCds /= nbCds
  n = Nexons($3)  
  if (n > 3) n = 3
  nbCdx[n]++
  next
}

/^FT   tRNA/ {
  nbTrna++
  next
}

/^FT   rRNA/ {
  nbRrna++
  next
}

/^FT   mRNA/ {
  next
}

/^FT   .*RNA/ {
  nbOrna++
  next
}

/^SQ / {
  inseq = 1
  seq = ""
  next
}

inseq && /^     / {
  s = $0
  gsub("[0-9]+", "", s)
  gsub(" ", "", s)
  seq = seq "" s
  next
}

/^\/\// {
  oklen = (len == length(seq) ? "ok" : "wrong")
  gc = GC(seq)
  print locus, orga, len, oklen, gc, nbCds+0, nbCdx[1]+0, \
         nbCdx[2]+0, nbCdx[3]+0, (nbCdx[1]+0)*100/Max(1, nbCds+0), \
         meanCds+0, nbTrna+0, nbRrna+0, nbOrna+0
  nbCds = nbTrna = nbRrna = nbOrna = len = inseq = meanCds = 0
  delete nbCdx
  orga = locus = "?"
}

  



