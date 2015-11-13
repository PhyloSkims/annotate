#
# get feature info from genbank
#

# @include libgbk.awk

function GC(s, _local_, i, len) {
  s = toupper(s)
  len = length(s)
  gsub("G|C", "", s)
  return ((len - length(s)) * 100 / len)
}

#
# rules
#

BEGIN {
  print "#locus orga len oklen gc nbCds nbCds_int0 nbCds_int1 nbCds_intsup1 perCds_noex meanCdsSize nbtRNA nbrRNA nboRNA"
}

/^LOCUS/ {
  locus = $2
  next
}

/^  ORGANISM/ {
  orga = substr($0, 13)
  split(orga, a, ";")
  orga = a[1]
  gsub(" ", "_", orga)
  next
}

/^     source/ {
  GetLoc($2, loc);
  len = loc[2];
  next
}

/^     CDS/ {
  meanCds = meanCds * nbCds + LenLocation($2)
  nbCds++
  meanCds /= nbCds
  n = Nexons($2)  
  if (n > 3) n = 3
  nbCdx[n]++
  next
}

/^     tRNA/ {
  nbTrna++
  next
}

/^     rRNA/ {
  nbRrna++
  next
}

/^     mRNA/ {
  next
}

/^     .*RNA/ {
  nbOrna++
  next
}

/^ORIGIN/ {
  inseq = 1
  seq = ""
  next
}

inseq && /^ +[1-9][0-9]*/ {
  s = substr($0, 11)
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

  



