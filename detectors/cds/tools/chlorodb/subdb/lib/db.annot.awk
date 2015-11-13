#

/^>/ {
  N++
  na = split($1, a, "@")
  if (a[na-1] > NEXMAX) NEXMAX = a[na-1]
  NEX[a[na-1]]++
  ANNOT[$NF]++
}

END {
  na = split(FILENAME, a, "/")
  na = split(a[na], a, "\\.")
  printf("%s %d ", a[1], N)
  s = ""
  for (i = 1 ; i <= NEXMAX ; i ++) {
    if (NEX[i] != 0)
      s = s "" i ":" NEX[i] "_"
  }
  gsub("_+$", "", s)
  printf("%s ", s)
  
  s = (NEXMAX == 1) ? "MONEX" : "POLYEX"
  printf("%s ", s)
  
  nmax = 0
  amax = "none"
  for (e in ANNOT) {
    if (ANNOT[e] > nmax) {
      nmax = ANNOT[e]
      amax = e
    }
  }
  print amax
  
}



