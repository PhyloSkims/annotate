#

{
  N[$NF]++
  E[$NF, N[$NF]] = $1
}

END {
  cmax = 1
  nmax = N[1]
  for (i in N) {
    if (N[i] > nmax) {
      nmax = N[i]
      cmax = i
    }
  }
  for (i = 1 ; i <= nmax ; i++)
    print E[cmax, i]
}
