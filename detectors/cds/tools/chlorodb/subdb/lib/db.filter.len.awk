#

function Abs(x) {
  return (x+0 >= 0) ? x+0 : 0-x
}

function Median(array, _local_, tmp, size, lhs, rhs) {
  size = asort(array, tmp)
  lhs = int((size - 1) / 2) + 1
  rhs = int(size / 2) + 1
  return ((tmp[lhs] + tmp[rhs]) / 2.0)
}

#

BEGIN {
  if (DELTA=="") DELTA=0.5
}

/^>/ {
  N++
  Id[N] = substr($1, 2)
  na = split($1, a, "@")
  Len[N] = a[na]
}

END {
  med = Median(Len)
  for (i = 1 ; i <= N ; i++) {
    if (Abs(Len[i]-med)/med <= DELTA)
      print Id[i]
  }
}
