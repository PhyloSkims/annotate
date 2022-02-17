#
#
# filter BlastX results for speedup
#

BEGIN {
  if (IDMIN == "")  IDMIN  = 50
  if (NBMIN == "")  NBMIN  = 50
  if (NBMAX == "")  NBMAX  = 200
  if (COMIN == "")  COMIN  = 0.8  
}

/^#/ { next }

{
  if ((($3+0 >= IDMIN) || (HitNum <= NBMIN)) && (HitIndex[$2] == 0)) {
    HitIndex[$2] = ++HitNum
    IndexHit[HitNum] = $2
  }
}

END {
  n = (HitNum > NBMAX) ? NBMAX : HitNum
  for (i = 1 ; i <= n ; i++)
    print IndexHit[i]
}


