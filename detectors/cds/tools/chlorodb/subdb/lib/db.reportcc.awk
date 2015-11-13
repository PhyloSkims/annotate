#
#

{
  cnt[$NF]++
}

END {
  n = asort(cnt)
  printf("cc_size %s", NAME)
  for (i = n ; i >= 1 ; i--)
    printf(" %d", cnt[i])
  print ""
}

