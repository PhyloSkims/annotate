#

{
  node[$1]++
  node[$2]++
  link[++M] = $1 " " $2
}


END {
 for (n in node)
   N++
 print "DL n=" N
 print "format = edgelist1"
 print "labels embedded:"
 print "data:"
 for (i = 1 ; i <= M ; i++)
   print link[i]
}


