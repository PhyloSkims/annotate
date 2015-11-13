#
BEGIN {
  print "id len"
}

/^>/ {
 na = split($1, a, "@")
 print substr($1, 2), a[na]
}

