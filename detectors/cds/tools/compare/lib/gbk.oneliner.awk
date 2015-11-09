#
# genbank oneLiner
#

/^            [^ ]/ {
  line = line "" substr($0, 12)
  next
}

/^                     [^\/]/ {
  line = line "" substr($0, 21)
  next
}

/^  ORGANISM/ {
  if (line != "") print line
  line = $0 ";"
  next
}

{
  if (line != "") print line
  line = $0
  next
}

END {
  print line
}
