#
# embl oneLiner
#

/^FT / {
  InFeat = 1
}

(InFeat == 0) && ($1 != "CC") && ($1 == pkey) && /^..   [^ ]/ {
  line = line " " substr($0, 6)
  next
}

(InFeat == 1) && /^FT                   [^\/]/ {
  line = line "" substr($0, 22)
  next
}

{
  if (line != "") print line
  line = $0
  pkey = $1
  next
}

END {
  if (line != "") print line
}
