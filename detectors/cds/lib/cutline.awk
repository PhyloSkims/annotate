#

{
  line = $0
  if (length(line) > 80) {
    print substr(line, 1, 80)
    rest = substr(line, 81)
    while (length(rest) > 59) {
      print "FT                   " substr(rest, 1, 59)
      rest = substr(rest, 60)
    }
    if (length(rest) > 0)
      print "FT                   " rest
  } 
  else {
    print line
  }
}
