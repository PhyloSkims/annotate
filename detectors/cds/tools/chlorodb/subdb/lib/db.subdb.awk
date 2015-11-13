#

BEGIN {
  if (FILE == "") FILE = "db.sel.txt"
  while (getline < FILE)
    INC[$1] = $1
  close(FILE)
}

/^>/ {
  name = substr($1, 2)
  ok = name in INC
}

ok {
  print $0
}
