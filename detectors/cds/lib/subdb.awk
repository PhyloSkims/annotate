#
# select subDB from fasta DB
#
# -v FILE

BEGIN {
  if (FILE == "") FILE = "dbsel.txt"
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
