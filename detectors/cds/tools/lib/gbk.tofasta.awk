#
# get fasta sequence from genbank
#

/^LOCUS/ {
  locus = $2
  next
}

/^ORIGIN/ {
  inseq = 1
  nln = 0
  delete seq
}

inseq && /^ +[1-9][0-9]*/ {
  s = substr($0, 11)
  gsub(" ", "", s)
  seq[++nln] = s
  next
}

/^\/\// {
  print ">" locus
  for (i = 1 ; i <= nln ; i++)
    print seq[i]
}

  



