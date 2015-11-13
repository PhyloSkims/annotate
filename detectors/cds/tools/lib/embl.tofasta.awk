#
# get fasta sequence from embl
#

/^ID / {
  locus = $2
  gsub(";", "", locus)
  next
}

/^SQ / {
  inseq = 1
  nln = 0
  delete seq
  next
}

/^\/\// {
  inseq = 0
  print ">" locus
  for (i = 1 ; i <= nln ; i++)
    print seq[i]
  next
}

inseq {
  s = $0
  gsub(" ", "", s)
  gsub("[0-9]+", "", s)
  seq[++nln] = s
  next
}
  



