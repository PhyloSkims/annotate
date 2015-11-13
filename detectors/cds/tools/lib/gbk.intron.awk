#
# get intron features from genbank
#

# @include libgbk.awk

BEGIN {
  print "#locus locustag genefam gene from to strand intron_num intron_nb acceptor-donor status"

  if (HEADONLY != "") exit(0)
  
  if (FASTA == "") Error("No FASTA file specified", 1)
  
  if (! TestPath(FASTA)) Error("Fasta file: '" FASTA "' not found", 1)

  Seq = tolower(ReadFasta(FASTA))
}

/^LOCUS/ {
  locus = $2
  next
}

/^     CDS/ {
  revstrand = match($2, "^complement")
  s = substr($0, 22)
  gsub("^complement", "", s)
  ok = ! match(s, "complement|order")
  if (! ok) next

  na = ParseLocation(s, locs)
  if (na < 2) next

  delete SINfo
  Ninfo = 0
  
  val = locs[1][1]
  for (i = 2 ; i <= na ; i++) {
    if (locs[i][1] < val) ok = 0
    val = locs[i][1]
  }
  if (! ok) next

  val = locs[1][2]
  for (i = 2 ; i <= na ; i++) {
    if (locs[i][2] < val) ok = 0
    val = locs[i][2]
  }
  if (! ok) next

  from = locs[1][2] + 1
  for (i = 2 ; i <= na ; i++) {
    to = locs[i][1] - 1
    inseq = SeqLocation(Seq, (from - 4) ".." (to + 4), revstrand)
    SINfo[++Ninfo] = from " " to " " (revstrand ? "R" : "D") " "\
                     (revstrand ? na-i+1 : i-1) " " na-1 " "\
                     substr(inseq, 1,4) "."\
                     substr(inseq, 5,6) "-"\
                     substr(inseq, length(inseq)-9, 6) "."\
                     substr(inseq, length(inseq)-3, 4) " "\
                     "ok"
    from = locs[i][2] + 1
  }
  
  gene = "none"
  locustag = "none"
  next
}

/^                     \/gene=/ {
  split($0, a, "=")
  gene = a[2]
  gsub("^[^a-z,A-Z]+", "", gene)
  gsub("\"", "", gene)
  gsub(" ", "_", gene)
  next
}

/^                     \/locus_tag=/ {
  split($0, a, "=")
  locustag = a[2]
  gsub("\"", "", locustag)
  gsub(" ", "_", locustag)
  next
}

/^                     \/translation=/ {
  for (i = 1 ; i <= Ninfo ; i++)
    print locus, locustag, GeneFamily(gene), gene, SINfo[i]
  Ninfo = 0
  next
}

/^\/\// {
  locus = "?"
}
