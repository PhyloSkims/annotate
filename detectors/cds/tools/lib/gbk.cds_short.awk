#
# get cds features from genbank (short version)
#

# @include lib.gbk.awk


BEGIN {
  if (MAXSPAN == "") MAXSPAN = 10000
  print "#genefam gene from to strand nexon length status protseq product"
}

/^     CDS/ {
  revstrand = match($2, "^complement")
  s = substr($0, 22)
  gsub("^complement", "", s)
  ok = ! match(s, "complement|order")
  nexon = Nexons(s)
  SpanLocation(s, sloc)
  spanlen = sloc[2] - sloc[1] + 1
  len = LenLocation(s)
  ok = ok && (len < MAXSPAN)
  cdsseq = ok ? SeqLocation(seq, s, revstrand) : "XXX"
  cstart = substr(cdsseq, 1,3)
  cstop  = substr(cdsseq, length(cdsseq)-2)

  gene = "none"
  locustag = "none"
  product = "none"
  translation = "X"
  incds = 1
  next
}

(incds && /^     [^ ]/) {
  print GeneFamily(gene), gene, sloc[1], sloc[2], (revstrand ? "R" : "D"), 
        nexon, len, (ok ? "Ok" : "Error"), translation, product
  incds = 0
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

/^                     \/product=/ {
  split($0, a, "=")
  product = a[2]
  gsub("\"", "", product)
  gsub(" ", "_", product)
  next
}

/^                     \/translation=/ {
  split($0, a, "=")
  translation = a[2]
  gsub("\"", "", translation)
  gsub(" ", "", translation)
  next
}

/^\/\// {
  locus = "?"
}

END {
  if (incds) {
    print GeneFamily(gene), gene, sloc[1], sloc[2], (revstrand ? "R" : "D"), 
          nexon, len, (ok ? "Ok" : "Error"), translation, product
  }
}
